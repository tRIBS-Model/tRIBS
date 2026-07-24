/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) 2025. tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/

/***************************************************************************
**
**  tPartition.cpp: In-process reach partitioning via libmetis.
**
**  This is a C++ reimplementation of the historical workflow:
**
**    connectivity2metis.pl  -> build a weighted graph of reaches
**    gpmetis                 -> METIS_PartGraphKway partitioning
**    metis2tribs.pl          -> emit the tRIBS .reach graphfile
**
**  The method -> (edges, weights) mapping is kept compatible
**  with the flags run_metis.zsh passed to connectivity2metis.pl:
**
**    SF    (1):  connectivity2metis ... 1 0 0 0 0  -> flow edges,  1 weight
**    SSF   (2):  connectivity2metis ... 1 0 0 0 2  -> flow+flux,   1 weight
**    SSF_H (3):  connectivity2metis ... 1 1 0 0 2  -> flow+flux,   2 weights
**                                       ^ ^         (pointCount, insideFlag)
**                                       | extra "inside" weight (SSF_H only)
**                                       node-count weight (always on)
**
***************************************************************************/

#include "src/tPartition/tPartition.h"

#include <metis.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace tPartition {

/*************************************************************************
**
**  tPartition::MethodToken
**
**  Returns the filename token for a partitioning method, matching the
**  naming used by the historical run_metis.zsh workflow.
**
**  Input:  method -- SF, SSF, or SSF_H
**  Output: "flow" (SF), "nconn" (SSF), or "upnconn" (SSF_H)
**
*************************************************************************/

const char* MethodToken(PartMethod method)
{
  switch (method) {
    case SF:    return "flow";
    case SSF:   return "nconn";
    case SSF_H: return "upnconn";
  }
  return "flow";
}

/*************************************************************************
**
**  tPartition::UsesFluxEdges
**
**  Whether a method includes subsurface flux edges in the reach graph.
**  SSF and SSF_H do; SF is surface-flow only.
**
**  Input:  method -- SF, SSF, or SSF_H
**  Output: true if subsurface flux edges should be added
**
*************************************************************************/

static bool UsesFluxEdges(PartMethod method)
{
  return method == SSF || method == SSF_H;
}

/*************************************************************************
**
**  tPartition::ReadConnectivity
**
**  Reads a connectivity.meshb file into a vector of reaches. Used by the
**  standalone RunPartition path; tRIBS itself builds the reach list from
**  its in-memory graph (tGraph::generatePartition) and does not need this.
**
**  Format (as written by MeshBuilder's bGraph::OutputConnectivity):
**    line 1: "#ReachID PointCount HeadNodeID OutletNodeID ..."  (header)
**    line 2: <numberOfReaches>
**    per reach: <id> <pointCount> <headID> <outletID>
**               <nDown> <d1..dNDown> <nFlux> <f1..fNFlux>
**
**  Input:  path -- path to a connectivity.meshb file
**  Output: reaches indexed by id (0..N-1)
**  Assumes: reach ids are contiguous 0..N-1
**
*************************************************************************/

std::vector<Reach> ReadConnectivity(const std::string& path)
{
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("tPartition: cannot open connectivity file '" +
                             path + "'");
  }

  std::string header;
  if (!std::getline(in, header)) {
    throw std::runtime_error("tPartition: '" + path + "' is empty");
  }

  int nReaches = 0;
  if (!(in >> nReaches) || nReaches < 0) {
    throw std::runtime_error("tPartition: bad reach count in '" + path + "'");
  }

  std::vector<Reach> reaches(nReaches);
  for (int i = 0; i < nReaches; ++i) {
    Reach r;
    int nDown = 0;
    if (!(in >> r.id >> r.pointCount >> r.headID >> r.outletID >> nDown)) {
      throw std::runtime_error("tPartition: truncated reach record in '" +
                               path + "' at line " + std::to_string(i));
    }
    if (r.id < 0 || r.id >= nReaches) {
      throw std::runtime_error("tPartition: reach id " + std::to_string(r.id) +
                               " out of range [0," + std::to_string(nReaches) +
                               ") -- reach ids must be contiguous");
    }
    r.downstream.resize(nDown);
    for (int j = 0; j < nDown; ++j) in >> r.downstream[j];

    int nFlux = 0;
    in >> nFlux;
    r.flux.resize(nFlux);
    for (int j = 0; j < nFlux; ++j) in >> r.flux[j];

    if (!in) {
      throw std::runtime_error("tPartition: parse error reading reach " +
                               std::to_string(r.id) + " in '" + path + "'");
    }
    reaches[r.id] = std::move(r);   // index by id, matching the perl scripts
  }
  return reaches;
}

/*************************************************************************
**
**  tPartition::WriteReachFile
**
**  Writes a tRIBS .reach graphfile from a partition array. Same layout as
**  the historical metis2tribs.pl: one "<partitionID> <reachID>" line per
**  reach, grouped by partition, reach ids ascending within a partition.
**
**  Input:  reaches -- the reach list (for point-count weights)
**          part    -- partition id per reach (as ComputePartition returns)
**          nParts  -- number of partitions
**          outPath -- path to write the .reach file
**  Output: writes outPath; prints per-partition node weights
**
*************************************************************************/

void WriteReachFile(const std::vector<Reach>& reaches,
                    const std::vector<int>& part,
                    int nParts,
                    const std::string& outPath)
{
  std::ofstream out(outPath);
  if (!out) {
    throw std::runtime_error("tPartition: cannot open output '" + outPath + "'");
  }

  const int effectiveParts = std::max(nParts, 1);
  for (int p = 0; p < effectiveParts; ++p) {
    long totalWeight = 0;
    for (std::size_t r = 0; r < part.size(); ++r) {
      if (part[r] == p) {
        out << p << " " << r << "\n";
        totalWeight += reaches[r].pointCount;
      }
    }
    std::cout << "Partition " << p << " weight = " << totalWeight << "\n";
  }
}

/*************************************************************************
**
**  tPartition::ComputePartition
**
**  Core partitioner. Builds a weighted, undirected reach graph in CSR form
**  and calls METIS_PartGraphKway. Flow edges (downstream links) are always
**  included; subsurface flux edges are added for SSF/SSF_H. Vertex weights
**  are reach node counts, with a second "inside/headwater" constraint for
**  SSF_H. This is the in-process replacement for connectivity2metis.pl +
**  gpmetis. nParts == 1 is handled trivially (METIS requires nParts >= 2).
**
**  Input:  reaches -- reach list with pointCount, downstream, flux
**          nParts  -- number of partitions (>= 1)
**          method  -- SF, SSF, or SSF_H
**  Output: partition id per reach (index == reach id)
**  Assumes: reach ids are contiguous 0..N-1
**
*************************************************************************/

std::vector<int> ComputePartition(const std::vector<Reach>& reaches,
                                  int nParts,
                                  PartMethod method)
{
  const idx_t nvtxs = static_cast<idx_t>(reaches.size());
  if (nvtxs == 0) {
    throw std::runtime_error("tPartition: no reaches to partition");
  }
  if (nParts < 1) {
    throw std::runtime_error("tPartition: nParts must be >= 1");
  }

  const bool useFlux = UsesFluxEdges(method);
  const idx_t ncon   = (method == SSF_H) ? 2 : 1;

  // Build undirected adjacency. std::set gives dedup + no duplicates,
  // matching the %seen collapse the perl did at output time.
  std::vector<std::set<idx_t>> adj(nvtxs);
  std::vector<idx_t> insideFlag(nvtxs, 0);   // old iweight (SSF_H 2nd weight)

  auto addEdge = [&](int a, int b) {
    if (a == b) return;                       // guard against self-loops
    if (a < 0 || a >= nvtxs || b < 0 || b >= nvtxs) {
      throw std::runtime_error("tPartition: neighbor id out of range");
    }
    adj[a].insert(b);
    adj[b].insert(a);
  };

  for (const Reach& r : reaches) {
    for (int d : r.downstream) {              // flow edges (all methods)
      addEdge(r.id, d);
      insideFlag[d] = 1;                      // d is downstream of something
    }
    if (useFlux) {
      for (int f : r.flux) addEdge(r.id, f);  // subsurface flux edges
    }
  }

  // Flatten to CSR
  std::vector<idx_t> xadj(nvtxs + 1);
  std::vector<idx_t> adjncy;
  adjncy.reserve(static_cast<std::size_t>(nvtxs) * 4);
  for (idx_t v = 0; v < nvtxs; ++v) {
    xadj[v] = static_cast<idx_t>(adjncy.size());
    for (idx_t nb : adj[v]) adjncy.push_back(nb);
  }
  xadj[nvtxs] = static_cast<idx_t>(adjncy.size());

  // Vertex weights, interleaved as METIS expects: vwgt[v*ncon + c]
  std::vector<idx_t> vwgt(static_cast<std::size_t>(nvtxs) * ncon);
  for (idx_t v = 0; v < nvtxs; ++v) {
    vwgt[v * ncon + 0] = static_cast<idx_t>(reaches[v].pointCount);
    if (ncon == 2) vwgt[v * ncon + 1] = insideFlag[v];
  }

  // Partition
  std::vector<idx_t> part(nvtxs, 0);
  idx_t objval = 0;

  if (nParts == 1) {
    // METIS requires nParts >= 2; the single-core case is trivial.
    std::fill(part.begin(), part.end(), 0);
  } else {
    idx_t ncon_i   = ncon;
    idx_t nparts_i = static_cast<idx_t>(nParts);

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;      // 0-based ids in/out

    idx_t nvtxs_i = nvtxs;
    int rc = METIS_PartGraphKway(
        &nvtxs_i, &ncon_i,
        xadj.data(), adjncy.data(),
        vwgt.data(), /*vsize*/ nullptr, /*adjwgt*/ nullptr,
        &nparts_i, /*tpwgts*/ nullptr, /*ubvec*/ nullptr,
        options, &objval, part.data());

    if (rc != METIS_OK) {
      throw std::runtime_error("tPartition: METIS_PartGraphKway failed (code " +
                               std::to_string(rc) + ")");
    }
  }

  std::cout << "tPartition: partitioned " << nvtxs << " reaches into "
            << nParts << " (edge-cut " << objval << ")\n";

  // Convert METIS idx_t -> int for the caller (tRIBS uses int partitions).
  return std::vector<int>(part.begin(), part.end());
}

/*************************************************************************
**
**  tPartition::PartitionReaches
**
**  Convenience wrapper: ComputePartition followed by WriteReachFile.
**
**  Input:  reaches -- reach list
**          nParts  -- number of partitions
**          method  -- SF, SSF, or SSF_H
**          reachOutPath -- path to write the .reach file
**  Output: writes the .reach file; returns 0
**
*************************************************************************/

int PartitionReaches(const std::vector<Reach>& reaches,
                     int nParts,
                     PartMethod method,
                     const std::string& reachOutPath)
{
  std::vector<int> part = ComputePartition(reaches, nParts, method);
  WriteReachFile(reaches, part, nParts, reachOutPath);
  std::cout << "tPartition: wrote " << reachOutPath << "\n";
  return 0;
}

/*************************************************************************
**
**  tPartition::RunPartition
**
**  Standalone driver: read a connectivity.meshb file, partition it, and
**  write the .reach file using the run_metis.zsh output-naming scheme
**  (<basename>_<token>_<nParts>nodes.reach).
**
**  Input:  connectivityPath -- path to connectivity.meshb
**          nParts   -- number of partitions
**          method   -- SF, SSF, or SSF_H
**          basename -- base name for the output file
**  Output: outPath -- (returned) path of the written .reach file
**  Returns: 0 on success (throws std::runtime_error on failure)
**
*************************************************************************/

int RunPartition(const std::string& connectivityPath,
                 int nParts,
                 PartMethod method,
                 const std::string& basename,
                 std::string& outPath)
{
  std::vector<Reach> reaches = ReadConnectivity(connectivityPath);

  std::ostringstream name;
  name << basename << "_" << MethodToken(method) << "_" << nParts << "nodes.reach";
  outPath = name.str();

  return PartitionReaches(reaches, nParts, method, outPath);
}

} // namespace tPartition

//=========================================================================
//
//                          End of tPartition.cpp
//
//=========================================================================
