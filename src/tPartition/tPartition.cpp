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
**    gpmetis                -> METIS_PartGraphKway partitioning
**    metis2tribs.pl         -> write the tRIBS .reach graphfile
**
**  The method -> (edges, weights) mapping is kept compatible
**  with the flags run_metis.zsh passed to connectivity2metis.pl:
**
**    SF    (0):  connectivity2metis ... 1 0 0 0 0  -> flow edges,  1 weight
**    SSF   (1):  connectivity2metis ... 1 0 0 0 2  -> flow+flux,   1 weight
**    SSF_H (2):  connectivity2metis ... 1 1 0 0 2  -> flow+flux,   2 weights
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
#include <stdexcept>
#include <string>
#include <vector>

namespace tPartition {

/*************************************************************************
**
**  tPartition::MethodToken
**
**  Returns the filename token for a partitioning method, used to build the
**  generated graphfile name <basename>_<token>_<nParts>nodes.reach.
**
**  Input:  method -- SF, SSF, or SSF_H
**  Output: "SF", "SSF", or "SSFH"
**
*************************************************************************/

const char* MethodToken(PartMethod method)
{
  switch (method) {
    case SF:    return "SF";
    case SSF:   return "SSF";
    case SSF_H: return "SSFH";
  }
  return "SF";
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
**  Output: writes outPath; prints per-partition node counts
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
    std::cout << "Partition " << p << ": " << totalWeight << " nodes\n";
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
**          edgeCut -- if non-null, receives the METIS edge-cut
**  Output: partition id per reach (index == reach id)
**  Assumes: reach ids are contiguous 0..N-1
**
**  Prints nothing: every MPI rank calls this, so reporting here would
**  duplicate the message once per rank. The caller reports via Cout.
**
*************************************************************************/

std::vector<int> ComputePartition(const std::vector<Reach>& reaches,
                                  int nParts,
                                  PartMethod method,
                                  int* edgeCut)
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

  if (edgeCut != 0) *edgeCut = static_cast<int>(objval);

  // Convert METIS idx_t -> int for the caller (tRIBS uses int partitions).
  return std::vector<int>(part.begin(), part.end());
}

} // namespace tPartition

//=========================================================================
//
//                          End of tPartition.cpp
//
//=========================================================================
