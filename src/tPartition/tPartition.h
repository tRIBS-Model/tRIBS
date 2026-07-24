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
**  tPartition.h: In-process reach partitioning via libmetis.
**
**  Replaces the external MeshBuilder->METIS workflow, i.e. the
**  connectivity2metis.pl + gpmetis + metis2tribs.pl chain driven by
**  run_metis.zsh, with a single call into the linked-in METIS library.
**
**  Input : connectivity.meshb  (written by bGraph::OutputConnectivity)
**  Output: <basename>_<method>_<nParts>nodes.reach  (tRIBS GRAPHFILE)
**
***************************************************************************/

#ifndef TPARTITION_H
#define TPARTITION_H

#include <string>
#include <vector>

namespace tPartition {

/// Partitioning strategy. Values match the OPT_Part argument historically
/// passed to run_metis.zsh so existing user muscle-memory carries over.
enum PartMethod {
  SF    = 1,   //!< Surface flow only                      (flow edges)
  SSF   = 2,   //!< Surface-subsurface flow                (flow + flux edges)
  SSF_H = 3    //!< Surface-subsurface flow with headwaters(flow + flux, 2 weights)
};

/// One reach as described by a line of connectivity.meshb.
struct Reach {
  int id         = 0;   //!< Reach id (assumed contiguous 0..N-1)
  int pointCount = 0;   //!< Number of nodes in the reach (vertex weight)
  int headID     = 0;   //!< Head node id
  int outletID   = 0;   //!< Outlet node id
  std::vector<int> downstream;  //!< Downstream reach ids (flow edges)
  std::vector<int> flux;        //!< Flux-coupled reach ids (subsurface edges)
};

/// Human-readable / filename token for a method: "flow", "nconn", "upnconn".
/// Mirrors the naming used by run_metis.zsh.
const char* MethodToken(PartMethod method);

/// Parse connectivity.meshb. Throws std::runtime_error on failure.
std::vector<Reach> ReadConnectivity(const std::string& path);

/// Core partitioner: build the weighted reach graph and call METIS.
/// Returns the partition id per reach (index == reach id). No file I/O.
/// This is the piece tRIBS calls directly on its in-memory reach graph.
/// Throws std::runtime_error on failure.
std::vector<int> ComputePartition(const std::vector<Reach>& reaches,
                                  int nParts,
                                  PartMethod method);

/// Write a tRIBS .reach graphfile from a partition array (as ComputePartition
/// returns). Format matches metis2tribs.pl: "<partitionID> <reachID>" per line.
void WriteReachFile(const std::vector<Reach>& reaches,
                    const std::vector<int>& part,
                    int nParts,
                    const std::string& reachOutPath);

/// Partition reaches across nParts and write the tRIBS .reach graphfile.
/// Thin wrapper: ComputePartition + WriteReachFile.
/// Returns the METIS objective (edge-cut). Throws std::runtime_error on failure.
int PartitionReaches(const std::vector<Reach>& reaches,
                     int nParts,
                     PartMethod method,
                     const std::string& reachOutPath);

/// Convenience driver: read connectivity file, partition, write .reach using
/// the run_metis.zsh output-naming scheme. Returns the output path via outPath.
/// Throws std::runtime_error on failure.
int RunPartition(const std::string& connectivityPath,
                 int nParts,
                 PartMethod method,
                 const std::string& basename,
                 std::string& outPath);

} // namespace tPartition

#endif // TPARTITION_H
