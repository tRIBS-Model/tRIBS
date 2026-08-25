/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) tRIBS Developers
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
**  tGraph::generatePartition builds the Reach list from the in-memory
**  mesh, calls ComputePartition, and writes the result via WriteReachFile
**  as <basename>_<method>_<nParts>nodes.reach (the tRIBS GRAPHFILE).
**
***************************************************************************/

#ifndef TPARTITION_H
#define TPARTITION_H

#include <string>
#include <vector>

namespace tPartition {

/// Partitioning strategy. Values match the GRAPHOPTION input keyword.
enum PartMethod {
  SF    = 0,   //!< Surface flow only                      (flow edges)
  SSF   = 1,   //!< Surface-subsurface flow                (flow + flux edges)
  SSF_H = 2    //!< Surface-subsurface flow with headwaters(flow + flux, 2 weights)
};

/// One reach of the stream network graph.
struct Reach {
  int id         = 0;   //!< Reach id (assumed contiguous 0..N-1)
  int pointCount = 0;   //!< Number of nodes in the reach (vertex weight)
  int headID     = 0;   //!< Head node id
  int outletID   = 0;   //!< Outlet node id
  std::vector<int> downstream;  //!< Downstream reach ids (flow edges)
  std::vector<int> flux;        //!< Flux-coupled reach ids (subsurface edges)
};

/// Filename token for a method: "SF", "SSF", or "SSFH". Used to build the
/// generated graphfile name <basename>_<token>_<nParts>nodes.reach.
const char* MethodToken(PartMethod method);

/// Core partitioner: build the weighted reach graph and call METIS.
/// Returns the partition id per reach (index == reach id). No file I/O and
/// no console output -- callers report results themselves, so under MPI the
/// message comes from the master rank only.
/// If edgeCut is non-null it receives the METIS edge-cut: the number of reach
/// couplings whose two reaches landed in different partitions.
/// This is the piece tRIBS calls directly on its in-memory reach graph.
/// Throws std::runtime_error on failure.
std::vector<int> ComputePartition(const std::vector<Reach>& reaches,
                                  int nParts,
                                  PartMethod method,
                                  int* edgeCut = 0);

/// Write a tRIBS .reach graphfile from a partition array (as ComputePartition
/// returns). Format matches metis2tribs.pl: "<partitionID> <reachID>" per line.
void WriteReachFile(const std::vector<Reach>& reaches,
                    const std::vector<int>& part,
                    int nParts,
                    const std::string& reachOutPath);

} // namespace tPartition

#endif // TPARTITION_H
