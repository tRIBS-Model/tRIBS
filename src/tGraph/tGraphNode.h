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
**  tGraphNode.h: Header for tBasin class and objects
**
**  tGraphNode Class used in tRIBS for the parallel version, providing 
**  information about upstream and downstream basins as graph nodes
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: tGraphNode Include and Define Statements
//
//
//=========================================================================

#ifndef TGRAPHNODE_H
#define TGRAPHNODE_H

#include <iostream>
#include <vector>
#include <set>

//=========================================================================
//
//
//                  Section 2: tGraphNode Class Definitions
//
//
//=========================================================================

class tGraphNode {

public:
  /// Default Constructor
  tGraphNode(); 
  /// Constructor that takes an ID
  tGraphNode(int nid);
  /// Constructor that takes an ID, upstream nodes, downstream nodes
  tGraphNode(int nid, std::vector<int>& nup, std::vector<int>& ndown);
  /// Destructor
  ~tGraphNode();

  /// Does another node flow into this node?
  bool hasUpstream() { if (upstream.size() > 0) return true;
                       else return false; }
  /// Does this node flow into another?
  bool hasDownstream() { if (downstream.size() > 0) return true;
                         else return false; }
  /// Does this node flow into another or get flowed into?
  bool hasFlux() { if (flux.size() > 0) return true;
                         else return false; }

  /// Return ID
  int getID() const { return id; }
  /// Return list of upstream nodes
  std::vector<int> getUpstream() const { return upstream; }
  /// Return list of downstream nodes
  std::vector<int> getDownstream() const { return downstream; }
  /// Return list of flux nodes
  std::vector<int> getFlux() const { return flux; }
  /// Return number of upstream nodes
  int getUpstreamCount() const { return upstream.size(); }
  // Return number of downstream nodes
  int getDownstreamCount() const { return downstream.size(); }
  /// Return number of flux nodes
  int getFluxCount() const { return flux.size(); }

  /// Add upstream node
  void addUpstream(int n) { upstream.push_back(n); }
  /// Add downstream node
  void addDownstream(int n) { downstream.push_back(n); }
  /// Add flux node
  void addFlux(int n) { flux.push_back(n); }
  /// Set ID
  void setId(int nid) { id = nid; }

private:

  int id;                        //!< Node ID
  std::vector<int> upstream;     //!< List of upstream nodes 
  std::vector<int> downstream;   //!< List of downstream nodes
  std::vector<int> flux;        //!< List of flux nodes
};

//! Write the tGraphNode to a given stream
std::ostream& operator<<(std::ostream& out, const tGraphNode& n);

#endif

//=========================================================================
//
//
//                          End of tGraphNode.h 
//
//
//=========================================================================
