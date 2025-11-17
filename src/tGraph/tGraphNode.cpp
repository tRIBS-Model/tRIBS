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
**  tGraphNode.cpp: Functions for class tGraphNode (see tGraphNode.h)
**
***************************************************************************/

#include "src/tGraph/tGraphNode.h"


/**************************************************************************
**
** Constructor that only creates object.
**
**************************************************************************/

tGraphNode::tGraphNode() : id(-1) {
}

/*************************************************************************
**
** Constructor that takes an id.
**
*************************************************************************/

tGraphNode::tGraphNode(int nid) : id(nid) {
}

/*************************************************************************
**
** Constructor that takes an id, upstream nodes, downstream nodes.
**
*************************************************************************/

tGraphNode::tGraphNode(int nid, std::vector<int>& nup, 
  std::vector<int>& ndown) : id(nid) {

  for (int i = 0; i < nup.size(); i++) {
    upstream.push_back(nup[i]);
  }
  for (int i = 0; i < ndown.size(); i++) {
    upstream.push_back(ndown[i]);
  }
}

/*************************************************************************
**
** Destructor
**
*************************************************************************/

tGraphNode::~tGraphNode(){
  id = -1;
  upstream.erase(upstream.begin(), upstream.end());
  downstream.erase(downstream.begin(), downstream.end());
  flux.erase(flux.begin(), flux.end());
} 

/*************************************************************************
**
** Write out tGraphNode
**
*************************************************************************/

std::ostream& operator<<(std::ostream& out, const tGraphNode& n) {
  out << "tGraphNode:" << std::endl;
  out << "  id = " << n.getID() << std::endl;
  std::vector<int> down = n.getDownstream();
  std::vector<int> up = n.getUpstream();
  std::vector<int> flux = n.getFlux();
  
  out << "  upstream = " << std::endl;
  for (int i = 0; i < up.size(); i++)
    out << up[i] << std::endl;
  out << "  downstream = " << std::endl;
  for (int i = 0; i < down.size(); i++) 
    out << down[i] << std::endl;
  out << "  flux = " << std::endl;
  for (int i = 0; i < flux.size(); i++) {
    out <<  flux[i] << std::endl;
  }
  out << std::endl;

  return out;
}

//=========================================================================
//
//
//                        End of tGraphNode.cpp
//
//
//=========================================================================
