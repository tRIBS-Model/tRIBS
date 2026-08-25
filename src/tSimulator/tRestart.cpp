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
**  tRestart.cpp: Functions for class tRestart (see tRestart.h)
**
***************************************************************************/

#include <cassert>
#include <cstdint>

#include "src/tSimulator/tRestart.h"

/*************************************************************************
**
** Constructor
**
*************************************************************************/

template< class tSubNode >
tRestart<tSubNode>::tRestart(
	tRunTimer* t,
	tMesh<tSubNode>* m,
	tKinemat* f,
	tSnowPack* s)
{
  timer = t;
  mesh = m;
  flow = f;
  snowpack = s;
}

template< class tSubNode >
int tRestart<tSubNode>::getSnowOpt() const
{
  return snowpack->getSnowOpt();
}

template< class tSubNode >
int tRestart<tSubNode>::getNumNodes() const
{
  // Must return the active (non-boundary) node count, not getSize(): the
  // restart write path (tMesh::writeRestart) serializes only IsActive() nodes,
  // so this count is both the header record count and the read-loop bound in
  // tMesh::readRestartGlobal. Using getSize() (which includes boundary nodes)
  // would advertise more records than are actually written.
  return mesh->getNodeList()->getActiveSize();
}

template< class tSubNode >
int tRestart<tSubNode>::getNumOutlets() const
{
  return flow->getReachOutletList().getSize();
}

template< class tSubNode >
void tRestart<tSubNode>::writeRestart(ostream & rStr)
{
  flow->writeRestart(rStr);
  mesh->writeRestart(rStr);
}

template< class tSubNode >
void tRestart<tSubNode>::writeRestartOutlets(ostream & rStr)
{
  flow->writeRestart(rStr);
}

template< class tSubNode >
void tRestart<tSubNode>::writeRestartNodes(ostream & rStr)
{
  mesh->writeRestart(rStr);
}

template< class tSubNode >
void tRestart<tSubNode>::readRestart(istream & rStr)
{
  flow->readRestart(rStr);
  mesh->readRestart(rStr);
}

template< class tSubNode >
void tRestart<tSubNode>::readRestartGlobal(istream & rStr, int32_t totalOutlets, int32_t totalNodes)
{
  flow->readRestartGlobal(rStr, static_cast<int>(totalOutlets));
  mesh->readRestartGlobal(rStr, static_cast<int>(totalNodes));
}

//=========================================================================
//
//
//                        End of tRestart.cpp
//
//
//=========================================================================
