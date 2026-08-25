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
**  tRestart.h: Header for tRestart class and objects
**
**  tRestart Class used in tRIBS for saving current state of a simulation
**  with the purpose of restarting at a later time
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: tRestart Include and Define Statements
//
//
//=========================================================================

#ifndef TRESTART_H
#define TRESTART_H

#include <iostream>
#include "src/tSimulator/tRunTimer.h"
#include "src/tFlowNet/tKinemat.h"
#include "src/tFlowNet/tReservoir.h" // JECR2015
#include "src/tFlowNet/tResData.h" // JECR2015
#include "src/tHydro/tSnowPack.h"

//=========================================================================
//
//
//                  Section 2: tRestart Class Definitions
//
//
//=========================================================================

template< class tSubNode >
class tRestart {

public:
  tRestart(
        tRunTimer* t,
        tMesh<tSubNode>* m,
        tKinemat* f,
        tSnowPack* s);

  ~tRestart() {}

  void writeRestart(ostream &);
  void writeRestartOutlets(ostream &);
  void writeRestartNodes(ostream &);
  void readRestart(istream &);
  void readRestartGlobal(istream &, int32_t, int32_t);

  int getSnowOpt() const;
  int getNumNodes() const;
  int getNumOutlets() const;

private:
  tRunTimer*         timer;           //!< Run timer
  tMesh<tSubNode>*   mesh;            //!< Mesh
  tKinemat*          flow;            //!< Kinematic flow
  tSnowPack*         snowpack;        //!< Snow pack structure
};

#endif

//=========================================================================
//
//
//                          End of tRestart.h 
//
//
//=========================================================================
