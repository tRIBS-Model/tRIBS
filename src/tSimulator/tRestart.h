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
  /// Constructor
  tRestart(
        tRunTimer* t,
        tMesh<tSubNode>* m,
        tKinemat* f,
        tWaterBalance* b,
        tHydroModel* h,
        tRainfall* r,
        tEvapoTrans* e,
        tIntercept* i,
        tSnowPack* s);

  /// Destructor
  ~tRestart() {}

  /// Write a restart dump
  void writeRestart(fstream &);
  /// Read a restart dump
  void readRestart(fstream &);

private:
  tRunTimer*         timer;           //!< Run timer
  tMesh<tSubNode>*   mesh;            //!< Mesh
  tKinemat*          flow;            //!< Kinematic flow
  tWaterBalance*     balance;         //!< Water balance
  tHydroModel*       hydro;           //!< Hydro model
  tRainfall*         rainfall;        //!< Rainfall
  tEvapoTrans*       evap;            //!< Evapotranspiration
  tIntercept*        intercept;       //!< Intercept structure
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
