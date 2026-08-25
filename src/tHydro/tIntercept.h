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
**  tIntercept.h: Header file for class tIntercept
**
**  This class encapsulates the rainfall interception routines
**  necessary for net precipitation computations. Interception is computed
**  using the Rutter et al (1971) method.
**
***************************************************************************/

#ifndef INTERCEPT_H
#define INTERCEPT_H

#include "src/Headers/Inclusions.h"

//=========================================================================
//
//
//            Section 1: tIntercept Class Declarations
//
//
//=========================================================================

class tIntercept
{
 public:
  tIntercept();
  tIntercept(SimulationControl*,tMesh<tCNode>*, tInputFile&,  tRunTimer*,
	     tResample*, tHydroModel*);
  ~tIntercept();

  void DeleteIntercept(); 
  int    getIoption();
  int    IsThereCanopy(tCNode *);
  void   callInterception(tCNode *, double);
  void   InterceptRutter(tCNode *, double);
  void   SetIntercpParameters(tCNode *);
  void   SetIntercpVariables(tInputFile&, tHydroModel*);
  double storageRungeKutta(double, double, double, double *);
  double RutterFn(double, double, double, double);
  double getCtoS(tCNode *);

  void readLUGrid(char*); // SKYnGM2008LU

  SimulationControl *simCtrl; 

 protected:
  tMesh<tCNode>   *gridPtr;
  GenericLandData *landPtr;
  tRunTimer       *timer;       
  int interceptOption;

  int luOption; // SKYnGM2008LU: added by AJR 2007
  
  // SKYnGM2008LU
  int nParmLU; 
  char **LUgridParamNames, **LUgridBaseNames, **LUgridExtNames; 
  char luFile[kName]; 

  int maxInterStormPeriod;
  double metTime;
  double coeffP, coeffS, coeffK, coeffb; //Rutter model
  double coeffV;
};

#endif

//=========================================================================
//
//
//                         End of tIntercept.h
//
//
//=========================================================================
