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
**  tWaterBalance.h:   Header file for tWaterBalance Class
**
**  Storage/computation class for the components of the water balance
**  in each node and for the entire watershed.
**
***************************************************************************/

#ifndef TWATERBALANCE_H 
#define TWATERBALANCE_H

//=========================================================================
//
//
//                  Section 1: tWaterBalance Include Statements
//
//
//=========================================================================

#include "src/Headers/Inclusions.h"

//=========================================================================
//
//
//                  Section 2: tWaterBalance Class Definitions
//
//
//=========================================================================

class tWaterBalance
{
public:
  tWaterBalance();
  tWaterBalance(SimulationControl*,tMesh<tCNode> *,tInputFile &);
  ~tWaterBalance();

  void initializeVariables();
  void SetWaterBalance(tInputFile &);
  void DeleteWaterBalance();
  void CanopyBalance();
  void UnSaturatedBalance();
  void SaturatedBalance();
  void BasinStorage(double);
  void Print(double *);
  void writeRestart(fstream &) const;
  void readRestart(fstream &);

protected:
  tMesh<tCNode> *gridPtr;      
  SimulationControl *simCtrl; 

  int finalTime;
  double metStep, unsStep, satStep; //WR 12192023: rounding errors converting from double in .in to int here
  double *BasinStorages;
};

#endif

//=========================================================================
//
//
//                       End of tWaterBalance.h
//
//
//=========================================================================
