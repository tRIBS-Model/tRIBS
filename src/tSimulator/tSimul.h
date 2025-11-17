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
**  tSimul.h: Header for tSimul.cpp 
**  	
**  tSimulator Class is a master class for the control functions of the
**  model simulation. Initialized and called in main.cpp. Used to control
**  calls to all tRIBS objects.
**
***************************************************************************/

#ifndef TSIMUL_H
#define TSIMUL_H

//=========================================================================
//
//
//                  Section 1: tSimul Include and Define Statements
//
//
//=========================================================================

#include "src/tSimulator/tControl.h"
#include "src/tSimulator/tRestart.h"
#include "src/tFlowNet/tKinemat.h"
#include "src/tFlowNet/tReservoir.h" // JECR2015
#include "src/tHydro/tSnowPack.h" // SKY2008Snow from AJR2007
#include "src/Headers/Inclusions.h"

//=========================================================================
//
//
//                  Section 2: Simulator Class Definition
//
//
//=========================================================================

class Simulator 
{
 public:
  Simulator(SimulationControl*, tRainfall*,
	    tRunTimer*, tCOutput <tCNode> *, 
       tRestart<tCNode> *);
  ~Simulator();

  SimulationControl *simCtrl;     // Pointer to simulation control 
  tRainfall         *rainIn;      // Pointer to rainfall input   
  tRunTimer         *timer;       // Pointer to timer 
  tCOutput<tCNode>  *outp;        // Pointer to output object
  tRestart<tCNode>  *restart;     // Pointer to restart object

  int count, fState;
  double dt_rain;                 // Time step of rain 
  double lfr_hour;                // Time tag of last forecasted rainfall 
  double lmr_hour;                // Time tag of last measured rainfall 
  double begin_hour;              // Time tag of initial time
  double met_hour;                // Time tag of last measured met
  double eti_hour;                // Time tag of last eti comp
  double GW_label;                // Label to check GW model run 
  
  int searchRain;                 // Search threshold (hours)


  int  check_mod_status();
  int  checkForecast();

  void get_next_mrain(int);
  void get_next_met();
  void get_next_gaugerain();
  void initialize_simulation(tEvapoTrans*, tSnowPack*, tInputFile&);
       //SMM 09252008 added parameters
  void end_simulation(tKinemat*);
  void simulation_loop(tHydroModel*, tKinemat*, tEvapoTrans*, 
		       tIntercept*, tWaterBalance*, tSnowPack*, // SKY2008Snow from AJR2007
		       tInputFile&); // SKY2008Snow
  void RunItAgain(tInputFile&, tHydroModel*, tKinemat*, 
		  tEvapoTrans*, tIntercept*, tWaterBalance*, tPreProcess*,  tSnowPack*); // SKY2008Snow from AJR2007
  void PrintRunTimeVars(tHydroModel *, int);
  void UpdatePrecipitationInput(int);
  void SurfaceHydroProcesses(tEvapoTrans *, tIntercept *, tSnowPack *); // SKY2008Snow from AJR2007
  void SubSurfaceHydroProcesses(tHydroModel *);
  void OutputSimulatedVars(tKinemat *);
  void UpdateWaterBalance(tWaterBalance *);
  void writeRestart(char*) const;
  void readRestart(tInputFile&);
};

#endif 

//=========================================================================
//
//
//                             End of tSimul.h 
//
//
//=========================================================================
