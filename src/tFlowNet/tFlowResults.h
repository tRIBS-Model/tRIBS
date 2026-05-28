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
**  tFlowResults.h: Functions for class tFlowResults (see tFlowResults.cpp)
**  
**  tFlowResults is a class for handling storage of the output discharge
**  and assigning values corresponding outputInterval in tRunTimer. It also
**  deals with runoff types.
**
***************************************************************************/

#ifndef TFLOWRESULTS_H
#define TFLOWRESULTS_H

//=========================================================================
//
//
//                  Section 1: tFlowResults Include & Define Statements
//
//
//=========================================================================

#include "src/tSimulator/tRunTimer.h"
#include "src/tSimulator/tControl.h"
#include "src/tInOut/tInputFile.h"
#include "src/Headers/Definitions.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>

//=========================================================================
//
//
//                  Section 2: tFlowResults Class Definitions
//
//
//=========================================================================

class tFlowResults 
{ 
public:
  tFlowResults( SimulationControl *, tInputFile &, tRunTimer *, double );
  ~tFlowResults();
    
  tRunTimer *timer;   		       // Time evolution 
  SimulationControl *simCtrl;          // Simulation Control

  char    baseHydroName[kMaxNameSize]; // Basename for output hydrograph 
  char    outlet[kMaxNameSize];        // Outlet name, root of the above
  char    Extension[kMaxExt];          // Extension for output hydrograph
  char    currHydroName[kMaxNameSize+kMaxExt]; // Current hydrograph name

  int     limit;   		       // Size of results array 
  int     iimax;   		       // Last non-zero of results array
  int     ribsOutput;                  // Compatibility with RIBS interphase
  int     writeFlag;                   // Flag for writing *.mrf header
  int     count;

  double *prr;     		// Previous rainfall  
  double *crr;     		// Rainfall in the time step 
  double *phydro;  		// Previous hydrograph (expected response)
  double *mhydro;  		// Response with measured rain 
  double *HsrfRout;             // Infiltration-excess
  double *SbsrfRout;            // Saturation-excess
  double *PsrfRout;             // Perched-return
  double *SatsrfRout;           // Groundwater exfiltration
  double *max;                  // Maximum rainfall in space
  double *min;                  // Minimum rainfall in space
  double *msm;                  // Mean soil moisture
  double *msmRt;                // Mean soil moisture in the root zone
  double *msmU;                 // Mean soil moisture in the unsaturated zone
  double *mgw;                  // Mean groundwater level in the basin
  double *met;                  // Mean evapotranspiration
  double *sat;                  // Basin saturated fraction
  double *frac;                 // Fractional rainfall coverage in basin

  // SKY2008Snow from AJR2007
  double *swe;			// Mean SWE in space
  double *melt;			// Mean melt in space
  double *snsub;		// Mean sn sublimation in space // CJC2020
  double *snevap;		// Mean sn evaporation in space // CJC2020
  double *stC;			// Mean snow temp (c) in space
  double *DUint;		// Mean change in energy in space
  double *slhf;			// Mean sn latent hf in space
  double *sshf;			// Mean sn sensibe hf in space
  double *sghf; 		// Mean sn ground hf in space
  double *sphf;			// Mean sn precipitation hf in space
  double *srli;			// Mean sn RLin in space
  double *srlo;			// Mean sn RLout in space
  double *srsi;			// Mean sn RSin in spce
  double *intsn;		// Mean int SWE in space
  double *intsub;		// Mean int sublimation in space
  double *intunl;		// Mean int unloading in space
  double *sca;			// Fraction snow covered area
  double *Perc; 	// Percolation from the channel bottom ASM percolation option
  double *qunsat;		// Mean net flow froom unsaturated zone CJC2025

  void  SetFlowResVariables(tInputFile &, double);
  void  writeAndUpdate(double);
  void  store_rain(double,double);
  void  store_maxminrain(double,double, int);
  void  store_saturation(double,double, int);
  void  free_results();
  void  write_inter_hyd(char *, char *);
  void  write_Runoff_Types(char *, char *);
  void  whenTimeIsOver( double );
  void  store_volume(double, double);          
  void  store_volume_Type(double, double, int); 
  void  add_m_volume_Type(double, int, int); 

  void  add_m_volume(double value, int iStep)   
    { mhydro[iStep]+= value/timer->getOutputIntervalSec(); } 

  double get_discharge(int ihour)
    {  return phydro[ihour]+mhydro[ihour]; }

  double get_discharge(double time) 
  { 
        int ihour;
        ihour=timer->getResStep(time); 
        return phydro[ihour]+mhydro[ihour]; 
  }

  void reset_meas_hyd()
    { for(int ii=0; ii < limit; ii++) mhydro[ii]=0.0; }
    

  void initialize_prev_hyd(double disch)
    { for(int ii=0; ii < limit; ii++) phydro[ii]=disch; }
   

  void update_prev_hyd()
    { for(int ii=0; ii < iimax; ii++) phydro[ii]+=mhydro[ii]; }

};

#endif

//=========================================================================
//
//
//                          End tFlowResults.h
//
//
//=========================================================================
