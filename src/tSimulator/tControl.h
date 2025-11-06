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
**  tControl.h: Header for tControl.cpp 
**  	
**  tControl Class is a class for the control functions of the
**  model simulation inputted by the user through the command line.
**
***************************************************************************/


#ifndef TCONTROL_H
#define TCONTROL_H

//=========================================================================
//
//
//                  Section 1: tControl Include and Define Statements
//
//
//=========================================================================


#ifdef ALPHA_64
  #include <iostream.h>
  #include <stdio.h>
  #include <stdlib.h>
#elif defined LINUX_32
  #include <iostream>
  #include <cstdio>
  #include <cstdlib>

#elif defined MAC
  #include <iostream>
  #include <cstdio>
  #include <cstdlib>

#elif defined WIN
  #include <iostream.h>
  #include <stdio.h>
  #include <stdlib.h>
#else 
  #include <iostream.h>
  #include <stdio.h>
  #include <stdlib.h>
#endif

using namespace std;

#define STD_INPUT  1
#define AUTO_INPUT 2

//=========================================================================
//
//
//                  Section 2: SimulationControl Class Definition
//
//
//=========================================================================

class SimulationControl {
 public:
  int  VerbID;            // ID of a verbose node 
  int  num_simul;         // # of simulation runs 
  char first_time;        // First computation loop Y or N 
  char mode;              // Mode of rainfall input
  bool inter_results;     // Write intermediate results Y or N,
  bool GW_model_label;    // Run groundwater model Y or N
  char Verbose_label;     // Verbose screen output Y or N
  char Check_label;       // Checking input file Y or N
  char *infile;           // Name of input file containing data
  char mod_is_on;         // The model stays on and waits for commands
  char hydro_visual;      // To turn on hydrograph visualization
  char Header_label;      // Suppress header information in outputs
  char fore_rain_label;   // Forecasted rain = Y or N 
  char smooth_weather;    // Special option: no randomness in climate
  char disp_time;
    // *Do NOT display it in help menu: confusing
  char debug;             // For debugging output for tGraph SMM

  SimulationControl(int, char **);
  ~SimulationControl();
};

#endif

//=========================================================================
//
//
//                          End of tControl.h
//
//
//=========================================================================
