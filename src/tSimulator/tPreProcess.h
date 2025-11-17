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
**  tPreProcess.h:  Header file for tPreProcess Class
**
**  This class is utilized for data pre-processing purposes which include
**  meteorological and rain gauge extraction from available sources. The
**  preprocessing class should be expandable to handle future data sources.
**
***************************************************************************/

#ifndef PREPROCESS_H
#define PREPROCESS_H

#include "src/Headers/Inclusions.h"

#ifdef ALPHA_64
#elif defined LINUX_32
  #ifndef _WIN32
  #include <unistd.h>
  #endif
#elif defined WIN
#else 
#endif

//=========================================================================
//
//
//            Section 1: tPreProcess Class Declaration
//
//
//=========================================================================

class tPreProcess
{
 public:
  tPreProcess();
  tPreProcess(SimulationControl*, tInputFile &);
  ~tPreProcess();
  
  SimulationControl *simCtrl;

  void CheckInputFile(tInputFile &);
  void CheckFileExists(tInputFile &, char *, const char*);
  void CheckPathNameCorrect(tInputFile &, char*, const char*);

  int    IterReadItem(tInputFile &,    int, const char *);
  void   IterReadItem(tInputFile &, char *, const char *);
  double IterReadItem(tInputFile &, double, const char *);

  int convertData;

};

#endif

//=========================================================================
//
//
//                    End of tPreProcess.h
//
//
//=========================================================================
