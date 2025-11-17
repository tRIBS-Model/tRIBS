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
**  tHydroMetConvert.h:   Header file for tHydroMetConvert.cpp 
**
**  Pre-processor program used to convert Operational RFC Point Informix Data
**  into tRIBS HMET_WES data format for hydrometeorological station format
**  (METDATAOPTION = 1). Class invoked only when CONVERTDATA Option = 1.
** 
***************************************************************************/

#ifndef THYDROMETCONVERT_H
#define THYDROMETCONVERT_H

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
//                  Section 1: tHydroMetConvert Class Declaration
//
//
//=========================================================================

class tHydroMetConvert
{
 public: 
  tHydroMetConvert();
  tHydroMetConvert(tInputFile &);
  ~tHydroMetConvert();
  void initialize();
  void callConvertRFC();
  void callConvertDMIP();
  void readMDI();
  void readLocData(int);
  void readPointData(int);
  void writeSDF(int);
  void writeMDF(int);
  void writeGaugeMDF(int);
  void writeGaugeSDF(int);
  void callMerge();
  void callGaugeMerge();
  void insertSort(double[], int);
  void readAndWriteDMIP();
  
 protected:
  char mdiFile[kName];
  char mdfFilebase[kName];
  char sdfFile[kName];
  char **fNameArray, **sNameArray, **pNameArray, **lNameArray, **param, **stat;
  double *date, *hour, *value, *elev, *latitude, *longitude, *dateArray;
  int *lookFor, *numTimes, *dNameArray;
  int numFiles, numStations, numParameters, Ncount, optMerge, convertData;
};

#endif


//=========================================================================
//
//
//                      End of tHydroMetConvert.h
//
//
//=========================================================================
