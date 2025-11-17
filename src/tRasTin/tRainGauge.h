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
**  tRainGauge.h:   Header file for class tRainGauge. Modeled after 
**                  tHydroMet
**
**  This class encapsulates the rainfall data from a raingauge station.
**  Functions provided to create and access measurements.
**  Data members include rainfall data and station information.
**
***************************************************************************/

#ifndef RAINGAUGE_H
#define RAINGAUGE_H

#include "src/Headers/Inclusions.h"

//=========================================================================
//
//
//                  Section 1: tRainGauge Class Declaration
//
//
//=========================================================================

class tRainGauge
{
 public:
  tRainGauge();
  ~tRainGauge();
  void setStation(int);
  void setLat(double);
  void setLong(double);

  // SKY2008Snow from AJR2007
  void setElev(double); //AJR @ NMT 2007

  void setTime(int);
  void setFileName(char*);
  void setYear(int *);
  void setMonth(int *);
  void setDay(int *);
  void setHour(int *);
  void setRain(double *);
  void setParm(int);
  
  int getStation();
  double getLat();
  double getLong();

  // SKY2008Snow from AJR2007
  double getElev(); //AJR @ NMT 2007

  int getTime();
  char* getFileName();
  int getYear(int);
  int getMonth(int);
  int getDay(int);
  int getHour(int);
  int getParm();
  double getRain(int);
  void writeRestart(fstream &) const;
  void readRestart(fstream &);
  
 protected:
  int numTimes, stationID, numParams;
  char *fileName;
  double basinLat, basinLong;

  // SKY2008Snow from AJR2007
  double elev;//AJR @ NMT 2007

  double *rain;
  int *year, *month, *day, *hour;
  
};

#endif


//=========================================================================
//
//
//                        End of tRainGauge.h
//
//
//=========================================================================
