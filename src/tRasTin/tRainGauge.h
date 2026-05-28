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
  void setFileName(char*);
  void setYear(const std::vector<int>&);
  void setMonth(const std::vector<int>&);
  void setDay(const std::vector<int>&);
  void setHour(const std::vector<int>&);
  void setRain(const std::vector<double>&);

  int getStation();
  double getLat();
  double getLong();

  // SKY2008Snow from AJR2007
  double getElev(); //AJR @ NMT 2007
  char* getFileName();
  int getYear(int);
  int getMonth(int);
  int getDay(int);
  int getHour(int);
  double getRain(int);

 protected:
  int stationID;
  char *fileName;
  double basinLat, basinLong;

  // SKY2008Snow from AJR2007
  double elev;//AJR @ NMT 2007
  std::vector<int>    year, month, day, hour;
  std::vector<double> rain;
};

#endif


//=========================================================================
//
//
//                        End of tRainGauge.h
//
//
//=========================================================================
