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
**  tHydroMet.h:   Header file for class tHydroMet
**
**  This class encapsulates the hydrometeorological from a weather station.
**  Functions provided to create and access micrometeorological measurements.
**  Data members include weather data and station information.
**
***************************************************************************/

#ifndef HYDROMET_H
#define HYDROMET_H

#include "src/Headers/Inclusions.h"

//=========================================================================
//
//
//                  Section 1: tHydroMet Class Declaration
//
//
//=========================================================================

class tHydroMet
{
 public:
  tHydroMet();
  ~tHydroMet();
  void setStation(int);
  void setLat(double,int);
  void setLong(double,int);
  void setOther(double);
  void setGmt(int);
  void setTime(int);
  void setParm(int);
  void setFileName(char*);
  void setYear(int *);
  void setMonth(int *);
  void setDay(int *);
  void setHour(int *);
  void setAirTemp(double *);
  void setDewTemp(double *);
  void setSurfTemp(double *);
  void setAtmPress(double *);
  void setSkyCover(double *);
  void setRHumidity(double *);
  void setWindSpeed(double *);
  void setNetRad(double *);
  void setPanEvap(double *);
  void setVaporPress(double *);
  void setRadGlobal(double *);
  void setRadDirect(double *);
  void setRadDiffuse(double *);
  void setRainMet(double *);

  int    getStation();
  double getLat(int);
  double getLong(int);
  double getOther();
  int    getGmt();
  int    getTime();
  int    getParm();
  char*  getFileName();
  int*   getYear();
  int    getYear(int);
  int*   getMonth(); 
  int    getMonth(int);
  int*   getDay();
  int    getDay(int);
  int*   getHour();
  int    getHour(int);
  double getAirTemp(int);
  double getMeanTemp();
  double getDewTemp(int);
  double getSurfTemp(int);
  double getAtmPress(int);
  double getSkyCover(int);
  double getRHumidity(int);
  double getWindSpeed(int);
  double getNetRad(int);
  double getPanEvap(int);
  double getVaporPress(int);
  double getRadGlobal(int);
  double getRadDirect(int);
  double getRadDiffuse(int);
  double getRainMet(int);

  void   writeRestart(fstream &) const;
  void   readRestart(fstream &);

 protected:
  int numTimes, numParams;
  int gmt, stationID;
  char fileName[kName];
  double latitude, longitude;
  double basinLat, basinLong;
  double otherVariable, meanTemp;
  double *airTemp, *dewTemp, *surfTemp;
  double *netRad, *panEvap, *vaporPress;
  double *rHumidity, *skyCover;
  double *windSpeed, *atmPress;
  double *RadGlobal, *RadDirect, *RadDiffuse, *RainMet;
  int *year, *month, *day, *hour;
  
};

#endif


//=========================================================================
//
//
//                        End of tHydroMet.h
//
//
//=========================================================================
