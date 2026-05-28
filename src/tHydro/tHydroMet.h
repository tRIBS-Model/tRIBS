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
#include <vector>

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
  void setLat(double);
  void setLong(double);
  void setOther(double);
  void setParm(int);
  void setFileName(char*);
  void setYear(const std::vector<int> &);
  void setMonth(const std::vector<int> &);
  void setDay(const std::vector<int> &);
  void setHour(const std::vector<int> &);
  void setAirTemp(const std::vector<double> &);
  void setSurfTemp(const std::vector<double> &);
  void setAtmPress(const std::vector<double> &);
  void setSkyCover(const std::vector<double> &);
  void setRHumidity(const std::vector<double> &);
  void setWindSpeed(const std::vector<double> &);
  void setPanEvap(const std::vector<double> &);
  void setRadGlobal(const std::vector<double> &);
  void setRadDirect(const std::vector<double> &);
  void setRadDiffuse(const std::vector<double> &);
  void setRainMet(const std::vector<double> &);

  int    getStation();
  double getLat();
  double getLong();
  double getOther();
  int    getTime() const;
  int    getParm();
  char*  getFileName();
  int    getYear(int) const;
  int    getMonth(int) const;
  int    getDay(int) const;
  int    getHour(int) const;
  double getAirTemp(int) const;
  double getMeanTemp();
  double getSurfTemp(int) const;
  double getAtmPress(int) const;
  double getSkyCover(int) const;
  double getRHumidity(int) const;
  double getWindSpeed(int) const;
  double getPanEvap(int) const;
  double getRadGlobal(int) const;
  double getRadDirect(int) const;
  double getRadDiffuse(int) const;
  double getRainMet(int) const;


 protected:
  int numParams;
  int stationID;
  char fileName[kName];
  double basinLat, basinLong;
  double otherVariable, meanTemp;
  std::vector<int>    year, month, day, hour;
  std::vector<double> airTemp, surfTemp;
  std::vector<double> panEvap;
  std::vector<double> rHumidity, skyCover;
  std::vector<double> windSpeed, atmPress;
  std::vector<double> RadGlobal, RadDirect, RadDiffuse, RainMet;
};

#endif


//=========================================================================
//
//
//                        End of tHydroMet.h
//
//
//=========================================================================
