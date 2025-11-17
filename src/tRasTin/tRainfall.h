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
**  tRainfall.h:   Header file for tRainfall classes (see tRainfall.h)
**
**  tRainfall class is used to read gridded radar rainfall input 
**
***************************************************************************/

#ifndef TRAINFALL_H
#define TRAINFALL_H

#include "src/Headers/Inclusions.h"
#include "src/tStorm/tStorm.h"

using namespace std;

class tRunTimer;

//=========================================================================
//
//
//                  Section 1: tRainfall Class Declarations
//
//
//=========================================================================

class tRainfall : public tStorm 
{
public:
  tRainfall();
  ~tRainfall();
  tRainfall(SimulationControl*, tMesh<tCNode> *, tInputFile &, tResample *);
  SimulationControl *simCtrl;    
  
  int  Compose_In_Mrain_Name(tRunTimer *);
  void SetRainVariables(tInputFile &);  
  void NewRain();
  void NewRain(double);

  void NewRain(tRunTimer *);
  void NewRainData(int);
  void InitializeGauge();
  void readGaugeStat(char *);
  void readGaugeData(int);
  void robustNess(double*, int);
  void assignStationToNode();
  void callRainGauge(tRunTimer *);
  void setToNode();
  void setfState(int);
  void writeRestart(fstream &) const;
  void readRestart(fstream &);
  
  char mrainfileIn[kMaxNameSize];
  int searchRain, rainfallType;
  double rainDt;
  tResample * getRsmplPtr() {return respPtr;}
  tMesh<tCNode> * getMeshPtr() {return gridPtr;}

protected:
  tMesh<tCNode> *gridPtr;
  tResample     *respPtr;
  tRainGauge * rainGauges;
  int numStations, arraySize, hourlyTimeStep, numRains;
  int *assignedRain, *currentTime;
  double *curRain, *latitude, *longitude, *gaugeRain; 

  // SKY2008Snow from AJR2007
  double precLapseRate;//AJR @ NMT 2007

  char inputname[kMaxNameSize];
  char forecastname[kMaxNameSize];
  char stationFile[kName];
  char extension[20]; 
  double aveMAP, cumMAP, climate;   
  int optForecast, fState, optMAP;
  ifstream infile; 
};

#endif

//=========================================================================
//
//
//                         End of tRainfall.h
//
//
//=========================================================================
