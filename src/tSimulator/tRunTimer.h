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
**  tRuntimer.h: Header for tRunTimer objects
**  	
**  tRunTimer objects are used to keep track of time in a time-evolving
**  simulation model. Their services include keeping track of when it's
**  time to write output, printing the current time to standard output if
**  desired, and writing the current time to a file every so often.
**
***************************************************************************/

#ifndef TRUNTIMER_H
#define TRUNTIMER_H

//=========================================================================
//
//
//                  Section 1: tRunTimer Include Statements
//
//
//=========================================================================

#include "src/tInOut/tInputFile.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <cassert>

//=========================================================================
//
//
//                  Section 2: tRunTimer Class Definition
//
//
//=========================================================================

class tRunTimer
{
public:
  tRunTimer( double duration, int optprint=1 );
  tRunTimer( tInputFile &infile );
  tRunTimer();
  double  getCurrentTime() const;
  double  getTimeStep() const;
  double  getMetStep() const;
  double  getEtIStep() const;
  double  getGWTimeStep() const;
  double  getRainDT() const;
  double  getEndTime() const;
  double  getRainTime() const;
  double  RemainingTime() const;
  double  getOutputInterval() const;
  double  getSpatialOutputInterval() const;
  double  RemainingTime(double);
  double  getOutputIntervalSec();
  double  getSpatialOutputIntervalSec();
  double  res_hour_mid(int);
  double  res_hour_begin(int);
  double  res_hour_end(int);
  double  get_abs_hour(double);
  double  getfTime();
  double  getfLength();
  double  getfLead();
  double  getMetTime(int);
  double  getStormTime();
  double  getPrevStormTime();	

  int     getStepForSpecifiedDT(double, double);
  int     getElapsedSteps(double);
  int     getElapsedETISteps(double);
  int     getElapsedMETSteps(double);
  int     getResStep(double);
  int     Advance( double );
  int     IsFinished();
  int     CheckOutputTime();
  int     CheckSpatialOutputTime();
  int     isGaugeTime(double);
  int     getoptForecast();

  void    InitializeTimer( tInputFile &infile );
  void    Start( double, double );
  void    correctCalendarTime(double, double, int *, 
			      int *, int *, int *, int *);
  void    addRainTime();
  void    addMetTime(int);
  void    res_time_begin(int, int *, int *);
  void    res_time_mid(int, int *, int *);
  void    res_time_end(int, int *, int *);
  void    UpdateStorm(double);
  void    writeRestart(fstream &) const;
  void    readRestart(fstream &);

  int days[13];
  int cumdays[13]; 

  int minute, hour, day, month, year; 
  int minuteRn, hourRn, dayRn, monthRn, yearRn;
  int minuteS, hourS, minS, dayS, monthS, yearS;
  double dtRain;    	        // time interval between rainfall files

private:
  double startHour;      	// Start hour of simulations
  double currentTime;
  double endTime;
  double nextOutputTime;
  double nextSPOutputTime;
  double tstep;                 // time step for UNsaturated zone comp.
  double gwatstep;       	// time step for Saturated zone comp.
  double metstep;        	// time step for Meteorological Data
  double outputInterval; 	// time step for output 
  double SPOutputInterval; 	// time step for spatial output
  double RainTime;  	        // keeps track of rainfall input time
  int    optForecast;           // Rainfall forecasting option
  double fTime, fLength, fLead; // Forecasting parameters
  double etistep;               // time step for ET and I comp.
  double MetTime;               // keeps track of met time
  double EtITime;               // keeps track of eti time
  double StormTime;             // keeps track of stochastic storms
  double StormTime_1;           // keeps track of previous storm time

  ofstream timeStatusFile;
};

inline double tRunTimer::getRainTime()    const { return RainTime; }
inline double tRunTimer::getRainDT()      const { return dtRain; }
inline double tRunTimer::getMetStep()     const { return metstep; }
inline double tRunTimer::getEtIStep()     const { return etistep; }
inline double tRunTimer::getCurrentTime() const { return currentTime; }
inline double tRunTimer::RemainingTime()  const { return (endTime - currentTime); }
inline double tRunTimer::getTimeStep()    const { return tstep; }
inline double tRunTimer::getGWTimeStep()  const { return gwatstep; }
inline double tRunTimer::getEndTime()     const { return endTime; }
inline double tRunTimer::getOutputInterval() const { return outputInterval; }
inline double tRunTimer::getSpatialOutputInterval() const { return SPOutputInterval; }

#endif
       
//=========================================================================
//
//
//                          End tRunTimer.h
//
//
//=========================================================================
