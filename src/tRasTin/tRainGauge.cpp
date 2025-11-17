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
**  tRainGauge.cpp: Functions for class tRainGauge (see tRainGauge.h)
**
***************************************************************************/

#include "src/tRasTin/tRainGauge.h"

//=========================================================================
//
//
//                  Section 1: tRainGauge Constructors and Destructors
//
//
//=========================================================================

tRainGauge::tRainGauge()
{
	numTimes = 0;
	stationID = 0;
	basinLong = 0.0;
	basinLat = 0.0;
	elev = 0.0;
}

tRainGauge::~tRainGauge()
{
	if (numTimes) {
		delete [] fileName;
		delete [] rain;
		delete [] year;
		delete [] month;
		delete [] day;
		delete [] hour;
	}
}

//=========================================================================
//
//
//                  Section 2: tRainGauge Functions
//
//
//=========================================================================

/***************************************************************************
**
** Set() and Get() Functions for Station Indentifiers
**
***************************************************************************/

void tRainGauge::setStation(int id){
	stationID = id;
}

int tRainGauge::getStation(){
	return stationID;
}

void tRainGauge::setLat(double lat){
	basinLat = lat;
}

double tRainGauge::getLat(){ 
	return basinLat;
}

void tRainGauge::setLong(double longit){
	basinLong = longit;
}

double tRainGauge::getLong(){
	return basinLong;
}

// SKY2008Snow from AJR2007 starts here
void tRainGauge::setElev(double value){
  elev = value;//AJR @ NMT 2007
}

double tRainGauge::getElev(){
  return elev;//AJR @ NMT 2007
}
// SKY2008Snow from AJR2007 ends here

void tRainGauge::setTime(int time){
	numTimes = time;
}

int tRainGauge::getTime(){
	return numTimes;
}


void tRainGauge::setFileName(char* file){
	fileName = new char[kName];
	for(int ct = 0;ct<kName;ct++){
		fileName[ct] = file[ct];
	}
}

char* tRainGauge::getFileName(){
	return fileName;
}

void tRainGauge::setParm(int par){
	numParams = par;
}

int tRainGauge::getParm(){
	return numParams;
}

/***************************************************************************
**
** Set() and Get() Functions for Rainfall Data
**
***************************************************************************/

void tRainGauge::setYear(int *syear){ 
	year = new int[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		year[ct] = syear[ct]; 
	}
}

int tRainGauge::getYear(int time){
	return year[time];
}

void tRainGauge::setMonth(int *smonth){
	month = new int[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		month[ct] = smonth[ct];
	}
}

int tRainGauge::getMonth(int time){
	return month[time];
}

void tRainGauge::setDay(int *sday){
	day = new int[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		day[ct] = sday[ct];
	}
}

int tRainGauge::getDay(int time){
	return day[time];
}

void tRainGauge::setHour(int *shour){
	hour = new int[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		hour[ct] = shour[ct];
	}
}


int tRainGauge::getHour(int time){
	return hour[time];
}

void tRainGauge::setRain(double *r){
	rain = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		rain[ct] = r[ct];
	}
}

double tRainGauge::getRain(int time){
	return rain[time];
}

/***************************************************************************
**
** tRainGauge::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/

void tRainGauge::writeRestart(fstream & rStr) const
{
  BinaryWrite(rStr, numTimes);
  BinaryWrite(rStr, stationID);
  BinaryWrite(rStr, numParams);
  for (int i = 0; i < numTimes; i++) {
    BinaryWrite(rStr, year[i]);
    BinaryWrite(rStr, month[i]);
    BinaryWrite(rStr, day[i]);
    BinaryWrite(rStr, hour[i]);
  }

  BinaryWrite(rStr, basinLat);
  BinaryWrite(rStr, basinLong);
  for (int i = 0; i < numTimes; i++)
    BinaryWrite(rStr, rain[i]);

  BinaryWrite(rStr, elev);      // snow
}


/***************************************************************************
**
** tRainGauge::readRestart() Function
**
***************************************************************************/

void tRainGauge::readRestart(fstream & rStr)
{
  BinaryRead(rStr, numTimes);
  BinaryRead(rStr, stationID);
  BinaryRead(rStr, numParams);
  for (int i = 0; i < numTimes; i++) {
    BinaryRead(rStr, year[i]);
    BinaryRead(rStr, month[i]);
    BinaryRead(rStr, day[i]);
    BinaryRead(rStr, hour[i]);
  }

  BinaryRead(rStr, basinLat);
  BinaryRead(rStr, basinLong);
  for (int i = 0; i < numTimes; i++)
    BinaryRead(rStr, rain[i]);

  BinaryRead(rStr, elev);      // snow
}

//=========================================================================
//
//
//                         End of tRainGauge.cpp
//
//
//=========================================================================
