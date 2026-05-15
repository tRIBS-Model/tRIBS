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
	stationID = 0;
	fileName = nullptr;
	basinLong = 0.0;
	basinLat = 0.0;
	elev = 0.0;
}

tRainGauge::~tRainGauge()
{
	delete[] fileName;
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
** Set() and Get() Functions for Station Identifiers
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

void tRainGauge::setElev(double value){
	elev = value;
}

double tRainGauge::getElev(){
	return elev;
}

void tRainGauge::setFileName(char* file){
	fileName = new char[kName];
	for (int ct = 0; ct < kName; ct++)
		fileName[ct] = file[ct];
}

char* tRainGauge::getFileName(){
	return fileName;
}

/***************************************************************************
**
** Set() and Get() Functions for Rainfall Data
**
***************************************************************************/

void tRainGauge::setYear(const std::vector<int>& v){
	year = v;
}

int tRainGauge::getYear(int time){
	return year[time];
}

void tRainGauge::setMonth(const std::vector<int>& v){
	month = v;
}

int tRainGauge::getMonth(int time){
	return month[time];
}

void tRainGauge::setDay(const std::vector<int>& v){
	day = v;
}

int tRainGauge::getDay(int time){
	return day[time];
}

void tRainGauge::setHour(const std::vector<int>& v){
	hour = v;
}

int tRainGauge::getHour(int time){
	return hour[time];
}

void tRainGauge::setRain(const std::vector<double>& v){
	rain = v;
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
	int sz = static_cast<int>(year.size());
	BinaryWrite(rStr, sz);
	BinaryWrite(rStr, stationID);
	for (int i = 0; i < sz; i++) {
		BinaryWrite(rStr, year[i]);
		BinaryWrite(rStr, month[i]);
		BinaryWrite(rStr, day[i]);
		BinaryWrite(rStr, hour[i]);
	}
	BinaryWrite(rStr, basinLat);
	BinaryWrite(rStr, basinLong);
	for (int i = 0; i < sz; i++)
		BinaryWrite(rStr, rain[i]);
	BinaryWrite(rStr, elev);
}


/***************************************************************************
**
** tRainGauge::readRestart() Function
**
***************************************************************************/

void tRainGauge::readRestart(fstream & rStr)
{
	int sz;
	BinaryRead(rStr, sz);
	BinaryRead(rStr, stationID);
	year.resize(sz);
	month.resize(sz);
	day.resize(sz);
	hour.resize(sz);
	for (int i = 0; i < sz; i++) {
		BinaryRead(rStr, year[i]);
		BinaryRead(rStr, month[i]);
		BinaryRead(rStr, day[i]);
		BinaryRead(rStr, hour[i]);
	}
	BinaryRead(rStr, basinLat);
	BinaryRead(rStr, basinLong);
	rain.resize(sz);
	for (int i = 0; i < sz; i++)
		BinaryRead(rStr, rain[i]);
	BinaryRead(rStr, elev);
}

//=========================================================================
//
//
//                         End of tRainGauge.cpp
//
//
//=========================================================================
