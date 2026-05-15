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
**  tHydroMet.cpp: Functions for class tHydroMet (see tHydroMet.h)
**
***************************************************************************/

#include "src/tHydro/tHydroMet.h"

//=========================================================================
//
//
//                  Section 1: tHydroMet Constructors and Destructors
//
//
//=========================================================================
tHydroMet::tHydroMet()
{
	numParams = 0;
	stationID = 0;
	basinLong = 0.0;
	basinLat = 0.0;
	otherVariable = 0.0;
	meanTemp = 0.0;
}

tHydroMet::~tHydroMet() = default;

//=========================================================================
//
//
//                  Section 2: tHydroMet Functions
//
//
//=========================================================================

/***************************************************************************
**
** Set() and Get() Functions for Station Identifiers
**
***************************************************************************/

void tHydroMet::setStation(int id){
	stationID = id;
}

int tHydroMet::getStation(){
	return stationID;
}

void tHydroMet::setLat(double lat){
	basinLat = lat;
}

double tHydroMet::getLat(){
	return basinLat;
}

void tHydroMet::setLong(double longit){
	basinLong = longit;
}

double tHydroMet::getLong(){
	return basinLong;
}

int tHydroMet::getTime() const {
	return static_cast<int>(year.size());
}

void tHydroMet::setParm(int parm){
	numParams = parm;
}

int tHydroMet::getParm(){
	return numParams;
}

void tHydroMet::setOther(double other){
	otherVariable = other;
}

double tHydroMet::getOther(){
	return otherVariable;
}

void tHydroMet::setFileName(char* file){
	snprintf(fileName, sizeof(fileName), "%s", file);
}

char* tHydroMet::getFileName(){
	return fileName;
}


/***************************************************************************
**
** Set() and Get() Functions for Hydrometeorologic Data
**
***************************************************************************/

void tHydroMet::setYear(const std::vector<int> &v){
	year = v;
}

int tHydroMet::getYear(int time) const {
	return year[time];
}

void tHydroMet::setMonth(const std::vector<int> &v){
	month = v;
}

int tHydroMet::getMonth(int time) const {
	return month[time];
}

void tHydroMet::setDay(const std::vector<int> &v){
	day = v;
}

int tHydroMet::getDay(int time) const {
	return day[time];
}

void tHydroMet::setHour(const std::vector<int> &v){
	hour = v;
}

int tHydroMet::getHour(int time) const {
	return hour[time];
}

void tHydroMet::setAirTemp(const std::vector<double> &v){
	int count = 0;
	meanTemp = 0.0;
	airTemp = v;
	for (double val : v) { meanTemp += val; count++; }
	if (count > 0) meanTemp /= count;
}

double tHydroMet::getAirTemp(int time) const {
	return airTemp[time];
}

double tHydroMet::getMeanTemp(){
	return meanTemp;
}

void tHydroMet::setSurfTemp(const std::vector<double> &v){
	surfTemp = v;
}

double tHydroMet::getSurfTemp(int time) const {
	return surfTemp[time];
}

void tHydroMet::setAtmPress(const std::vector<double> &v){
	atmPress = v;
}

double tHydroMet::getAtmPress(int time) const {
	return atmPress[time];
}

void tHydroMet::setSkyCover(const std::vector<double> &v){
	skyCover = v;
}

double tHydroMet::getSkyCover(int time) const {
	return skyCover[time];
}

void tHydroMet::setRHumidity(const std::vector<double> &v){
	rHumidity = v;
}

double tHydroMet::getRHumidity(int time) const {
	return rHumidity[time];
}

void tHydroMet::setWindSpeed(const std::vector<double> &v){
	windSpeed = v;
}

double tHydroMet::getWindSpeed(int time) const {
	return windSpeed[time];
}

void tHydroMet::setPanEvap(const std::vector<double> &v){
	panEvap = v;
}

double tHydroMet::getPanEvap(int time) const {
	return panEvap[time];
}

void tHydroMet::setRadGlobal(const std::vector<double> &v){
	RadGlobal = v;
}

double tHydroMet::getRadGlobal(int time) const {
	return RadGlobal[time];
}

void tHydroMet::setRadDirect(const std::vector<double> &v){
	RadDirect = v;
}

double tHydroMet::getRadDirect(int time) const {
	return RadDirect[time];
}

void tHydroMet::setRadDiffuse(const std::vector<double> &v){
	RadDiffuse = v;
}

double tHydroMet::getRadDiffuse(int time) const {
	return RadDiffuse[time];
}

void tHydroMet::setRainMet(const std::vector<double> &v){
	RainMet = v;
}

double tHydroMet::getRainMet(int time) const {
	return RainMet[time];
}

/***************************************************************************
**
** tHydroMet::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/

void tHydroMet::writeRestart(fstream & rStr) const
{
  int sz = static_cast<int>(year.size());
  BinaryWrite(rStr, sz);
  BinaryWrite(rStr, numParams);
  BinaryWrite(rStr, stationID);
  BinaryWrite(rStr, basinLat);
  BinaryWrite(rStr, basinLong);
  BinaryWrite(rStr, otherVariable);
  BinaryWrite(rStr, meanTemp);
  for (int i = 0; i < sz; i++) {
    BinaryWrite(rStr, year[i]);
    BinaryWrite(rStr, month[i]);
    BinaryWrite(rStr, day[i]);
    BinaryWrite(rStr, hour[i]);
  }
	// Giuseppe DEBUG Restart 2012 - START 
	// I have introduced an IF and ELSE to deal with the case
	// of PAN ET data
  if (numParams == 5) {
    for (int i = 0; i < sz; i++)
      BinaryWrite(rStr, panEvap[i]);
  } else {// Giuseppe DEBUG Restart 2012 - END	
    for (int i = 0; i < sz; i++) {
      BinaryWrite(rStr, airTemp[i]);
      BinaryWrite(rStr, surfTemp[i]);
      BinaryWrite(rStr, rHumidity[i]);
      BinaryWrite(rStr, skyCover[i]);
      BinaryWrite(rStr, windSpeed[i]);
      BinaryWrite(rStr, atmPress[i]);
      BinaryWrite(rStr, RadGlobal[i]);
    }
  }
}

/***************************************************************************
**
** tHydroMet::readRestart() Function
**
***************************************************************************/

void tHydroMet::readRestart(fstream & rStr)
{
  int sz = 0;
  BinaryRead(rStr, sz);
  BinaryRead(rStr, numParams);
  BinaryRead(rStr, stationID);
  BinaryRead(rStr, basinLat);
  BinaryRead(rStr, basinLong);
  BinaryRead(rStr, otherVariable);
  BinaryRead(rStr, meanTemp);

  year.resize(sz); month.resize(sz); day.resize(sz); hour.resize(sz);
  for (int i = 0; i < sz; i++) {
    BinaryRead(rStr, year[i]);
    BinaryRead(rStr, month[i]);
    BinaryRead(rStr, day[i]);
    BinaryRead(rStr, hour[i]);
  }
	// Giuseppe DEBUG Restart 2012 - START 
	// I have introduced an IF and ELSE to deal with the case
	// of PAN ET data
  if (numParams == 5) {
    panEvap.resize(sz);
    for (int i = 0; i < sz; i++)
      BinaryRead(rStr, panEvap[i]);
  } else {// Giuseppe DEBUG Restart 2012 - END	
    airTemp.resize(sz); surfTemp.resize(sz); rHumidity.resize(sz);
    skyCover.resize(sz); windSpeed.resize(sz); atmPress.resize(sz);
    RadGlobal.resize(sz);
    for (int i = 0; i < sz; i++) {
      BinaryRead(rStr, airTemp[i]);
      BinaryRead(rStr, surfTemp[i]);
      BinaryRead(rStr, rHumidity[i]);
      BinaryRead(rStr, skyCover[i]);
      BinaryRead(rStr, windSpeed[i]);
      BinaryRead(rStr, atmPress[i]);
      BinaryRead(rStr, RadGlobal[i]);
    }
  }
}

//=========================================================================
//
//
//                         End of tHydroMet.cpp
//
//
//=========================================================================
