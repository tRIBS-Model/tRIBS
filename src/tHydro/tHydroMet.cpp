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
	numTimes = 0;
	numParams = 0;
	gmt = 0;
	stationID = 0;
	latitude = 0.0;
	longitude = 0.0;
	basinLong = 0.0;
	basinLat = 0.0;
	otherVariable = 0.0;
	meanTemp = 0.0;
}

tHydroMet::~tHydroMet()
{
	if (numTimes > 0 ) { 
		delete [] year;
		delete [] month;
		delete [] day;
		delete [] hour;
		
		// BUG fixed by Giuseppe Mascaro to allow safe deleting when
		// using PAN stations - June 2012
		if (numParams == 5) {
			
			if (panEvap){
				delete [] panEvap;
			}
		}
		else {
			if (airTemp){
				delete [] airTemp;
				delete [] dewTemp;
				delete [] surfTemp;
				delete [] rHumidity;
				delete [] skyCover;
				delete [] windSpeed;
				delete [] atmPress;
				delete [] vaporPress;
				delete [] RadGlobal;
			}
			if (numParams >= 11) {
				delete [] netRad;
			}
		}
		// END of BUG Fixed by Giuseppe Mascaro - June 2012
	}
}

//=========================================================================
//
//
//                  Section 2: tHydroMet Functions
//
//
//=========================================================================

/***************************************************************************
**
** Set() and Get() Functions for Station Indentifiers
**
***************************************************************************/

void tHydroMet::setStation(int id){
	stationID = id;
}

int tHydroMet::getStation(){
	return stationID;
}

void tHydroMet::setLat(double lat, int option){
	if(option == 1)
		latitude = lat;
	else
		basinLat = lat;
}

double tHydroMet::getLat(int option){
	if(option==1)
		return latitude;
	else 
		return basinLat;
}

void tHydroMet::setLong(double longit, int option){
	if(option == 1)
		longitude = longit;
	else
		basinLong = longit;
}

double tHydroMet::getLong(int option){
	if(option == 1)
		return longitude;
	else
		return basinLong;
}

void tHydroMet::setGmt(int gmtime){
	gmt = gmtime;
}

int tHydroMet::getGmt(){
	return gmt;
}

void tHydroMet::setTime(int time){
	numTimes = time;
}

int tHydroMet::getTime(){
	return numTimes;
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
	snprintf(fileName,sizeof(fileName),"%s", file);//WR--09192023: 'sprintf' is deprecated: This function is provided for compatibility reasons only.
}

char* tHydroMet::getFileName(){
	return fileName;
}


/***************************************************************************
**
** Set() and Get() Functions for Hydrometeorologic Data
**
***************************************************************************/

void tHydroMet::setYear(int *syear){ 
	year = new int[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		year[ct] = syear[ct]; 
	}
}

int* tHydroMet::getYear(){
	return year;
}

int tHydroMet::getYear(int time){
	return year[time];
}

void tHydroMet::setMonth(int *smonth){
	month = new int[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		month[ct] = smonth[ct];
	}
}

int* tHydroMet::getMonth(){
	return month;
}

int tHydroMet::getMonth(int time){
	return month[time];
}

void tHydroMet::setDay(int *sday){
	day = new int[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		day[ct] = sday[ct];
	}
}

int* tHydroMet::getDay(){
	return day;
}

int tHydroMet::getDay(int time){
	return day[time];
}


void tHydroMet::setHour(int *shour){
	hour = new int[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		hour[ct] = shour[ct];
	}
}

int* tHydroMet::getHour(){
	return hour;
}

int tHydroMet::getHour(int time){
	return hour[time];
}

void tHydroMet::setAirTemp(double *air){
	int count = 0; 
	airTemp = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		meanTemp += air[ct];
		airTemp[ct] = air[ct];
		count++;
	}
	meanTemp = meanTemp/count;
}

double tHydroMet::getAirTemp(int time){
    return airTemp[time];
}

double tHydroMet::getMeanTemp(){
	return meanTemp;
}

void tHydroMet::setDewTemp(double *dew){
	dewTemp = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		dewTemp[ct] = dew[ct];
	}
}

double tHydroMet::getDewTemp(int time){
	return dewTemp[time];
}


void tHydroMet::setSurfTemp(double *surf){
	surfTemp = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		surfTemp[ct] = surf[ct];
	}
}

double tHydroMet::getSurfTemp(int time){
	return surfTemp[time];
}

void tHydroMet::setAtmPress(double *atm){
	atmPress = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		atmPress[ct] = atm[ct];
	}
}

double tHydroMet::getAtmPress(int time){
	return atmPress[time];
}

void tHydroMet::setSkyCover(double *sky){
	skyCover = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		skyCover[ct]= sky[ct];
	}
}

double tHydroMet::getSkyCover(int time){
	return skyCover[time];
}

void tHydroMet::setRHumidity(double *rh){
	rHumidity = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		rHumidity[ct] = rh[ct];
	}
}

double tHydroMet::getRHumidity(int time){
	return rHumidity[time];
}


void tHydroMet::setVaporPress(double *rh){
	vaporPress= new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		vaporPress[ct] = rh[ct];
	}
}

double tHydroMet::getVaporPress(int time){
	return vaporPress[time];
}

void tHydroMet::setWindSpeed(double *wind){
	windSpeed = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		windSpeed[ct] = wind[ct];
	}
}

double tHydroMet::getWindSpeed(int time){
	return windSpeed[time];
}

void tHydroMet::setNetRad(double *nrad){
	netRad = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		netRad[ct] = nrad[ct];
	}
}

double tHydroMet::getNetRad(int time){
	return netRad[time];
}

void tHydroMet::setPanEvap(double *panE){
	panEvap = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		panEvap[ct] = panE[ct];
	}
}

double tHydroMet::getPanEvap(int time){
	return panEvap[time];
}

void tHydroMet::setRadGlobal(double *rad){
	RadGlobal = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		RadGlobal[ct] = rad[ct];
	}
}

double tHydroMet::getRadGlobal(int time){
	return RadGlobal[time];
}

void tHydroMet::setRadDirect(double *rad){
	RadDirect = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		RadDirect[ct] = rad[ct];
	}
}

double tHydroMet::getRadDirect(int time){
	return RadDirect[time];
}

void tHydroMet::setRadDiffuse(double *rad){
	RadDiffuse = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		RadDiffuse[ct] = rad[ct];
	}
}

double tHydroMet::getRadDiffuse(int time){
	return RadDiffuse[time];
}

void tHydroMet::setRainMet(double *rain){
	RainMet = new double[numTimes];
	for(int ct=0;ct<numTimes;ct++){
		RainMet[ct] = rain[ct];
	}
}

double tHydroMet::getRainMet(int time){
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
  BinaryWrite(rStr, numTimes);
  BinaryWrite(rStr, numParams);
  BinaryWrite(rStr, gmt);
  BinaryWrite(rStr, stationID);
  BinaryWrite(rStr, latitude);
  BinaryWrite(rStr, longitude);
  BinaryWrite(rStr, basinLat);
  BinaryWrite(rStr, basinLong);
  BinaryWrite(rStr, otherVariable);
  BinaryWrite(rStr, meanTemp);
  for (int i = 0; i < numTimes; i++) {
    BinaryWrite(rStr, year[i]); 
    BinaryWrite(rStr, month[i]);
    BinaryWrite(rStr, day[i]);
    BinaryWrite(rStr, hour[i]);
  }
	// Giuseppe DEBUG Restart 2012 - START 
	// I have introduced an IF and ELSE to deal with the case
	// of PAN ET data
	if (numParams == 5){		
		for (int i = 0; i < numTimes; i++) {
			BinaryWrite(rStr, panEvap[i]);		
		}
	}
	else {// Giuseppe DEBUG Restart 2012 - END	
		for (int i = 0; i < numTimes; i++) {
			BinaryWrite(rStr, airTemp[i]);
			BinaryWrite(rStr, dewTemp[i]);
			BinaryWrite(rStr, surfTemp[i]);
			BinaryWrite(rStr, rHumidity[i]);
			BinaryWrite(rStr, skyCover[i]);
			BinaryWrite(rStr, windSpeed[i]);
			BinaryWrite(rStr, atmPress[i]);
			BinaryWrite(rStr, vaporPress[i]);
			BinaryWrite(rStr, RadGlobal[i]);
			//BinaryWrite(rStr, RadDirect[i]);
			//BinaryWrite(rStr, RadDiffuse[i]);
			//BinaryWrite(rStr, RainMet[i]);
			//BinaryWrite(rStr, panEvap[i]);
			if (numParams >= 11)
				BinaryWrite(rStr, netRad[i]);
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
  BinaryRead(rStr, numTimes);
  BinaryRead(rStr, numParams);
  BinaryRead(rStr, gmt);
  BinaryRead(rStr, stationID);
  BinaryRead(rStr, latitude);
  BinaryRead(rStr, longitude);
  BinaryRead(rStr, basinLat);
  BinaryRead(rStr, basinLong);
  BinaryRead(rStr, otherVariable);
  BinaryRead(rStr, meanTemp);
  for (int i = 0; i < numTimes; i++) {
    BinaryRead(rStr, year[i]);
    BinaryRead(rStr, month[i]);
    BinaryRead(rStr, day[i]);
    BinaryRead(rStr, hour[i]);
  }
	// Giuseppe DEBUG Restart 2012 - START 
	// I have introduced an IF and ELSE to deal with the case
	// of PAN ET data
	if (numParams == 5){		
	for (int i = 0; i < numTimes; i++) {
		BinaryRead(rStr, panEvap[i]);		
	}
  }
  else {// Giuseppe DEBUG Restart 2012 - END	
	  for (int i = 0; i < numTimes; i++) {
		  BinaryRead(rStr, airTemp[i]);
		  BinaryRead(rStr, dewTemp[i]);
		  BinaryRead(rStr, surfTemp[i]);
		  BinaryRead(rStr, rHumidity[i]);
		  BinaryRead(rStr, skyCover[i]);
		  BinaryRead(rStr, windSpeed[i]);
		  BinaryRead(rStr, atmPress[i]);
		  BinaryRead(rStr, vaporPress[i]);
		  BinaryRead(rStr, RadGlobal[i]);
		  //BinaryRead(rStr, RadDirect[i]);
		  //BinaryRead(rStr, RadDiffuse[i]);
		  //BinaryRead(rStr, RainMet[i]);
		  //BinaryRead(rStr, panEvap[i]);
		  if (numParams >= 11)
			  BinaryRead(rStr, netRad[i]);
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
