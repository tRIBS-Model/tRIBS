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
**  tResData.cpp: Functions for class tResData (see tResData.h)
**
***************************************************************************/

#include "src/tFlowNet/tResData.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//                  Section 1: tResData Constructors and Destructors
//
//
//=========================================================================

tResData::tResData()
{
    std::cerr << "Constructing tResData object \n"; //-WR debug
	const int rSIZE = 500;
	ResType = new int[rSIZE];
	rElev = new double[rSIZE];
	rDischarge = new double[rSIZE];
	rStorage = new double[rSIZE];
	rInflow = new double[1000]; //1000 is just in the initialization (See: 'setResArraySize()')
	rEDS = new double[rSIZE];
	int setNum = 0;
	int EDScount = 0;
	int ResLines = 0;
	int ResNodeID = 0;
	int ResIDtype = 0;
	double ResInElev = 0;
	int rStep = 0;

	rSTQnext = new double[1000]; //1000 is just in the initialization (See: 'setResArraySize()')

}

tResData::~tResData()
{
        std::cerr << "Destructing tResData object \n"; //-WR debug
        delete ResType;
        delete rElev;
        delete rDischarge;
        delete rStorage;
        delete rInflow;
        delete rEDS;
        delete rSTQnext;
}

//=========================================================================
//
//
//                  Section 2: tResData Functions
//
//
//=========================================================================

/***************************************************************************
**
** Set() and Get() Functions for Reservoir Data
**
***************************************************************************/
int tResData::getnumLines(char *resfile){
	int numLines = 0;
 	string line;
	ifstream readFile(resfile);
		while (getline(readFile, line)){
       			++numLines;
		}
	return numLines;
}

void tResData::setRNum(int rN){
	setNum = rN;
	rStep = rN;
}

int tResData::getRNum(){
	return setNum;
}

void tResData::setResType(int rt){
	ResType[setNum] = rt;
}

int tResData::getResType(int type){
	return ResType[type];
}

void tResData::setResElev(double re){
	rElev[setNum] = re;
}

double tResData::getResElev(int elev){
	return rElev[elev];
}

void tResData::setResDischarge(double rd){
	rDischarge[setNum] = rd;
}

double tResData::getResDischarge(int dis){
	return rDischarge[dis];
}

void tResData::setResStorage(double rs){
	rStorage[setNum] = rs;
	setNum++; // Counter set here since Storage is the last value to be read.
}

double tResData::getResStorage(int stor){
	return rStorage[stor];
}

void tResData::setResEDS(double reds, int EDScount){
	rEDS[EDScount] = reds;
}

double tResData::getResEDS(int eldis){
	return rEDS[eldis];
}

void tResData::setResLines(int rl){
	ResLines = rl;
}

int tResData::getResLines(){
	return ResLines;
}

/****************************************************/
/******************RES*Polygon*ID********************/
void tResData::setResNodeID(int rID){
	ResNodeID = rID;
}

int tResData::getResNodeID(){
	return ResNodeID;
}

void tResData::setResNodeType(int rIDtype){
	ResIDtype = rIDtype;
}

int tResData::getResNodeType(){
	return ResIDtype;
}

void tResData::setInitial_H(double rH){
	ResInElev = rH;
}

double tResData::getInitial_H(){
	return ResInElev;
}

void tResData::setResArraySize(int arrSize){
	rInflow = new double[arrSize];
	rSTQnext = new double[arrSize]; 
// Array size is estimated based on RUNTIME and TIMESTEP
// See tKinemat.cpp Line 736
}

/******************END*Polygon*ID********************/
/****************************************************/


void tResData::setSTQnext(double rSTQ, int STQstep)
{
	rSTQnext[STQstep] = rSTQ;
}

double tResData::getSTQnext(int timeSTQ)
{
	return rSTQnext[timeSTQ];
}

int tResData::getRoutingStep(){
	rStep++;
	return rStep;
}

void tResData::setInflow(double rIn){
	if (rStep == 1){
		rInflow[0] = 0.0;
		rInflow[rStep] = rIn;
	}
	else {
		rInflow[rStep] = rIn;
	}
}

double tResData::getInflow(int ResIn){
	return rInflow[ResIn];
}

//=========================================================================
//
//
//                          End tResData.cpp
//
//
//=========================================================================
