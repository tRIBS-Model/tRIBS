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
**  tReservoir.cpp: Functions for class tReservoir (see tReservoir.h)
**           A Finite-Element Kinematic Wave Routing Algorithm
**
***************************************************************************/

#include "src/tFlowNet/tReservoir.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//                  Section 1: tReservoir Constructors/Destructors
//
//
//=========================================================================

/****************************************************************************
**  
**  tReservoir::tReservoir()
**
**  Constructor for tRIBS model use
**
*****************************************************************************/
tReservoir::tReservoir()
{

}


tReservoir::tReservoir(tInputFile &inFile)
{
	nReservoirs = 0;
	SetResVariables(inFile);
	SetResNodes(inFile);
}

/****************************************************************************
**  
**  tReservoir::~tReservoir()
**
**  Destructor
**
*****************************************************************************/
tReservoir::~tReservoir()
{
    /*TODO: this code was causing a runtime error [EXC_BAD_ACCESS(Code=1)].
    I think this was because it was using delete [], which is used for arrays of objects
    also tReservoir class creates a pointer to tResData that should be initialized as a null pointer in the case there is no reservoir*/

    delete reservoirNodes;
    delete reservoirTypes;

}

//=========================================================================
//
//
//                  Section 2: tReservoir Functions
//
//
//=========================================================================

/*****************************************************************************\
**  
**  tReservoir::RunLevelPoolRouting()
**
**  Runs the Reservoir component using the Level Pool Routing method.
**  
\*****************************************************************************/
void tReservoir::RunLevelPoolRouting(double Qin)
{	
	/* Defines from which tResData class TYPE will the code read from */
	rType = getCurrResType();

	/* Defines from which tResData class NODE will the code read from */
	rNode = getCurrResNode();	

	ComputeInflow(Qin);
	ComputeSTQnext();
			 
	return;
}

/*****************************************************************************\
**  
**  tReservoir::ComputeInflow()
**
**  Computes the Inflow into the reservoir from the upstream boundary node
**  in the current time step and the previous one.
**  
\*****************************************************************************/
void tReservoir::ComputeInflow(double Qinflow)
{
	// Call a function to Get current time step
	RStep = reservoirNodes[rNode].getRoutingStep();
	reservoirNodes[rNode].setInflow(Qinflow);

	if (RStep == 1) {
		ResQinflow2 = Qinflow;
		resInflow = 0 + ResQinflow2; // At initial time step Inflow = 0
	}

	else {
		ResQinflow2 = Qinflow; // Get inflow at the upper BND node
		ResQinflow = reservoirNodes[rNode].getInflow(RStep-1);
		resInflow = ResQinflow + ResQinflow2; // = Ij + Ij+1
	}

	return;
}

/***************************************************************************
**
** tReservoir::ComputeInitialSTQ()
**
** Function to compute the inital [2Sj/dt - Qj] based on the initial water 
** level at the Reservoir specified by the user.
**
** where S = Storage; dt = timestep; Q = discharge. 
** (STQ = Storage/Time/Discharge)
**
***************************************************************************/
void tReservoir::ComputeInitialSTQ()
{
	lengthH = reservoirTypes[rType].getResLines();

	//Read provided table for Elevation-Storage-Discharge data
	initialH = reservoirNodes[rNode].getInitial_H();
	
	for (int h=0; h!=lengthH; h++) {
		elevData = reservoirTypes[rType].getResElev(h);
		if (initialH > elevData)
			continue;
		else
			interNum = h;
			break;
	}
	
	//Get Discharge Q and Storage S at elevation H by linear interpolation.
	resQ = reservoirTypes[rType].getResDischarge(interNum-1);
	resH = reservoirTypes[rType].getResElev(interNum-1);
	resS = reservoirTypes[rType].getResStorage(interNum-1);
	resQ2 = reservoirTypes[rType].getResDischarge(interNum);
	resH2 = reservoirTypes[rType].getResElev(interNum);
	resS2 = reservoirTypes[rType].getResStorage(interNum);

	initialS = resS + (resS2-resS)*((initialH-resH)/(resH2-resH));
	initialQ = resQ + (resQ2-resQ)*((initialH-resH)/(resH2-resH)); 
	
	double timestepUsed = getModelTimeStep();

	//Obtain the Storage-Discharge function (2S/dt + Q) for the first time step.
	STQ_0 = resInflow + ((2.0*initialS/timestepUsed)-initialQ); //Read timestep used in the model
	
	return;
}

/***************************************************************************
**
** tReservoir::ComputeSTQnext()
**
** Function to compute the [2Sj+1/dt - Qj+1] based on the Outflow.
**
***************************************************************************/
void tReservoir::ComputeSTQnext() 
{
	// Computes discharge for initial conditions or subsequent time steps
	if (RStep == 1){
		ComputeInitialSTQ();
		ComputeResQ();

		STQnext = STQ_0 - 2.0*Q_0;
		reservoirNodes[rNode].setSTQnext(STQnext, RStep);
	}
	
	else {
		ComputeSTQ();
		ComputeResQ();

		STQnext = STQ - 2.0*Q_0;
		reservoirNodes[rNode].setSTQnext(STQnext, RStep);
	}
		
	return;
}

/***************************************************************************
**
** tReservoir::ComputeSTQ()
**
** Function to compute the Storage-Discharge relation [2Sj+1/dt + Qj+1]
** of the Reservoir at the next time step.
**
***************************************************************************/
void tReservoir::ComputeSTQ()
{
	STQprev = reservoirNodes[rNode].getSTQnext(RStep-1); //Read STQ from previous time step
	STQ = resInflow + STQprev;

	return;
}

/***************************************************************************
**
** tReservoir::ComputeResQ()
**
** Function to compute the Outflow [Qj+1] from the Reservoir using linear
** interpolation from the provided data table.
**
***************************************************************************/
void tReservoir::ComputeResQ()
{
	if (RStep == 1){
		Q_0 = initialQ;
		setResDischargeOut(Q_0);
		setResElevOut(initialH);
	}

	else {
		for (int x=0; x!=lengthH; x++) {
			EDSdata = reservoirTypes[rType].getResEDS(x);
			if (STQ > EDSdata)
				continue;
			else
				interNum2 = x;
				break;
		}

	resQ = reservoirTypes[rType].getResDischarge(interNum2-1);
	resQ2 = reservoirTypes[rType].getResDischarge(interNum2);
	EDSdata = reservoirTypes[rType].getResEDS(interNum2-1);
	EDSdata2 = reservoirTypes[rType].getResEDS(interNum2);
	resH = reservoirTypes[rType].getResElev(interNum2-1);
	resH2 = reservoirTypes[rType].getResElev(interNum2);

	Q_0 = resQ + (resQ2 - resQ)*((STQ-EDSdata)/(EDSdata2-EDSdata));
	H_0 = resH + (resH2 - resH)*((STQ-EDSdata)/(EDSdata2-EDSdata));

	setResDischargeOut(Q_0);
	setResElevOut(H_0);
	}

	return;
}

//=========================================================================
//
//
//                  Section 3: tReservoir Read File Functions
//
//
//=========================================================================

/***************************************************************************
**
**  tReservoir::SetResVariables(tInputFile &inFile)
**
**  Initializes tResData object
**
***************************************************************************/
void tReservoir::SetResVariables(tInputFile &inFile)
{
	inFile.ReadItem(resfile, "RESDATA");
	readReservoirFile(resfile);

	return;
}

/***************************************************************************
**
**  tReservoir::SetResNodes(tInputFile &inFile)
**
**  Initializes tResData object
**
***************************************************************************/
void tReservoir::SetResNodes(tInputFile &inFile)
{
	inFile.ReadItem(resNodeFile, "RESPOLYGONID");
	readResNodeFile(resNodeFile);

	return;
}

/***************************************************************************
**
** tReservoir::readResNodeFile() Function
**
**
** Reads the Reservoir Polygon ID File which provides information concerning
** the selected nodes to be represented as Reservoirs.
**
** Format for the Reservoir Polygon ID File:
**
** Header:
** nReservoirs	nNodeParams (3)
**
** Body:
** NodeID	ResNodeType	Initial_H
**
** Description:
** nReservoirs = # of Reservoirs; nNodeParams = # of parameters (3)
** NodeID       (int)	 Node selected by the user to be a Reservoir
** ResNodeType  (int)    Type of Reservoir associated to the Node
** Initial_H   (double)  Initial Water Surface Elevation at the Reservoir [m]
**
***************************************************************************/
void tReservoir::readResNodeFile(char *resNodeFile)
{
	int nNodeParams;
	int NodeID, ResNodeType;
	double Initial_H;

	Cout<<"\nReading Reservoir Polygon ID File '";
	Cout<< resNodeFile<<"'..."<<endl<<flush;

	ifstream readFile(resNodeFile); 
	if (!readFile) {
		cout << "File "<<resNodeFile<<" not found." << endl;
		cout<<"Exiting Program...\n\n"<<endl;
		exit(1);
	}

	readFile >> nReservoirs;
	setNReservoirs(nReservoirs);
	readFile >> nNodeParams;
	reservoirNodes = new tResData[nReservoirs];
	assert(reservoirNodes != 0);

	for (int count=0;count<nReservoirs;count++) {
		reservoirNodes[count].setResArraySize(getResArraySize());
		reservoirNodes[count].setRNum(0);
		for (int ct=0;ct<nNodeParams;ct++) {
			if (ct==0) {
				readFile >> NodeID;
				reservoirNodes[count].setResNodeID(NodeID);
			}
			if (ct==1) {
				readFile >> ResNodeType;
				reservoirNodes[count].setResNodeType(ResNodeType);
			}
			if (ct==2) {
				readFile >> Initial_H;
				reservoirNodes[count].setInitial_H(Initial_H);
			}
		}
	}
	readFile.close();
}

/***************************************************************************
**
** tReservoir::readReservoirFile() Function
**
**
** Reads the Reservoir File which provides information concerning
** the different types of Reservoir specified by the user. Creates an
** array of tResData objects for storing data. (see tResData.h)
**
** Format for the Reservoir Data File:
**
** Header:
** nTypes nResParams (4)
**
** Body:
** Type# Elevation Discharge Storage
**
** Description:
** nTypes = # of different types of reservoirs; nResParams = # of parameters (4)
** Type#       	     (int)	1->N
** Elevation   	     (double)  	Stage or water elevation at the Reservoir [m]
** Discharge   	     (double)  	Discharge for each elevation  [m^3/s]
** Storage     	     (double)  	Storage for each elevation  [1000 m^3]
**
***************************************************************************/
void tReservoir::readReservoirFile(char *resfile)
{
	int nTypes, nResParams, nLines;
	int ResType;
	double rElev, rDischarge, rStorage;
	double EDSdt, EDSstorage, EDSdischarge, EDSvalue;
	
	Cout<<"\nReading Reservoir Data File '";
	Cout<< resfile<<"'..."<<endl<<flush;
	
	ifstream readFile(resfile); 
	if (!readFile) {
		cout << "File "<<resfile<<" not found." << endl;
		cout<<"Exiting Program...\n\n"<<endl;
		exit(1);
	}

	readFile >> nTypes;
	readFile >> nResParams;

	reservoirTypes = new tResData[nTypes];
	nLines = reservoirTypes[0].getnumLines(resfile);
	assert(reservoirTypes != 0);

	int currType = 0; //Initializes Reservoir Type
	reservoirTypes[currType].setRNum(0);

	for (int countLine=0;countLine < (nLines-1);countLine++) { //Reads all lines from the file (excludes header)
		for (int ct=0;ct < nResParams;ct++) { // Reads parameters from each line (4)
			if (ct==0) {
				readFile >> ResType;
				if (ResType != currType) {
					currType++;
					reservoirTypes[currType].setResType(ResType);
				}
				else {
					reservoirTypes[currType].setResType(ResType);
					int NumType = reservoirTypes[currType].getRNum();
				}
			}
			if (ct==1) {
				readFile >> rElev;
				reservoirTypes[currType].setResElev(rElev);
			}
			if (ct==2) {
				readFile >> rDischarge;
				reservoirTypes[currType].setResDischarge(rDischarge);
				int NumDis = reservoirTypes[currType].getRNum();
			}
			if (ct==3) {
				readFile >> rStorage;
				reservoirTypes[currType].setResLines(reservoirTypes[currType].getRNum());
				reservoirTypes[currType].setResStorage(rStorage);
			}
		}
	}
	readFile.close();
 
	EDSdt = getModelTimeStep(); // Units in seconds
	for (int iter=0;iter < nTypes;iter++) {
		for (int param=0;param < reservoirTypes[iter].getResLines();param++) {
			EDSstorage = reservoirTypes[iter].getResStorage(param);
			EDSdischarge = reservoirTypes[iter].getResDischarge(param);
			EDSvalue = ((2*EDSstorage)/EDSdt) + EDSdischarge; // 2S/dt + Q
			reservoirTypes[iter].setResEDS(EDSvalue, param);
		}
	}
}
//=========================================================================
//
//
//                  Section 4: tReservoir Set/Get Functions
//
//
//=========================================================================

void tReservoir::setNReservoirs(int resNum){NumRes = resNum;}
int  tReservoir::getNReservoirs(){return NumRes;}

void tReservoir::setCurrResNode(int resN){ResNID = resN;}
int  tReservoir::getCurrResNode(){return ResNID;}

void tReservoir::setCurrResType(int resT){ResTyp = resT;}
int  tReservoir::getCurrResType(){return ResTyp;}

void tReservoir::setResDischargeOut(double resQout){ResOutflow = resQout;}
double  tReservoir::getResDischargeOut(){return ResOutflow;}

void tReservoir::setResElevOut(double resHout){ResOutElev = resHout;}
double  tReservoir::getResElevOut(){return ResOutElev;}

void tReservoir::setModelTimeStep(double resTimeStep){ResTimeDt = resTimeStep;}
double  tReservoir::getModelTimeStep(){return ResTimeDt;}

void tReservoir::setResArraySize(int resArray){ArraySize = resArray;}
int  tReservoir::getResArraySize(){return ArraySize;}

// Get function for NODE ID and Type access from tKinemat
int tReservoir::getNodetKinemat(int tKNode)
{
	return reservoirNodes[tKNode].getResNodeID();
}

int tReservoir::getTypetKinemat(int tKType)
{
	return reservoirNodes[tKType].getResNodeType();
}

//=========================================================================
//
//
//                          End tReservoir.cpp
//
//
//=========================================================================
