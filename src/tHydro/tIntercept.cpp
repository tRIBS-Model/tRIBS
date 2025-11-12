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
**  tIntercept.h: Header file for class tIntercept 
**
***************************************************************************/


#include "src/tHydro/tIntercept.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//            Section 1: Constructor and Initialization Functions
//
//
//=========================================================================

tIntercept::tIntercept(){
	gridPtr = 0;
}

tIntercept::tIntercept(SimulationControl *simCtrPtr, tMesh<tCNode> *gridRef, 
					   tInputFile &inFile, tRunTimer *tmrptr, tResample *resamp,
					   tHydroModel *hydro)
{
	gridPtr = gridRef;
	timer   = tmrptr;
	SetIntercpVariables( inFile, hydro );
}

void tIntercept::SetIntercpVariables(tInputFile &inFile, tHydroModel *hydro)
{
	interceptOption = inFile.ReadItem(interceptOption,"OPTINTERCEPT");
	maxInterStormPeriod = inFile.ReadItem(maxInterStormPeriod, "INTSTORMMAX");
	maxInterStormPeriod = maxInterStormPeriod/2;
	metTime = inFile.ReadItem(metTime, "METSTEP");

	luOption = inFile.ReadItem(luOption, "OPTLANDUSE"); // SKYnGM2008LU: added by SY, TM: 11/19/07
	nParmLU = 0; 

	if(interceptOption != 0){
		landPtr = hydro->landPtr;  

		// SKYnGM2008LU 
		if (luOption == 1 || luOption == 2) {
			inFile.ReadItem(luFile, "LUGRID"); 
			readLUGrid(luFile);
		}

	}
}

tIntercept::~tIntercept()
{
    // Only attempt to delete grid info if interception was active in the first place
	if (interceptOption != 0) { 
		if (luOption == 1 || luOption == 2) { 
			DeleteIntercept();
		}
	}

	Cout<<"tIntercept Object has been destroyed..."<<endl;
}

/***************************************************************************
**
** tIntercept::DeleteIntercept()
**
** Auxiliary function used by the destructor
**  
***************************************************************************/
void tIntercept::DeleteIntercept()
{

	for (int sz=0;sz<nParmLU;sz++) {
		delete [] LUgridParamNames[sz];
		delete [] LUgridBaseNames[sz];
		delete [] LUgridExtNames[sz];
	}
	delete [] LUgridParamNames;
	delete [] LUgridBaseNames;
	delete [] LUgridExtNames;

	return;
}

// SKYnGM2008LU
/***************************************************************************
**
** tIntercept::readLUGrid() Function
**
** Reads a file (*.gdf) from LUGRID keyword containing the base names of 
** the various input land use parameter grids along with the extension
** used for the filename. These follow a string that identifies the line
** with the parameters (ie. AL,TF,VH,SR,VF,CS,IC,CC,DC,DE,OT,LA). If no data 
** available for any parameters, the string NO_DATA should be input 
** instead of the path name and extension name. In this version, a single
** value of latitude, longitude and GMT is used for entire grids. The
** impact of this assumption should be small for small grids.
**
** Number of parameters 
** Latitude Longitude GMT
** ParamName Base Name1 Extension1 (../PATH/parameterName1)
** ParamName Base Name2  Extension2
** ...
** Example:
** 3
** 32.5 -100.3 -6
** AL ../PATH/ALbase txt
** TF ../PATH/TFbase txt
** VH NO_DATA NO_DATA
**
***************************************************************************/
void tIntercept::readLUGrid(char *gridFile) 
{
	int numParameters;
	double LUgridlat, LUgridlong;
	int LUgridgmt;
	
	cout<<"\nReading Land-Use Data Grid File in tIntercept: "; 
	cout<< gridFile<<"..."<<endl<<flush;
	
	ifstream readFile(gridFile); 
	if (!readFile) {
		cout << "\nFile "<<gridFile<<" not found!" << endl;
		cout << "Exiting Program...\n\n"<<endl;
		exit(1);}
	
	readFile >> numParameters; 
	nParmLU = numParameters;
	readFile >> LUgridlat;
	readFile >> LUgridlong;
	readFile >> LUgridgmt;

	LUgridBaseNames = new char*[numParameters];
	LUgridExtNames = new char*[numParameters];
	LUgridParamNames = new char*[numParameters];
	
	for (int ct=0;ct<numParameters;ct++) {
		LUgridParamNames[ct] = new char[10];
		LUgridBaseNames[ct] = new char[kName];
		LUgridExtNames[ct] = new char[10];
		readFile >> LUgridParamNames[ct];
		readFile >> LUgridBaseNames[ct];
		readFile >> LUgridExtNames[ct];
	}
	return;
}

//=========================================================================
//
//
//               Section 2: External Calling Function
//
//
//=========================================================================

/***************************************************************************
**
** tIntercept::callInterception() Function
**
** Sets the Interception Option in the Constructor
** Error message for non-valid options in InputFile
** Calls the appropriate Method InterceptGray() or InterceptRutter()
**
***************************************************************************/
void tIntercept::callInterception(tCNode *cNode, double Ep)
{
	if(interceptOption == 1){
		InterceptGray(cNode);
	}
	else if(interceptOption == 2){
		InterceptRutter(cNode, Ep);
	}
	else{
		cout << "\nInterception Option "<<interceptOption<<" not valid."<<endl;
		cout << "\tPlease use :" << endl;
		cout << "\t\t(1) for Gray (1970) Method: Two Parameter Model"<<endl;
		cout << "\t\t(2) for Rutter (1971) Method: Four Parameter Model"<<endl;
		cout << "Exiting Program...\n\n"<<endl;
		exit(1);
	}

	return;
}

//=========================================================================
//
//
//               Section 3: Interception Functions
//
//
//=========================================================================


/***************************************************************************
**
** tIntercept::InterceptGray() Function based on Gray(1970) 
** 	(Also see Bras,1993, p233)
**
** Uses coeffA, coeffB as parameters read from Land Use Table
**
**    I(t) = R(t)          while cumI <= coeffA    (mm/hr)
**    I(t) = coeffB*R(t)   while cumI >  coeff
**
***************************************************************************/
void tIntercept::InterceptGray(tCNode *cNode)
{ 
	double cumIntercept, intercept, rainfall, interStormLength;
	double netPrecipitation;
	double minRainAmount = 1.0;
	
	SetIntercpParameters( cNode );
	
	
	
	rainfall = cNode->getRain();
	
	if (rainfall <= minRainAmount)                //Modify InterStorm Length
		cNode->setStormLength(1, timer->getEtIStep());
	else
		cNode->setStormLength(0, 0);
	
	cumIntercept = cNode->getCumIntercept();      //Calculate Interception
	
	if (cumIntercept <= coeffA)
		intercept = rainfall;   
	else
		intercept = rainfall*coeffB;
	
	cNode->setInterceptLoss(intercept);           //For the _VEGETATED_ fract
	
	interStormLength = cNode->getStormLength();   //Set Cumulative Interception
	
	if (interStormLength < maxInterStormPeriod)   //For the _VEGETATED_ fract
		cNode->setCumIntercept(cumIntercept+intercept*timer->getEtIStep());
	else
		cNode->setCumIntercept(0.0);
	
	// Rate for the _ENTIRE_ cell:
	netPrecipitation = coeffV*(rainfall-intercept) + (1-coeffV)*rainfall;
	
	cNode->setNetPrecipitation(netPrecipitation); //Assign Net Precipitation
	
	return;
}

/***************************************************************************
**
** tIntercept::InterceptRutter() Function based on Rutter et al(1971,1975) 
**                            Also see Shuttleworth(1979,1988)
**
** Uses coeffP, coeffS, coeffK, coeffb
**
**     dC(t)/dt = (1-coeffP)*R(t) - D(t) - (C(t)/coeffS)*(Ep(t))  
**     Ep(t)  (mm/hr)
**
**     D(t) = 60*coeffK*exp(coeffb*(C(t)-coeffS))   (mm/hr)
**     T(t) = coeffP*R(t)
**
**     netR(t) = D(t)+T(t)                          (mm/hr)
**
**     I(t) = R(t) - netR(t)                        (mm/hr)
**       
**     Note: Interception in this context is a mixture of residual canopy 
**           storage and evaporation loss. Need to account for difference
**           when considering evapotranspiration scheme. Will leave as
**           lumped quantity for now. 
**
** According to Shuttleworth: coeffP ~ 0.25 (canopy), 0.02 - 0.2 (trunk)
**                            coeffS ~ 1 mm (canopy), 0.1 mm (trunk)
**                            coeffK ~ 0.002mm/min
**                            coeffb ~ 4 mm-1
**
***************************************************************************/
void tIntercept::InterceptRutter(tCNode *cNode, double evaporation)
{
	double rainfall, prevStorage, currentStorage, can_flx;
	double throughfall, drainage, netPrecipitation, interception;
	double interStormLength, cumIntercept, ctos;
	double minRainAmount = 1.0;                   //Minimum rainfall rate (mm/hr)
	double dt = timer->getEtIStep(); 			  // Timestep in hours
	
	rainfall     = cNode->getRain();

	prevStorage  = cNode->getCanStorage(); //Refers only to the VEGETATED frct
	ctos         = getCtoS(cNode);
	if (ctos > 1)
		ctos = 0;
	
	if (rainfall <= minRainAmount)                 //Modify InterStorm Length
		cNode->setStormLength(1, timer->getEtIStep());
	else
		cNode->setStormLength(0, 0);

	//If simulation can be skipped - just skip it
	if (prevStorage > evaporation*ctos*timer->getEtIStep() || rainfall > 0) { // dt here?

		cumIntercept = cNode->getCumIntercept();
		SetIntercpParameters( cNode );
		
		// Runge-Kutta Solution C(t) when C(t-1) is <=> S
		drainage = 0.0;                       
		currentStorage = storageRungeKutta(prevStorage, rainfall, 
										   evaporation, &drainage);

		double rain_on_canopy_vol = (1 - coeffP) * rainfall * dt;

		// The total water available to evaporate is the initial storage plus rainfall.
        double total_available_water = prevStorage + rain_on_canopy_vol;
        if (currentStorage > total_available_water) {
             // This can happen due to RK4 overshoot. If so, cap it.
             currentStorage = total_available_water;
        }

		// The RK4 solver can undershoot and produce a non-physical negative storage.
		// Correct this at the source before it propagates.
		if (currentStorage < 0.0) {
			currentStorage = 0.0;
		}
		
		// Calculate the physically possible evaporation RATE
		double avg_storage = (prevStorage + currentStorage) / 2.0;
		double avg_drain_rate_prop = (avg_storage > 0) ? coeffK * exp(coeffb * (avg_storage - coeffS)) : 0.0;
		double evapWetCanopy_rate; // mm/hr;
		if (avg_storage < coeffS) {
			evapWetCanopy_rate = (avg_storage / coeffS) * evaporation; // 'evaporation' is Ep
		} else {
			evapWetCanopy_rate = evaporation;
		}
		if (evapWetCanopy_rate < 0) evapWetCanopy_rate = 0.0;

		// Calculate the total water lost from the canopy.
        double storage_change_vol = currentStorage - prevStorage;      // [mm]
		double total_loss_vol = rain_on_canopy_vol - storage_change_vol; // [mm]
        if (total_loss_vol < 0.0) total_loss_vol = 0.0; // Loss cannot be negative

		// Calculate the evaporation VOLUME for the timestep
		double evapWetCanopy_vol = evapWetCanopy_rate * dt; // [mm]

		// The actual evaporated volume cannot be more than the total water lost.
        if (evapWetCanopy_vol > total_loss_vol) evapWetCanopy_vol = total_loss_vol;
        if (evapWetCanopy_vol < 0.0) evapWetCanopy_vol = 0.0;

		// Calculate drainage volume as the residual of the water balance
		double drainage_vol = total_loss_vol - evapWetCanopy_vol; // [mm]
        if (drainage_vol < 0.0) drainage_vol = 0.0;

		double residual = 0.0;

		// To avoid Inf time to C = 0 -> dump it if lower threshold
		if (currentStorage * 100 / coeffS <= 3) { // dump if less than 3% of coeffS
			// Check if adding the residual to evaporation would violate the energy budget.
			double potential_evap_vol = evaporation * dt;
			if ((evapWetCanopy_vol + currentStorage) <= potential_evap_vol) {
				// Energetically possible to evaporate the residual 
				evapWetCanopy_vol += currentStorage;
			} else {
				// Not enough energy. Send residual to drainage instead. 
				drainage_vol += currentStorage;
			}
			currentStorage = 0.0; // Empty the canopy in either case.
		}
		
		// Rainfall is constant during the dt interval
		throughfall = coeffP*rainfall; 
		
		// Drainage represents an average value over the interval
		// For the _VEGETATED_ fraction of a cell:
		netPrecipitation = (drainage_vol/dt) + throughfall;
		
		if (netPrecipitation <= rainfall)
			interception = rainfall - netPrecipitation;
		else
			interception = 0.0;
		
		// Rate for the _ENTIRE_ cell:
		netPrecipitation = coeffV*netPrecipitation + (1-coeffV)*rainfall;
		
		if(currentStorage < 0.0)
			currentStorage = 0.0;
			
		// Set the dynamic variables to tCNode
		cNode->setNetPrecipitation(netPrecipitation); //For the _ENTIRE_ fract
		cNode->setInterceptLoss(interception);        //For the _VEGETATED_ fract
		cNode->setCanStorage(currentStorage);         //For the _VEGETATED_ fract
		cNode->setEvapWetCanopy(evapWetCanopy_vol/dt);    //For the _VEGETATED_ fract

		interStormLength = cNode->getStormLength();
		if(interStormLength < maxInterStormPeriod)
			cNode->setCumIntercept(cumIntercept + interception*timer->getEtIStep());
		else
			cNode->setCumIntercept(0.0);

	}
	else {
		cNode->setNetPrecipitation(0.0); 
		cNode->setInterceptLoss(0.0);
		cNode->setCanStorage(0.0);
		cNode->setCumIntercept(0.0);
		cNode->setEvapWetCanopy(0.0);

	}
	return;
}

/***************************************************************************
** 
** tIntercept::storageRungeKutta() and RutterFn() Functions
**
** Calculates the storage C(t) when C(t-1) is known < S using
** the Rutter Formula given in Rutter(1971), Eq 13a.
** 
** Numerical Method given non-linear ODE in C(t) adapted from
** Lerman (1993) p 117. 
** 
***************************************************************************/
double tIntercept::RutterFn(double t, double C, double R, double Ep) 
{
	double diffC;
	// If _less_ than the canopy field capacity S:
	if (C < coeffS)
		diffC = ((1-coeffP)*R - (C/coeffS)*Ep - coeffK*exp(coeffb*(C-coeffS)));
	
	// Else -> the surface is wet -> evaporation rate = potential:
	else 
		diffC = ((1-coeffP)*R - Ep - coeffK*exp(coeffb*(C-coeffS)));
	return diffC;
}
double tIntercept::storageRungeKutta(double prevStore, double R, double Ep,
									 double *drainage)
{
	int cnt = 0;
	double ta, tb, tc, td;
	double h, t0, tlast, c0, c00, numiter;
	
	numiter = ((floor)(R*3.0));
	if (numiter < 50.0) numiter = 50.0;

	h = timer->getEtIStep()/numiter; 
	t0 = 0.0;                     
	tlast = timer->getEtIStep();
	c0 = prevStore;
	
	for (t0=0.0; t0 < tlast; t0+=h) {
		ta = h*RutterFn(t0, c0, R, Ep);
		tb = h*RutterFn(t0+h/2.0, c0+ta/2.0, R, Ep);
		tc = h*RutterFn(t0+h/2.0, c0+tb/2.0, R, Ep);
		td = h*RutterFn(t0+h, c0+tc, R, Ep);
		c00 = c0;
		c0 += (ta/6.0 + tb/3.0 + tc/3.0 + td/6.0);
		
		*drainage = *drainage + coeffK*exp(coeffb*( (c00+c0)/2 - coeffS ));
		cnt++;
	}
	*drainage = *drainage/cnt;
	return c0;
}

/***************************************************************************
**
** tIntercept::SetIntercpParameters(tCNode *cNode) Function
**
** Assigns values of parameters for the current Voronoi cell
**
***************************************************************************/
void tIntercept::SetIntercpParameters(tCNode *cNode)
{
	landPtr->setLandPtr(cNode->getLandUse());
	
	if (interceptOption == 1) {
		coeffA = landPtr->getLandProp(1);  //Storage Capacity
		coeffB = landPtr->getLandProp(2);  //Interception Coefficient
		coeffV = landPtr->getLandProp(11); //Vegetation Fraction 
	}
	else if (interceptOption == 2) {
		coeffP = landPtr->getLandProp(3);          //Free Throughfall Coefficient
		coeffS = landPtr->getLandProp(4);          //Canopy Storage Capacity
		coeffK = landPtr->getLandProp(5);          //Drainage Coefficient
		coeffb = landPtr->getLandProp(6);          //Drainage Exponential Parameter
		coeffV = landPtr->getLandProp(11);         //Vegetation Fraction 
	}

	if (luOption == 1 || luOption == 2) {
		for (int ct=0;ct<nParmLU;ct++) { 
			if (strcmp(LUgridParamNames[ct],"CS")==0) {
				if (interceptOption == 1) {
					coeffA = cNode->getCanStorParam();  //Canopy Storage Parameter
				}
			}
			if (strcmp(LUgridParamNames[ct],"IC")==0) {
				if (interceptOption == 1) {
					coeffB = cNode->getIntercepCoeff();  //Interception Coefficient
				}
			}
			if (strcmp(LUgridParamNames[ct],"TF")==0) {
				if (interceptOption == 2) {
					coeffP = cNode->getThroughFall();  //Free Throughfall Coefficient
				}
			}
			if (strcmp(LUgridParamNames[ct],"CC")==0) {
				if (interceptOption == 2) {
					coeffS = cNode->getCanFieldCap();  //Canopy Storage Capacity
				}
			}
			if (strcmp(LUgridParamNames[ct],"DC")==0) {
				if (interceptOption == 2) {
					coeffK = cNode->getDrainCoeff();  //Drainage Coefficient
				}
			}
			if (strcmp(LUgridParamNames[ct],"DE")==0) {
				if (interceptOption == 2) {
					coeffb = cNode->getDrainExpPar();  //Drainage Exponential Parameter
				}
			}
			if (strcmp(LUgridParamNames[ct],"VF")==0) {
				if ( (interceptOption == 1) || (interceptOption == 2) ) {
					coeffV = cNode->getVegFraction();  //Vegetation Fraction 
				}
			}
		}
	}
	
	return;
}

/***************************************************************************
**
** tIntercept::getIoption() Function
**
***************************************************************************/
int tIntercept::getIoption()
{
	return interceptOption;
}

/***************************************************************************
**
** tIntercept::getCtoS() Function
** 
** Returns relation of current storage to canopy field capacity
**
***************************************************************************/
double tIntercept::getCtoS(tCNode *cNode)
{
	double ctos {};
	if (interceptOption == 1)
		ctos = 1;
	else if (interceptOption == 2) {
		landPtr->setLandPtr(cNode->getLandUse());
		coeffS = landPtr->getLandProp(4);          //Canopy Field Capacity
		ctos = cNode->getCanStorage()/coeffS;
	}
	return ctos;
}

/***************************************************************************
**
** tIntercept::IsThereCanopy() Function
** 
** Returns '1' if there is canopy interception, '0'- rainfall passes through
**
***************************************************************************/
int tIntercept::IsThereCanopy(tCNode *cNode)
{
	int answer {};
	if (interceptOption == 1)
		answer = 1;
	else if (interceptOption == 2) {
		landPtr->setLandPtr(cNode->getLandUse());
		coeffP = landPtr->getLandProp(3);       //Free Throughfall Coefficient
		if (coeffP < 0.999)
			answer = 1;
		else 
			answer = 0;
	}
	return answer;
}

/***************************************************************************
**
** tIntercept::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/
void tIntercept::writeRestart(fstream & rStr) const
{
  BinaryWrite(rStr, interceptOption);
  BinaryWrite(rStr, maxInterStormPeriod);
  BinaryWrite(rStr, metTime);
  BinaryWrite(rStr, coeffA);
  BinaryWrite(rStr, coeffB);
  BinaryWrite(rStr, coeffP);
  BinaryWrite(rStr, coeffS);
  BinaryWrite(rStr, coeffK);
  BinaryWrite(rStr, coeffb);
  BinaryWrite(rStr, coeffV);
}

/***************************************************************************
**
** tIntercept::readRestart() Function
**
***************************************************************************/
void tIntercept::readRestart(fstream & rStr)
{
  BinaryRead(rStr, interceptOption);
  BinaryRead(rStr, maxInterStormPeriod);
  BinaryRead(rStr, metTime);
  BinaryRead(rStr, coeffA);
  BinaryRead(rStr, coeffB);
  BinaryRead(rStr, coeffP);
  BinaryRead(rStr, coeffS);
  BinaryRead(rStr, coeffK);
  BinaryRead(rStr, coeffb);
  BinaryRead(rStr, coeffV);
}

//=========================================================================
//
//
//                      End of tIntercept.cpp
//
//
//=========================================================================
