/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 * VERSION 5.2
 *
 * Copyright (c) 2024. tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/

/***************************************************************************
**
**  tEvapoTrans.cpp:   Function file for tEvapoTran class (see tEvapoTrans.h)
**
***************************************************************************/

#include "src/tHydro/tEvapoTrans.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//            Section 1: Constructor and Initialization Routines
//
//
//=========================================================================

/***************************************************************************
**
** tEvapoTrans() Constructor and Destructor
** 
** Point Data from HydroMeteorological Station
** Grid Data from HydroMetGrids 
**
***************************************************************************/
tEvapoTrans::tEvapoTrans()
{
    elevation =  0.0;
	simCtrl = nullptr;
    weatherSimul = nullptr;
	respPtr = nullptr;
	timer = nullptr;
	timeCount = 0.0;
	rainPtr = nullptr;
	weatherStations = nullptr;
	hydrPtr = nullptr;
	landPtr = nullptr;
	soilPtr = nullptr;
	currentTime = nullptr;
	assignedStation = nullptr;
	gridParamNames = nullptr;
	gridBaseNames = nullptr;
	gridExtNames = nullptr;
	LUgridParamNames = nullptr;
	LUgridBaseNames = nullptr;
	LUgridExtNames = nullptr;
	airpressure = nullptr;
	dewtemperature = nullptr;
	skycover = nullptr;
	windspeed = nullptr;
	airtemperature = nullptr;
	surftemperature = nullptr;
	netradiation = nullptr;
	incomingsolar = nullptr; //E.R.V. 3/6/2012
	evapotranspiration = nullptr;
	relhumidity = nullptr;
	vaporpressure= nullptr;
	LandUseAlbGrid = nullptr;
	ThroughFallGrid = nullptr;
	VegHeightGrid = nullptr;
	StomResGrid = nullptr;
	VegFractGrid = nullptr;
	CanStorParamGrid = nullptr;
	IntercepCoeffGrid = nullptr;
	CanFieldCapGrid = nullptr;
	DrainCoeffGrid = nullptr;
	DrainExpParGrid = nullptr;
	OptTransmCoeffGrid = nullptr;
	LeafAIGrid = nullptr;
	ALgridhours = nullptr;
	TFgridhours = nullptr;
	VHgridhours = nullptr;
	SRgridhours = nullptr;
	VFgridhours = nullptr;
	CSgridhours = nullptr;
	ICgridhours = nullptr;
	CCgridhours = nullptr;
	DCgridhours = nullptr;
	DEgridhours = nullptr;
	OTgridhours = nullptr;
	LAgridhours = nullptr;
	ALgridFileNames = nullptr;
	TFgridFileNames = nullptr;
	VHgridFileNames = nullptr;
	SRgridFileNames = nullptr;
	VFgridFileNames = nullptr;
	CSgridFileNames = nullptr;
	ICgridFileNames = nullptr;
	CCgridFileNames = nullptr;
	DCgridFileNames = nullptr;
	DEgridFileNames = nullptr;
	OTgridFileNames = nullptr;
	LAgridFileNames = nullptr;

	gridPtr = nullptr; 	nParmLU = 0;
}

tEvapoTrans::tEvapoTrans(SimulationControl *simCtrPtr, tMesh<tCNode> *gridRef,
						 tInputFile &infile, tRunTimer *t, tResample *resamp, 
						 tHydroModel *hydro, tRainfall *storm)
{
	timeCount = 0;	
	simCtrl = nullptr;
        weatherSimul = nullptr;
	gridPtr = nullptr;
	respPtr = nullptr;
	timer = nullptr;
	rainPtr = nullptr;
	weatherStations = nullptr;
	hydrPtr = nullptr;
	landPtr = nullptr;
	soilPtr = nullptr;
	currentTime = nullptr;
	assignedStation = nullptr;
	gridParamNames = nullptr;
	gridBaseNames = nullptr;
	gridExtNames = nullptr;
	LUgridParamNames = nullptr;
	LUgridBaseNames = nullptr;
	LUgridExtNames = nullptr;
	airpressure = nullptr;
	dewtemperature = nullptr;
	skycover = nullptr;
	windspeed = nullptr;
	airtemperature = nullptr;
	surftemperature = nullptr;
	netradiation = nullptr;
	incomingsolar = nullptr; //E.R.V. 3/6/2012
	evapotranspiration = nullptr;
	relhumidity = nullptr;
	vaporpressure = nullptr;
	LandUseAlbGrid = nullptr;
	ThroughFallGrid = nullptr;
	VegHeightGrid = nullptr;
	StomResGrid = nullptr;
	VegFractGrid = nullptr;
	CanStorParamGrid = nullptr;
	IntercepCoeffGrid = nullptr;
	CanFieldCapGrid = nullptr;
	DrainCoeffGrid = nullptr;
	DrainExpParGrid = nullptr;
	OptTransmCoeffGrid = nullptr;
	LeafAIGrid = nullptr;
	ALgridhours = nullptr;
	TFgridhours = nullptr;
	VHgridhours = nullptr;
	SRgridhours = nullptr;
	VFgridhours = nullptr;
	CSgridhours = nullptr;
	ICgridhours = nullptr;
	CCgridhours = nullptr;
	DCgridhours = nullptr;
	DEgridhours = nullptr;
	OTgridhours = nullptr;
	LAgridhours = nullptr;
	ALgridFileNames = nullptr;
	TFgridFileNames = nullptr;
	VHgridFileNames = nullptr;
	SRgridFileNames = nullptr;
	VFgridFileNames = nullptr;
	CSgridFileNames = nullptr;
	ICgridFileNames = nullptr;
	CCgridFileNames = nullptr;
	DCgridFileNames = nullptr;
	DEgridFileNames = nullptr;
	OTgridFileNames = nullptr;
	LAgridFileNames = nullptr; 	nParmLU = 0;

	gridPtr = gridRef;
	respPtr = resamp;
	rainPtr = storm;
	hydrPtr = hydro;
	timer = t;
	simCtrl = simCtrPtr;
	DeriveAspect();
	SetEvapTVariables(infile, hydro);
}

tEvapoTrans::~tEvapoTrans()
{
	DeleteEvapoTrans();
	Cout<<"tEvapoTrans Object has been destroyed..."<<endl;
}

/***************************************************************************
**
** tEvapoTrans::SetEvapTVariables() 
**
** Auxiliary function used by the constructor
**  
***************************************************************************/
void tEvapoTrans::SetEvapTVariables(tInputFile &infile, tHydroModel *hydro)
{
	evapotransOption = infile.ReadItem(evapotransOption,"OPTEVAPOTRANS");
	Ioption = infile.ReadItem(Ioption,"OPTINTERCEPT");
	metdataOption = infile.ReadItem(metdataOption, "METDATAOPTION");
	gFluxOption = infile.ReadItem(gFluxOption, "GFLUXOPTION");
	timeStep = infile.ReadItem(timeStep, "METSTEP");
	Tlinke = infile.ReadItem(Tlinke, "TLINKE");

	luOption = infile.ReadItem(luOption, "OPTLANDUSE"); // SKYnGM2008LU: added by SY, TM: 11/19/07
	if (luOption == 1) {
		luInterpOption = infile.ReadItem(luInterpOption, "OPTLUINTERP"); // SKYnGM2008LU
	}

	snowOption = infile.ReadItem(snowOption, "OPTSNOW"); //read in snow options -- AJR 2007 @ NMT
	shelterOption = infile.ReadItem(shelterOption, "OPTRADSHELT"); //read in shading option -- AJR 2007 @ NMT
	tempLapseRate = infile.ReadItem(tempLapseRate,"TEMPLAPSE");//K/m -- RINEHART 2007 @ NMT
  
	rainInt = infile.ReadItem(rainInt,  "RAININTRVL");
	
	if (evapotransOption != 0) {
		landPtr = hydro->landPtr;
		soilPtr = hydro->soilPtr;
		initializeVariables();
		SetSunVariables();

		
		if (gFluxOption != 1 && gFluxOption != 2) {
			Cout<<"\nGround Heat Flux Option "<< gFluxOption;
			Cout<<" not valid if evapotranspiration routine active."<<endl;
			Cout<<"\tPlease use :"<<endl;
			Cout << "\t\t(1) for Temperature gradient method" << endl;
			Cout << "\t\t(2) for Force-restore method"<< endl; 
			Cout << "\nExiting Program..."<<endl<<endl;
			exit(1);
		}

	}
	}

/***************************************************************************
**
** tEvapoTrans::CreateHydroMetAndLU() 
**
** Auxiliary function used by the constructor
**  
***************************************************************************/
void tEvapoTrans::CreateHydroMetAndLU(tInputFile &infile)
{
	if (evapotransOption != 0) {
				
		// If stochastic rainfall is used - simulate meterologic variables
		if (rainPtr->getoptStorm()) {
			weatherSimul = new tHydroMetStoch(gridPtr, timer, infile, this, rainPtr);
			newHydroMetStochData(0);
		}
		// Read observed data from files otherwise
		else if (metdataOption == 1) {
			infile.ReadItem(stationFile, "HYDROMETSTATIONS");
			readHydroMetStat(stationFile);
			for (int ct=0;ct<numStations;ct++) {
				readHydroMetData(ct);
			}
			assignStationToNode();
		}
		else if (metdataOption == 2) {
			infile.ReadItem(stationFile, "HYDROMETGRID");
			readHydroMetGrid(stationFile);
			createVariant(); 
			//timeCount = 0;  // bug fixed by Pat - June 2009
		}
		else {
			Cout <<"\nMeteorological Data Option "<< metdataOption;
			Cout <<" not valid."<<endl;
			Cout << "\tPlease use :" << endl;
			Cout << "\t\t(1) for Weather Station Point Data" << endl;
			Cout << "\t\t(2) for Gridded Meteorological Data"<< endl;  
			Cout << "\nExiting Program..."<<endl<<endl;
			exit(1);
		}

		// Read observed data from files otherwise
		if (luOption == 0) {
			Cout << "\nUsing ID Land-Use Map" << endl;
		}
		else if (luOption == 1) {
			cout << "\nUsing dynamic Land-Use Data Grids" << endl;
			infile.ReadItem(luFile, "LUGRID"); // corrected by SY, TM: 11/19/07
			readLUGrid(luFile);
			createVariantLU(); 
			AtFirstTimeStepLUFlag=1;
			if  ((luInterpOption != 0) && (luInterpOption != 1) ) {
				Cout <<"\nLand Use Grid Interpolation Data Option "<< luOption <<" not valid."<<endl;
				Cout << "\tPlease use :" << endl;
				Cout << "\t\t(0) for using current/past values till next incoming grid" << endl;
				Cout << "\t\t(1) for interpolating between curent/past and next incoming grids"<< endl;  
				Cout << "\nExiting Program..."<<endl<<endl;
				exit(1);
			}	
		}
		else {
			Cout <<"\nLand Use Data Option "<< luOption<<" not valid."<<endl;
			Cout << "\tPlease use :" << endl;
			Cout << "\t\t(0) for ID Base Map" << endl;
			Cout << "\t\t(1) for Gridded Land Use Data"<< endl;  
			Cout << "\nExiting Program..."<<endl<<endl;
			exit(1);
		}
	}
}

/***************************************************************************
**
** tEvapoTrans::DeleteEvapoTrans()
**
** Auxiliary function used by the destructor
**  
***************************************************************************/
void tEvapoTrans::DeleteEvapoTrans()
{
	if (evapotransOption != 0) { 
		delete [] currentTime;
		delete [] assignedStation;
		delete [] weatherStations;
		
		if (rainPtr->getoptStorm())
			delete weatherSimul;
		
		if (metdataOption == 2) {
			if (evapotransOption != 4) {
				delete airpressure;
				delete dewtemperature;
				delete skycover;
				delete windspeed;
				delete airtemperature;
				delete surftemperature;
				delete netradiation;
				delete incomingsolar; //E.R.V. 3/6/2012
				delete relhumidity;
				delete vaporpressure;
				
				for (int sz=0;sz<nParm;sz++) {
					delete [] gridParamNames[sz]; //TODO: EXC_BAD_ACCESS (code=1, address=0x0) -WR
					delete [] gridBaseNames[sz];
					delete [] gridExtNames[sz];
				}
				delete [] gridParamNames;
				delete [] gridBaseNames;
				delete [] gridExtNames;
			}
			else
				delete evapotranspiration;
		}

		if (luOption ==1) {
		  deleteLUGrids();
		}	
	}
	}

/***************************************************************************
**
** tEvapoTrans::initializeVariables() Function
**
**
***************************************************************************/
void tEvapoTrans::initializeVariables()
{
	arraySize = gridPtr->getNodeList()->getActiveSize();
	
	assignedStation = new int[arraySize];
	
	VerbID = simCtrl->VerbID;
	
	airTemp = 0.0; dewTemp = 0.0; surfTemp = 0.0;
	atmPress = 0.0; windSpeed = 0.0; rHumidity = 0.0;
	skyCover = 0.0; netRad = 0.0; vPress = 0.0;
	latitude = 0.0; longitude = 0.0; gmt = 0;
	nodeHour = 0; Tso = 0.0; Tlo = 0.0; Gso = 0.0;
	inLongR = 0.0; outLongR = 0.0; inShortR = 0.0;  
	dewTempC = 0.0; surfTempC = 0.0; atmPressC = 0.0;
	rHumidityC = 0.0; skyCoverC = 0.0; netRadC = 0.0;

    skycover_flag = 0;

	Is = Ic = Ics = Id = Ids = 0.0;
	gFlux = 0.0; hFlux = 0.0; lFlux = 0.0;
	potEvap = 0.0; actEvap = 0.0;
	rain = 0.0; betaS = 0.0; betaT = 0.0; windSpeedC = 0.0;
	slope = aspect = 0.0;
	Io=alphaD=sinAlpha=del=phi=tau=circ=0.0; 
	SunRisHrLoc=SunSetHrLoc=deltaT=0.0;
	coeffAl = 0.0; coeffH = 0.0; coeffKt = 0.0;
	coeffRs = 0.0; coeffKs = 0.0; coeffCs = 0.0; coeffV = 0.0;
	panEvap = 0.0; coeffPan = 0.0; // Giuseppe June 2012	

	coeffLAI = 0.0; //RINEHART 2007 @ NMT
	//Horizon angles
	ha0000 = 0.0; ha0225 = 0.0;
	ha0450 = 0.0; ha0675 = 0.0;
	ha0900 = 0.0; ha1125 = 0.0;
	ha1350 = 0.0; ha1575 = 0.0;
	ha1800 = 0.0; ha2025 = 0.0;
	ha2250 = 0.0; ha2475 = 0.0;
	ha2700 = 0.0; ha2925 = 0.0;
	ha3150 = 0.0; ha3375 = 0.0;

	hourlyTimeStep = 0; thisStation = 0; oldTimeStep = 0;
	vapOption = tsOption = nrOption = 0;
	
	currentTime = new int[4];
	for (int count=0;count<4;count++) {
		currentTime[count] = 0;
	}
}

/***************************************************************************
**
** tEvapoTrans::assignStationToNode() Function
**
** Obtains the reference latitude and longitude for each station in the grid
** projection specified for the basin as well as the id code for each 
** tHydroMet station. Calls tResample for obtaining the Thiessen polygons
** assignments of each node to a station.
**
***************************************************************************/
void tEvapoTrans::assignStationToNode()
{
	double *stationLong, *stationLat;
	int *stationID;
	int *id_st_tmp;
	
	stationID   = new int[numStations];
	stationLong = new double[numStations];
	stationLat  = new double[numStations];
	
	for (int ct=0;ct<numStations;ct++) {
		stationID[ct]   = weatherStations[ct].getStation();
		stationLong[ct] = weatherStations[ct].getLong(2);
		stationLat[ct]  = weatherStations[ct].getLat(2);
	}
	
	id_st_tmp = respPtr->doIt(stationID, stationLong, stationLat, numStations); 
	
	for (int ct=0; ct < arraySize; ct++)
		assignedStation[ct] = id_st_tmp[ct];
	
	delete [] stationID; 
	delete [] stationLong; 
	delete [] stationLat;
}

/***************************************************************************
**
** tEvapoTrans::setTime() Function
**
***************************************************************************/
void tEvapoTrans::setTime(int time)
{
	// Get current calendar time
	if (rainPtr->getoptStorm()) {
		currentTime[0] = timer->year;
		currentTime[1] = timer->month;
		currentTime[2] = timer->day;
		currentTime[3] = timer->hour;
		
		// Obtain current hour
		nodeHour = timer->hour;
	}
	else if (metdataOption == 1) {
		currentTime[0] = weatherStations[0].getYear(time); 
		currentTime[1] = weatherStations[0].getMonth(time);
		currentTime[2] = weatherStations[0].getDay(time);
		currentTime[3] = weatherStations[0].getHour(time);
	}
	else if (metdataOption == 2){
		currentTime[0] = timer->year;
		currentTime[1] = timer->month;
		currentTime[2] = timer->day;
		currentTime[3] = timer->hour;
	}
	}


//=========================================================================
//
//
//               Section 2: External Calling Routine
//
//
//=========================================================================
/***************************************************************************
**
** tEvapoTrans::SetEnvironment()
**
***************************************************************************/
void tEvapoTrans::SetEnvironment() {

    // Set time variables ('hourlyTimeStep' is a 'tEvapoTrans' variable)
    setTime(hourlyTimeStep);

    // Compute Sun variables for given hour (basin average,
    // although spatially distributed variables can be easily
    // obtained by re-defining lat/long values for each node)
    SetSunVariables();

    // If stochastic rainfall is used -- use simulated
    // hydrometeorological variables (spatially uniform)
    if (rainPtr->getoptStorm()) {
        tCNode *cNode;
        tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());

        newHydroMetStochData(hourlyTimeStep);

        cNode = nodeIter.FirstP();
        while (nodeIter.IsActive()) {
            // Set simulated values to the node

            cNode->setAirTemp(airTemp);
            cNode->setDewTemp(dewTemp);
            cNode->setRelHumid(rHumidity);
            cNode->setVapPressure(vPress);
            cNode->setSkyCover(skyCover);
            cNode->setWindSpeed(windSpeed);
            cNode->setAirPressure(atmPress);
            if (!hourlyTimeStep) {
                cNode->setSoilTemp(Tlo - 273.15);
                cNode->setSurfTemp(Tso - 273.15);
            }
            cNode = nodeIter.NextP();
        }
    }
}

/***************************************************************************
**
** tEvapoTrans::callEvapoPotential() Function
**
** Called from tSimulator during simulation loop
** 
***************************************************************************/
void tEvapoTrans::callEvapoPotential()
{

	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	int cnt {};
	int count {};
	int tempIndex {};
	double EP {};
	double SkyC {};
	
	if (simCtrl->Verbose_label == 'Y')
		cout << "\nPotential Evaporation Routine Call..."<<endl;
	
	// Set time, sun, and meteorological variables
	SetEnvironment();
	
	// Resample Meteorological Grids, if option
	if (metdataOption == 2){
		resampleGrids(timer);
	}

	if (luOption == 1) { // resampling Land Use grids done here, i.e., dynamic case
	  if (AtFirstTimeStepLUFlag) {			
	    initialLUGridAssignment();
	    AtFirstTimeStepLUFlag=0;
	  }
	  else {
	    LUGridAssignment();
	  }	
	} 

	// Loop through all nodes for this time period
	cNode = nodeIter.FirstP();
	while (nodeIter.IsActive()) {
	
	  if (luOption == 1) { 
	    if ( luInterpOption == 1) { // LU values linearly interpolated between 'previous' and 'until' values

	   //   Cout <<"\nLand Use Grid Interpolation Data Option "<< luOption <<" does not currently work."<<endl;
	   //   Cout << "\tPlease use :" << endl;
	   //   Cout << "\t\t(0) for using current/past values till next incoming grid" << endl;
	   //   Cout << "\nExiting Program..."<<endl<<endl;
	   //   exit(1);

	      interpolateLUGrids(cNode);
	    }
	  }
	  
	  // Elapsed MET steps from the beginning, used for averaging dynamic LU grid values below over time for integ. output
	  auto te = (double)timer->getElapsedMETSteps(timer->getCurrentTime());
	  integratedLUVars(cNode, te);

	  // Use ID for debugging purposes 
	  ID = cNode->getID();
	  
	  // Get rainfall
	  rain = cNode->getRain();
	  
	  // Set Elevation, Slope and Aspect
	  elevation = cNode->getZ();
	  slope = fabs(atan(cNode->getFlowEdg()->getSlope()));
	  aspect = cNode->getAspect();
	  
	  //Develop sky view factors and horizon angles, if necessary
	  if ((shelterOption > 0)&&(shelterOption < 4)) { //CHANGED IN 2008
	    for (tempIndex = 0; tempIndex < 16; tempIndex++) {
	      switch ( tempIndex ) {
	      case 0:
		ha2700 = cNode->getHorAngle2700();
		break;
	      case 1:
		ha2925 = cNode->getHorAngle2925();
		break;
	      case 2:
		ha3150 = cNode->getHorAngle3150();
		break;
	      case 3:
		ha3375 = cNode->getHorAngle3375();
		break;
	      case 4:
		ha0000 = cNode->getHorAngle0000();
		break;
	      case 5:
		ha0225 = cNode->getHorAngle0225();
		break;
	      case 6:
		ha0450 = cNode->getHorAngle0450();
		break;
	      case 7:
		ha0675 = cNode->getHorAngle0675();
		break;
	      case 8:
		ha0900 = cNode->getHorAngle0900();
		break;
	      case 9:
		ha1125 = cNode->getHorAngle1125();
		break;
	      case 10:
		ha1350 = cNode->getHorAngle1350();
		break;
	      case 11:
		ha1575 = cNode->getHorAngle1575();
		break;
	      case 12:
		ha1800 = cNode->getHorAngle1800();
		break;
	      case 13:
		ha2025 = cNode->getHorAngle2025();
		break;
	      case 14:
		ha2250 = cNode->getHorAngle2250();
		break;
	      case 15:
		ha2475 = cNode->getHorAngle2475();
		break;
	      default:
		cout << "\nCheck tempInd -- did not exist or assign" << endl;
	      }//end-switch
	    }//end-for
	    shelterFactorGlobal = cNode->getSheltFact();
	  }
	  else if (shelterOption == 0) {
	    shelterFactorGlobal = 0.5*(1 + cos(slope));
	    cNode->setSheltFact(shelterFactorGlobal);
	  }
	  else {
	    shelterFactorGlobal = 1;
	  }
	  
	  // Set Coefficients - override if dynamic land use
	  if (luOption == 1) {
	    newLUGridData(cNode);
	    if (gFluxOption == 1 || gFluxOption == 2) {
			// Giuseppe 2016 - Begin changes to allow reading soil properties from grids
            //	      coeffKs = soilPtr->getSoilProp(10);
            //	      coeffCs = soilPtr->getSoilProp(11);
            coeffKs = cNode->getVolHeatCond();
            coeffCs = cNode->getSoilHeatCap();
			// Giuseppe 2016 - End changes to allow reading soil properties from grids
	    }
	  }
	  else{
	    setCoeffs(cNode);
	  }


      //updates meteorological variables if not in stochastic mode
      if (!rainPtr->getoptStorm()) {
          if (metdataOption == 1) {
              thisStation = assignedStation[count];
              newHydroMetData(hourlyTimeStep); //read in met data from station file -- inherited function

              if (fabs(skyCover-9999.99)<1.0E-3){ // work around since nodata from grids only set to 9999.99 once in tvariannt
                  skycover_flag =1;
              }

          } else if (metdataOption == 2) {
              //resampleGrids(timerET); // read in met grid data -- inherited function
              newHydroMetGridData(cNode); // set up and get appropriate data -- inherited function

              if (fabs(skyCover-9999.99)<1.0E-3){ // work around since nodata from grids only set to 9999.99 once in tvariannt
                  skycover_flag =1;
              }
          }

          // Set the observed values to the node:
          // they will be required by other function calls
          vPress = vaporPress(); //-- ADDED IN ORDER TO SET RH... CORRECTLY FOR SNOW
          cNode->setAirTemp(airTemp); // celsius
          cNode->setDewTemp(dewTemp);
          cNode->setRelHumid(rHumidity);
          cNode->setVapPressure(vPress);

          // Check/modify cloud cover values
          if (skycover_flag == 1) {
              skyCover = compSkyCover();
          }

          cNode->setSkyCover(skyCover);
          cNode->setWindSpeed(windSpeed);
          cNode->setAirPressure(atmPress);
          cNode->setShortRadIn(inShortR);

          //Set Soil/Surface Temperature
          if (hourlyTimeStep == 0) {
              cNode->setSoilTemp(Tlo - 273.15);
              cNode->setSurfTemp(Tso - 273.15);
          }

      }

	  
	  if (Ioption == 0) {
	    cNode->setNetPrecipitation(rain);
	  }
	  
	  // Call Beta functions
	  betaFunc(cNode); 
	  betaFuncT(cNode);
	  
	  // Get Soil/Surface Temperature
	  Tso = cNode->getSurfTemp() + 273.15;
	  Tlo = cNode->getSoilTemp() + 273.15;
	  
	  // Calculate the Potential and Actual Evaporation
	  if (evapotransOption == 1) {   
	    EvapPenmanMonteith(cNode);
	  }
	  else if (evapotransOption == 2) {
	    EvapDeardorff(cNode);
	  }
	  else if (evapotransOption == 3) {
	    EvapPriestlyTaylor(cNode);
	  }
	  else if (evapotransOption == 4) {
	    EvapPan();
	  }
	  else {
	    Cout << "\nEvapotranspiration Option " << evapotransOption;
	    Cout <<" not valid." << endl;
	    Cout << "\tPlease use :" << endl;
	    Cout << "\t\t(1) for Penman-Monteith Method" << endl;
	    Cout << "\t\t(2) for Deardorff Method"<< endl;
	    Cout << "\t\t(3) for Priestly-Taylor Method" << endl;
	    Cout << "\t\t(4) for Pan Evaporation Measurements" << endl;
	    Cout << "Exiting Program...\n\n"<<endl;
	    exit(1);
	  }
	  // Set computed values to the node variables
	  setToNode(cNode);
	  
	  // Estimate average Ep and cloudiness
	  if (rainPtr->getoptStorm() && Io > 0.0) {
	    potEvap = cNode->getPotEvap();
	    SkyC += skyCover;
	    cnt++;
	  }
	  
	  cNode = nodeIter.NextP();
	  count++;
	}

	timeCount++; // bug fixed by Pat - June 2009

	// Assign old time
	oldTimeStep = hourlyTimeStep;
	
	// Update hourly time 
	hourlyTimeStep++;
	
	// Submit the basin average value
	if (rainPtr->getoptStorm()) {
		// SKY2008Snow -- Following bug corrected to account for SkyC division by cnt again in the next ComputeDailyEpCld call
		//cnt ? SkyC /= ((double)cnt) : SkyC = skyCover;
		cnt ? SkyC = SkyC : SkyC = skyCover;
		
		// Get approximate EP from Priestley-Taylor method
		EP = ApproximateEP();
		
		// Submit values to climate simulator
		weatherSimul->ComputeDailyEpCld(EP, SkyC/cnt);
		
		// Assign the radiation variables to the 'tHydrometStoch'
		if (!count) {
			weatherSimul->setSunH(alphaD);
			weatherSimul->setSinH(sinAlpha);
			weatherSimul->setIo(Io);
			weatherSimul->setIdir(Ics);
			weatherSimul->setIdif(Ids);
			weatherSimul->OutputHydrometVars();
		}
	}
	}

/***************************************************************************
**
**  tEvapoTrans::ApproximateEP()
**
**  Estimate Ep from Priestley-Taylor method assuming Rn(AirTemp).
**  That way, simulation of the dew point temperature will not depend
**  neither on vegetation state nor on soil moisture conditions.
**
***************************************************************************/
double tEvapoTrans::ApproximateEP()
{
	double cc, psy, denom, lam, sigma, emiss_soi;
	double Rn, G, EP;
	double v1;

	if (shelterOption < 3)
		v1 = shelterFactorGlobal;
	else
		v1 = 1;

	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	sigma = 5.67E-8;                      // [W m^-2 K^-4]
	emiss_soi = 0.96;                     // [-]
	
	// The following variables depend on 'airTemp' and 'dewTemp'
	// (these need to be assigned before the function call)
	// Right now, those corresponding to the last node are used
	lam = latentHeat();
	cc  = clausClap();
	psy = 100.0*psychoMetric();
	denom = ((cc/psy) + 1.0);
	
	// Get longwave from the atmosphere using first node
	inLongR  = inLongWave( nodeIter.LastP() );
	
	// Assume air temperature for estimation of L  // [W m^-2]
	outLongR = v1*emiss_soi*sigma*pow((airTemp+273.15),4.0);
	
	// Approximate Rn (term on the right accounts for reflection)  // [W m^-2]
	// Note: could use actual estimated albedo (not 0.17)
	if (sinAlpha > 0)
		Rn = Is*(1.0-0.17) + (inLongR - outLongR);
	else 
		Rn = (inLongR - outLongR);
	
	// Assume ground heat flux 10% of Rn
	// Multiply by 3600 to convert kg/m2/s to mm/hr due to 1kg/m2 = 1mm water
	if (Rn > 0.0) {
		G   = 0.1*Rn;
		EP  = 3600.*(1.26/lam)*(Rn-G)*((cc/psy)/denom);
	}
	else 
		EP = 0.0;
	return EP;
}

/***************************************************************************
**
** tEvapoTrans::setCoeffs() Function
**
** P1 reserved for Land Use ID integer from Anderson (1976) classification
** P2-P7 used for interception parameters
** 
** Properties in Combination Method assigned to:
** 
** P8     CoeffAl  (Surface Albedo)
** P9     CoeffH   (Vegetation Height)
** P10    CoeffKt  (Optical Transmission Coefficient for Vegetation)
** P11    CoeffRs  (Stomatal Resistance)
** P12    CoeffV   (Vegetation Fraction)
**
** Properties in Deardorff Method assigned to:
** 
** P8    CoeffAl  (Surface Albedo)
** P12   CoeffV   (Vegetation Fraction)
**
** Properties in Priestly-Taylor Method assigned to:
** 
** P8    CoeffAl  (Surface Albedo)
** P12   CoeffV   (Vegetation Fraction)
**
***************************************************************************/
void tEvapoTrans::setCoeffs(tCNode* cNode) 
{ 
	landPtr->setLandPtr(cNode->getLandUse());

	if (snowOption == 1)
		coeffLAI = landPtr->getLandProp(12);

	if (evapotransOption == 1) {
		coeffAl = landPtr->getLandProp(7); 
		coeffH  = landPtr->getLandProp(8); 
		coeffKt = landPtr->getLandProp(9);
		coeffRs = landPtr->getLandProp(10);
		coeffV  = landPtr->getLandProp(11);
	}
	else if (evapotransOption == 2 || evapotransOption == 3) {
		coeffAl = landPtr->getLandProp(7); 
		coeffV  = landPtr->getLandProp(11);
	}
	else if (evapotransOption == 4)
		coeffV  = landPtr->getLandProp(11);

    if (coeffV >= 1.0) //prevents loss of snow when unloaded from canopy WR 05/12/2024
        coeffV = 0.99;

	if (gFluxOption == 1 || gFluxOption == 2) {
        coeffKs = cNode->getVolHeatCond();
        coeffCs = cNode->getSoilHeatCap();
	}
	
	}

/***************************************************************************
**
** tEvapoTrans::callEvapoTrans() Function
**
** Called from tSimulator during simulation loop
** 
** Calculate various proportions of evapotranspiration due to the interception
** loss, transpiration and bare soil evaporation. Following Wigmosta(1994).
**
** Evaporation From Wet Canopy:
**
**    Ewc(t) = min(Ep(t),I(t))     Minimum of Interception and PotEvaporation
**  
**    Potential evaporation is backed out from the evaporation calculated
**    in the surface energy balance (includes surface resistances) and called
**    potEvaporation and the transpiration factor.
**
** Evaporation From Dry Canopy (Transpiration):
**
**    Edc(t) = (Ep(t)-Ewc(t))*{CC(t)+Psy(t)}/{CC(t) + Psy(t)*(1-rs(t)/ra(t))}
**
** Evaporation From Soil:
**
**    Es(t) = B*Ep(t)              Surface control on Potential Rate
** 
** Evapotranspiration (Total):
**
**    ET(t) = Ewc(t) + Edc(t) + Ebs(t)    (mm/hr)
**
***************************************************************************/
void tEvapoTrans::callEvapoTrans(tIntercept *Intercept, int flag) 
{

	int count = 0;
	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	
	if (simCtrl->Verbose_label == 'Y') {
		cout<<"EvapoTranspiration Routine Call..."<<endl<<endl;
	}

	// SKYnGM2008LU: Following handles the 'Interception ON' case in SurfaceHydroProcesses
	if (getEToption() == 0 && Intercept->getIoption() == 1) {
	  if (luOption == 1) { // resampling Land Use grids done here, i.e., dynamic case
	    if (AtFirstTimeStepLUFlag) {
	      initialLUGridAssignment();
	      AtFirstTimeStepLUFlag=0;
	    }
	    else {
	      LUGridAssignment();
	    }
	  }
	}

	cNode = nodeIter.FirstP();
	while ( nodeIter.IsActive() ) {

	  if (getEToption() == 0 && Intercept->getIoption() == 1) {
	    if (luOption == 1) { 
	      if ( luInterpOption == 1) { // LU values linearly interpolated between 'previous' and 'until' values
		
		//Cout <<"\nLand Use Grid Interpolation Data Option "<< luOption <<" not valid."<<endl;
		//Cout << "\tPlease use :" << endl;
		//Cout << "\t\t(0) for using current/past values till next incoming grid" << endl;
		//Cout << "\t\t(1) for interpolating between current/past and next incoming grids"<< endl;  
		//Cout << "\nExiting Program..."<<endl<<endl;
		//exit(1);

		interpolateLUGrids(cNode);
	      }
	    }
	  }
		
	  // Elapsed MET steps from the beginning, used for averaging dynamic LU grid values below over time for integ. output
	  auto te = (double)timer->getElapsedMETSteps(timer->getCurrentTime());
	  integratedLUVars(cNode, te);			

	  ID = cNode->getID();
	  elevation = cNode->getZ();
	  ComputeETComponents(Intercept, cNode, count, flag);
	  cNode = nodeIter.NextP();
	  count++;
	}
	timeCount++;
}

/***************************************************************************
**
** tEvapoTrans::ComputeETComponents() Function
**
** Implements the above function at the scale of a Voronoi cell
** 'flag' indicates if the model has to be run with the interception scheme
** 
***************************************************************************/
void tEvapoTrans::ComputeETComponents(tIntercept *Intercept, tCNode *cNode, 
									  int count, int flag)
{
	double potEvaporation, evapoTranspiration;
	double evapWetCanopy, evapDryCanopy, evapSoil, CanStorage;
	double cc, ra, rs, psy, transFactor, actEvaporation, ctos;
	evapWetCanopy = evapDryCanopy = evapSoil = CanStorage = ctos = 0.0;
        potEvaporation = evapoTranspiration = 0.0; 
        cc = ra = rs = psy = transFactor = actEvaporation = 0.0;

	if ( evapotransOption ) {
		
		// Set Coefficients
		if (luOption == 1) {
			newLUGridData(cNode);
			if (gFluxOption == 1 || gFluxOption == 2) {;
                coeffKs = cNode->getVolHeatCond();
                coeffCs = cNode->getSoilHeatCap();
            }
		}
		else{
		  setCoeffs(cNode);
		}
		
		// Call Beta function for transpiration
		betaFuncT(cNode);
	
		// Assign hydromet vars only if real rainfall data are used
		// Assign (as spatially uniform) otherwise
		if (!rainPtr->getoptStorm()) {
			// Get Met Data
			if (metdataOption == 1) {
				thisStation = assignedStation[count];
				newHydroMetData(hourlyTimeStep);//AJR 2008 -- CHANGED FROM OLDTIMESTEP TO HOURLYTIMESTEP
			}
			else if (metdataOption == 2) {
				newHydroMetGridData(cNode);
			}
		}

		// The computed Latent Heat LE -
		// - "pseudo-resistance" evaporation: transFactor*Ep
		potEvaporation = cNode->getPotEvap();
		
		// The soil-moisture controlled evaporation: beta*transFactor*Ep
		actEvaporation = cNode->getActEvap();
		
		if (potEvaporation < 0) {
			potEvaporation = 0.0;
			cNode->setPotEvap(0.0);
		}
		if (actEvaporation < 0.0) {
			actEvaporation = 0.0;
			cNode->setActEvap(0.0);
		}
		
		// Transpiration
		cc = clausClap();
		psy = psychoMetric();
		ra = aeroResist();
		rs = stomResist();
		transFactor = (cc + psy)/(cc + psy*(1+rs/ra));
		
		// Check if the interception scheme is turned on
		if ( flag && (coeffV > 0) && Intercept->IsThereCanopy( cNode )) {
			// Get quantities and make checks depending on Interception model
			if (Ioption == 1) {
				CanStorage = cNode->getCumIntercept();
				ctos = 1.0;
			}
			else if (Ioption == 2) {
				CanStorage = cNode->getCanStorage();
				ctos = Intercept->getCtoS( cNode ); // To scale Ep
				if (ctos > 1)
					ctos = 1.0;
			}
			
			// Evaporation from Wet Canopy
			// 'ctos' is C/S - that gives the term in the Rutter equation
			if (CanStorage >= ctos*potEvaporation*timer->getEtIStep())
				evapWetCanopy = potEvaporation;
			else {
				evapWetCanopy = CanStorage/timer->getEtIStep(); 
				ctos = 1;
			}
			// Call to Interception  Model
			Intercept->callInterception(cNode, potEvaporation);
		}
		else {
			cNode->setNetPrecipitation(cNode->getRain());
		}
		
		// The actual amount extracted from the canopy storage
		evapWetCanopy *= ctos;
		
		// Evaporation from Dry Canopy:
		//  If canopy is wet (C>=S) - transpiration does not occur
		//  If canopy is semi-wet (C<S) - transpiration uses the rest of the energy
		//  Following Elathir and Bras (1993)    
		evapDryCanopy = betaT*((potEvaporation - evapWetCanopy)*transFactor);
		
		// Account for the vegetation fraction
		evapWetCanopy *= coeffV;
		evapDryCanopy *= coeffV;
		
		// Evaporation from Bare Soil

		if(cNode->getBedrockDepth() <= 0) //
        {
            evapSoil = 0;

        }
        else{
		evapSoil = (1-coeffV)*(actEvaporation);
        }

		// Total Evapotranspiration
		evapoTranspiration = evapWetCanopy + evapDryCanopy + evapSoil;
		
		// Assignments
		if (Ioption == 1) {
			CanStorage = cNode->getCumIntercept();
			cNode->setCumIntercept(CanStorage - evapWetCanopy/coeffV*timer->getEtIStep()); 
		}
		
		cNode->setEvapWetCanopy(evapWetCanopy);
		cNode->setEvapDryCanopy(evapDryCanopy);
		cNode->setEvapSoil(evapSoil);
		cNode->setEvapoTrans(evapoTranspiration);
		cNode->addTotEvap(evapoTranspiration); // add to cumulative totals CJC2020
		cNode->addBarEvap(evapSoil); // add to cumulative totals CJC2020
		
		// Update average ET rate from an element
		auto te = (double)timer->getElapsedETISteps(timer->getCurrentTime());
		if (fabs(te - 1.0) < 1.0E-6)
			cNode->setAvET(evapoTranspiration);
		else if (te > 1.0) 
			cNode->setAvET((cNode->getAvET()*(te-1.0) + evapoTranspiration)/te);
	}
	else {
		if ( flag && Intercept->IsThereCanopy( cNode ))
			Intercept->callInterception(cNode, 0.0);
	}
	}

//=========================================================================
//
//
//                  Section 3: Meteorologic Equations
//
//
//=========================================================================


/***************************************************************************
**
** tEvapoTrans::clausClap() Function
**
** Calculates the Clausius-Clayperon relationship in the form:
**                                See Rogers and Yau(1989) p14
**
**        CC(t) = (L(t)*es(t))/(Rv*T(t)^2)  (N/m2*K) = Pa/K
**  
**        where L(t)  latent heat of vaporization   (J/kg)
**              es(t)  vapor pressure at saturation  (mb)
**              Rv   vapor gas constant = 461.5 J/kg/K 
**              T(t)  air temperature  (C)
**              Conversion between mb and Pa = 100
**
***************************************************************************/
double tEvapoTrans::clausClap() 
{
	double latHeat, esat, rv, airTempK, cc;
	latHeat = latentHeat();
	esat = satVaporPress();
	rv = 461.5;
	airTempK = airTemp + 273.15;
	cc = 100.0*(latHeat*esat)/(rv*pow(airTempK,2.0));
	return cc;
}

/***************************************************************************
**
** tEvapoTrans::latentHeat() Function
**
** Calculates the Latent Heat of Vaporization as a function of temperature
** using an expression in Bras(1990) p84:
**
**         L(t) = (597.3 - 0.57*T(t))*4186.8   (J/kg)
**
**         where T(t) = air temperature in degrees C
**               Conversion from cal/g to J/kg = 4186.8
**
***************************************************************************/
double tEvapoTrans::latentHeat() 
{
	return ((597.3-airTemp*0.57)*(4186.8));
}

/***************************************************************************
**
** tEvapoTrans::satVaporPress() Function
**
** Calculates the Saturation Vapor Pressure es(t) using an empirical
** equation by Bolton(1980) shown to be valid over -30C<T<35C range
** See Rogers and Yau(1989) p16.
**
**        es(t) = 6.112*exp((17.67*T(t))/(T+243.5))
** 
**        T(t) in degrees Celsius
**        es(t) in millibar (mb)
**
***************************************************************************/
double tEvapoTrans::satVaporPress() 
{
	return(6.112*exp((17.67*airTemp)/(airTemp+243.5)));
}

/***************************************************************************
**
** tEvapoTrans::satVaporPress()
**
** Calculates the Saturation Vapor Pressure es(t) using an empirical
** equation by Bolton(1980) shown to be valid over -30C<T<35C range
** See Rogers and Yau(1989) p16.
**
**        es(t) = 6.112*exp((17.67*T(t))/(T+243.5))
** 
**        T(t) in degrees Celsius
**        es(t) in millibar (mb)
**
***************************************************************************/
double tEvapoTrans::satVaporPress(double toC)
{
	return (6.112*exp((17.67*toC)/(toC+243.5)));
}

/***************************************************************************
**
** tEvapoTrans::vaporPress() Function
**
** Calculates vapor pressure from relative humidity (rh as percentage)
** or from dew point temperature (Td in C) depending on option set. 
** Returns vapor pressure in millibar. See Entekhabi(1997) p6
** 
** Option 1:    e(t) = es(t)*rh(t)/100;                        in mb
** Option 2:    e(t) = eo*exp((L(t)/Rv)*(1/To - 1/Td(t)))      in mb
**
** Computes dewpoint temperature from relative humidity using
** and vice versa, assigns to corrected arrays.
**
***************************************************************************/
double tEvapoTrans::vaporPress() 
{
	double dewTempK;
	double rv = 461.5;
	double eo = 6.112;
	double to = 273.15;
	
	if (dewHumFlag == 0) {
		vPressC = satVaporPress()*rHumidity/100.0;
		dewTempK = 1.0/((1.0/to) - (log(vPressC/eo)*rv/latentHeat()));
		dewTempC = dewTempK - 273.15;
		rHumidityC = rHumidity;

		// AJR2008, SKY2008Snow
		dewTemp = dewTempC;
		vPress = vPressC;

	}
	else if (dewHumFlag == 1) {
		dewTempK = dewTemp + 273.15;
		dewTempC = dewTemp;
		vPressC = eo*exp((latentHeat()/rv)*((1.0/to)-(1.0/dewTempK)));
		rHumidityC = 100.0*vPressC/satVaporPress();
		
		// AJR2008, SKY2008Snow
		rHumidity = rHumidityC;
		vPress = vPressC;

	}
	else if (dewHumFlag == 2) {
		vPressC = vPress;
		rHumidityC = 100.0*vPressC/satVaporPress();		
		dewTempK = 1.0/((1.0/to) - (log(vPressC/eo)*rv/latentHeat()));
		dewTempC = dewTempK - 273.15;

		// AJR2008, SKY2008Snow
		rHumidity = rHumidityC;
		dewTemp = dewTempC;

	}

	return vPressC;
}

/***************************************************************************
** 
** tEvapoTrans::psychoMetric() Function
**
** Calculates the psychometric constant using Maidment et al(1993) p4.13
**
**        Psy(t) = cp*P(t)/(eps*L(t))    (mb/K)
**
**        cp  Specific Heat of Moist Air = 1013 J/K*kg
**        P(t) = Patm(t) + e(t)
**        eps = Rd/Rv = 0.622
**
***************************************************************************/
double tEvapoTrans::psychoMetric() 
{
	double eps = 0.622;
	double cp = 1013.0;
	return ((cp*totalPress())/(eps*latentHeat()));
}

/***************************************************************************
** 
** tEvapoTrans::totalPress() Function
**
** Calculates the total Pressure P(t) 
**
**         P(t) = Patm(t) + e(t)  (mb)
**
** Assigns corrected Atm Pressure.
**
***************************************************************************/
double tEvapoTrans::totalPress() 
{
	double totPressure, atmospress;
	
	if (fabs(atmPress-9999.99) < 1.0E-3)
		atmospress = 1000.0;
	else
		atmospress = atmPress;
	atmPressC = atmospress;
	totPressure = atmospress + vaporPress();
	return totPressure;
}

/***************************************************************************
** 
** tEvapoTrans::densityMoist() Function
**
** Calculates the moist air density rho(t) as in Bras(1990) p85.
**
**         rho(t) = 100*(P(t)/R*T(t))*(1 - 0.378*(e(t)/P(t))  (kg/m3)
**  
**         R  Dry-Air Gas Constant  287 J/K*kg
**
***************************************************************************/
double tEvapoTrans::densityMoist() 
{
	double rho, Rconst, airTempK;
	Rconst = 287.6;
	airTempK = airTemp + 273.15;
	rho = (100.0*totalPress()/(Rconst*airTempK))*
		(1.0 - 0.378*vaporPress()/totalPress());
	return rho;
}

/***************************************************************************
**
** tEvapoTrans::aeroResist() Function
**
** Estimates the aerodynamic resistance coefficient ra(t) from wind velocity
** data and vegetation height estimates based on land-use reclassification.
** The wind velocity data assumed taken from a 2-meter meteorological tower.
** Vegetation height read from Landuse reclassification table into coeffH.
** 
** Approach is modeled as Shuttleworth in Maidment et al. (1993) p4.12
**
**        ra(t) = ln((zm-d)/zom)*ln((zm-d)/zov)/(u(t)*k^2)  (s/m)
**
**        zm Measurement height above vegetation zm = 2 + coeffH (m)
**        d  Displacement heigth  d = 0.67*coeffH  (m)
**        zom Momentum roughness height  zom = 0.123*coeffH  (m)
**        zov Vapor roughness height  zov = 0.0123*coeffH    (m)
**        k  Von Karman constant 0.41 
**        u(t) Wind speed (m/s)
**
** Consideration for the vegetation fraction is made such that the 
** aerodynamic resistance is a weighted fraction of bare soil and 
** vegetation ra, each calculated with a specific height (vegHeight)
** and vegBare (minimum = 10 cm). The vegetation fraction is an 
** input to the model as coeffV.
**
***************************************************************************/
double tEvapoTrans::aeroResist() 
{
	double vonKarm = 0.41;
	double ra, vegHeight, vegFrac, vegBare;
	double zm, zom, zov, d, rav, ras;
	
	if (coeffH == 0)
		vegHeight = 0.1;
	else 
		vegHeight = coeffH;
	
	vegBare = 1.0;   //ERV 04/2009 used to be 0.1
	vegFrac = coeffV;
	
	if (windSpeed == 0.0 || fabs(windSpeed-9999.99)<1.0E-3)
		windSpeedC = 0.1;    //Minimum wind speed (m/s)
	else
		windSpeedC = windSpeed;
	
	zm = 2.0 + vegHeight;
	zom = 0.123*vegHeight;
	zov = 0.0123*vegHeight;
	d = 0.67*vegHeight;
	rav = log((zm-d)/zom)*log((zm-d)/zov)/(windSpeedC*pow(vonKarm,2.0));
	
	zm = 2.0 + vegBare;
	zom = 0.123*vegBare;
	zov = 0.0123*vegBare;
	d = 0.67*vegBare;
	ras = log((zm-d)/zom)*log((zm-d)/zov)/(windSpeedC*pow(vonKarm,2.0));
	
	ra = (1-vegFrac)*ras + vegFrac*rav;
	
	return ra;
}

/***************************************************************************
**
** tEvapoTrans::stomResist() Function
**
** Estimates the stomatal resistance based on the parameter value for 
** coeffRs give as a reclassification of landuse parameters. A modification
** to account for the diurnal cycle of rs is made by using a time-dependent
** coefficient rsRatio. Depends on having set currentTime.
**
***************************************************************************/
double tEvapoTrans::stomResist() 
{
	double rs;
	int currenthour;
	double rsRatio[24] = {3.837, 3.589, 3.21, 2.43, 1.617, 1.196, 1.067, 1.014, 
		0.995, 0.976, 0.976, 1.0, 1.053, 1.167, 1.354, 1.637,
		2.043, 2.66, 3.215, 3.507, 3.689, 3.818, 3.923, 4.024};
	
	currenthour = currentTime[3];
	rs = coeffRs*rsRatio[currenthour];
	
	// A simple way to constrain transpiration during hours 
	// when there is no incoming solar radiation
	if (alphaD < 0.0)
		rs *= 1000.0;
	
	return rs;
}

//=========================================================================
//
//
//                  Section 4: Radiation Equations
//
//
//=========================================================================


/***************************************************************************
**
** tEvapoTrans::julianDay() Function
**
** Calculates Julian Day (Jday). Requires having set the currentTime. Uses 
** an array daysInMonth to specify the number of days in each month.
**
***************************************************************************/
int tEvapoTrans::julianDay() 
{
	int JDay = 0;
	int dayInMonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	int cmonth, cday, cyear;
	
	cyear  = currentTime[0];
	cmonth = currentTime[1];
	cday   = currentTime[2];
	
	if (cmonth==1)  
		JDay = cday;
	else if (cmonth==2)
		JDay = cday+31;
	else if (cmonth>2) {
		JDay += cday;
		for (int count=0;count<cmonth-1;count++) {
			JDay += dayInMonth[count];
		}
		if ((cyear%4) == 0) {
			JDay += 1;
		}
	}
	return JDay;
}

/***************************************************************************
**
** tEvapoTrans::SetSunVariables() Function
**
** Calculate the Incoming Solar Radiation Intensity at Top of Atmosphere Io
**       Io = (Wo/r^2)sin(alpha)          (J/m2*s)            Bras(1990) 2.9
**     Solar Constant Wo
**         Wo = 1353 W/m2
**     Ratio of Actual Earth-Sun to Mean Earth-Sun Distance (r)          2.10
**         r = 1.0 + 0.017*cos((2*pi/365)*(186-D))  
**     Julian Day D
**     Solar Altitude sin(alpha)         Bras(1990) from Eagleson(1970)  2.4
**         sin(alpha) = sin(del)*sin(phi) + cos(del)*cos(phi)*cos(tau)  
**     Declination of sun  del (radians)                                 2.5
**              del = (23.45*pi/180)*cos((2*pi/365)*(172-D))   
**     Local Latitude  phi (radians)
**              phi = latitude*pi/180  
**
***************************************************************************/
void tEvapoTrans::SetSunVariables()
{
	int JDay, month;
	double r, dgmt, alphaR;
	double Ts, Tsp1, longSM, longM;
	double tau1, tau2;
	double Wo = 1367.0;
	double pi = 4.0*atan(1.0);
	double CIRC[12] = {.07, .10, .12, .14, .16, .19, .23, .20, .16, .12, .08, .07};


	dgmt  = gmt;
	longSM = 15.0*abs(gmt);
	longM = fabs(longitude);
	Ts    =  nodeHour;      // Current Hour
	JDay  = julianDay();    // Current Julian day
	month = currentTime[1]; // Current Month
							//circ  = CIRC[month-1]; // disabled
	
	r = 1.0 + 0.017*cos((2.0*pi/365.0)*(186.0-JDay));
	del = (23.45*pi/180.0)*cos((2.0*pi/365.0)*(172.0-JDay));
	phi = latitude*pi/180.0;
	
	if (gmt) 
		deltaT = (1.0/15.0)*(longSM - longM)*(dgmt/abs(dgmt));
	else 
		// If location is in the "0"th zone, 
		// the following won't work for west longitude
		deltaT = -longM*15.0; 
	
	// Compute hour angle for the current and following hour
	tau1 = ComputeHourAngle( Ts+1.0, deltaT );
	if (Ts < 23.0)
		Tsp1 = Ts+timer->getEtIStep(); 
	else
		Tsp1 = 0.0;
	tau2 = ComputeHourAngle( Tsp1+1.0, deltaT );
	
	// Take an average 'tau' representative for the current hour
	// Around local noon, dicrepancies may arise due to non-
	// synchronization with true solar noon: so, adjust it
	if (fabs(tau1-tau2) >= pi) {
		if (tau2 < tau1)
			tau2 += (2.0*pi);
		else 
			tau1 += (2.0*pi);
	}
	tau = (tau1 + tau2)/2.0;
	if (tau > 2.0*pi)
		tau -= 2.0*pi;
	
	// Take an average of 'sinAlpha'
	//sinAlpha = ((timer->getEtIStep())*sin(del)*sin(phi)+12.0/pi*cos(del)*cos(phi)*
	//(sin(tau2)-sin(tau1))); REMOVED ERV 02/2009 (NOT WORKING)

	// Compute sinAlpha using Bras (1990), eq 2.9
	sinAlpha = (sin(del)*sin(phi)+cos(del)*cos(phi)*cos(tau));

	// Compute an angle
	alphaR = asin(sinAlpha);
	alphaD = alphaR*180.0/pi;

	
	// Compute sun's azimuth using the "hour angle method" [rad from North]
	sunaz = atan(-sin(tau)/(tan(del)*cos(phi)-sin(phi)*cos(tau)));
	if (tau > 0.0 && tau <= pi) {
		if (sunaz > 0.0) 
			sunaz += pi;
		else
			sunaz += (2.0*pi);
	}
	else if (tau >= pi && tau <= 2.0*pi) {
		if (sunaz < 0.0) 
			sunaz += pi;
	}
	
	// Compute extraterrestrial radiation [W m^-2]
	// at the top of the atmosphere
	Io = (Wo/pow(r,2.0));
	
	// These are calculated in LOCAL time. To obtain values in
	// standard meridian time, add 'deltaT' to both
	if (!nodeHour) {
		SunRisHrLoc = (2.*pi-acos(-tan(del)*tan(phi)))*180./pi/15. - 12.;
		SunSetHrLoc = acos(-tan(del)*tan(phi))*180./pi/15. + 12.;
	}
	
	// The total day length
	DayLength = 2.*acos(-tan(del)*tan(phi))*180./pi/15.;
}

/***************************************************************************
**
** tEvapoTrans::ComputeHourAngle() Function
**
**     Hour angle of Sun  tau  (radians)
**              tau = 15*(Ts + 12 - deltaT)  if Ts <= 12                 2.6
**              tau = 15*(Ts - 12 - deltaT)  if Ts  > 12                 2.7
**                   Ts = hour at each location
**                   deltaT = (i/15)*(longSM-longM)                      2.8
**                        longSM = 15*fabs(gmt)         
**                        longM = longitude
**                        i = -1 for West , 1 for East    
**
***************************************************************************/
double tEvapoTrans::ComputeHourAngle(double TT, double delta)
{
	double dtau;
	double pi = 4.0*atan(1.0);
	if (TT < 12.0 + delta)
		dtau = (TT + 12.0 - delta)*15.0;
	else if (fabs(TT-12.0-delta) < 1.0E-3) {
		if (delta > 0.0)
			dtau = (TT + 12.0 - delta)*15.0;
		else
			dtau = (TT - 12.0 - delta)*15.0;
	}
	else
		dtau = (TT - 12.0 - delta)*15.0;
	
	return(dtau*pi/180.0);
}

/***************************************************************************
**
** tEvapoTrans::inShortWave() Function
**
**
** Calculate the Clear Sky Incoming Direct Solar Radiation Intensity Ic 
**                                               for horizontal surface
**       Ic = Io*t^(1/sin(alpha))           Wilson and Gallant (Eq 4.7)
**     Transmission coefficient (t)
**        t = 0.65 + 0.00008*Z      Z = node elevation (meters) (Eq 4.8)
**     
**     Adjustment for Circumsolar Diffuse Radiation
**       Ic = Ic + Id*CIRC                    W&G Eq 4.17
**
** Calculate the Clear Sky Incoming Diffusive Radiation Intensity Id   
**                                            for horizontal surface 
**       Id = (0.271 - 0.294*t^(1/sin(alpha)))*Io
**     Adjustment for Circumsolar Diffuse Radiation
**       Id = Id - Id*CIRC                    W&G Eq 4.18
**       CIRC = {.07, .10, .12, .14, .16, .19, .23, .20, .16, .12, .08, .07}
**           for each month Jan - Dec.
** 
** Calculate the Direct Solar Radiation on Sloping Surface
**
**       Ics = Ic*cosi                        W&G Eq 4.20 - 4.24
**       cosi = A + B*cos(tau) + C*sin(tau)
**       A = sin(del)*sin(phi)*cos(beta) + sin(beta)*cos(asp)*cos(phi)
**       B = cos(del)*(cos(phi)*cos(beta) - sin(phi)*cos(asp)*sin(beta))
**       C = sin(beta)*cos(del)*sin(asp)
**
**       where beta = slope angle (radians), asp = aspect angle (radians)
**
** Calculate the Diffuse Solar Radiation on Sloping Surface
**
**       Ids = Id * v
**       v = sky view fraction ~ cos2(beta/2)         Moore et al (1993)
**
** Calculate Reflection Radiation Component
**
**       Ir = (Ic + Id)*(1 - v)*A      A = albedo     W&G Eq 4.26
**
** Total Incoming Solar Radiation on Sloping Surface (Clear Sky)
**
**       Ith = (Ics + Ids) + Ir
**
** Calculate the Cloudy Sky Incoming Solar Radiation Intensity Is
**       Is = (1-0.65*N^2)*(Ics + Ids) + Ir                              2.29
**     Fractional Sky Cover N 
**          N = skycover/10;
**
** Calculate the Effect of Vegetation Absorption on Incoming Solar Radiation
**       Iv = Kt*Is;
**     Optical Transmission Coefficient of Surface/Vegetation Kt
**     Landuse-derived coefficient coeffKt    Range (0-1)  bareground = 1;
**
** Calculate the Albedo Effect on Incoming Solar Radiation
**       Isw = Iv(1-A)                                                   2.32
**     Albedo Coefficient of Surface/Vegetation
**     Landuse-derived coefficient coeffA     Range (0-1) See Bras(1990) p.37
**
***************************************************************************/
double tEvapoTrans::inShortWave(tCNode *cNode)
{
	double  N, Iv, Isw, Ir; //
	double v, cosi, scover;
	double RadGlobClr;

	Ic=Is=Id=Ir=Ids=Ics=Isw=Iv=0.0;

	// Elevation, slope and aspect have been set before

	// SKY2008Snow from AJR2007
	SunHour = 0.0;

	if (alphaD > 0.0) {

		// Atmospheric turbidity Tlinke: 2.0-3.0 - for rural sites: CALIBRATE

		// Estimate direct beam and diffuse fluxes for horizontal surface
		// Raditaion fluxes 'Ic' and 'Id' will be estimated
		DirectDiffuse(elevation);

        // because inShortR is incoming solar radiation--sky cover is already factored in--WR

        //		// Cloud cover information
        //		if (fabs(skyCover-9999.99) < 1.0E-3) {
        //			skyCover = compSkyCover();//ADDED BY RINEHART 2007 @ NMT
        //			scover = skyCover;
        //		}
        //		else
        //			scover = skyCover;
        //        skyCoverC = scover;


        N = 0; //10.0;

		// If observations (for a horizontal surface) exist -
		// use them, at least in an approximate manner
		if (tsOption > 1 && !rainPtr->getoptStorm()) {
			RadGlobClr = (inShortR / (1.0 - 0.65 * pow(N, 2.0)));
			Ic = Ic/(Ic*sinAlpha + Id)*RadGlobClr;
			Id = RadGlobClr - Ic*sinAlpha;
		}

		// 1) Slope aspect for direct beam radiation
		// Account for the aspect and slope of the element
		// Estimate 'cosi' and compare it with the Sun position
		//  'cosi' = cos(i), where 'i' is the angle between
		//  the sun beam and the normal to the slope surface

		// SKY2008Snow from AJR2007
		//  RINEHART 2007 @ NMT
		//  Only do this computation if sheltering is turned on.
		if (shelterOption < 4) { //CHANGED IN 2008

			cosi = cos(slope)*sinAlpha + sin(slope)*cos(asin(sinAlpha))*cos(sunaz-aspect);

			// SKY2008Snow, AJR2008
			//if (cosi >= 0.0)
			if (cosi >= 0) {
				Ics = Ic*cosi;
				SunHour = 1.0; //YES SUN
			}
			else {
				Ics = 0.0;
				SunHour = 0.0; //NO SUN
			}

		}
		else {
			Ics = Ic;
			SunHour = 1.0;
		}
		if ( (shelterOption == 2) || (shelterOption == 1) ) { //CHANGED IN 2008
			//RINEHART 2007 @ NMT
			//  Finds out whether sun is visible or not.
			Ics *= aboveHorizon(ID);
			SunHour *= aboveHorizon(ID);
		}

		// Components of diffuse
		// Adjust circumsolar factor (disabled)
		//Ic = Ic + Id*circ;
		//Id = Id - Id*circ;

		// 2) Horizon factor for diffuse radiation
		// SKY2008Snow from AJR2007
		//	RINEHART 2007 @ NMT
		//	Uses the different options to decide whether or not to
		//	account for no, local or remote shading.
		if ((shelterOption > 0)&&(shelterOption < 3)) {
			v = shelterFactorGlobal;//include remote
		}
		else if ( (shelterOption == 0) || (shelterOption == 3) ) {//CHANGED IN 2008

		v = 0.5*(1.0 + cos(slope)); //only local

		// SKY2008Snow from AJR2007
		}
		else {
			v = 1.0; //no shading
		}

		Ids = Id*v;

		// 3) Account for cloud cover - the result is
		// the Global Shortwave Irradiance [W m^-2]
		Is = (1.0-0.65*pow(N,2.0))*(Ics + Ids);

		// 4) Reflected radiation from surrounded sites

		// SKY2008Snow from AJR2007
		if (shelterOption == 0) {

			Ir = coeffAl*Is*(1-cos(slope))*0.5;
			Is += Ir;

		// SKY2008Snow from AJR2007
		}
		else if ((shelterOption > 1)&&(shelterOption < 4)) { //CHANGED IN 2008
			Ir = coeffAl*Is*( 0.5*(1 + cos(slope)) - shelterFactorGlobal);
			landRefGlobal = 0.5*(1 + cos(slope)) - shelterFactorGlobal;
			Is += Ir;
		}
		else {
			Ir = 0.0;//AJR2008, SKY2008Snow
		}

		// 5) Account for vegetation
		if (evapotransOption == 1)
			Iv = Is*coeffKt*coeffV + Is*(1.0-coeffV);
		else
			Iv = Is;

		// Account for albedo
		Isw = Iv*(1.0-coeffAl);
  }
	else
		Ic=Is=N=Iv=Isw=Id=Ids=Ics=Ir=0.0;

	// Assign the radiation variables to the 'tHydrometStoch' for ID = 0
	if (rainPtr->getoptStorm() && (ID == 0)) {
		weatherSimul->setSunH(alphaD);
		weatherSimul->setIdir(Ics);
		weatherSimul->setIdir_vis(0.5*Ics);
		weatherSimul->setIdir_nir(0.5*Ics);
		weatherSimul->setIdif(Ids+Ir);//AJR2008, SKY2008Snow
		weatherSimul->setIdif_vis(0.5*(Ids+Ir));//AJR2008, SKY2008Snow
		weatherSimul->setIdif_nir(0.5*(Ids+Ir));//AJR2008, SKY2008Snow
		weatherSimul->OutputHydrometVars();
	}

	// Set shortwave variables to the node (partition is approximate)
	if (tsOption > 1 && !rainPtr->getoptStorm()) {
        cNode->setShortRadIn(inShortR); //or set(Is), they must be equal
    }
	else {
        cNode->setShortRadIn(Isw);
    }
	cNode->setShortRadIn_dir(Ics*(1.0-0.65*pow(N,2.0)));
	cNode->setShortRadIn_dif((Ids+Ir)*(1.0-0.65*pow(N,2.0)));//AJR2008, SKY2008Snow

	return Isw;
}

/***************************************************************************
**
** tEvapoTrans::DirectDiffuse() Function
**
** Estimates clear sky direct beam and diffuse radiative fluxes
**
***************************************************************************/
void tEvapoTrans::DirectDiffuse( double elev)
{
	double h0, m, pp0, Dh0ref, h0ref, drm;
	double TnTLK, Fdh0, A1p, A1, A2, A3;
	double pi = 4.0*atan(1.0);
	
	// 1.) ------------------- Direct beam -------------------
	// The following as according to (Kasten, 1996; Kasten and Young, 1989)
	h0 = alphaD;  //known from SetSunVariables()
	
	// Correction for a given elevation
	pp0 = exp(-elev/8434.5);
	
	// Atmospheric refraction component
	Dh0ref = 0.061359*(0.1594+1.123*h0+0.065656*h0*h0)/(1+28.9344*h0+277.3971*h0*h0);
	h0ref = h0 + Dh0ref;
	
	// Relative optical air mass [-]
	m = pp0/(sin(h0ref*pi/180.)+0.50572*pow((h0ref+6.07995),-1.6364));
	
	// The Rayleigh optical thickness at air mass
	if (m <= 20.0)
		drm = 1/(6.6296+1.7513*m-0.1202*m*m+0.0065*pow(m,3.0)-0.00013*pow(m,4.0)); 
	else 
		drm = 1/(10.4 + 0.718*m);
	
	// Finally, account for radiation attenuation:
	//   obtain the beam irradiance normal to the solar beam [W m^-2]
	//   for cloudless conditions of the atmosphere
	Ic = Io*exp(-0.8662*Tlinke*m*drm);
	
	
	// 2.) ------------------- Diffuse -------------------
	TnTLK = -0.015843 + 0.030543*Tlinke + 0.0003797*pow(Tlinke,2.0);
	A1p = 0.26463 - 0.061581*Tlinke + 0.0031408*pow(Tlinke,2.0);
	if (A1p*TnTLK < 0.0022)
		A1 = 0.0022/TnTLK;
	else
		A1 = A1p;
	A2 = 2.04020 + 0.018945*Tlinke - 0.011161*pow(Tlinke,2.0);
	A3 = -1.3025 + 0.039231*Tlinke + 0.0085079*pow(Tlinke,2.0);
	
	// The solar altitude function
	Fdh0 = A1 + A2*sinAlpha + A3*pow(sinAlpha,2.0);
	
	// Estimation of diffuse radiation on horizontal surface [W m^-2]
	Id = Io*TnTLK*Fdh0;
}

/***************************************************************************
**
** tEvapoTrans::inLongWave() Function based on gray-body (Bras(1990) p44)
**
**        Rlin(t) = K(t)*Ea(t)*S*T(t)^4   (J/(m^2*s)) 
**
**        S      Stefan-Boltzmann Constant = 5.67e-8  J/(s*m^2*K^4)
**        K(t)   Cloud cover coefficient K(t) = (1+0.17*N^2) []
**        N(t)   Sky Cover/10
**        Ea(t)  Atmospheric Thermal Emissivity  Idso(1981)
**               Ea(t) = 0.74 + 0.0049*e(t)  []
**        T(t)   Air Temperature (K)  
**
** Assigns corrected SkyCover 
**
***************************************************************************/
double tEvapoTrans::inLongWave(tCNode *cNode) 
{
	double sigma, kCloud, airTempK, Ea, Rlin, scover, N,v0;
	sigma = 5.67E-8;
	
	if (shelterOption < 3)
		v0 = shelterFactorGlobal;
	else
		v0 = 1;

	// Check/modify cloud cover values
	if (skycover_flag == 1) {

		// SKY2008Snow from AJR2007
		skyCover = compSkyCover();//added by RINEHART 2007 @ NMT
						//  uses relative humidity and rainfall to compute
						//  sky cover.
		scover = skyCover;

	}
	else 
		scover = skyCover;
	skyCoverC = scover;
	
	N = scover/10.0;
	kCloud = 1.0 + 0.17*pow(N,2.0);
	airTempK = airTemp+273.15;
	Ea = 0.74 + 0.0049*vaporPress();
	Rlin = v0*kCloud*Ea*sigma*pow(airTempK,4.0);
	return Rlin;
}

// SKY2008Snow from AJR2007
/***************************************************************************
**
**  tEvapoTrans::compSkyCover()
**
**    RINEHART 2007 @ NMT
**
**    Computes sky cover as a function of relative humidity and rainfall.
**    Taken from Benjamin and Carlson (1986)
**
***************************************************************************/

double tEvapoTrans::compSkyCover() {
  
  double sc;

  if (rain > 0) {
    sc = 10.0; // if raining out, then assume very cloudy
  }
  else {
    sc = round( 10.0*(3.2*rHumidityC/100.0 - 2.4)/0.8 );
  }
  
  //esnure sc to be within limits
  if (sc < 0) 
    sc = 0;
  else if (sc > 10)
    sc = 10.0;

  return sc;
}

/***************************************************************************
**
**  tEvapoTrans::aboveHorizon()
**
**  RINEHART 2007 @ NEW MEXICO TECH
**  
**    This function (a) finds the 22.5 degree azimuth closest to the sun 
**    azimuth, (b) compares the sun altitude to the horizon angle in this
**    direction, (c) set yesOrNo appropriately, and (d) returns whether or
**    not the point can see the sun.
**
***************************************************************************/
double tEvapoTrans::aboveHorizon(int IDinside) {
  double yesOrNo;
  int count(0);
  int dummy(0);

  for(count = 0; count < 16; count++ ) {

    dummy = 0;
    yesOrNo = 1.0;
    
    
    if ( ((-count*22.5 + 180*sunaz/3.1416) < 12.25) || 
		    ((-count*22.5 + 180*sunaz/3.1416) >= -12.25) ) {

	 switch ( count ) {
	    case 0:
		    
	      dummy++;
	      if ( (ha0000)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;
	      
	      
	    case 1:
	      
	      dummy++;
	      if ( (ha0225)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 2:

	      dummy++;
	      if ( (ha0450)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 3:

	      dummy++;
	      if ( (ha0675)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 4:
	      
	      dummy++;
	      if ( (ha0900)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 5:

	      dummy++;
	      if ( (ha1125)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 6:

	      dummy++;
	      if ( (ha1350)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 7:

	      dummy++;
	      if ( (ha1575)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 8:

	      dummy++;
	      if ( (ha1800)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 9:

	      dummy++;
	      if ( (ha2025)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 10:

	      dummy++;
	      if ( (ha2250)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 11:

	      dummy++;
	      if ( (ha2475)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 12:

	      dummy++;
	      if ( (ha2700)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 13:

	      dummy++;
	      if ( (ha2925)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 14:

	      dummy++;
	      if ( (ha3150)*(180/3.1416) < alphaD)
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;

	      
	    case 15:

	      dummy++;
	      if ( (ha3375)*(180/3.1416) < alphaD )
		yesOrNo = 1.0;
	      else
	        yesOrNo = 0.0;
	      break;


	    default:
	      cout << "\ncheck tempInd in aboveHorizon()" << endl;
	 }//switch
    }//if (sun at correct angle)
	
   
    if (dummy)
      break;

  }//for

  if (dummy == 0)
    cout << "Did not find horizon angle azimuth" << endl;
  
  return yesOrNo;
}  

#define TOLF   1.0e-5
/***************************************************************************
**
**  energyBalance() Function
**
**  The function estimates the ground 'Tg' temperature that leads to
**  the radiation balance at the ground in the element. 
**
***************************************************************************/
double tEvapoTrans::energyBalance(tCNode* cNode)
{
	int cnt;
	double f, df;
	double Tg, eps;
	double cosi,v;
	double Isw; // Following were removed as shadows WR: Ic,Ics,Id,Ids,Is

    SunHour=0;
	Ic=Ics=Id=Ids=Is=Isw=0.0;
    elevation = cNode->getZ(); //SMM 10172008
	HeatTransferProperties( cNode );

	
	// Compute INcoming longwave radiation from the atmosphere
	inLongR = inLongWave( cNode );
	cNode->setLongRadIn(inLongR);
	
	// Compute the amount of absorbed shortwave radation
	// and retrieve the computed values (it is assumed that
	// Sun variables have been set for the current hour)
	
	//E.R.V. 3/6/2012 Note: This works for gridded IS input, need to test robustly with other options.
	if(metdataOption == 2){
	  if(fabs(inShortR-0.0) < 1.0E-3){
	    inShortR = inShortWave( cNode );
	    cNode->setShortAbsbSoi(inShortR);
	  }
	  else {
	    //tiantian 05/15/2012: incorporate adjustments for sloping  surface, vegetation and albedo
		if(alphaD > 0.0)
		{
			DirectDiffuse(elevation);//estimate direct and diffuse radiation
			Ic=inShortR;//replace direct radiation with grid NLDAS input

		//  Only do this computation if sheltering is turned on.
			if (shelterOption < 4) { //CHANGED IN 2008
		
				cosi = cos(slope)*sinAlpha + sin(slope)*cos(asin(sinAlpha))*cos(sunaz-aspect);

				// SKY2008Snow, AJR2008
				//if (cosi >= 0.0)
			if (cosi >= 0) {
				Ics = Ic*cosi;
				SunHour = 1.0; //YES SUN
			}
			else {
				Ics = 0.0;
				SunHour = 0.0; //NO SUN
			}	

		}
		else {
			Ics = Ic;
			SunHour = 1.0;
		}
		if ( (shelterOption == 2) || (shelterOption == 1) ) { //CHANGED IN 2008
			//RINEHART 2007 @ NMT
			//  Finds out whether sun is visible or not.
			Ics *= aboveHorizon(ID);
			SunHour *= aboveHorizon(ID);
		}
		v=0.5*(1.0 + cos(slope));
		Ids=Id*v;
		Is=(Ics+Ids)+(Ics+Ids)*coeffAl*0.5*(1.0-cos(slope)); //plus reflection radiation component-> total incoming solar on sloping surface
		if (evapotransOption == 1)
			Isw=(Is*coeffKt*coeffV + Is*(1.0-coeffV))*(1-coeffAl); //acount for vegetation and albedo
		inShortR=Isw;
		cNode->setShortAbsbSoi(inShortR);
		}
	  	}
	}
	else{
	    inShortR = inShortWave( cNode );
	    cNode->setShortAbsbSoi(inShortR);
	}

	// Compute resistances for turbulent fluxes
	Rah  = aeroResist();
	Rstm = stomResist();

	Tg = Tso;  // Surface toK from the previous time step
	
	// Iteratively find the ground temperature
	Tg = rtsafe_mod_energy(cNode, 223.15, 373.15, 
						   TOLF, Tg, inShortR, &eps, &cnt);
	
	// To make sure the fluxes correspond to the obtained Tg
	FunctionAndDerivative(cNode, Tg, f, df, inShortR);
	
	if (simCtrl->Verbose_label == 'Y' && ID == VerbID) {
		cout<<"\n\t******** GROUND ENERGY BALANCE: ********"<<endl;
		cout<<"\tRabsb_soi = "<<inShortR<<";   NetLongRad = "<<(outLongR-inLongR)<<endl;
		cout<<"\t    lEsoi = "<<lFlux
			<<";\t Hsoi = "<<hFlux<<";\t G = "<<gFlux<<endl<<flush;
		cout<<"\n\t---> Tg = "<<Tg-273.15
			<<";\t Tair = "<<airTemp<<endl<<flush;
		cout<<"\tImbalance (fT) = "<<f<<endl;
		cout<<"\t****************************************"<<endl<<endl;
		}

	// Set the temperatures values to the node
	Tso = Tg;
	Tlo += (gFlux*timeStep*60.0)/
		(coeffCs*33.862683*pow((coeffKs/coeffCs/3.6361E-05),0.5));  
	cNode->setSurfTemp(Tg  - 273.15);
	cNode->setSoilTemp(Tlo - 273.15);
	
	// Set element-scale fluxes
	cNode->setNetRad(netRadC);
	cNode->setGFlux(gFlux);
	cNode->setHFlux(hFlux);
	cNode->setLFlux(lFlux);
	cNode->setLongRadOut(outLongR);

	// Set the snow fluxes to zero
	cNode->setSnLHF(0.0);
	cNode->setSnSHF(0.0);
	cNode->setSnSub(0.0); // CJC2020
	cNode->setSnEvap(0.0); // CJC2020
	cNode->setSnPHF(0.0);
	cNode->setSnGHF(0.0);
	cNode->setSnRLin(0.0);
	cNode->setSnRLout(0.0);
	cNode->setSnRSin(0.0);  
	
	return potEvap;
}
#undef TOLF

/*****************************************************************************\
**  
** HeatTransferProperties()
**
** Set element's soil heat properties: conductivity, capacity, diffusivity
**
\*****************************************************************************/
void tEvapoTrans::HeatTransferProperties(tCNode* cNode)
{
    // Giuseppe 2016 - Begin changes to allow reading soil properties from grids
	//	soilPtr->setSoilPtr(cNode->getSoilID());
    //	SoilHeatCondTh = soilPtr->getSoilProp(10);      // Conductivity [W m^-1 K^-1]
    //	SoilHeatCpctTh = soilPtr->getSoilProp(11);      // Capacity     [J m^-3 K^-1]
    SoilHeatCondTh = cNode->getVolHeatCond(); // Conductivity [W m^-1 K^-1]
    SoilHeatCpctTh = cNode->getSoilHeatCap(); // Capacity     [J m^-3 K^-1]
	// Giuseppe 2016 - End changes to allow reading soil properties from grids
	
	SoilHeatDiffTh = SoilHeatCondTh/SoilHeatCpctTh; }

#define MAXITER 100
/*****************************************************************************\
**  
**  rtsafe_mod_energy()
**
**  Finds a root of the polynomial which lies in the range [x1 and x2] 
**  starting from initial guess - xguess. Accuracy of estimation is xacc. 
**  c1, c2, c3 are the polynomial coefficients
**
**  rtsafe_mod(C1, C2, C3, 0.0, 1.0, DX_TOL, VegT->getAllocationLeaf());
** 
\*****************************************************************************/
double tEvapoTrans::rtsafe_mod_energy(tCNode* cNode, double x1, double x2,
									  double xacc, double xguess,
									  double Rabsb_soi, double *eps, int *cnt)
{
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;
	fh=fl=df=dx=0.0;
	
	(*eps) = 0.0;
	(*cnt) = 0;
	
	// -- 'fl' & 'fh' below are the evaluation function values --
	// -- corresponding to arguments 'x1' and 'x2' --
	FunctionAndDerivative(cNode, x1, fl, df, Rabsb_soi);
	if (fl == 0.0) return x1;
	
	FunctionAndDerivative(cNode, x2, fh, df, Rabsb_soi);
	if (fh == 0.0) return x2;
	
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		//cerr<<"tEvapotrans: Root must be bracketed by negative "
		//    <<"and positive f_n values!"<<endl;
		//cerr<<"\tfl = "<<fl<<"; fh = "<<fh<<"; xguess = "<<xguess<<endl<<flush;
		//cerr<<"\tID = "<<ID<<"; wind = "<<windSpeed<<"; inShortR = "<<inShortR
		//	 <<"; airTemp = "<<airTemp<<endl<<endl<<flush;
		//double tt = 1.0;
		//cout<<1./(tt-1.0)<<endl;
	}
	
	(*eps) = 9999.0;
	(*cnt) = 0;
	
	if (fl < 0.0) {  // Orient the search so that f(xl) < 0
		xl = x1;
		xh = x2;
	}
	else {
		xh = x1;
		xl = x2;
	}
	
	rts = xguess;         // A better guess than central value
	dxold = fabs(x2-x1);  // the "stepsize before last",
	dx = dxold;           // and the last step
	
	FunctionAndDerivative(cNode, rts, f, df, Rabsb_soi);
	
	for (j=1; j <= MAXITER; j++) {  // Loop over allowed iterations
		(*eps) = f;
		(*cnt)=j;
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) // Bisect if Newton is out of range
			|| (fabs(2.0*f) > fabs(dxold*df))) {    // or not decreasing fast enough
			dxold = dx;
			dx = 0.5*(xh-xl);
			rts = xl+dx;
			(*eps) = f;
			if (xl == rts) return rts; // Change in root is
		}                              // negligible, take it
		else {
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			(*eps) = f;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts; // Convergence criterion
		
		FunctionAndDerivative(cNode, rts, f, df, Rabsb_soi);
		(*eps) = f;
		
		if (f < 0.0) // <-- Maintain the bracket on the root
			xl=rts;
		else
			xh=rts;
	}
	cerr<<"\n\t\ttEvapotrans: Energy balance: NO convergence in "<<MAXITER<<"\n";
	cerr<<"\t ERROR = "<<f
		<<";  Initial = "<<xguess
		<<";  Last estimate = "<<rts<<";  ID = "<<ID
		<<"\n\t##### Initial value is kept."<<endl<<endl<<flush;
	return xguess;
}
#undef MAXITER

/***************************************************************************
**
**  FunctionAndDerivative()
** 
**  Returns value of the function in non-closed form (fv) as well as 
**  the value of its derivative (dv). Used in Newton iterative procedure  
**  to find root of an equation.
**
***************************************************************************/
void tEvapoTrans::FunctionAndDerivative(tCNode* cNode,
										double Tg, double& fv, double& dv, 
										double Rabsb_soi)
{
	double Rn, Lsoi, Hsoi, lEsoi, G(0), dQ, num;
	double Ep(0), Eps, LE(0), H(0), ccTs;
	double alpha = 1.26;
	double sigma = 5.67E-8;
	double Es = 0.98;   //Calibrated for site (may need to be land parameter)
	double Cp = 1013.0;
	double Ch = 0.0025;
	double rv = 461.5;
	double rho    = densityMoist();
	double P      = totalPress();
	double satVap = satVaporPress();
	double es     = vaporPress();
	double lam    = latentHeat();
	double cc     = clausClap();
	double psy    = 100.0*psychoMetric();
	double denom   = ((cc/psy) + 1.0);
	double denomrs = ((cc/psy) + 1.0 + Rstm/Rah);
	double pi = 4.0*atan(1.0);
	double DTime  = timeStep*60.0;  //Seconds
	
	double dRndTg(0), dLdTg, dHdTg(0), dlEdTg(0), dGdTg(0);
	double soiFct, vegFct, esat, qhSatTs, qhTa;
	double v1;

	if (shelterOption < 3)
		v1 = shelterFactorGlobal;
	else
		v1 = 1;

	// Set temperature to the node
	cNode->setSurfTemp(Tg - 273.15);
	
	// 1.) == NET longwave radiation ==
	inLongR  = cNode->getLongRadIn();
	outLongR = v1*Es*sigma*pow(Tg,4.0);
	Lsoi = outLongR - inLongR;
	
	// 2.) == Ground heat flux ==
	if (gFluxOption == 1) { // Surface toC from the previous time step Tso is used
        G = pow((4.0 * coeffKs * coeffCs / (pi * DTime)), 0.5) * (Tg - Tso);
    }
	else if (gFluxOption == 2) {
        G = ForceRestore(Tg, 1);
    }
	
	// 3.) == NET radiation at the surface ==
	Rn = Rabsb_soi - Lsoi;
	
	// 4.) == Sensible heat flux ==
	if (evapotransOption == 1 || evapotransOption == 3) {
        H = rho * Cp * (Tg - (airTemp + 273.15)) / Rah;
    }
	else if (evapotransOption == 2) {
        H = rho * Cp * (Tg - (airTemp + 273.15)) * Ch * windSpeedC;
    }
	Hsoi = H;
	
	// 5.) == Latent heat flux ==
	if (evapotransOption == 1) {
		dQ = (0.622/P)*(satVap-es)*rho*lam/Rah;
		soiFct = (1.0-coeffV)*betaS;
		vegFct = coeffV*betaT;
		num = (cc/psy)*(Rn-G) + dQ;
		
		Eps = num/denom/lam;
		Ep  = num/denomrs/lam;
		LE = soiFct*lam*Eps + vegFct*lam*Ep;
		Ep *= (denomrs/denom);
	}
	else if (evapotransOption == 2) {
		esat = satVaporPress(Tg-273.15);
		ccTs = (lam*esat)/(rv*pow(Tg,2.0));  
		qhSatTs = (0.622/P)*esat;
		qhTa    = (0.622/P)*es;
		(qhSatTs > qhTa ? Ep = rho*Ch*(qhSatTs-qhTa)*windSpeedC : Ep = 0.0);
		LE = Ep*lam*betaS;
	}
	else if (evapotransOption == 3) {
		(Rn > G ? Ep = (alpha/lam)*(Rn-G)*((cc/psy)/denom) : Ep = 0.0);
		LE = Ep*lam*betaS;
	}
	lEsoi = LE;
	
	
	// ----------------------------------------------
	// Energy Balance as function of soil temperature
	fv = -Rabsb_soi + Lsoi + Hsoi + lEsoi + G;
	
	
	// Compute partial derivatives
 	dLdTg = 4.0*Es*sigma*pow(Tg,3.0);
	dRndTg = -dLdTg;
	
	if (gFluxOption == 1)
		dGdTg = pow((4.0*coeffKs*coeffCs/(pi*DTime)),0.5);
	else if (gFluxOption == 2)
		dGdTg = ForceRestore(Tg, 2); 
	
	if (evapotransOption == 1) {
		dHdTg = rho*Cp/Rah;
		dlEdTg = soiFct*(cc/psy)*(dRndTg-dGdTg)/denom
			+ vegFct*(cc/psy)*(dRndTg-dGdTg)/denomrs;
	}
	else if (evapotransOption == 2) {
		dHdTg = rho*Cp*Ch*windSpeedC;
		dlEdTg = lam*rho*Ch*windSpeedC*(0.622/P)*ccTs*betaS;
	}
	else if (evapotransOption == 3) {
		dHdTg = rho*Cp/Rah;
		dlEdTg = betaS*(alpha)*((cc/psy)/denom)*(dRndTg - dGdTg);
	}
	
	// ----------------------------------------------
	// Derivative with respect to the ground temperature
	dv = dLdTg + dHdTg + dlEdTg + dGdTg;        
	
	// Assign values to the 
	hFlux   = Hsoi;
	lFlux   = lEsoi;
	netRadC = Rn;
	gFlux   = G;
	potEvap = Ep;
	surfTempC = surfTemp = Tg - 273.15;
	
	if (simCtrl->Verbose_label == 'Y' && ID == VerbID && false) {
		cout<<"\n ---------> Tg = "<<Tg-273.15<<endl<<flush;
		cout<<"\t---> Tair = "<<airTemp<<"; dewTemp = "<<dewTemp
			<<"; dewTempC = "<<dewTempC<<";"<<endl;
		cout<<"\trHumidity = "<<rHumidity<<"; vapPress = "<<vPress<<";"<<endl;
		cout<<"\tatmPress = "<<atmPress<<"; windSpeed = "<<windSpeed
			<<"; skyCover = "<<skyCover<<";"<<endl;
		cout<<"== GROUND ENERGY BALANCE: =="<<endl;
		cout<<"\tInShortR (Rabsb_soi) = "<<Rabsb_soi<<endl
			<<"\tInLongR = "<<inLongR
			<<";\t OutLongR = "<<(0.94*5.67E-8*pow(Tg,4.0))
			<<"\t NetLongRsoi = "<<Lsoi<<";"<<endl;
		cout<<"\tlEsoi = "<<lEsoi
			<<";\t Hsoi = "<<Hsoi<<";\t G = "<<G<<endl<<flush;
		cout<<"\tdLdTg = "<<dLdTg<<";\t dHdTg = "<<dHdTg
			<<";\t dlEdTg = "<<dlEdTg<<";\t dGdTg = "<<dGdTg<<endl<<flush;
		cout<<">>>>>>> fT = "<<fv<<";\t dfdT = "<<dv<<endl<<flush;
	}
	}

/***************************************************************************
**
** tEvapoTrans::ForceRestore()
** TODO see if the computations here can be optimized, notably pow is expensive -WR
** Calculates the ground heat flux and the deep soil temperature using
** the Force-restore method. Depending on the option used (1 or 2) this
** function, can return G or dG/dTg. See Lin (1980) and Hu and Islam (1995).
**
** The force-restore expression for G(0,t) is:
** 
** G(0,t) = 0.5*cs*d1( alpha* dTg/dt + w1*(Tg - Tl))     (1)
**
** where w1 = daily frequency [sec-1]
**       k = heat diffusivity = ks/cs  [m2/s]
**       d1 = damping depth [m] = (2*k/w1)^0.5
**       cs = heat capacity [J/m3K]
**       del = depth of upper soil layer = 0.1 m
**       Tg = soil heat in upper del layer [K]
**       Tl = deep soil heat [K]
**
** dTl/dt = G(0,t)/(cs*(365pi)^0.5*d1         (2)
**
** Derivatives approximated at (T(t) - T(t-1))/dt
** Equation (2) solved for Tl(t) and substituted in (1)
** Equation (1) solved for G(0,t)
**
** dG/dTg = (0.5*cs*d1)/(1+(w1*dt)/(2*pow((365*pi),0.5)))*[alpha/dt + w]
**
***************************************************************************/
double tEvapoTrans::ForceRestore(double Ts, int option) 
{
	double w1, d1, k, cs, ddel, pi, dt;    // input
	double alpha, nd, G, Tg, Tl, dg, dTgdt, dGdTg, dTK, tempo = {};
	
	dTK = 0.001;                  // Delta T 
	pi = 4.0*atan(1.0);
	w1 = 2.*pi/86400.;            // Daily frequency [s^-1]
	dt = timeStep*60.0;           // Time step [s]
	
	cs = coeffCs;                 // Heat capacity     [J m^-3 K^-1]
	k  = coeffKs/coeffCs;         // Heat diffusivity  [m^2 s^-1]
	d1 = pow((2.0*k/w1),0.5);     // Damping depth of the diurnal temperature [m] (not Tlo)
				      // d1*sqrt(365)               // Penetration depth of the annual temperature wave  [m]
	ddel = 0.1;                   // Soil layer thickness  [m]

	//ERV 04/2009 soil layer thickness was 0.01 m (1 cm)
	
	nd = ddel/d1;                  // Normalized depth
	
	if (nd >= 0.0 && nd <= 5.0) {
		alpha = 1+0.943*nd+0.223*pow(nd,2.0)+0.0168*pow(nd,3.0)-0.00527*pow(nd,4.0);
	}
	else {
		cout<<"\nWarning: Normalized depth for Force-Restore Equation out of ";
		cout<<"\nvalid range (0 <= ddel/d1 <= 5). Please modify the soil heat ";
		cout<<"\nconductivy (ks) and heat capacity (cs) accordingly.";
		cout<<"\n\nExiting Program..."<<endl<<endl;
		exit(1);
	}
	
	Tg = Ts;     // degree Kelvin
	Tl = Tlo;    // degree Kelvin  
    
	dg = (0.5*cs*d1)/(1.0 + (w1*dt)/(2.0*pow(365.,0.5)));
	
	if (option == 1) {
		dTgdt = (Tg - Tso)/dt;
		G = dg*(alpha*dTgdt + w1*(Tg - Tl)); 
		Gso = G;   //Computed value of G at i-th iteration
		tempo = G;
	}
	// Numerical approximation of the gradient 
	else if (option == 2) {
		dTgdt = (Tg+dTK - Tso)/dt;
		dGdTg = (dg*(alpha*dTgdt + w1*(Tg - Tl)) - Gso)/dTK;
		tempo = dGdTg;
	}
	return tempo;
}

/***************************************************************************
**
** tEvapoTrans::DeriveAspect() Function
**
** Derives aspect of a Voronoi cell
** 
** Called if (evapotransOption != 0)
** 
***************************************************************************/
void tEvapoTrans::DeriveAspect() 
{
	double alpha, d1, x1, y1, x2, y2;
	tCNode *cNode;
	tEdge  *flowedg;
	tArray<double> xy(2), xy1(2);
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	
	cNode = nodeIter.FirstP();
	while (nodeIter.IsActive()) {
		
		flowedg = cNode->getFlowEdg();
		d1  = flowedg->getLength();
		xy  = flowedg->getOriginPtrNC()->get2DCoords();
		xy1 = flowedg->getDestinationPtrNC()->get2DCoords();
		
		x1  = xy1[0]-xy[0]; 
		y1  = xy1[1]-xy[1]; 
		x2  = 0.0; 
		y2  = d1;  // or (d1,0.0) if aspect from the east
		
		// Compute the angle (radians) between the 
		// reference vector (North) and flow edge
		alpha = acos((x1*x2 + y1*y2)/(d1*d1)); 
		
		// Adjustment of angle
		if (xy1[0] < xy[0])
			alpha = 8.0*atan(1.0)-alpha;  
		
		cNode->setAspect(alpha);
		cNode = nodeIter.NextP();
	}
	}

//=========================================================================
//
//
//                  Section 5: Evaporation Equations
//
//
//=========================================================================


/***************************************************************************
**
** tEvapoTrans::EvapPenmanMonteith() Function based on Penman-Monteith Method
**             See Shuttleworth(1979) or Wigmosta et al(1994)
**
**   Actual Evaporation 
**
**        Ep(t) = CC(t)*(Rn(t)-G(t))/Psy(t)+rho(t)*L(t)*dQ(t)/ra(t)  (kg/m2/s))
**                --------------------------------
**                       L(t)*(1+(CC(t)/Psy(t))+rs/ra)
**
**        CC(t) Clausius-Clayperon Relationship in clausClap()  N/m2/K
**        dQ(t) Vapor Pressure Deficit  (qs(t) - q(t)) mb
**        Psy(t) Psychometric constant in psychoMetric() mb/K
**        ra(t) Aerodynamic resistance in aeroResist() s/m
**        cp  Specific Heat of water  J/K/kg
**        rho(t) Moist Air density in densityMoist()  kg/m3
**        L(t) Latent Heat of Vaporization in latentHeat()  J/kg
**        Rn(t) Net radiation in netRadiation()  J/m2/s
**        rs(t) Stomatal Resistance  s/m
**
** Multiply by 3600 to convert kg/m2/s to mm/hr due to 1kg/m2 = 1mm water
**
***************************************************************************/
void tEvapoTrans::EvapPenmanMonteith(tCNode* cNode) 
{
	potEvap = 3600.0*energyBalance(cNode);   // Actual rate, including resistances
	actEvap = 3600.0*(lFlux/(latentHeat()));  
}

/***************************************************************************
**
** tEvapoTrans::EvapDeardorff() Function  
**
** Potential Evaporation computed from Deardorff Equation
**       
**        Ep(t) = rho(t)*Ch*u(t)*(qhSatTs(t)-qhTa(t))   (kg/m2/s)
**     
**        rho(t) Moist Air Density (kg/m3) using densityMoist()
**        Ch     Heat and Moisture Coefficient Ch = 0.0025
**        u(t)   Wind Speed (m/s)
**        qhSatTs(t)  Specific humidity at surface Temperature Ts  []
**            qhSatTs(t) = (0.622*esat(Ts))/(atmPress)          3.28
**               esat(Ts) as in vaporPress() (mb)
**               atmPress from totalPress() (mb)  
**        qhTa   Specific Humidity at Air Temperature []
**            qhTa = 0.622*vaporPress/atmPress                     
**        Ts     Surface Temperature in K
**
** EvapDeardorff() requires:
**       airTemp, windSpeed, dewTemp/rHumidity, atmPress, skyCover
**       coeffA  for inShortWave()
**
** EvapDeardorff() does not require:
**       netRad, surfTemp (optional)
**       coeffH  for aeroResist()
**       coeffRs for stomResist()
**       coeffKt for inShortWave()  bare-soil evaporation (=1)
**
** Multiply by 3600 to convert kg/m2/s to mm/hr due to 1kg/m2 = 1mm water
**
***************************************************************************/
void tEvapoTrans::EvapDeardorff(tCNode* cNode) 
{
	potEvap = 3600.0*energyBalance(cNode);
	actEvap = potEvap*betaS;
}

/***************************************************************************
**
** tEvapoTrans::EvapPriestlyTaylor() Function
**
** Calculates evaporation using the Priestly-Taylor formula which is useful
** under conditions of minimum advection (energy-dominated).
**
**       Ep(t) = (alpha/L(t))*(CC(t)/(CC(t)+Psy(t)))*(Rn(t)-G(t))
**
**       alpha  Priestly Taylor coefficient = 1.26 Semi-humid, Humid regions
**
** EvapPriestlyTaylor() requires (same as Deardorff):
**       airTemp, windSpeed, dewTemp/rHumidity, atmPress, skyCover
**       coeffA  for inShortWave()
**       netRad, surfTemp (calculated as in PM)
**
** EvapPriestlyTaylor() does not require:
**       coeffH  for aeroResist()
**       coeffRs for stomResist()
**       coeffKt for inShortWave()  bare-soil evaporation (=1)
**
** Multiply by 3600 to convert kg/m2/s to mm/hr due to 1kg/m2 = 1mm water
**
***************************************************************************/
void tEvapoTrans::EvapPriestlyTaylor(tCNode* cNode) 
{
	potEvap = 3600.0*energyBalance(cNode);
	actEvap = potEvap*betaS;
}

/***************************************************************************
**
** EvapPan() Function
**
** Calculates actual pan evaporation given the measurements from point
** locations (treated as Meteorological Stations). Pan coefficient used
** to correct for artifacts of pan geometry. This coefficient is inputed
** through the station file otherVariable column. User must use mm/hr
** for pan evaporation in MDF format.
**
***************************************************************************/
void tEvapoTrans::EvapPan() 
{
	if (metdataOption == 1) {
		potEvap = coeffPan*panEvap;
	}
	else if (metdataOption == 2) {
		potEvap = panEvap;
	}
	
	// Very approximate: we would need stomatal resistance 
	// to obtain the transpiration component right
	actEvap = potEvap*(betaS*(1.0-coeffV) + coeffV*betaT);
}

/***************************************************************************
**
** betaFunc() Function
**
** This function calculates the beta parameter in the relationship
** between the potential evaporation Ep (potEvap) and the actual 
** evaporation Ea (actEvap). The B-method options are well described in
** Mahfouf and Noilhan (1991). Here, we choose the simple parameterization
** by Deardorff (1978):
**
**   beta(theta) = min(1,theta/thetaS*0.75)
**
**   LEa = beta(theta)*LEp
**
** where theta = soil moisture in upper 10 cm at time t
**       thetaS = saturation soil moisture parameter
**       LEa = actual latent heat flux
**       LEp = potential latent heat flux
**       beta(theta) = beta as function of soil moisture
**       thetaS*0.75 = field capacity soil moisture
**
***************************************************************************/
void tEvapoTrans::betaFunc(tCNode* cNode) 
{
	double beta, Ths, Thr, Th, ratio, Th_star, soilzone_cutoff;

    Ths = cNode->getThetaS(); // Saturation moisture content
    Thr = cNode->getThetaR(); // Residual moisture content
	// Giuseppe 2016 - End changes to allow reading soil properties from grids

    soilzone_cutoff = cNode->getSoilCutoff();

	Th_star = landPtr->getLandProp(13);

	// Check that Thw <= Th_star <= Ths
	if ((Th_star > Ths) || (Th_star < Thr))
		{ 
		printf("Th_star %f out of the range residual Th_r-Th_s in land cover class %d\n", Th_star, cNode->getLandUse());
		printf("Modify the value of the land cover class\n");
		exit(1);
		}
	// End of modifications by Luis Mendez and Giuseppe Mascaro (April 2013)

	Th = cNode->getSoilMoisture(); // avg volumetric soil moisture content down to 100 mm
	ratio = (Th - Thr)/(Th_star - Thr);

    if (ratio > 1.0)
		beta = 1.0;
	else if (ratio < 0.0)
		beta = 0.0;
	else
		beta = ratio;

    if(Th<=soilzone_cutoff)
        beta = 0.0;

	betaS = beta;
}

/***************************************************************************
**
** betaFuncT() Function
**
** This function calculates the beta parameter in the relationship
** between the potential transpiration Etp and the actual transpiration
** Eta. The B-method options are well described in Mahfouf et al (1996)
** and Feddes et al (2001). Here, we choose a simple parameterization
** based on the CSIRO9 model that resembles the betaFunc() for bare soil.
**
**   betaT(theta) = min(1,(theta - thetaW)/(thetaS*0.75 - thetaW)
**
**   Eta = beta(theta)*Etp
**
** where theta = soil moisture in top 1 meter of soil
**       thetaS = saturation soil moisture parameter
**       thetaW = wilting point moisture = residual
**       Eta = actual transpiration
**       Etp = potential transpiration
**       betaT(theta) = beta as function of soil moisture
**       thetaS*0.75 = field capacity soil moisture
**
***************************************************************************/
void tEvapoTrans::betaFuncT(tCNode* cNode) 
{

	double beta, Ths, Thw, Th, ratio, Th_star, rootzone_cutoff;

    Ths = cNode->getThetaS(); // Saturation moisture content
    Thw = cNode->getThetaR(); // Residual moisture content
	// Giuseppe 2016 - End changes to allow reading soil properties from grids

    rootzone_cutoff = cNode->getRootCutoff();

	Th_star = landPtr->getLandProp(14); 
	// Check that Thw <= Th_star <= Ths
	if ((Th_star > Ths) || (Th_star < Thw))
		{ 
		printf("Th_star out of the range residual Th_r-Th_s in land cover class %d\n", cNode->getLandUse());
		printf("Modify the value of the land cover class");
		exit(1);
		}
	// End of modifications by Luis Mendez and Giuseppe Mascaro (April 2013)


	Th = cNode->getRootMoisture();// avg volumetric soil moisture content down to 1000 mm
	ratio = (Th - Thw)/(Th_star - Thw);
    
	if (ratio < 1.0)
		beta = ratio;
    else
		beta = 1.0;

    if(Th<=rootzone_cutoff)
        beta = 0.0;

	betaT = beta;
}

//=========================================================================
//
//
//                  Section 6: Set and Get Functions
//
//
//=========================================================================


/***************************************************************************
**
** tEvapoTrans::getEToption() Functions
**
***************************************************************************/
int tEvapoTrans::getEToption() 
{
	return evapotransOption;
}

/***************************************************************************
**
** tEvapoTrans::setToNode() Function
**
** Functions used to set the values of tCNode for each time step and for
** each voronoi cell. 
**
***************************************************************************/
void tEvapoTrans::setToNode(tCNode* cNode) 
{ 
	if (dewHumFlag == 0) {
		cNode->setDewTemp(dewTempC);
		cNode->setVapPressure(vPressC);
	}
	else if (dewHumFlag == 1) {
		cNode->setVapPressure(vPressC);
		cNode->setRelHumid(rHumidityC);
	}
	else if (dewHumFlag == 2) {
		cNode->setDewTemp(dewTempC);
		cNode->setRelHumid(rHumidityC);
	}
	
	cNode->setPotEvap(potEvap);
	cNode->setActEvap(actEvap);
	cNode->setAirTemp(airTemp);
	cNode->setAirPressure(atmPressC);
	cNode->setSkyCover(skyCoverC);
	cNode->setWindSpeed(windSpeedC);
	cNode->addCumHrsSun(SunHour);

	// Elapsed MET steps from the beginning
	auto te = (double)timer->getElapsedMETSteps(timer->getCurrentTime());
	// Estimated evaporative fraction
	double tmp = 0;
	if ((fabs(hFlux)+fabs(lFlux)) > 1.0)
		tmp = fabs(lFlux)/(fabs(hFlux)+fabs(lFlux));
	
	if (fabs(te - 1.0) <= 1.0E-6)
		cNode->setAvEvapFract(tmp);
	else if (te > 1.0)
		cNode->setAvEvapFract((cNode->getAvEvapFract()*(te-1.0) + tmp)/te);

	cNode->setSheltFact(shelterFactorGlobal);
	cNode->setLandFact(landRefGlobal);	
	cNode->addRSin(inShortR*3600.0); }

//=========================================================================
//
//
//                       Section 7: Input Routines
//
//
//=========================================================================


/***************************************************************************
**
** tEvapoTrans::readHydroMetStat() Function
**
**
** Reads the HydroMetStation File which provides information concerning
** the weather stations used for hydrometeorological data. Creates an
** array of tHydroMet objects for storing data. (see tHydroMet.h)
**
** Format for HydroMetStation File:
**
** Header:
** nStations nParams (10)
**
** Body:
**
** StationID# FilePath Latitude Longitude GMT RecordLength NumParameters Other
**
** The file should have N rows corresponding to the N number of stations.
** and M columns corresponding to the data for each station.
**
** StationID  (int)         1->N
** FilePath   (string)      File name and path for the station datafile
** Absolute Latitude   (double)  Station latitude  (decimal degree)
** Reference Latitude  (double)  Station latitude  (basin projection)
** Absolute Longitude  (double)  Station longitude (decimal degree)
** Reference Longitude (double)  Station longitude (basin projection)
** GMT        (int)         Offset from Greenwich Mean Time (hrs)
** RecordLength  (int)      Number records for each station
** NumParameters (int)      Number of parameters for each station (11)
** Other (double)           Other variable (ex. elevation)
**                          For pan evap, other used for pan coefficient
**                          Set to Zero if not used.
**
***************************************************************************/
void tEvapoTrans::readHydroMetStat(char *stationfile) 
{
	int nStations, nParams;
	int Gmt, stationID, numTimes, numParams;
	double alat, along, rlat, rlong, otherVar;
	char fileName[kName];
	//assert(fileName != 0);//WR--09192023: comparison of array 'fileName' not equal to a null pointer is always true
	
	Cout<<"\nReading HydroMeteorological Station File '";
	Cout<< stationfile<<"'..."<<endl<<flush;
	
	ifstream readFile(stationfile); 
	if (!readFile) {
		cout << "File "<<stationfile<<" not found." << endl;
		cout<<"Exiting Program...\n\n"<<endl;
		exit(1);
	}
	
	readFile >> nStations;
	readFile >> nParams;
	
	weatherStations = new tHydroMet[nStations];
	assert(weatherStations != nullptr);
	
	numStations = nStations;
	
	for (int count=0;count < nStations;count++) {
		
		for (int ct=0;ct < nParams;ct++) {
			if (ct==0) {
				readFile >> stationID;
				weatherStations[count].setStation(stationID);
			}
			if (ct==1) {
				readFile >> fileName;
				weatherStations[count].setFileName(fileName);
			}
			if (ct==2) {
				readFile >> alat;
				weatherStations[count].setLat(alat,1);
			}
			if (ct==3) {
				readFile >> rlat;
				if (!readFile) {
					cout<<"\nError in the SDF File "<<fileName<< " !!"<<endl;
					cout<<"Please replace BasinLat with Reference Latitude ..."<<endl;
					cout<<"Exiting Program...\n\n"<<endl;
					exit(2);
				}
				weatherStations[count].setLat(rlat,2);
			}
			if (ct==4) {
				readFile >> along;
				weatherStations[count].setLong(along,1);
			}
			if (ct==5) {
				readFile >> rlong;
				if (!readFile) {
					cout<<"\nError in the SDF File "<<fileName<< " !!"<<endl;
					cout<<"Please replace BasinLong with Reference Longitude..."<<endl;
					cout<<"Exiting Program...\n\n"<<endl;
					exit(3);
				} 
				weatherStations[count].setLong(rlong,2);
			}
			if (ct==6) {
				readFile >> Gmt;
				if (!readFile) {
					cout<<"\nError in the SDF File "<<fileName<< " !!"<<endl;
					cout<<"Please replace GMT Tag with GMT Value..."<<endl;
					cout<<"Exiting Program...\n\n"<<endl;
					exit(4);
				}
				weatherStations[count].setGmt(Gmt);
			}
			if (ct==7) {
				readFile >> numTimes;
				weatherStations[count].setTime(numTimes);
			}
			if (ct==8) {
				readFile >> numParams;
				weatherStations[count].setParm(numParams);
			}
			if (ct==9) {
				readFile >> otherVar;
				weatherStations[count].setOther(otherVar);
			}
		}
	}
	readFile.close();
}

/***************************************************************************
**
** tEvapoTrans::readHydroMetData() Function
**
**
** Reads and assigns data values to tHydroMet objects. 
**
** Format is based on HMET_WES, as described by Ogden (1998). Minor 
** modifications issued in terms of units, numerical format and parameters. 
** Hydrometeorological data in other formats must be converted prior to use.
**
** File Format:
** 
** Description Line: 
**      Abbreviations of Parameters as character strings (abbreviations 
**      Ex. Y M D H PA TD XC US TA TS NR     
**
** Body Lines:
**      Values for each parameters. Read in as ints and doubles
**      Ex. Year (4 digit number), Month, Day, Hour (int)
**          AtmPressure (mb)  double     
**          RelativeHumidity (%) double   | either
**          DewTemperature (C)  double    |   or
**          SkyCover (tenths) double
**          WindSpeed (m/s) double
**          AirTemperature (C) double
**          SurfTemperature (C) double
**          NetRadiation (W/m2) double
**      No Data Flag as 9999.99 for doubles
**      These variables are time-series. Not to be confused with protected
**      class member functions which are based on Node geometry.
**
**      Can be used for Pan Evaporation Measurements through option 4
**      Only stores Pan Evaporation Data, no Meteorological Data
**          PanEvaporation (mm/hour) double
**
***************************************************************************/
void tEvapoTrans::readHydroMetData(int num) 
{
	int numParams, numTimes;
	char fileName[kName];
	char paramNames[10] = {};
	char notUsed[10];
	char paramNames2[10] = {};
	char paramNames3[10] = {};
	char *tmpstr;
	
	int *year, *month, *day, *hour;
	double *AtmPressure, *DewTemperature, *AirTemperature;
	double *SkyCover, *WindSpeed, *RelativeHumidity, *VaporPressure;
	double *NetRadiation, *SurfTemperature, *PanEvap;
	double *GlobRadiation;
	double tempo;
	
	tmpstr = weatherStations[num].getFileName();
    snprintf(fileName,sizeof(fileName),"%s", tmpstr);//WR--09192023: 'sprintf' is deprecated: This function is provided for compatibility reasons only.
	numParams = weatherStations[num].getParm();
	numTimes  = weatherStations[num].getTime();
	
	Cout<<"\nReading HydroMeteorological Data File '";
	Cout<<fileName<<"'..."<<endl<<flush;
	
	ifstream readDataFile(fileName);
	if (!readDataFile) {
		cout << "\nFile " <<fileName<<" not found!" << endl;
		cout << "Exiting Program...\n\n"<<endl;
		exit(2);}
	
	for (int cnt = 0; cnt<numParams; cnt++) {
		if (cnt==5)
			readDataFile >> paramNames; 
		else if (cnt==9) 
			readDataFile >> paramNames2;
		else if (cnt==10)
			readDataFile >> paramNames3;
		else
			readDataFile >> notUsed;
	}
	
	// "6th" data element in the input data (air humidity in some form)
	if (strcmp(paramNames,"RH")==0)
		vapOption = 1;
	else if (strcmp(paramNames,"TD")==0)
		vapOption = 2;
	else if (strcmp(paramNames,"VP")==0)
		vapOption = 3;
	
	if (strcmp(paramNames2,"TS")==0)
		tsOption = 1;
	else
		tsOption = 2;

	if (strcmp(paramNames3,"NR")==0)
		nrOption = 1;
	else
		nrOption = 2;
	
	year  = new int[numTimes];
	month = new int[numTimes];
	day   = new int[numTimes];
	hour  = new int[numTimes];
	
	if (evapotransOption != 4) {
		AtmPressure = new double[numTimes];
		DewTemperature= new double[numTimes];
		AirTemperature = new double[numTimes];
		SkyCover = new double[numTimes];
		WindSpeed = new double[numTimes];
		RelativeHumidity = new double[numTimes];
		NetRadiation= new double[numTimes];
		SurfTemperature = new double[numTimes];
		VaporPressure = new double[numTimes];
		GlobRadiation = new double[numTimes];
		
		for (int count = 0; count < numTimes; count++) {
			for (int ct = 0; ct < numParams ;ct++) {
				if (ct==0) {
					readDataFile >> year[count];}
				else if (ct==1) {
					readDataFile >> month[count];}
				else if (ct==2) {
					readDataFile >> day[count];}
				else if (ct==3) {
					readDataFile >> hour[count];}
				
				else if (ct==4) {
					readDataFile >> tempo;
					
					// SKY2008Snow, AJR2008	
					// if (tempo < 700 || tempo > 1200)
					if (tempo < 600 || tempo > 1200)
					
						AtmPressure[count] = 9999.99;
					else 
						AtmPressure[count] = tempo;
				}
				else if (ct==5) {
					if (vapOption == 1) {
						readDataFile >> tempo;
						if (tempo < 0 || tempo > 100)
							RelativeHumidity[count]= 9999.99;
						else 
							RelativeHumidity[count] = tempo;
						DewTemperature[count] = 9999.99;
						VaporPressure[count] = 9999.99;
					}
					else if (vapOption == 2) {
						readDataFile >> tempo;
						if (tempo < -50 || tempo > 60)
							DewTemperature[count] = 9999.99;
						else 
							DewTemperature[count] = tempo;
						RelativeHumidity[count] = 9999.99;
						VaporPressure[count] = 9999.99;
					}
					else if (vapOption == 3) {
						readDataFile >> tempo;
						if (tempo < 0  || tempo > 200)
							VaporPressure[count] = 9999.99;
						else 
							VaporPressure[count] = tempo;
						RelativeHumidity[count] = 9999.99;
						DewTemperature[count] = 9999.99;
					}
				}
				else if (ct==6) {
					readDataFile >> tempo;
					if (tempo < 0 || tempo > 10)
						SkyCover[count]= 9999.99;
					else 
						SkyCover[count] = tempo;
				}
				else if (ct==7) {
					readDataFile >> tempo;
					if (tempo < 0 || tempo > 30)
						WindSpeed[count] = 9999.99;
					else 
						WindSpeed[count] = tempo;
				}
				else if (ct==8) {
					readDataFile >> tempo;
					if (tempo < -50 || tempo > 60)
						AirTemperature[count] = 9999.99;
					else 
						AirTemperature[count] = tempo;
				}
				else if (ct==9) {
					readDataFile >> tempo;

					// SKY2008Snow, AJR2008
					if (fabs(tempo-9999.99) > 1.0E-3) { 

						if (tsOption == 1) {
							if (tempo < -60 || tempo > 70)
								SurfTemperature[count]= 9999.99;
							else
								SurfTemperature[count] = tempo;
							GlobRadiation[count] = 9999.99;
						}
						else {
							GlobRadiation[count] = tempo;
							SurfTemperature[count]= 9999.99;
						}
					
					// SKY2008Snow, AJR2008
					}
					else {
						tsOption = 0;
						GlobRadiation[count] = 9999.99;
						SurfTemperature[count]= 9999.99;
					}


				}
				else if (ct==10) {
					readDataFile >> tempo;
					if (nrOption == 1) {
						if (tempo < -1000 || tempo > 1000)
							NetRadiation[count]= 9999.99;
						else
							NetRadiation[count] = tempo;
					}
					else
						NetRadiation[count] = 9999.99;
				}
				// If 'numParams' exceeds 11 - for compatability 
				// with other implementations of tRIBS
				else 
					readDataFile >> tempo;
				
			} // loop through 'numParams' 
		}   // loop through 'numTimes'
	}
	else {
		PanEvap = new double[numTimes]; 
		assert(PanEvap != nullptr);
		
		for (int sount = 0;sount<numTimes;sount++) {
			for (int sct = 0;sct<numParams;sct++) {
				if (sct==0) {
					readDataFile >> year[sount];}
				else if (sct==1) {
					readDataFile >> month[sount];}
				else if (sct==2) {
					readDataFile >> day[sount];}
				else if (sct==3) {
					readDataFile >> hour[sount];}
				else if (sct==4) {
					readDataFile >> PanEvap[sount];
				}
			}
		}
	}
	
	readDataFile.close();
	
	weatherStations[num].setYear(year);
	weatherStations[num].setMonth(month);
	weatherStations[num].setDay(day);
	weatherStations[num].setHour(hour);
	
	if (evapotransOption != 4) {
		robustNess(AirTemperature, numTimes);
		robustNess(DewTemperature, numTimes);
		robustNess(AtmPressure, numTimes);
		robustNess(SkyCover, numTimes);
		robustNess(RelativeHumidity, numTimes);
		robustNess(WindSpeed, numTimes);
		robustNess(SurfTemperature, numTimes);
		robustNess(VaporPressure, numTimes);
		robustNess(GlobRadiation, numTimes);
		
		weatherStations[num].setAirTemp(AirTemperature);
		weatherStations[num].setDewTemp(DewTemperature);
		weatherStations[num].setAtmPress(AtmPressure);
		weatherStations[num].setSkyCover(SkyCover);
		weatherStations[num].setRHumidity(RelativeHumidity);
		weatherStations[num].setWindSpeed(WindSpeed);
		weatherStations[num].setSurfTemp(SurfTemperature);
		weatherStations[num].setVaporPress(VaporPressure);
		weatherStations[num].setRadGlobal(GlobRadiation);
		
		//cout << "\treadHydroMetData GlobRadiation: " << GlobRadiation<<endl;
		
		if (numParams >= 11) {
			robustNess(NetRadiation, numTimes);
			weatherStations[num].setNetRad(NetRadiation);
		}
		
		delete [] AtmPressure; 
		delete [] DewTemperature; 
		delete [] AirTemperature;
		delete [] RelativeHumidity; 
		delete [] SkyCover;
		delete [] WindSpeed;
		delete [] NetRadiation; 
		delete [] SurfTemperature;
		delete [] VaporPressure;
		delete [] GlobRadiation;
	}
	else {
		robustNess(PanEvap, numTimes);
		weatherStations[num].setPanEvap(PanEvap);
		delete[] PanEvap;
	}
	delete [] year; 
	delete [] month; 
	delete [] day; 
	delete [] hour;
}

/***************************************************************************
**
** tEvapoTrans::robustNess() Function
**
** This function is used to check variables for NO_DATA flag = 9999.99
** It searches the double array forward and backwards, substituting the
** NO_DATA value with the previously read valid entry.
**
***************************************************************************/
void tEvapoTrans::robustNess(double *variable, int size)
{
	double lastStored = 9999.99;
	double firstStored = 9999.99;
	
	for (int ct=0;ct<size;ct++) {
		if (fabs(variable[ct]-9999.99) > 1.0E-3)
			lastStored = variable[ct];   
		if (fabs(variable[ct]-9999.99) < 1.0E-3)
			variable[ct] = lastStored; 
	}   
	
	for (int dt=size-1;dt>=0;dt--) {
		if (fabs(variable[dt]-9999.99) > 1.0E-3)
			firstStored = variable[dt];
		if (fabs(variable[dt]-9999.99) < 1.0E-3)
			variable[dt] = firstStored;
	}
	}

/***************************************************************************
**
** tEvapoTrans::readHydroMetGrid() Function
**
** Reads a file (*.gdf) from HYDROMETGRID keyword containing the base names of 
** the various input meteorologic parameter grids along with the extension
** used for the filename. These follow a string that identifies the line
** with the parameters (ie. PA, TD, XC, US, TA, TS, NR). If no data is
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
** PA ../PATH/PAbase txt
** TD ../PATH/TDbase txt
** NR NO_DATA NO_DATA
**
***************************************************************************/
void tEvapoTrans::readHydroMetGrid(char *gridFile) 
{
	int numParameters;
	
	Cout<<"\nReading HydroMeteorological Grid File: "; 
	Cout<< gridFile<<"..."<<endl<<flush;
	
	ifstream readFile(gridFile); 
	if (!readFile) {
		cout << "\nFile "<<gridFile<<" not found!" << endl;
		cout << "Exiting Program...\n\n"<<endl;
		exit(1);}
	
	readFile >> numParameters; 
	nParm = numParameters;
	readFile >> gridlat;
	readFile >> gridlong;
	readFile >> gridgmt;
	
	gridBaseNames = new char*[numParameters];
	gridExtNames = new char*[numParameters];
	gridParamNames = new char*[numParameters];
	
	for (int ct=0;ct<numParameters;ct++) {
		gridParamNames[ct] = new char[10];
		gridBaseNames[ct] = new char[kName];
		gridExtNames[ct] = new char[kMaxExt];
		readFile >> gridParamNames[ct];
		readFile >> gridBaseNames[ct];
		readFile >> gridExtNames[ct];
	}
	}

/***************************************************************************
**
** tEvapoTrans::readLUGrid() Function
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
**
**  Remember that Base Name and Extension for LU Grids cannot be NO_DATA now
**  (though they can be in the HydroMetGrids)
**
** Example:
** 3
** 32.5 -100.3 -6
** AL ../PATH/ALbase txt
** TF ../PATH/TFbase txt
** VH ../PATH/VHbase txt
**
***************************************************************************/
void tEvapoTrans::readLUGrid(char *gridFile) 
{
	int numParameters;
	
	Cout<<"\nReading Land-Use Data Grid File: "; 
	Cout<< gridFile<<"..."<<endl<<flush;
	
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
		LUgridParamNames[ct] = new char[kMaxExt];
		LUgridBaseNames[ct] = new char[kName];
		LUgridExtNames[ct] = new char[kMaxExt];
		readFile >> LUgridParamNames[ct];

		if ( (strcmp(LUgridParamNames[ct],"AL")!=0) &&
				(strcmp(LUgridParamNames[ct],"TF")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"VH")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"SR")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"VF")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"CS")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"IC")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"CC")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"DC")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"DE")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"OT")!=0) &&
		  		(strcmp(LUgridParamNames[ct],"LA")!=0) ) {
			
			Cout << "\nA land use parameter name in the LU gdf file is an unexpected one."<<endl;
			Cout << "\nExpected variables: AL,TF,VH,SR,VF,CS,IC,CC,DC,DE,OT or LA" << endl;
			Cout << "\tCheck and re-run the program" << endl;
			Cout << "\nExiting Program..."<<endl<<endl;
			exit(1);

		}

		readFile >> LUgridBaseNames[ct];

		if (strcmp(LUgridBaseNames[ct],"NO_DATA")==0) {
			Cout << "\nCannot use NO_DATA for LU Grids"<<endl;
			Cout << "\nExiting Program..."<<endl<<endl;
			exit(1);
		}

		readFile >> LUgridExtNames[ct];
	}
	}

/***************************************************************************
**
** tEvapoTrans::createVariant() Function
**
** Initializes the tVariant pointer objects for each of the time-varying
** meteorologic parameters. Assigns zero using newVariable(char *) to
** the tCNode property. If NO_DATA encountered, assigns 9999.999 to 
** variable in tCNode by calling noData member function in tVariant.
**
***************************************************************************/
void tEvapoTrans::createVariant() 
{

	if (evapotransOption != 4) {
		for (int ct=0;ct<nParm;ct++) { 
			if (strcmp(gridParamNames[ct],"PA")==0) {
				airpressure = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					airpressure->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					airpressure->newVariable(gridParamNames[ct]);}
				else
					airpressure->noData(gridParamNames[ct]);
			}
			if (strcmp(gridParamNames[ct],"TD")==0) {
				dewtemperature = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					dewtemperature->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					dewtemperature->newVariable(gridParamNames[ct]);}
				else
					dewtemperature->noData(gridParamNames[ct]);
			}
			if (strcmp(gridParamNames[ct],"XC")==0) {
				skycover = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					skycover->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					skycover->newVariable(gridParamNames[ct]);}
				else
					skycover->noData(gridParamNames[ct]);
			}
			if (strcmp(gridParamNames[ct],"US")==0) {
				windspeed = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					windspeed->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					windspeed->newVariable(gridParamNames[ct]);}
				else
					windspeed->noData(gridParamNames[ct]);
			}
			if (strcmp(gridParamNames[ct],"TA")==0) {
				airtemperature = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					airtemperature->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					airtemperature->newVariable(gridParamNames[ct]);}
				else
					airtemperature->noData(gridParamNames[ct]);
			}
			if (strcmp(gridParamNames[ct],"TS")==0) {
				surftemperature = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					surftemperature->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					surftemperature->newVariable(gridParamNames[ct]);}
				else
					surftemperature->noData(gridParamNames[ct]);
			}
			if (strcmp(gridParamNames[ct],"NR")==0) {
				netradiation = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					netradiation->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					netradiation->newVariable(gridParamNames[ct]);}
				else
					netradiation->noData(gridParamNames[ct]);
			}
			if (strcmp(gridParamNames[ct],"RH")==0) {
				relhumidity = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					relhumidity->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					relhumidity->newVariable(gridParamNames[ct]);}
				else
					relhumidity->noData(gridParamNames[ct]);
			}
			if (strcmp(gridParamNames[ct],"VP")==0) {
				vaporpressure = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					vaporpressure->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					vaporpressure->newVariable(gridParamNames[ct]);}
				else
					vaporpressure->noData(gridParamNames[ct]);
			}
			if (strcmp(gridParamNames[ct],"IS")==0) {  //E.R.V. 3/6/2012
				incomingsolar = new tVariant(gridPtr,respPtr);
				if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
					incomingsolar->setFileNames(gridBaseNames[ct], gridExtNames[ct]);
					incomingsolar->newVariable(gridParamNames[ct]);}
				else
					incomingsolar->noData(gridParamNames[ct]);
			}

		}
	}   
	else {
		if (strcmp(gridParamNames[0],"ET")==0) {
			evapotranspiration = new tVariant(gridPtr,respPtr);
			evapotranspiration->setFileNames(gridBaseNames[0], gridExtNames[0]);
			evapotranspiration->newVariable(gridParamNames[0]);}
		else {
			cout<<"\nError in ET Grid Input...."<<endl;
			evapotranspiration->noData(gridParamNames[0]);
			exit(1);}
	
	}
}

/***************************************************************************
**
** tEvapoTrans::createVariantLU() Function
**
** Initializes the tVariant pointer objects for each of the time-varying
** LAND-USE parameters. Assigns zero using newVariable(char *) to
** the tCNode property. 
** --REMEMBER-- that NO_DATA should not be encountered now, hence no call 
** to noData member function in tVariant!!!
**
***************************************************************************/
void tEvapoTrans::createVariantLU() 
{
	for (int ct=0;ct<nParmLU;ct++) { 
		if (strcmp(LUgridParamNames[ct],"AL")==0) {
			LandUseAlbGrid = new tVariant(gridPtr,respPtr);
			LandUseAlbGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(LandUseAlbGrid, LUgridParamNames[ct]);
			LandUseAlbGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"TF")==0) {
			ThroughFallGrid = new tVariant(gridPtr,respPtr);
			ThroughFallGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(ThroughFallGrid, LUgridParamNames[ct]);
			ThroughFallGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"VH")==0) {
			VegHeightGrid = new tVariant(gridPtr,respPtr);
			VegHeightGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(VegHeightGrid, LUgridParamNames[ct]);
			VegHeightGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"SR")==0) {
			StomResGrid = new tVariant(gridPtr,respPtr);
			StomResGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(StomResGrid, LUgridParamNames[ct]);
			StomResGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"VF")==0) {
			VegFractGrid = new tVariant(gridPtr,respPtr);
			VegFractGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(VegFractGrid, LUgridParamNames[ct]);
			VegFractGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"CS")==0) {
			CanStorParamGrid = new tVariant(gridPtr,respPtr);
			CanStorParamGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(CanStorParamGrid, LUgridParamNames[ct]);
			CanStorParamGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"IC")==0) {
			IntercepCoeffGrid = new tVariant(gridPtr,respPtr);
			IntercepCoeffGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(IntercepCoeffGrid, LUgridParamNames[ct]);
			IntercepCoeffGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"CC")==0) {
			CanFieldCapGrid = new tVariant(gridPtr,respPtr);
			CanFieldCapGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(CanFieldCapGrid, LUgridParamNames[ct]);
			CanFieldCapGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"DC")==0) {
			DrainCoeffGrid = new tVariant(gridPtr,respPtr);
			DrainCoeffGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(DrainCoeffGrid, LUgridParamNames[ct]);
			DrainCoeffGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"DE")==0) {
			DrainExpParGrid = new tVariant(gridPtr,respPtr);
			DrainExpParGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(DrainExpParGrid, LUgridParamNames[ct]);
			DrainExpParGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"OT")==0) {
			OptTransmCoeffGrid = new tVariant(gridPtr,respPtr);
			OptTransmCoeffGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(OptTransmCoeffGrid, LUgridParamNames[ct]);
			OptTransmCoeffGrid->newVariable(LUgridParamNames[ct]);
		}
		if (strcmp(LUgridParamNames[ct],"LA")==0) {
			LeafAIGrid = new tVariant(gridPtr,respPtr);
			LeafAIGrid->setFileNames(LUgridBaseNames[ct], LUgridExtNames[ct]);
			SetGridTimeInfoVariables(LeafAIGrid,LUgridParamNames[ct]);
			LeafAIGrid->newVariable(LUgridParamNames[ct]);
		}
	}
}

/***************************************************************************
**
** tEvapoTrans::newHydroMetStochData() Function
**
** Assigns the values of the current meteorological parameters that are
** simulated by functions of the class 'tHydroMetStoch'
**
***************************************************************************/
void tEvapoTrans::newHydroMetStochData(int time) 
{
	weatherSimul->SimulateHydrometVars();
	
	// Obtain simulated values from tHydroMetStoch
	airTemp   = weatherSimul->getAirTemp();
	windSpeed = weatherSimul->getWindSpeed();
	skyCover  = weatherSimul->getSkyCover();
	dewTemp   = weatherSimul->getDewTemp();
	rHumidity = weatherSimul->getRHumidity();
	vPress    = weatherSimul->getVaporPress();
	atmPress  = weatherSimul->getAtmPress();
	
	// These are currently assigned 9999.99
	surfTemp  = weatherSimul->getSurfTemp();
	netRad    = weatherSimul->getNetRad();
	
	if (!time) {
		latitude =  weatherSimul->getLat(1);
		longitude = weatherSimul->getLong(1);
		gmt = weatherSimul->getGmt();
		Tso = weatherSimul->getAirTemp() + 273.15;
		Tlo = Tso; //Assumption
	}
	}

/***************************************************************************
**
** tEvapoTrans::newHydroMetData() Function
**
** Assigns the values of the current meteorological parameters to the 
** nodes based on the results of the tResample Thiessen polygon routine.
**
***************************************************************************/
void tEvapoTrans::newHydroMetData(int time) 
{

	// Obtain values from tHydroMet
	for (int i=0; i<numStations;i++) {
		if (thisStation == weatherStations[i].getStation()) {
			if (evapotransOption != 4) {
				airTemp = weatherStations[i].getAirTemp(time);

				// SKY2008Snow from AJR2007
				airTemp += tempLapseRate*(elevation - weatherStations[i].getOther()); //lapse rate added by AJR 2007 @ NMT

				dewTemp = weatherStations[i].getDewTemp(time);
				surfTemp = weatherStations[i].getSurfTemp(time);
				rHumidity = weatherStations[i].getRHumidity(time);
				vPress = weatherStations[i].getVaporPress(time);
				atmPress = weatherStations[i].getAtmPress(time);
				windSpeed = weatherStations[i].getWindSpeed(time);
				skyCover = weatherStations[i].getSkyCover(time);
				nodeHour = weatherStations[i].getHour(time);
                inShortR = weatherStations[i].getRadGlobal(time);

				// For run-time checks of input data
				if (false) {
					cout<<"\t---> Time = "<<time<<"; Station ID = "<<i<<endl;
					cout<<"\t---> Tair = "<<airTemp<<"; dewTemp = "<<dewTemp<<";"<<endl;
					cout<<"\trHumidity = "<<rHumidity<<"; vapPress = "<<vPress<<";"<<endl;
					cout<<"\tatmPress = "<<atmPress<<"; windSpeed = "<<windSpeed
						<<"; skyCover = "<<skyCover<<";"<<endl;
				}
				
				if (weatherStations[i].getParm() >= 11)
					netRad = weatherStations[i].getNetRad(time);
				
				if (time == 0) {
					latitude = weatherStations[i].getLat(1);
					longitude = weatherStations[i].getLong(1);
					gmt = weatherStations[i].getGmt();
					Tso = weatherStations[i].getAirTemp(time) + 273.15;
					Tlo = Tso;
					
					//Find the Available Humidity Data
					if (fabs(dewTemp-9999.99)<1.0E-3 && fabs(vPress-9999.99)<1.0E-3){
						dewHumFlag = 0;}
					else if (fabs(rHumidity-9999.99)<1.0E-3 && fabs(vPress-9999.99)<1.0E-3){
						dewHumFlag = 1;}
					else if (fabs(rHumidity-9999.99)<1.0E-3 && fabs(dewTemp-9999.99)<1.0E-3){
						dewHumFlag = 2;} 
				}
			}
			else {
				panEvap = weatherStations[i].getPanEvap(time);
				if (time == 0) {
					coeffPan = weatherStations[i].getOther();
				}
			}
		}
	}
	}

/***************************************************************************
**
** tEvapoTrans::newHydroMetGridData() Function
**
** Assigns the values of the current meteorological parameters to the 
** nodes based on the results of the tResample grid input.
**
***************************************************************************/
void tEvapoTrans::newHydroMetGridData(tCNode * cNode) {   
	if (evapotransOption != 4) {
		airTemp = cNode->getAirTemp();
		dewTemp = cNode->getDewTemp();
		surfTemp = cNode->getSurfTemp();
		rHumidity = cNode->getRelHumid();
		atmPress = cNode->getAirPressure();
		windSpeed = cNode->getWindSpeed();
		skyCover = cNode->getSkyCover();
		netRad = cNode->getNetRad();
        inShortR = cNode->getShortRadIn(); //E.R.V 3/6/2012
		vPress = cNode->getVapPressure();
		nodeHour = timer->hour;

		if (timeCount == 0) {
			latitude = gridlat;
			longitude = gridlong;
			gmt = gridgmt;
			Tso = cNode->getAirTemp() + 273.15; 
			Tlo = Tso; 
			
			//Find the Available Humidity Data
			if (fabs(dewTemp-9999.99)<1.0E-3 && fabs(vPress-9999.99)<1.0E-3){
				dewHumFlag = 0;}
			else if (fabs(rHumidity-9999.99)<1.0E-3 && fabs(vPress-9999.99)<1.0E-3){
				dewHumFlag = 1;}
			else if (fabs(rHumidity-9999.99)<1.0E-3 && fabs(dewTemp-9999.99)<1.0E-3) {
				dewHumFlag = 2;}
		}
	}
	else {
		panEvap = cNode->getGridET();
	}	
	}

/***************************************************************************
**
** tEvapoTrans::newLUGridData() Function
**
** Assigns values of current land use parameters to nodes based on results of 
** tResample grid input from updateLUVarOfBothGrids &/or updateLUVarOfPrevGrid.
** TODO check that same landcover parameters are being updated as provided from a table
***************************************************************************/
void tEvapoTrans::newLUGridData(tCNode * cNode) 
{ 
	for (int ct=0;ct<nParmLU;ct++) { 
		if (strcmp(LUgridParamNames[ct],"AL")==0) {
			if ( (evapotransOption == 1) ||
					(evapotransOption == 2) ||
					(evapotransOption == 3) ){
				coeffAl = cNode->getLandUseAlb();
			}
		}
		if (strcmp(LUgridParamNames[ct],"VH")==0) {
			if (evapotransOption == 1) {
				coeffH = cNode->getVegHeight();
			}
		}
		if (strcmp(LUgridParamNames[ct],"OT")==0) {
			if (evapotransOption == 1) {
				coeffKt = cNode->getOptTransmCoeff();
			}
		}
		if (strcmp(LUgridParamNames[ct],"SR")==0) {
			if (evapotransOption == 1) {
				coeffRs = cNode->getStomRes();	
			}
		}		
		if (strcmp(LUgridParamNames[ct],"VF")==0) {
			if ( (evapotransOption == 1) ||
					(evapotransOption == 2) ||
					(evapotransOption == 3) ||
					(evapotransOption == 4) ){
				coeffV = cNode->getVegFraction();

                if (coeffV >= 1.0) //prevents loss of snow when unloaded from canopy WR 05/12/2024
                    coeffV = 0.99;
			}
		}
		if (strcmp(LUgridParamNames[ct],"LA")==0) {			
			coeffLAI = cNode->getLeafAI(); // SKY2008Snow
		}
	}

	if (IfNotFirstTStepLU == 0) { // modified from correction by SY, TM: 11/19/07
		latitude = LUgridlat;
		longitude = LUgridlong;
		gmt = LUgridgmt;
		IfNotFirstTStepLU = 1;
	}

	}

/***************************************************************************
**
** tEvapoTrans::resampleGrids() Function
**
**
***************************************************************************/
void tEvapoTrans::resampleGrids(tRunTimer *t) 
{
	if (evapotransOption!=4) {
		for (int ct=0;ct<nParm;ct++) { 
			if (strcmp(gridBaseNames[ct],"NO_DATA")!=0) {
				if (strcmp(gridParamNames[ct],"TA")==0) {
					airtemperature->composeFileName(t);
					airtemperature->updateVariable(gridParamNames[ct]);}
				if (strcmp(gridParamNames[ct],"PA")==0) {
					airpressure->composeFileName(t);
					airpressure->updateVariable(gridParamNames[ct]);}
				if (strcmp(gridParamNames[ct],"XC")==0) {
					skycover->composeFileName(t);
					skycover->updateVariable(gridParamNames[ct]);}
				if (strcmp(gridParamNames[ct],"US")==0) {
					windspeed->composeFileName(t);
					windspeed->updateVariable(gridParamNames[ct]);}
				if (strcmp(gridParamNames[ct],"TD")==0) {
					dewtemperature->composeFileName(t);
					dewtemperature->updateVariable(gridParamNames[ct]);}
				if (strcmp(gridParamNames[ct],"RH")==0) {
					relhumidity->composeFileName(t);
					relhumidity->updateVariable(gridParamNames[ct]);}
				if (strcmp(gridParamNames[ct],"TS")==0) {
					surftemperature->composeFileName(t);
					surftemperature->updateVariable(gridParamNames[ct]);}
				if (strcmp(gridParamNames[ct],"NR")==0) {
					netradiation->composeFileName(t);
					netradiation->updateVariable(gridParamNames[ct]);}
				if (strcmp(gridParamNames[ct],"VP")==0) {
					vaporpressure->composeFileName(t);
					vaporpressure->updateVariable(gridParamNames[ct]);}
				if (strcmp(gridParamNames[ct],"IS")==0) {   //E.R.V 3/6/2012
					incomingsolar->composeFileName(t);
					incomingsolar->updateVariable(gridParamNames[ct]);}
			}  
		}
	}
	else {
		if (strcmp(gridBaseNames[0],"NO_DATA")!=0) {
			if (strcmp(gridParamNames[0],"ET")==0) {
				evapotranspiration->composeFileName(t);
				evapotranspiration->updateVariable(gridParamNames[0]);}
		}
	}
	}


/***************************************************************************
**
** tEvapoTrans::initialLUGridAssignment() Function
**
**
***************************************************************************/
void tEvapoTrans::initialLUGridAssignment()
{

  for (int ct=0;ct<nParmLU;ct++) {
    
    if (strcmp(LUgridParamNames[ct],"AL")==0) {
      if ( (timer->getCurrentTime())>(double(ALgridhours[NowTillWhichALgrid]) && numALfiles >1) ) {
	while ( (timer->getCurrentTime())>(double(ALgridhours[NowTillWhichALgrid])) ) {
	  NowTillWhichALgrid++;}
	LandUseAlbGrid->updateLUVarOfBothGrids("AL", ALgridFileNames[NowTillWhichALgrid]);
	LandUseAlbGrid->updateLUVarOfPrevGrid("AL", ALgridFileNames[NowTillWhichALgrid-1]);
      }
      else {
	LandUseAlbGrid->updateLUVarOfBothGrids("AL", ALgridFileNames[1]);
	LandUseAlbGrid->updateLUVarOfPrevGrid("AL", ALgridFileNames[1]);}
      
    }
    if (strcmp(LUgridParamNames[ct],"TF")==0) {
      if ( (timer->getCurrentTime())>(double(TFgridhours[NowTillWhichTFgrid])) && numTFfiles > 1) {
	while ( (timer->getCurrentTime())>(double(TFgridhours[NowTillWhichTFgrid])) ) {
	  NowTillWhichTFgrid++;}
	ThroughFallGrid->updateLUVarOfBothGrids("TF", TFgridFileNames[NowTillWhichTFgrid]);
	ThroughFallGrid->updateLUVarOfPrevGrid("TF", TFgridFileNames[NowTillWhichTFgrid-1]);
      }
      else {
	ThroughFallGrid->updateLUVarOfBothGrids("TF", TFgridFileNames[1]);
	ThroughFallGrid->updateLUVarOfPrevGrid("TF", TFgridFileNames[1]);}
    }
    if (strcmp(LUgridParamNames[ct],"VH")==0) {
      if ( (timer->getCurrentTime())>(double(VHgridhours[NowTillWhichVHgrid])) && numVHfiles > 1) {
	while ( (timer->getCurrentTime())>(double(VHgridhours[NowTillWhichVHgrid])) ) {
	  NowTillWhichVHgrid++;
	}
	VegHeightGrid->updateLUVarOfBothGrids("VH", VHgridFileNames[NowTillWhichVHgrid]);
	VegHeightGrid->updateLUVarOfPrevGrid("VH", VHgridFileNames[NowTillWhichVHgrid-1]);
      }
      else {
	VegHeightGrid->updateLUVarOfBothGrids("VH", VHgridFileNames[1]);
	VegHeightGrid->updateLUVarOfPrevGrid("VH", VHgridFileNames[1]);
      }
    }
    if (strcmp(LUgridParamNames[ct],"SR")==0) {
      if ( (timer->getCurrentTime())>(double(SRgridhours[NowTillWhichSRgrid])) && numSRfiles >1) {
	while ( (timer->getCurrentTime())>(double(SRgridhours[NowTillWhichSRgrid])) ) {
	  NowTillWhichSRgrid++;
	}
	StomResGrid->updateLUVarOfBothGrids("SR", SRgridFileNames[NowTillWhichSRgrid]);
	StomResGrid->updateLUVarOfPrevGrid("SR", SRgridFileNames[NowTillWhichSRgrid-1]);
      }
      else {
	StomResGrid->updateLUVarOfBothGrids("SR", SRgridFileNames[1]);
	StomResGrid->updateLUVarOfPrevGrid("SR", SRgridFileNames[1]);
      }
    }
    if (strcmp(LUgridParamNames[ct],"VF")==0) {
      if ( (timer->getCurrentTime())>(double(VFgridhours[NowTillWhichVFgrid])) && numVFfiles > 1 ) {
	while ( (timer->getCurrentTime())>(double(VFgridhours[NowTillWhichVFgrid])) ) {
	  NowTillWhichVFgrid++;
	}
	VegFractGrid->updateLUVarOfBothGrids("VF", VFgridFileNames[NowTillWhichVFgrid]);
	VegFractGrid->updateLUVarOfPrevGrid("VF", VFgridFileNames[NowTillWhichVFgrid-1]);
      }
      else {
	VegFractGrid->updateLUVarOfBothGrids("VF", VFgridFileNames[1]);
	VegFractGrid->updateLUVarOfPrevGrid("VF", VFgridFileNames[1]);
      }
    }
    if (strcmp(LUgridParamNames[ct],"CS")==0) {
      if ( (timer->getCurrentTime())>(double(CSgridhours[NowTillWhichCSgrid])) && numCSfiles > 1) {
	while ( (timer->getCurrentTime())>(double(CSgridhours[NowTillWhichCSgrid])) ) {
	  NowTillWhichCSgrid++;
	}
	CanStorParamGrid->updateLUVarOfBothGrids("CS", CSgridFileNames[NowTillWhichCSgrid]);
	CanStorParamGrid->updateLUVarOfPrevGrid("CS", CSgridFileNames[NowTillWhichCSgrid-1]);
      }
      else {
	CanStorParamGrid->updateLUVarOfBothGrids("CS", CSgridFileNames[1]);
	CanStorParamGrid->updateLUVarOfPrevGrid("CS", CSgridFileNames[1]);
      }
    }
    if (strcmp(LUgridParamNames[ct],"IC")==0) {
      if ( (timer->getCurrentTime())>(double(ICgridhours[NowTillWhichICgrid])) && numICfiles > 1) {
	while ( (timer->getCurrentTime())>(double(ICgridhours[NowTillWhichICgrid])) ) {
	  NowTillWhichICgrid++;
	}
	IntercepCoeffGrid->updateLUVarOfBothGrids("IC", ICgridFileNames[NowTillWhichICgrid]);
	IntercepCoeffGrid->updateLUVarOfPrevGrid("IC", ICgridFileNames[NowTillWhichICgrid-1]);
      }
      else {
	IntercepCoeffGrid->updateLUVarOfBothGrids("IC", ICgridFileNames[1]);
	IntercepCoeffGrid->updateLUVarOfPrevGrid("IC", ICgridFileNames[1]);
      }
    }
    if (strcmp(LUgridParamNames[ct],"CC")==0) {
      if ( (timer->getCurrentTime())>(double(CCgridhours[NowTillWhichCCgrid])) && numCCfiles > 1 ) {
	while ( (timer->getCurrentTime())>(double(CCgridhours[NowTillWhichCCgrid])) ) {
	  NowTillWhichCCgrid++;
	}
	CanFieldCapGrid->updateLUVarOfBothGrids("CC", CCgridFileNames[NowTillWhichCCgrid]);
	CanFieldCapGrid->updateLUVarOfPrevGrid("CC", CCgridFileNames[NowTillWhichCCgrid-1]);
      }
      else {
	CanFieldCapGrid->updateLUVarOfBothGrids("CC", CCgridFileNames[1]);
	CanFieldCapGrid->updateLUVarOfPrevGrid("CC", CCgridFileNames[1]);
      }
    }
    if (strcmp(LUgridParamNames[ct],"DC")==0) {
      if ( (timer->getCurrentTime())>(double(DCgridhours[NowTillWhichDCgrid])) && numDCfiles > 1) {
	while ( (timer->getCurrentTime())>(double(DCgridhours[NowTillWhichDCgrid])) ) {
	  NowTillWhichDCgrid++;
	}
	DrainCoeffGrid->updateLUVarOfBothGrids("DC", DCgridFileNames[NowTillWhichDCgrid]);
	DrainCoeffGrid->updateLUVarOfPrevGrid("DC", DCgridFileNames[NowTillWhichDCgrid-1]);
      }
      else {
	DrainCoeffGrid->updateLUVarOfBothGrids("DC", DCgridFileNames[1]);
	DrainCoeffGrid->updateLUVarOfPrevGrid("DC", DCgridFileNames[1]);
      }
    }
    if (strcmp(LUgridParamNames[ct],"DE")==0) {
      if ( (timer->getCurrentTime())>(double(DEgridhours[NowTillWhichDEgrid])) && numDEfiles > 1) {
	while ( (timer->getCurrentTime())>(double(DEgridhours[NowTillWhichDEgrid])) ) {
	  NowTillWhichDEgrid++;
	}
	DrainExpParGrid->updateLUVarOfBothGrids("DE", DEgridFileNames[NowTillWhichDEgrid]);
	DrainExpParGrid->updateLUVarOfPrevGrid("DE", DEgridFileNames[NowTillWhichDEgrid-1]);
      }
      else {
	DrainExpParGrid->updateLUVarOfBothGrids("DE", DEgridFileNames[1]);
	DrainExpParGrid->updateLUVarOfPrevGrid("DE", DEgridFileNames[1]);
      }
    }
    if (strcmp(LUgridParamNames[ct],"OT")==0) {
      if ( (timer->getCurrentTime())>(double(OTgridhours[NowTillWhichOTgrid])) && numOTfiles > 1) {
	while ( (timer->getCurrentTime())>(double(OTgridhours[NowTillWhichOTgrid])) ) {
	  NowTillWhichOTgrid++;
	}
	OptTransmCoeffGrid->updateLUVarOfBothGrids("OT", OTgridFileNames[NowTillWhichOTgrid]);
	OptTransmCoeffGrid->updateLUVarOfPrevGrid("OT", OTgridFileNames[NowTillWhichOTgrid-1]);
      }
      else {
	OptTransmCoeffGrid->updateLUVarOfBothGrids("OT", OTgridFileNames[1]);
	OptTransmCoeffGrid->updateLUVarOfPrevGrid("OT", OTgridFileNames[1]);
      }
    }
    if (strcmp(LUgridParamNames[ct],"LA")==0) {
      if ( (timer->getCurrentTime())>(double(LAgridhours[NowTillWhichLAgrid])) && numLAfiles > 1 ) {
	while ( (timer->getCurrentTime())>(double(LAgridhours[NowTillWhichLAgrid])) ) {
	  NowTillWhichLAgrid++;
	}
	LeafAIGrid->updateLUVarOfBothGrids("LA", LAgridFileNames[NowTillWhichLAgrid]);
	LeafAIGrid->updateLUVarOfPrevGrid("LA", LAgridFileNames[NowTillWhichLAgrid-1]);
      }
      else {
	LeafAIGrid->updateLUVarOfBothGrids("LA", LAgridFileNames[1]);
	LeafAIGrid->updateLUVarOfPrevGrid("LA", LAgridFileNames[1]);
      }
    }
  } // end for loop

  }

/***************************************************************************
**
** tEvapoTrans::LUGridAssignment() Function
**
**
***************************************************************************/

void tEvapoTrans::LUGridAssignment()
{
  for (int ct=0;ct<nParmLU;ct++) { 
    if (strcmp(LUgridParamNames[ct],"AL")==0) { //if parameter exists
      if (NowTillWhichALgrid <= numALfiles) { // if you have not exceeded the file amount
	if ((timer->getCurrentTime())>(double(ALgridhours[NowTillWhichALgrid]))) {// the time is greater then the ALgrid hours
	  NowTillWhichALgrid++;
	  if ( (NowTillWhichALgrid-1)<numALfiles) {
	    LandUseAlbGrid->updateLUVarOfBothGrids("AL", ALgridFileNames[NowTillWhichALgrid]);
	  }					
	  else {
	    LandUseAlbGrid->updateLUVarOfPrevGrid("AL", ALgridFileNames[numALfiles]);
	  }
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"TF")==0) {
      if (NowTillWhichTFgrid <= numTFfiles) {
	if ((timer->getCurrentTime())>(double(TFgridhours[NowTillWhichTFgrid]))) { 
	  NowTillWhichTFgrid++;
	  if ((NowTillWhichTFgrid-1)<numTFfiles) {							
	    ThroughFallGrid->updateLUVarOfBothGrids("TF", TFgridFileNames[NowTillWhichTFgrid]);
	  }
	  else {
	    ThroughFallGrid->updateLUVarOfPrevGrid("TF", TFgridFileNames[numTFfiles]);
	  }
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"VH")==0) {
      if (NowTillWhichVHgrid<=numVHfiles) {		
	if ((timer->getCurrentTime())>(double(VHgridhours[NowTillWhichVHgrid]))) {
	  NowTillWhichVHgrid++;
	  if ((NowTillWhichVHgrid-1)<numVHfiles) {
	    VegHeightGrid->updateLUVarOfBothGrids("VH", VHgridFileNames[NowTillWhichVHgrid]);
	  }
	  else {
	    VegHeightGrid->updateLUVarOfPrevGrid("VH", VHgridFileNames[numVHfiles]);
	  }
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"SR")==0) {
      if (NowTillWhichSRgrid<=numSRfiles) {
	if ((timer->getCurrentTime())>(double(SRgridhours[NowTillWhichSRgrid]))) { 
	  NowTillWhichSRgrid++;
	  if ((NowTillWhichSRgrid-1)<numSRfiles) {
	    StomResGrid->updateLUVarOfBothGrids("SR", SRgridFileNames[NowTillWhichSRgrid]);
	  }
	  else {
	    StomResGrid->updateLUVarOfPrevGrid("SR", SRgridFileNames[numSRfiles]);
	  }	
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"VF")==0) {
      if (NowTillWhichVFgrid<=numVFfiles) {
	if ((timer->getCurrentTime())>(double(VFgridhours[NowTillWhichVFgrid]))) { 
	  NowTillWhichVFgrid++;
	  if ((NowTillWhichVFgrid-1)<numVFfiles) {
	    VegFractGrid->updateLUVarOfBothGrids("VF", VFgridFileNames[NowTillWhichVFgrid]);
	  }
	  else {
	    VegFractGrid->updateLUVarOfPrevGrid("VF", VFgridFileNames[numVFfiles]);
	  }	
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"CS")==0) {
      if (NowTillWhichCSgrid<=numCSfiles) {
	if ((timer->getCurrentTime())>(double(CSgridhours[NowTillWhichCSgrid]))) { 
	  NowTillWhichCSgrid++;
	  if ((NowTillWhichCSgrid-1)<numCSfiles) {
	    CanStorParamGrid->updateLUVarOfBothGrids("CS", CSgridFileNames[NowTillWhichCSgrid]);
	  }
	  else {
	    CanStorParamGrid->updateLUVarOfPrevGrid("CS", CSgridFileNames[numCSfiles]);
	  }	
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"IC")==0) {
      if (NowTillWhichICgrid<=numICfiles) {
	if ((timer->getCurrentTime())>(double(ICgridhours[NowTillWhichICgrid]))) { 
	  NowTillWhichICgrid++;
	  if ((NowTillWhichICgrid-1)<numICfiles) {
	    IntercepCoeffGrid->updateLUVarOfBothGrids("IC", ICgridFileNames[NowTillWhichICgrid]);
	  }
	  else {
	    IntercepCoeffGrid->updateLUVarOfPrevGrid("IC", ICgridFileNames[numICfiles]);
	  }
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"CC")==0) {
      if (NowTillWhichCCgrid<=numCCfiles) {
	if ((timer->getCurrentTime())>(double(CCgridhours[NowTillWhichCCgrid]))) { 
	  NowTillWhichCCgrid++;
	  if ((NowTillWhichCCgrid-1)<numCCfiles) {							
	    CanFieldCapGrid->updateLUVarOfBothGrids("CC", CCgridFileNames[NowTillWhichCCgrid]);
	  }
	  else {
	    CanFieldCapGrid->updateLUVarOfPrevGrid("CC", CCgridFileNames[numCCfiles]);
	  }
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"DC")==0) {
      if (NowTillWhichDCgrid<=numDCfiles) {
	if ((timer->getCurrentTime())>(double(DCgridhours[NowTillWhichDCgrid]))) { 
	  NowTillWhichDCgrid++;
	  if ((NowTillWhichDCgrid-1)<numDCfiles) {
	    DrainCoeffGrid->updateLUVarOfBothGrids("DC", DCgridFileNames[NowTillWhichDCgrid]);
	  }
	  else {
	    DrainCoeffGrid->updateLUVarOfPrevGrid("DC", DCgridFileNames[numDCfiles]);
	  }
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"DE")==0) {
      if (NowTillWhichDEgrid<=numDEfiles) {
	if ((timer->getCurrentTime())>(double(DEgridhours[NowTillWhichDEgrid]))) {
	  NowTillWhichDEgrid++;
	  if ((NowTillWhichDEgrid-1)<numDEfiles) {	
	    DrainExpParGrid->updateLUVarOfBothGrids("DE", DEgridFileNames[NowTillWhichDEgrid]);
	  }
	  else {
	    DrainExpParGrid->updateLUVarOfPrevGrid("DE", DEgridFileNames[numDEfiles]);
	  }
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"OT")==0) {
      if (NowTillWhichOTgrid<=numOTfiles) {
	if ((timer->getCurrentTime())>(double(OTgridhours[NowTillWhichOTgrid]))) {
	  NowTillWhichOTgrid++;
	  if ((NowTillWhichOTgrid-1)<numOTfiles) {	
	    OptTransmCoeffGrid->updateLUVarOfBothGrids("OT", OTgridFileNames[NowTillWhichOTgrid]);
	  }
	  else {
	    OptTransmCoeffGrid->updateLUVarOfPrevGrid("OT", OTgridFileNames[numOTfiles]);
	  }
	}
      }
    }
    if (strcmp(LUgridParamNames[ct],"LA")==0) {
      if (NowTillWhichLAgrid<=numLAfiles) {
	if ((timer->getCurrentTime())>(double(LAgridhours[NowTillWhichLAgrid]))) { 
	  NowTillWhichLAgrid++;
	  if ((NowTillWhichLAgrid-1)<numLAfiles) {
	    LeafAIGrid->updateLUVarOfBothGrids("LA", LAgridFileNames[NowTillWhichLAgrid]);
	  }
	  else {
	    LeafAIGrid->updateLUVarOfPrevGrid("LA", LAgridFileNames[numLAfiles]);
	  }
	}
      }
    }
  }
}


/***************************************************************************
**
** interpolateLUGrids() Function
**
***************************************************************************/
void tEvapoTrans::interpolateLUGrids(tCNode* cNode)
{
  for (int ct=0;ct<nParmLU;ct++) { 
    if ( (strcmp(LUgridParamNames[ct],"AL")==0) && (NowTillWhichALgrid > 1) &&
	 ( NowTillWhichALgrid < (numALfiles+1) ) ) 
      {
	cNode->setLandUseAlb( cNode->getLandUseAlbInPrevGrid()+
			      ( cNode->getLandUseAlbInUntilGrid() - cNode->getLandUseAlbInPrevGrid() )*
			      ( timer->getCurrentTime() - double(ALgridhours[NowTillWhichALgrid-1]) )/
			       ( double(ALgridhours[NowTillWhichALgrid])-double(ALgridhours[NowTillWhichALgrid-1]) ) ) ;	
      }
    if ( (strcmp(LUgridParamNames[ct],"TF")==0) && (NowTillWhichTFgrid > 1) &&
	 ( NowTillWhichTFgrid < (numTFfiles+1) ) ) 
      {
	cNode->setThroughFall( cNode->getThroughFallInPrevGrid()+
			       ( cNode->getThroughFallInUntilGrid() - cNode->getThroughFallInPrevGrid() )*
			       ( timer->getCurrentTime() - double(TFgridhours[NowTillWhichTFgrid-1]) )/
				(double(TFgridhours[NowTillWhichTFgrid])-double(TFgridhours[NowTillWhichTFgrid-1])) );
      }
    if ( (strcmp(LUgridParamNames[ct],"VH")==0) && (NowTillWhichVHgrid > 1) &&
	 ( NowTillWhichVHgrid < (numVHfiles+1) ) ) 
      {
	cNode->setVegHeight( cNode->getVegHeightInPrevGrid()+
			     ( cNode->getVegHeightInUntilGrid() - cNode->getVegHeightInPrevGrid() )*
			     ( timer->getCurrentTime() - double(VHgridhours[NowTillWhichVHgrid-1]) )/
			      (double(VHgridhours[NowTillWhichVHgrid])-double(VHgridhours[NowTillWhichVHgrid-1])) );
	
      }
    if ( (strcmp(LUgridParamNames[ct],"SR")==0) && (NowTillWhichSRgrid > 1) &&
	 ( NowTillWhichSRgrid < (numSRfiles+1) ) ) 
      {
	cNode->setStomRes( cNode->getStomResInPrevGrid()+
			   ( cNode->getStomResInUntilGrid() - cNode->getStomResInPrevGrid() )*
			   ( timer->getCurrentTime() - double(SRgridhours[NowTillWhichSRgrid-1]) )/
			    (double(SRgridhours[NowTillWhichSRgrid])-double(SRgridhours[NowTillWhichSRgrid-1])) );
      }
    if ( (strcmp(LUgridParamNames[ct],"VF")==0) && (NowTillWhichVFgrid > 1) &&
	 ( NowTillWhichVFgrid < (numVFfiles+1) ) ) 
      {
	cNode->setVegFraction( cNode->getVegFractionInPrevGrid()+
			       (cNode->getVegFractionInUntilGrid() - cNode->getVegFractionInPrevGrid() )*
			       ( timer->getCurrentTime() - double(VFgridhours[NowTillWhichVFgrid-1]) )/
				( double(VFgridhours[NowTillWhichVFgrid]) - double(VFgridhours[NowTillWhichVFgrid-1]) ) );
      }
    if ( (strcmp(LUgridParamNames[ct],"CS")==0) && (NowTillWhichCSgrid > 1) &&
	 ( NowTillWhichCSgrid < (numCSfiles+1) ) ) 
      {
	cNode->setCanStorParam( cNode->getCanStorParamInPrevGrid()+
				( cNode->getCanStorParamInUntilGrid() - cNode->getCanStorParamInPrevGrid() )*
				( timer->getCurrentTime() - double(CSgridhours[NowTillWhichCSgrid-1]) )/
				 ( double(CSgridhours[NowTillWhichCSgrid])-double(CSgridhours[NowTillWhichCSgrid-1])) );
      }
    if ( (strcmp(LUgridParamNames[ct],"IC")==0) && (NowTillWhichICgrid > 1) &&
	 ( NowTillWhichICgrid < (numICfiles+1) ) ) 
      {
	cNode->setIntercepCoeff( cNode->getIntercepCoeffInPrevGrid()+
				 ( cNode->getIntercepCoeffInUntilGrid() - cNode->getIntercepCoeffInPrevGrid() )*
				 ( timer->getCurrentTime() - double(ICgridhours[NowTillWhichICgrid-1]) )/
				  (double(ICgridhours[NowTillWhichICgrid])-double(ICgridhours[NowTillWhichICgrid-1])) );
      }
    if ( (strcmp(LUgridParamNames[ct],"CC")==0) && (NowTillWhichCCgrid > 1) &&
	 ( NowTillWhichCCgrid < (numCCfiles+1) ) ) 
      {
	cNode->setCanFieldCap( cNode->getCanFieldCapInPrevGrid()+
			       ( cNode->getCanFieldCapInUntilGrid() - cNode->getCanFieldCapInPrevGrid() )*
			       ( timer->getCurrentTime() - double(CCgridhours[NowTillWhichCCgrid-1]) )/
				( double(CCgridhours[NowTillWhichCCgrid])-double(CCgridhours[NowTillWhichCCgrid-1])) );
      }
    if ( (strcmp(LUgridParamNames[ct],"DC")==0) && (NowTillWhichDCgrid > 1) &&
	 ( NowTillWhichDCgrid < (numDCfiles+1) ) ) 
      {
	cNode->setDrainCoeff( cNode->getDrainCoeffInPrevGrid()+
			      ( cNode->getDrainCoeffInUntilGrid() - cNode->getDrainCoeffInPrevGrid() )*
			      ( timer->getCurrentTime() - double(DCgridhours[NowTillWhichDCgrid-1]) )/
			       ( double(DCgridhours[NowTillWhichDCgrid])-double(DCgridhours[NowTillWhichDCgrid-1]) ) );
      }
    if ( (strcmp(LUgridParamNames[ct],"DE")==0) && (NowTillWhichDEgrid > 1) &&
	 ( NowTillWhichDEgrid < (numDEfiles+1) ) ) 
      {
	cNode->setDrainExpPar( cNode->getDrainExpParInPrevGrid()+
			       ( cNode->getDrainExpParInUntilGrid() - cNode->getDrainExpParInPrevGrid() )*
			       ( timer->getCurrentTime() - double(DEgridhours[NowTillWhichDEgrid-1]) )/
				(double(DEgridhours[NowTillWhichDEgrid])-double(DEgridhours[NowTillWhichDEgrid-1])) );
      }
    if ( (strcmp(LUgridParamNames[ct],"OT")==0) && (NowTillWhichOTgrid > 1) &&
	 ( NowTillWhichOTgrid < (numOTfiles+1) ) ) 
      {
	cNode->setOptTransmCoeff( cNode->getOptTransmCoeffInPrevGrid()+
				  ( cNode->getOptTransmCoeffInUntilGrid() - cNode->getOptTransmCoeffInPrevGrid() )*
				  ( timer->getCurrentTime() - double(OTgridhours[NowTillWhichOTgrid-1]) )/
				   (double(OTgridhours[NowTillWhichOTgrid])-double(OTgridhours[NowTillWhichOTgrid-1])) );
      }
    if ( (strcmp(LUgridParamNames[ct],"LA")==0) && (NowTillWhichLAgrid > 1) &&
	 ( NowTillWhichLAgrid < (numLAfiles+1) ) ) 
      {
	cNode->setLeafAI( cNode->getLeafAIInPrevGrid()+
			  ( cNode->getLeafAIInUntilGrid() - cNode->getLeafAIInPrevGrid() )*
			  ( timer->getCurrentTime() - double(LAgridhours[NowTillWhichLAgrid-1]) )/
			   (double(LAgridhours[NowTillWhichLAgrid])-double(LAgridhours[NowTillWhichLAgrid-1])) );
      }
  } // end for loop
  
  return;

}

/***************************************************************************
**
** constantLUGrids() Function
**
***************************************************************************/
void tEvapoTrans::constantLUGrids(tCNode* cNode)
{
    for (int ct=0;ct<nParmLU;ct++) {
        if ( (strcmp(LUgridParamNames[ct],"AL")==0))
        {
            cNode->setLandUseAlb( cNode->getLandUseAlbInPrevGrid()) ;
        }
        if ( (strcmp(LUgridParamNames[ct],"TF")==0))
        {
            cNode->setThroughFall(cNode->getThroughFallInPrevGrid());
        }
        if ( (strcmp(LUgridParamNames[ct],"VH")==0))
        {
            cNode->setVegHeight( cNode->getVegHeightInPrevGrid() );
        }
        if ( (strcmp(LUgridParamNames[ct],"SR")==0))
        {
            cNode->setStomRes( cNode->getStomResInPrevGrid());
        }
        if ( (strcmp(LUgridParamNames[ct],"VF")==0))
        {
            cNode->setVegFraction( cNode->getVegFractionInPrevGrid() );
        }
        if ( (strcmp(LUgridParamNames[ct],"CS")==0))
        {
            cNode->setCanStorParam( cNode->getCanStorParamInPrevGrid());
        }
        if ( (strcmp(LUgridParamNames[ct],"IC")==0))
        {
            cNode->setIntercepCoeff( cNode->getIntercepCoeffInPrevGrid());
        }
        if ( (strcmp(LUgridParamNames[ct],"CC")==0))
        {
            cNode->setCanFieldCap( cNode->getCanFieldCapInPrevGrid());
        }
        if ( (strcmp(LUgridParamNames[ct],"DC")==0) )
        {
            cNode->setDrainCoeff( cNode->getDrainCoeffInPrevGrid());
        }
        if ( (strcmp(LUgridParamNames[ct],"DE")==0) )
        {
            cNode->setDrainExpPar( cNode->getDrainExpParInPrevGrid() );
        }
        if ( (strcmp(LUgridParamNames[ct],"OT")==0))
        {
            cNode->setOptTransmCoeff( cNode->getOptTransmCoeffInPrevGrid() );
        }
        if ( (strcmp(LUgridParamNames[ct],"LA")==0))
        {
            cNode->setLeafAI( cNode->getLeafAIInPrevGrid() );
        }
    } // end for loop

    return;

}

/***************************************************************************
**
** integratedLUVars() Function
**
***************************************************************************/
void tEvapoTrans::integratedLUVars(tCNode* cNode, double te){

  for (int ct=0;ct<nParmLU;ct++) { 
    if (strcmp(LUgridParamNames[ct],"AL")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvLandUseAlb(cNode->getLandUseAlb());
      else if (te > 1.0) 
	cNode->setAvLandUseAlb((cNode->getAvLandUseAlb()*(te-1.0) + cNode->getLandUseAlb())/te);
    }
    if (strcmp(LUgridParamNames[ct],"TF")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvThroughFall(cNode->getThroughFall());
      else if (te > 1.0) 
	cNode->setAvThroughFall((cNode->getAvThroughFall()*(te-1.0) + cNode->getThroughFall())/te);
    }
    if (strcmp(LUgridParamNames[ct],"VH")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvVegHeight(cNode->getVegHeight());
      else if (te > 1.0) 
	cNode->setAvVegHeight((cNode->getAvVegHeight()*(te-1.0) + cNode->getVegHeight())/te);
    }
    if (strcmp(LUgridParamNames[ct],"SR")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvStomRes(cNode->getStomRes());
      else if (te > 1.0) 
	cNode->setAvStomRes((cNode->getAvStomRes()*(te-1.0) + cNode->getStomRes())/te);
    }
    if (strcmp(LUgridParamNames[ct],"VF")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvVegFraction(cNode->getVegFraction());
      else if (te > 1.0) 
	cNode->setAvVegFraction((cNode->getAvVegFraction()*(te-1.0) + cNode->getVegFraction())/te);
    }
    if (strcmp(LUgridParamNames[ct],"CS")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvCanStorParam(cNode->getCanStorParam());
      else if (te > 1.0) 
	cNode->setAvCanStorParam((cNode->getAvCanStorParam()*(te-1.0) + cNode->getCanStorParam())/te);
    }
    if (strcmp(LUgridParamNames[ct],"IC")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvIntercepCoeff(cNode->getIntercepCoeff());
      else if (te > 1.0) 
	cNode->setAvIntercepCoeff((cNode->getAvIntercepCoeff()*(te-1.0) + cNode->getIntercepCoeff())/te);
    }
    if (strcmp(LUgridParamNames[ct],"CC")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvCanFieldCap(cNode->getCanFieldCap());
      else if (te > 1.0) 
	cNode->setAvCanFieldCap((cNode->getAvCanFieldCap()*(te-1.0) + cNode->getCanFieldCap())/te);
    }
    if (strcmp(LUgridParamNames[ct],"DC")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvDrainCoeff(cNode->getDrainCoeff());
      else if (te > 1.0) 
	cNode->setAvDrainCoeff((cNode->getAvDrainCoeff()*(te-1.0) + cNode->getDrainCoeff())/te);
    }
    if (strcmp(LUgridParamNames[ct],"DE")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvDrainExpPar(cNode->getDrainExpPar());
      else if (te > 1.0) 
	cNode->setAvDrainExpPar((cNode->getAvDrainExpPar()*(te-1.0) + cNode->getDrainExpPar())/te);
    }
    if (strcmp(LUgridParamNames[ct],"OT")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvOptTransmCoeff(cNode->getOptTransmCoeff());
      else if (te > 1.0) 
	cNode->setAvOptTransmCoeff((cNode->getAvOptTransmCoeff()*(te-1.0) + cNode->getOptTransmCoeff())/te);
    }
    if (strcmp(LUgridParamNames[ct],"LA")==0) {  
      if (fabs(te - 1.0) < 1.0E-6)
	cNode->setAvLeafAI(cNode->getLeafAI());
      else if (te > 1.0) 
	cNode->setAvLeafAI((cNode->getAvLeafAI()*(te-1.0) + cNode->getLeafAI())/te);
    }
  }		
  
  return; 
}

/***************************************************************************
**
** SetGridTimeInfoVariables() Function
**
***************************************************************************/
void tEvapoTrans::SetGridTimeInfoVariables(tVariant *VariantLU, char *LUgridParamName)
{
	int checkhour, checkday, checkmonth, checkyear, checkminute;
	int currentTimeLU, GridHourCounter, numFilesCounter;
	int *tempgridhours;
	int DiffBetweenGridHours;

	checkhour = timer->hourS; 
	checkday = timer->dayS;
	checkmonth = timer->monthS;
	checkyear = timer->yearS;
	checkminute = 0;

	if (strcmp(LUgridParamName,"AL")==0) {numALfiles = 0;}
	else if (strcmp(LUgridParamName,"TF")==0) {numTFfiles = 0;}
	else if (strcmp(LUgridParamName,"VH")==0) {numVHfiles = 0;}
	else if (strcmp(LUgridParamName,"SR")==0) {numSRfiles = 0;}
	else if (strcmp(LUgridParamName,"VF")==0) {numVFfiles = 0;}
	else if (strcmp(LUgridParamName,"CS")==0) {numCSfiles = 0;}
	else if (strcmp(LUgridParamName,"IC")==0) {numICfiles = 0;}
	else if (strcmp(LUgridParamName,"CC")==0) {numCCfiles = 0;}
	else if (strcmp(LUgridParamName,"DC")==0) {numDCfiles = 0;}
	else if (strcmp(LUgridParamName,"DE")==0) {numDEfiles = 0;}
	else if (strcmp(LUgridParamName,"OT")==0) {numOTfiles = 0;}
	else if (strcmp(LUgridParamName,"LA")==0) {numLAfiles = 0;}

	numFilesCounter = 0;
	currentTimeLU = 0;

	while( (double)(currentTimeLU) <= timer->getEndTime() ) { 

		if ( VariantLU->infile2.is_open() )
			VariantLU->infile2.close();

		snprintf(VariantLU->fileIn,sizeof(VariantLU->fileIn),"%s%02d%02d%04d%02d.%s", VariantLU->getInputName(),
			checkmonth, checkday, checkyear, checkhour, VariantLU->getExtension() );//WR--09192023: 'sprintf' is deprecated: This function is provided for compatibility reasons only.

		VariantLU->infile2.open(VariantLU->fileIn);
	
		if ( VariantLU->infile2.is_open() ) {

			if (strcmp(LUgridParamName,"AL")==0) {numALfiles++;}
			else if (strcmp(LUgridParamName,"TF")==0) {numTFfiles++;}
			else if (strcmp(LUgridParamName,"VH")==0) {numVHfiles++;}
			else if (strcmp(LUgridParamName,"SR")==0) {numSRfiles++;}
			else if (strcmp(LUgridParamName,"VF")==0) {numVFfiles++;}
			else if (strcmp(LUgridParamName,"CS")==0) {numCSfiles++;}
			else if (strcmp(LUgridParamName,"IC")==0) {numICfiles++;}
			else if (strcmp(LUgridParamName,"CC")==0) {numCCfiles++;}
			else if (strcmp(LUgridParamName,"DC")==0) {numDCfiles++;}
			else if (strcmp(LUgridParamName,"DE")==0) {numDEfiles++;}
			else if (strcmp(LUgridParamName,"OT")==0) {numOTfiles++;}
			else if (strcmp(LUgridParamName,"LA")==0) {numLAfiles++;}

			numFilesCounter++;
		}

		timer->correctCalendarTime((double)(currentTimeLU), (double)(1), &checkminute, 
						&checkhour, &checkday, &checkmonth, &checkyear);

		currentTimeLU++;
	}

	if (numFilesCounter == 0) {
		Cout << "\tWhere are the varying Land Use grids"<<endl;
		Cout << "\tfor parameter " << LUgridParamName << " ?" << endl;
		Cout << "\nExiting Program..."<<endl<<endl;
		exit(1);
	}

	if (numFilesCounter == 1) {
		Cout << "\tThe dynamic grid option for parameter " << LUgridParamName <<endl;
		Cout << "\thas only one grid, assumed over the entire simulation time " << endl;
		Cout << "\nContinuing execution anyway..."<<endl<<endl;
	}

	checkhour = timer->hourS; 
	checkday = timer->dayS;
	checkmonth = timer->monthS;
	checkyear = timer->yearS;
	checkminute = 0;

	if (strcmp(LUgridParamName,"AL")==0) {
		ALgridhours = new int [numALfiles+1]; // The +1 means that the 0th value is not used
		ALgridFileNames = new char*[numALfiles+1];
		for (int ct=0;ct<numALfiles+1;ct++) {
			ALgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"TF")==0) {
		TFgridhours = new int [numTFfiles+1];
		TFgridFileNames = new char*[numTFfiles+1];
		for (int ct=0;ct<numTFfiles+1;ct++) {
			TFgridFileNames[ct]=new char[kName];
		}
	} 
	else if (strcmp(LUgridParamName,"VH")==0) {
		VHgridhours = new int [numVHfiles+1];
		VHgridFileNames = new char*[numVHfiles+1];
		for (int ct=0;ct<numVHfiles+1;ct++) {
			VHgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"SR")==0) {
		SRgridhours = new int [numSRfiles+1];
		SRgridFileNames = new char*[numSRfiles+1];
		for (int ct=0;ct<numSRfiles+1;ct++) {
			SRgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"VF")==0) {
		VFgridhours = new int [numVFfiles+1];
		VFgridFileNames = new char*[numVFfiles+1];
		for (int ct=0;ct<numVFfiles+1;ct++) {
			VFgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"CS")==0) {
		CSgridhours = new int [numCSfiles+1];
		CSgridFileNames = new char*[numCSfiles+1];
		for (int ct=0;ct<numCSfiles+1;ct++) {
			CSgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"IC")==0) {
		ICgridhours = new int [numICfiles+1];
		ICgridFileNames = new char*[numICfiles+1];
		for (int ct=0;ct<numICfiles+1;ct++) {
			ICgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"CC")==0) {
		CCgridhours = new int [numCCfiles+1];
		CCgridFileNames = new char*[numCCfiles+1];
		for (int ct=0;ct<numCCfiles+1;ct++) {
			CCgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"DC")==0) {
		DCgridhours = new int [numDCfiles+1];
		DCgridFileNames = new char*[numDCfiles+1];
		for (int ct=0;ct<numDCfiles+1;ct++) {
			DCgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"DE")==0) {
		DEgridhours = new int [numDEfiles+1];
		DEgridFileNames = new char*[numDEfiles+1];
		for (int ct=0;ct<numDEfiles+1;ct++) {
			DEgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"OT")==0) {
		OTgridhours = new int [numOTfiles+1];
		OTgridFileNames = new char*[numOTfiles+1];
		for (int ct=0;ct<numOTfiles+1;ct++) {
			OTgridFileNames[ct]=new char[kName];
		}
	}
	else if (strcmp(LUgridParamName,"LA")==0) {
		LAgridhours = new int [numLAfiles+1];
		LAgridFileNames = new char*[numLAfiles+1];
		for (int ct=0;ct<numLAfiles+1;ct++) {
			LAgridFileNames[ct]=new char[kName];
		}
	}

	tempgridhours = new int [numFilesCounter+1];
	tempgridhours[0]=0;

	currentTimeLU = 0;
	GridHourCounter = 0;

	while( (double)(currentTimeLU) <= timer->getEndTime() ) { 

		if (VariantLU->infile2.is_open() ) {
			VariantLU->infile2.close();
		}

		snprintf(VariantLU->fileIn,sizeof(VariantLU->fileIn),"%s%02d%02d%04d%02d.%s", VariantLU->getInputName(),
			checkmonth, checkday, checkyear, checkhour, VariantLU->getExtension() );//WR--09192023: 'sprintf' is deprecated: This function is provided for compatibility reasons only.

		VariantLU->infile2.open(VariantLU->fileIn);
	
		if (VariantLU->infile2.is_open()) {

			GridHourCounter++;

			if (strcmp(LUgridParamName,"AL")==0) {
				ALgridhours[GridHourCounter]=currentTimeLU;
				strcpy(ALgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"TF")==0) {
				TFgridhours[GridHourCounter]=currentTimeLU;
				strcpy(TFgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"VH")==0) {
				VHgridhours[GridHourCounter]=currentTimeLU;
				strcpy(VHgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"SR")==0) {
				SRgridhours[GridHourCounter]=currentTimeLU;
				strcpy(SRgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"VF")==0) {
				VFgridhours[GridHourCounter]=currentTimeLU;
				strcpy(VFgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"CS")==0) {
				CSgridhours[GridHourCounter]=currentTimeLU;
				strcpy(CSgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"IC")==0) {
				ICgridhours[GridHourCounter]=currentTimeLU;
				strcpy(ICgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"CC")==0) {
				CCgridhours[GridHourCounter]=currentTimeLU;
				strcpy(CCgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"DC")==0) {
				DCgridhours[GridHourCounter]=currentTimeLU;
				strcpy(DCgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"DE")==0) {
				DEgridhours[GridHourCounter]=currentTimeLU;
				strcpy(DEgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"OT")==0) {
				OTgridhours[GridHourCounter]=currentTimeLU;
				strcpy(OTgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			else if (strcmp(LUgridParamName,"LA")==0) {
				LAgridhours[GridHourCounter]=currentTimeLU;
				strcpy(LAgridFileNames[GridHourCounter],VariantLU->fileIn);
			}
			
			tempgridhours[GridHourCounter]=currentTimeLU;

		}

		timer->correctCalendarTime((double)(currentTimeLU), (double)(1), &checkminute, 
						&checkhour, &checkday, &checkmonth, &checkyear);

		currentTimeLU++;
	}

	for (int ii=1; ii<=numFilesCounter;ii++) {
		DiffBetweenGridHours=tempgridhours[ii]-tempgridhours[ii-1];
		
		if (DiffBetweenGridHours > 2928) { //2928 hours is approximately four months
			Cout << "\tYou have a large time interval of " << DiffBetweenGridHours <<"hours"<<endl;
			Cout << "\tbetween grid " << ii << "and its preceding grid"<<endl;
			Cout << "\nContinuing execution anyway..."<<endl<<endl;	
		}
	}
	
	delete [] tempgridhours;

	if (strcmp(LUgridParamName,"AL")==0) {NowTillWhichALgrid = 1;}
	else if (strcmp(LUgridParamName,"TF")==0) {NowTillWhichTFgrid = 1;}
	else if (strcmp(LUgridParamName,"VH")==0) {NowTillWhichVHgrid = 1;}
	else if (strcmp(LUgridParamName,"SR")==0) {NowTillWhichSRgrid = 1;}
	else if (strcmp(LUgridParamName,"VF")==0) {NowTillWhichVFgrid = 1;}
	else if (strcmp(LUgridParamName,"CS")==0) {NowTillWhichCSgrid = 1;}
	else if (strcmp(LUgridParamName,"IC")==0) {NowTillWhichICgrid = 1;}
	else if (strcmp(LUgridParamName,"CC")==0) {NowTillWhichCCgrid = 1;}
	else if (strcmp(LUgridParamName,"DC")==0) {NowTillWhichDCgrid = 1;}
	else if (strcmp(LUgridParamName,"DE")==0) {NowTillWhichDEgrid = 1;}
	else if (strcmp(LUgridParamName,"OT")==0) {NowTillWhichOTgrid = 1;}
	else if (strcmp(LUgridParamName,"LA")==0) {NowTillWhichLAgrid = 1;}

	return;

}

/***************************************************************************
**
** tEvapoTrans::deleteLUGrids() Function
**
**
***************************************************************************/

void tEvapoTrans::deleteLUGrids() 
{
  for (int ct=0;ct<nParmLU;ct++) { 
	if (strcmp(LUgridParamNames[ct],"AL")==0) {
	  delete LandUseAlbGrid;
	  delete [] ALgridhours;
	  for (int sz=0;sz<numALfiles+1;sz++) {
	    delete [] ALgridFileNames[sz];
	  }
	  delete [] ALgridFileNames;		
	}
	if (strcmp(LUgridParamNames[ct],"TF")==0) {
	  delete ThroughFallGrid;
	  delete [] TFgridhours; 
	  for (int sz=0;sz<numTFfiles+1;sz++) {
	    delete [] TFgridFileNames[sz];
	  }
	  delete [] TFgridFileNames;
	}
	if (strcmp(LUgridParamNames[ct],"VH")==0) {
	  delete VegHeightGrid;
	  delete [] VHgridhours;
	  for (int sz=0;sz<numVHfiles+1;sz++) {
	    delete [] VHgridFileNames[sz];
	  }
	  delete [] VHgridFileNames;
	}
	if (strcmp(LUgridParamNames[ct],"SR")==0) {
	  delete StomResGrid;
	  delete [] SRgridhours;
	  for (int sz=0;sz<numSRfiles+1;sz++) {
	    delete [] SRgridFileNames[sz];
	  }
	  delete [] SRgridFileNames;
	}
	if (strcmp(LUgridParamNames[ct],"VF")==0) {
	  delete VegFractGrid;
	  delete [] VFgridhours;
	  for (int sz=0;sz<numVFfiles+1;sz++) {
	    delete [] VFgridFileNames[sz];
	  }
	  delete [] VFgridFileNames;
	}
	if (strcmp(LUgridParamNames[ct],"CS")==0) {
	  delete CanStorParamGrid;
	  delete [] CSgridhours;
	  for (int sz=0;sz<numCSfiles+1;sz++) {
	    delete [] CSgridFileNames[sz];
	  }
	  delete [] CSgridFileNames;
	}
	if (strcmp(LUgridParamNames[ct],"IC")==0) {
	  delete IntercepCoeffGrid;
	  delete [] ICgridhours;
	  for (int sz=0;sz<numICfiles+1;sz++) {
	    delete [] ICgridFileNames[sz];
	  }
	  delete [] ICgridFileNames;					
	}
	if (strcmp(LUgridParamNames[ct],"CC")==0) {
	  delete CanFieldCapGrid;
	  delete [] CCgridhours;
	  for (int sz=0;sz<numCCfiles+1;sz++) {
	    delete [] CCgridFileNames[sz];
	  }
	  delete [] CCgridFileNames;
	}
	if (strcmp(LUgridParamNames[ct],"DC")==0) {
	  delete DrainCoeffGrid;
	  delete [] DCgridhours;
	  for (int sz=0;sz<numDCfiles+1;sz++) {
	    delete [] DCgridFileNames[sz];
	  }
	  delete [] DCgridFileNames;
	}
	if (strcmp(LUgridParamNames[ct],"DE")==0) {
	  delete DrainExpParGrid;
	  delete [] DEgridhours;
	  for (int sz=0;sz<numDEfiles+1;sz++) {
	    delete [] DEgridFileNames[sz];
	  }
	  delete [] DEgridFileNames;
	}
	if (strcmp(LUgridParamNames[ct],"OT")==0) {
	  delete OptTransmCoeffGrid;
	  delete [] OTgridhours;
	  for (int sz=0;sz<numOTfiles+1;sz++) {
	    delete [] OTgridFileNames[sz];
	  }
	  delete [] OTgridFileNames;
	}
	if (strcmp(LUgridParamNames[ct],"LA")==0) {
	  delete LeafAIGrid;
	  delete [] LAgridhours;
	  for (int sz=0;sz<numLAfiles+1;sz++) {
	    delete [] LAgridFileNames[sz];
	  }
	  delete [] LAgridFileNames;					
	}
  }
  
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


//=========================================================================
//
//
//                  Section 8: Debug Routines
//
//
//=========================================================================
void tEvapoTrans::Debug(int time, int flag) 
{
	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	
	if (flag == 1) {
		cNode = nodeIter.FirstP();
		while ( nodeIter.IsActive() ) {
			//Print out Statements Here
			cNode = nodeIter.NextP();
		}
	}
	}

/***************************************************************************
**
** tEvapoTrans::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/
void tEvapoTrans::writeRestart(fstream & rStr) const
{ 
  BinaryWrite(rStr, VerbID);
  BinaryWrite(rStr, vapOption);
  BinaryWrite(rStr, tsOption);
  BinaryWrite(rStr, nrOption);
  BinaryWrite(rStr, Rah);
  BinaryWrite(rStr, Rstm);
  BinaryWrite(rStr, SoilHeatCondTh);
  BinaryWrite(rStr, SoilHeatCpctTh);
  BinaryWrite(rStr, SoilHeatDiffTh);
  BinaryWrite(rStr, Tlinke);
  BinaryWrite(rStr, Is);
  BinaryWrite(rStr, Ic);
  BinaryWrite(rStr, Ics);
  BinaryWrite(rStr, Id);
  BinaryWrite(rStr, Ids);
  BinaryWrite(rStr, vPressC);
  BinaryWrite(rStr, Epot);;

  BinaryWrite(rStr, numStations);
  BinaryWrite(rStr, arraySize);
  BinaryWrite(rStr, hourlyTimeStep);
  BinaryWrite(rStr, nParm);
  BinaryWrite(rStr, gridgmt);
  BinaryWrite(rStr, metdataOption);
  BinaryWrite(rStr, Ioption);
  BinaryWrite(rStr, gFluxOption);
  BinaryWrite(rStr, dewHumFlag);
  BinaryWrite(rStr, ID);
  BinaryWrite(rStr, gmt);
  BinaryWrite(rStr, nodeHour);
  BinaryWrite(rStr, thisStation);
  BinaryWrite(rStr, oldTimeStep);
  for (int i = 0; i < arraySize; i++)
    BinaryWrite(rStr, assignedStation[i]);

  BinaryWrite(rStr, timeStep);
  BinaryWrite(rStr, timeCount);
  BinaryWrite(rStr, gridlat);
  BinaryWrite(rStr, gridlong);
  BinaryWrite(rStr, coeffH);
  BinaryWrite(rStr, coeffKt);
  BinaryWrite(rStr, coeffAl);
  BinaryWrite(rStr, coeffRs);
  BinaryWrite(rStr, coeffV);
  BinaryWrite(rStr, coeffKs);
  BinaryWrite(rStr, coeffCs);
  BinaryWrite(rStr, coeffPan);
  BinaryWrite(rStr, potEvap);
  BinaryWrite(rStr, actEvap);
  BinaryWrite(rStr, panEvap);
  BinaryWrite(rStr, betaS);
  BinaryWrite(rStr, betaT);
  BinaryWrite(rStr, airTemp);
  BinaryWrite(rStr, dewTemp);
  BinaryWrite(rStr, surfTemp);
  BinaryWrite(rStr, Tso);
  BinaryWrite(rStr, Tlo);
  BinaryWrite(rStr, rHumidity);
  BinaryWrite(rStr, atmPress);
  BinaryWrite(rStr, windSpeed);
  BinaryWrite(rStr, skyCover);
  BinaryWrite(rStr, netRad);
  BinaryWrite(rStr, latitude);
  BinaryWrite(rStr, longitude);
  BinaryWrite(rStr, vPress);
  BinaryWrite(rStr, inLongR);
  BinaryWrite(rStr, inShortR);
  BinaryWrite(rStr, outLongR);
  BinaryWrite(rStr, elevation);
  BinaryWrite(rStr, slope);
  BinaryWrite(rStr, aspect);
  BinaryWrite(rStr, atmPressC);
  BinaryWrite(rStr, surfTempC);
  BinaryWrite(rStr, skyCoverC);
  BinaryWrite(rStr, rHumidityC);
  BinaryWrite(rStr, dewTempC);
  BinaryWrite(rStr, windSpeedC);
  BinaryWrite(rStr, netRadC);
  BinaryWrite(rStr, gFlux);
  BinaryWrite(rStr, hFlux);
  BinaryWrite(rStr, lFlux);
  BinaryWrite(rStr, rain);
  BinaryWrite(rStr, Gso);
  BinaryWrite(rStr, Io);
  BinaryWrite(rStr, alphaD);
  BinaryWrite(rStr, sinAlpha);
  BinaryWrite(rStr, del);
  BinaryWrite(rStr, phi);
  BinaryWrite(rStr, tau);
  BinaryWrite(rStr, circ);
  BinaryWrite(rStr, sunaz);
  BinaryWrite(rStr, SunRisHrLoc);
  BinaryWrite(rStr, SunSetHrLoc);
  BinaryWrite(rStr, DayLength);
  BinaryWrite(rStr, deltaT);
  BinaryWrite(rStr, RadDirObs);
  BinaryWrite(rStr, RadDifObs);

  BinaryWrite(rStr, snowOption); // Snow and Shelter
  BinaryWrite(rStr, shelterOption);
  BinaryWrite(rStr, LUgridgmt);
  BinaryWrite(rStr, luOption);
  BinaryWrite(rStr, nParmLU);
  BinaryWrite(rStr, luInterpOption);
  BinaryWrite(rStr, IfNotFirstTStepLU);
  BinaryWrite(rStr, metHour); 
  BinaryWrite(rStr, etHour); 
  BinaryWrite(rStr, rainInt);
  BinaryWrite(rStr, LUgridlat); 
  BinaryWrite(rStr, LUgridlong);
  BinaryWrite(rStr, coeffLAI);
  BinaryWrite(rStr, shelterFactorGlobal); 
  BinaryWrite(rStr, landRefGlobal);
  BinaryWrite(rStr, horizonAngle);
  BinaryWrite(rStr, ha0000); 
  BinaryWrite(rStr, ha0225); 
  BinaryWrite(rStr, ha0450); 
  BinaryWrite(rStr, ha0675); 
  BinaryWrite(rStr, ha0900); 
  BinaryWrite(rStr, ha1125);
  BinaryWrite(rStr, ha1350); 
  BinaryWrite(rStr, ha1575); 
  BinaryWrite(rStr, ha1800); 
  BinaryWrite(rStr, ha2025); 
  BinaryWrite(rStr, ha2250);
  BinaryWrite(rStr, ha2475);
  BinaryWrite(rStr, ha2700);
  BinaryWrite(rStr, ha2925);
  BinaryWrite(rStr, ha3150);
  BinaryWrite(rStr, ha3375);
  BinaryWrite(rStr, tempLapseRate);
  BinaryWrite(rStr, SunHour);
  BinaryWrite(rStr, AtFirstTimeStepLUFlag);
    
    // to get right time vegetation parameters after reading restart files. Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichALgrid); 
    BinaryWrite(rStr, NowTillWhichTFgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichVHgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichSRgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichVFgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichCSgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichICgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichCCgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichDCgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichDEgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichOTgrid); // Ara Ko 2017
    BinaryWrite(rStr, NowTillWhichLAgrid); // Ara Ko 2017

  if (evapotransOption != 0) {
    for (int i = 0; i < 3; i++) 
       BinaryWrite(rStr, currentTime[i]);
    if (rainPtr->getoptStorm())
      weatherSimul->writeRestart(rStr);
    if (metdataOption == 1)
      for (int i = 0; i < numStations; i++)
        weatherStations[i].writeRestart(rStr);
  }
}


/***************************************************************************
**
** tEvapoTrans::readRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/

void tEvapoTrans::readRestart(fstream & rStr)
{
  BinaryRead(rStr, VerbID);
  BinaryRead(rStr, vapOption);
  BinaryRead(rStr, tsOption);
  BinaryRead(rStr, nrOption);
  BinaryRead(rStr, Rah);
  BinaryRead(rStr, Rstm);
  BinaryRead(rStr, SoilHeatCondTh);
  BinaryRead(rStr, SoilHeatCpctTh);
  BinaryRead(rStr, SoilHeatDiffTh);
  BinaryRead(rStr, Tlinke);
  BinaryRead(rStr, Is);
  BinaryRead(rStr, Ic);
  BinaryRead(rStr, Ics);
  BinaryRead(rStr, Id);
  BinaryRead(rStr, Ids);
  BinaryRead(rStr, vPressC);
  BinaryRead(rStr, Epot);

  BinaryRead(rStr, numStations);
  BinaryRead(rStr, arraySize);
  BinaryRead(rStr, hourlyTimeStep);
  BinaryRead(rStr, nParm);
  BinaryRead(rStr, gridgmt);
  BinaryRead(rStr, metdataOption);
  BinaryRead(rStr, Ioption);
  BinaryRead(rStr, gFluxOption);
  BinaryRead(rStr, dewHumFlag);
  BinaryRead(rStr, ID);
  BinaryRead(rStr, gmt);
  BinaryRead(rStr, nodeHour);
  BinaryRead(rStr, thisStation);
  BinaryRead(rStr, oldTimeStep);
  for (int i = 0; i < arraySize; i++)
    BinaryRead(rStr, assignedStation[i]);

  BinaryRead(rStr, timeStep);
  BinaryRead(rStr, timeCount);
  BinaryRead(rStr, gridlat);
  BinaryRead(rStr, gridlong);
  BinaryRead(rStr, coeffH);
  BinaryRead(rStr, coeffKt);
  BinaryRead(rStr, coeffAl);
  BinaryRead(rStr, coeffRs);
  BinaryRead(rStr, coeffV);
  BinaryRead(rStr, coeffKs);
  BinaryRead(rStr, coeffCs);
  BinaryRead(rStr, coeffPan);
  BinaryRead(rStr, potEvap);
  BinaryRead(rStr, actEvap);
  BinaryRead(rStr, panEvap);
  BinaryRead(rStr, betaS);
  BinaryRead(rStr, betaT);
  BinaryRead(rStr, airTemp);
  BinaryRead(rStr, dewTemp);
  BinaryRead(rStr, surfTemp);
  BinaryRead(rStr, Tso);
  BinaryRead(rStr, Tlo);
  BinaryRead(rStr, rHumidity);
  BinaryRead(rStr, atmPress);
  BinaryRead(rStr, windSpeed);
  BinaryRead(rStr, skyCover);
  BinaryRead(rStr, netRad);
  BinaryRead(rStr, latitude);
  BinaryRead(rStr, longitude);
  BinaryRead(rStr, vPress);
  BinaryRead(rStr, inLongR);
  BinaryRead(rStr, inShortR);
  BinaryRead(rStr, outLongR);
  BinaryRead(rStr, elevation);
  BinaryRead(rStr, slope);
  BinaryRead(rStr, aspect);
  BinaryRead(rStr, atmPressC);
  BinaryRead(rStr, surfTempC);
  BinaryRead(rStr, skyCoverC);
  BinaryRead(rStr, rHumidityC);
  BinaryRead(rStr, dewTempC);
  BinaryRead(rStr, windSpeedC);
  BinaryRead(rStr, netRadC);
  BinaryRead(rStr, gFlux);
  BinaryRead(rStr, hFlux);
  BinaryRead(rStr, lFlux);
  BinaryRead(rStr, rain);
  BinaryRead(rStr, Gso);
  BinaryRead(rStr, Io);
  BinaryRead(rStr, alphaD);
  BinaryRead(rStr, sinAlpha);
  BinaryRead(rStr, del);
  BinaryRead(rStr, phi);
  BinaryRead(rStr, tau);
  BinaryRead(rStr, circ);
  BinaryRead(rStr, sunaz);
  BinaryRead(rStr, SunRisHrLoc);
  BinaryRead(rStr, SunSetHrLoc);
  BinaryRead(rStr, DayLength);
  BinaryRead(rStr, deltaT);
  BinaryRead(rStr, RadDirObs);
  BinaryRead(rStr, RadDifObs);

  BinaryRead(rStr, snowOption); // Snow and Shelter
  BinaryRead(rStr, shelterOption);
  BinaryRead(rStr, LUgridgmt);
  BinaryRead(rStr, luOption);
  BinaryRead(rStr, nParmLU);
  BinaryRead(rStr, luInterpOption);
  BinaryRead(rStr, IfNotFirstTStepLU);
  BinaryRead(rStr, metHour);
  BinaryRead(rStr, etHour);
  BinaryRead(rStr, rainInt);
  BinaryRead(rStr, LUgridlat);
  BinaryRead(rStr, LUgridlong);
  BinaryRead(rStr, coeffLAI);
  BinaryRead(rStr, shelterFactorGlobal);
  BinaryRead(rStr, landRefGlobal);
  BinaryRead(rStr, horizonAngle);
  BinaryRead(rStr, ha0000);
  BinaryRead(rStr, ha0225);
  BinaryRead(rStr, ha0450);
  BinaryRead(rStr, ha0675);
  BinaryRead(rStr, ha0900);
  BinaryRead(rStr, ha1125);
  BinaryRead(rStr, ha1350);
  BinaryRead(rStr, ha1575);
  BinaryRead(rStr, ha1800);
  BinaryRead(rStr, ha2025);
  BinaryRead(rStr, ha2250);
  BinaryRead(rStr, ha2475);
  BinaryRead(rStr, ha2700);
  BinaryRead(rStr, ha2925);
  BinaryRead(rStr, ha3150);
  BinaryRead(rStr, ha3375);
  BinaryRead(rStr, tempLapseRate);
  BinaryRead(rStr, SunHour);
  BinaryRead(rStr, AtFirstTimeStepLUFlag);
  BinaryRead(rStr, NowTillWhichALgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichTFgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichVHgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichSRgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichVFgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichCSgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichICgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichCCgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichDCgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichDEgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichOTgrid); // Ara Ko 2017
  BinaryRead(rStr, NowTillWhichLAgrid); // Ara Ko 2017


  if (evapotransOption != 0) {
    for (int i = 0; i < 3; i++) 
       BinaryRead(rStr, currentTime[i]);
    if (rainPtr->getoptStorm())
      weatherSimul->readRestart(rStr);
    if (metdataOption == 1)
      for (int i = 0; i < numStations; i++)
        weatherStations[i].readRestart(rStr);
  }
}


//=========================================================================
//
//
//			End of tEvapoTrans.cpp                   
//
//
//=========================================================================

