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
**  tRainfall.cpp:   Functions for tRainfall classes (see tRainfall.h)
**
***************************************************************************/

#include "src/tRasTin/tRainfall.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//                  Section 1: tRainfall Constructors and Destructors
//
//
//=========================================================================

// Default Constructor
tRainfall::tRainfall() 
{
    gridPtr = 0;
}

// Constructor
tRainfall::tRainfall(SimulationControl *simCtrPtr, tMesh<tCNode> *gridRef, 
					 tInputFile &inFile, tResample *resamp) 
: tStorm( simCtrPtr, inFile )
{
	gridPtr = gridRef;
	respPtr = resamp;
	simCtrl = simCtrPtr;
	
	SetRainVariables( inFile );

}

/***************************************************************************
**
**  tRainfall::SetRainVariables(tInputFile &inFile)
**
**  Initializes tRainfall object
**
***************************************************************************/
void tRainfall::SetRainVariables(tInputFile &inFile)
{ 
	rainfallType = inFile.ReadItem(rainfallType, "RAINSOURCE");
	rainDt = inFile.ReadItem(rainDt, "RAININTRVL");
	fState = 0;     //Default values
	optMAP = 0;  
	climate = 0.0; 
	optForecast = 0;  
	
	// If data are used as the model forcing
	if (!getoptStorm()) {
		
		if (rainfallType == 1)
			Cout<<"Rainfall Source: \t\tNEXRAD StageIII/P1 Product [cm/hr]\n";
		else if (rainfallType == 2)
			Cout<<"Rainfall Source: \t\tWSI Product [mm/hr]\n";
		
		if (rainfallType == 1 || rainfallType == 2) { 
			inFile.ReadItem(inputname, "RAINFILE");
			inFile.ReadItem(extension, "RAINEXTENSION");
			
			// To make compatible with existing model setups
			if (inFile.IsItemIn( "RAINSEARCH" ))
				searchRain = inFile.ReadItem(searchRain, "RAINSEARCH");
			else 
				searchRain = 24; //Default value
			optMAP = inFile.ReadItem(optMAP, "RAINDISTRIBUTION");
			Cout<<"Rainfall Input Path: \t\t'"<<inputname<<"'"<<endl;
			Cout<<"Rainfall File Extension: \t"<<extension<<endl;
			NewRain();  
		}
		else if (rainfallType == 3) {
			Cout<<"Rainfall Source: \t\tRaingauge Stations"<<endl;
			inFile.ReadItem(stationFile, "GAUGESTATIONS");

			// SKY2008Snow from AJR2007
			if (inFile.IsItemIn("PRECLAPSE"))
				precLapseRate = inFile.ReadItem(precLapseRate, "PRECLAPSE"); //AJR @ NMT 2007
			else
				precLapseRate = 0.0;

			readGaugeStat(stationFile);
			for (int ct=0;ct<numStations;ct++) {
				readGaugeData(ct);
			}

			assignStationToNode();
			InitializeGauge();

		}
		else {
			Cout<<"\nRainfall Source Option " << rainfallType;
			Cout<<" not valid." <<endl;
			Cout<<"\tPlease use: "<<endl;
			Cout<<"\t\t(1) NEXRAD Stage III Radar"<<endl;
			Cout<<"\t\t(2) WSI Radar Rainfall"<<endl;
			Cout<<"\t\t(3) Rain Gauge Station Rainfall"<<endl;
			Cout << "Exiting Program...\n\n"<<endl;
			exit(1);
		}
	}
	
	// Get Forecast File Directory, use same extension
	optForecast = inFile.ReadItem(optForecast, "FORECASTMODE"); 

	if (optForecast == 1) {
		inFile.ReadItem(forecastname, "FORECASTFILE");
	}
	else if (optForecast == 3) {
		climate = inFile.ReadItem(climate, "CLIMATOLOGY");
	}
	
	// Re-initialize average, cumulative MAP and # rainfall files
	numRains = 0;
	aveMAP = 0.0;
	cumMAP = 0.0;
}

// Destructor
tRainfall::~tRainfall() 
{
	gridPtr = nullptr;
	respPtr = nullptr;
	simCtrl = nullptr;
	
	if (rainfallType == 1 || rainfallType == 2) {
#ifdef ALPHA_64
		if ( infile )
			infile.close();
#elif defined LINUX_32
		if ( infile.is_open() )
			infile.close();
#elif defined WIN
		if ( infile )
			infile.close();
#else 
		if ( infile )
			infile.close();
#endif
	}
	if (rainfallType == 3) {
		delete [] currentTime; 
		delete [] latitude;
		delete [] longitude; 
		delete [] gaugeRain;
		delete [] rainGauges;
		// delete [] assignedRain;  -- GMnSKY2008MLE
	}
	Cout<<"tRainfall Object has been destroyed..."<<endl<<flush;
}

//=========================================================================
//
//
//                  Section 2: tRainfall Functions
//
//
//=========================================================================

/***************************************************************************
**
**  tRainfall::Compose_In_Mrain_Name(tRunTime *t)
**
**  Composes filename to rainfall file and opens it.
**
***************************************************************************/
int tRainfall::Compose_In_Mrain_Name(tRunTimer *t) 
{ 
#ifdef ALPHA_64
    if ( infile )
		infile.close();
#elif defined LINUX_32
    if ( infile.is_open() )
		infile.close();
#elif defined WIN
    if ( infile )
		infile.close();
#else 
    if ( infile )
		infile.close();
#endif
	
	// Read rainfall file depending on forecast state
	if (fState == 0) {
		
		if (t->minuteRn || t->dtRain < 1) //If 'minute' is NOT equal to '0'  
			snprintf(mrainfileIn,sizeof(mrainfileIn), "%s%02d%02d%04d%02d%02d.%s", inputname,
					t->monthRn, t->dayRn, t->yearRn, t->hourRn, t->minuteRn, extension);//WR--09192023:'sprintf' is deprecated: This function is provided for compatibility reasons only.
		else {           //If 'minute' IS equal to '0'
			snprintf(mrainfileIn,sizeof(mrainfileIn),"%s%02d%02d%04d%02d.%s", inputname,
					t->monthRn, t->dayRn, t->yearRn, t->hourRn, extension);//WR--09192023:'sprintf' is deprecated: This function is provided for compatibility reasons only.
		}
		infile.open(mrainfileIn);
	}
	
	else if (fState == 1) {
		
		if (t->getoptForecast() == 1) {
			if (t->minuteRn || t->dtRain < 1) //If 'minute' is NOT equal to '0'  
				snprintf(mrainfileIn,sizeof(mrainfileIn), "%s%02d%02d%04d%02d%02d.%s", forecastname,
						t->monthRn, t->dayRn, t->yearRn, t->hourRn, t->minuteRn, extension); //WR--09192023:'sprintf' is deprecated: This function is provided for compatibility reasons only.
			else {           //If 'minute' IS equal to '0'
				snprintf(mrainfileIn,sizeof(mrainfileIn),"%s%02d%02d%04d%02d.%s", forecastname,
						t->monthRn, t->dayRn, t->yearRn, t->hourRn, extension);//WR--09192023:'sprintf' is deprecated: This function is provided for compatibility reasons only.
			}
			infile.open(mrainfileIn);
		}
		else if (t->getoptForecast() == 2) {   //Persistence
			infile.open(mrainfileIn);
		}
		else if (t->getoptForecast() == 3) {  //Climatology
			return 1;
		}
	}
	
	else if (fState == 2)
		return 1;
	
	// Check if file opened
#ifdef ALPHA_64
    if ( !infile )
		return 0;
#elif defined LINUX_32
    if ( !(infile.is_open()) )
		return 0;
#elif defined WIN
    if ( !infile )
		return 0;
#else 
    if ( !infile )
		return 0;
#endif
	else {
		infile.close();
		return 1;
	}
}

/***************************************************************************
**
**  tRainfall::NewRain()
**
**  Initializes object to zero
**
***************************************************************************/
void tRainfall::NewRain() 
{
	int id;
	tCNode * cn;
	tMeshListIter<tCNode> nodeIter( gridPtr->getNodeList() );
	
	id = 0;
	cn = nodeIter.FirstP();
	while( nodeIter.IsActive() ) {
		cn->setRain( 0.0 );
		cn = nodeIter.NextP();
		id++;
	}
	return;
}

/***************************************************************************
**
**  tRainfall::NewRain(double)
**
**  Initializes object to uniform rate
**
***************************************************************************/
void tRainfall::NewRain(double rain) 
{
	int id;
	tCNode * cn;
	tMeshListIter<tCNode> nodeIter( gridPtr->getNodeList() );
	
	id = 0;
	cn = nodeIter.FirstP();
	while( nodeIter.IsActive() ) {
		cn->setRain( rain );
		cn = nodeIter.NextP();
		id++;
	}
	return;
}

/***************************************************************************
**
**  tRainfall::NewRain(tRunTime *t)
**
**  Reads in new rain corresponding to a current time tag, uses tResample
**  and assigns values to the Voronoi nodes
**
***************************************************************************/
void tRainfall::NewRain(tRunTimer *t) 
{
	int id;
	tCNode * cn;
	tMeshListIter<tCNode> nodeIter( gridPtr->getNodeList() );
	double sumRain = 0.0;
	double sumArea = 0.0;
	double maxRain = 250.0;  //Maximum valid rainfall (mm/hr)
	
	id = 0;
	cn = nodeIter.FirstP();
	
	// FSTATE == 0 OR 1
	if (fState == 0 || fState == 1) {
		
		curRain = respPtr->doIt(mrainfileIn, 1);
		
		// Check Valid Rainfall and Compute MAP 
		while( nodeIter.IsActive() ) {
			if (curRain[id] < 0.0 || curRain[id] > maxRain*t->getRainDT())
				curRain[id] = 0.0;
			sumRain = sumRain + cn->getVArea()*curRain[id];   
			sumArea = sumArea + cn->getVArea();
			cn = nodeIter.NextP();
			id++; 
		}
		
		// Compute cum and ave MAPs over time prior to forecast
		// conditioned on rain occuring (sumRain > 0)
		if (fState == 0 && sumRain > 0) {
			numRains++;
			cumMAP += sumRain/sumArea;
			aveMAP = cumMAP/numRains;
		}
		
		// Assign Rainfall Values
		id = 0;
		cn = nodeIter.FirstP();
		while( nodeIter.IsActive() ) { 
			// Assign MAP to curRain for optMAP = 1
			if (optMAP == 1)
				curRain[id]=sumRain/sumArea;
			
			if (rainfallType == 1)
				cn->setRain( curRain[id]*10.0/t->getRainDT()); //NEXRAD - cm/hour
			else if (rainfallType == 2)
				cn->setRain( curRain[id]/t->getRainDT());      //WSI - mm/hour
			cn = nodeIter.NextP();
			id++;
		}
		
		// Climatological forecast for FState == 1
		cn = nodeIter.FirstP();
		if (fState == 1 && t->getoptForecast() == 3) {
			while( nodeIter.IsActive() ) {
				cn->setRain( climate ); 
				cn = nodeIter.NextP();
			}
		}
	}
	
	// FSTATE == 2
	else if (fState == 2) {
		if (t->getoptForecast() != 3) {
			while( nodeIter.IsActive() ) {    
				if (rainfallType == 1)
					cn->setRain( aveMAP*10.0 ); //NEXRAD - cm/hour
				else if (rainfallType == 2)
					cn->setRain( aveMAP );      //WSI - mm/hour
				cn = nodeIter.NextP();
			}
		}
		else {    //Climatological forecast
			while( nodeIter.IsActive() ) {    
				cn->setRain( climate ); 
				cn = nodeIter.NextP();
			}
		}
	}

	return;
}

//=========================================================================
//
//
//                  Section 2: tRainfall Functions for Gauges
//
//
//=========================================================================


/***************************************************************************
**
** tRainfall::InitializeGauge() Function
**
** Initialize Variables
**
***************************************************************************/
void tRainfall::InitializeGauge() 
{
	arraySize = gridPtr->getNodeList()->getActiveSize();
	hourlyTimeStep = 0;
	
	gaugeRain = new double[arraySize];
	latitude  = new double[arraySize];
	longitude = new double[arraySize];
	// assignedRain = new int[arraySize]; -- GMnSKY2008MLE
	
	for (int ct=0;ct<arraySize;ct++) {
		gaugeRain[ct] = 0.0;
		latitude[ct]  = 0.0;
		longitude[ct] = 0.0;
	}
	
	currentTime = new int[4];
	for (int count=0;count<4;count++) {
		currentTime[count] = 0;
	}
	return;
}

/***************************************************************************
**
** tRainfall::callRainGauge() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/
void tRainfall::callRainGauge(tRunTimer *t) 
{
	NewRainData(hourlyTimeStep);


	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	cNode = nodeIter.FirstP();
	while(nodeIter.IsActive()){

		nodeIter.NextP();

	}

	setToNode();

	hourlyTimeStep++;
	return;
}

/***************************************************************************
**
** tRainfall::NewRainData() Function
**
** Assigns the values of the current raingauge values to the 
** nodes based on the results of the tResample Thiessen polygon routine.
** Assign the values from tRainGauge to tCNode.
**
***************************************************************************/
void tRainfall::NewRainData(int time) 
{  
	currentTime[0] = rainGauges[0].getYear(time);
	currentTime[1] = rainGauges[0].getMonth(time);
	currentTime[2] = rainGauges[0].getDay(time);
	currentTime[3] = rainGauges[0].getHour(time);
	
	assignStationToNode();
	
	// SKY2008Snow from AJR2007
	int ct = 0;
	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	cNode = nodeIter.FirstP();
	while(nodeIter.IsActive()){ 
		for(int i=0; i<numStations;i++){
			if(assignedRain[ct] == rainGauges[i].getStation()){



				if ( (rainGauges[i].getRain(time) >= 1e-5) &&
						( rainGauges[i].getRain(time) + 
						  precLapseRate*(cNode->getZ() - rainGauges[i].getElev()) >= 1e-5) ) {

					gaugeRain[ct] = rainGauges[i].getRain(time) +
						precLapseRate*(cNode->getZ() - rainGauges[i].getElev());// AJR @ NMT 2007

				}
				else {
					gaugeRain[ct] = 0.0;
				}


				if(time == 0){
					latitude[ct] = rainGauges[i].getLat();
					longitude[ct] = rainGauges[i].getLong();
				}
			}
		}
		cNode->setRain(gaugeRain[ct]);
		cNode = nodeIter.NextP();
		ct++;
	}

	/*
	// Obtain values from tRainGauge
	for (int ct=0; ct<arraySize;ct++) {
		for (int i=0; i<numStations;i++) {
			if (assignedRain[ct] == rainGauges[i].getStation()) {
				gaugeRain[ct] = rainGauges[i].getRain(time);
				if (time == 0) {
					latitude[ct] = rainGauges[i].getLat();
					longitude[ct] = rainGauges[i].getLong();
				}
			}
		}
	}
	*/
}

/***************************************************************************
**
** tRainfall::readGaugeStat() Function
**
**
** Reads the Rain Gauge Station File which provides information concerning
** the raingauges used for rainfall input. Creates an array of tRainGauge
** for storing data. (see tRainGauge.h)
**
** Format for RainGaugeStation File:
**
** Header:
** nStations nParams (6)
**
** Body:
** StationID# FilePath Latitude Longitude RecordLength NumParm
**
** The file should have N rows corresponding to the N number of stations.
** and M columns corresponding to the data for each station.
**
** StationID  (int)         1->N
** FilePath   (string)      File name and path for the station datafile
** Reference Latitude  (double)  Station latitude  (basin projection)
** Reference Longitude (double)  Station longitude (basin projection)
** RecordLength  (int)      Number records for each station
** NumParm (int)            Number of parameter in each station record
**
***************************************************************************/
void tRainfall::readGaugeStat(char *stationfile) 
{
	int nStations, nParams;
	int stationID, numTimes, numParams;
	double rlat, rlong;

	// SKY2008Snow from AJR2007
	double elevation; // AJR @ NMT 2007

	char fileName[kName];

	// SKY2008Snow from AJR2007
	//assert(fileName != 0); //WR--09192023: comparison of array 'fileName' not equal to a null pointer is always true
	
	Cout<<"\nReading Rain Gauge Station File '";
	Cout<< stationfile<<"'..."<<endl<<flush;
	
	ifstream readFile(stationfile); 
	if (!readFile) {
		cout << "File "<<stationfile<<" not found." << endl;
		cout<<"Exiting Program...\n\n"<<endl;
		exit(1);
	}
	
	readFile >> nStations;
	readFile >> nParams;

	rainGauges = new tRainGauge[nStations];
	assert(rainGauges != 0);
	numStations = nStations;
	
	for (int count=0;count<nStations;count++) {
		for (int ct=0;ct<nParams;ct++) {
			if (ct==0) {
				readFile >> stationID;
				rainGauges[count].setStation(stationID);
			}
			if (ct==1) {
				readFile >> fileName;
				rainGauges[count].setFileName(fileName);
			}
			if (ct==2) {
				readFile >> rlat;
				rainGauges[count].setLat(rlat);
			}
			if (ct==3) {
				readFile >> rlong;
				rainGauges[count].setLong(rlong);
			}
			if (ct==4) {
				readFile >> numTimes;
				rainGauges[count].setTime(numTimes);
			}
			if (ct==5) {
				readFile >> numParams;
				rainGauges[count].setParm(numParams);
			}

			// SKY2008Snow from AJR2007
			if(ct==6){
				readFile >> elevation;
				rainGauges[count].setElev(elevation);
			}

		}
	}
	readFile.close();
}

/***************************************************************************
**
** tRainfall::readGaugeData() Function
**
**
** Reads and assigns data values to tRainGauge objects. 
** File Format:
** 
** Description Line: 
**      Abbreviations of Parameters as character strings
**      Ex. Y M D H R     
**
** Body Lines:
**      Values for each parameters. Read in as ints and doubles
**      Ex. Year (4 digit number), Month, Day, Hour (int)
**          Rain (double) in mm/hr
**
**      No Data Flag as 9999.99 for doubles
**
***************************************************************************/
void tRainfall::readGaugeData(int num) 
{
	int numParams, numTimes;
	char fileName[kMaxNameSize];
	int *year, *month, *day, *hour;
	double *Rain;
	double tempo;
	char *tmpstr;
	
	tmpstr = rainGauges[num].getFileName();
    snprintf(fileName,sizeof(fileName),"%s", tmpstr);//WR--09192023: 'sprintf' is deprecated: This function is provided for compatibility reasons only.
	numParams = rainGauges[num].getParm();
	numTimes  = rainGauges[num].getTime();
	
	cout<<"\nReading RainGauge Data File '";
	cout<< fileName<<"'..."<<endl<<flush;

	ifstream readDataFile(fileName);
	if (!readDataFile) {
		cout << "\nFile " <<fileName<<" not found!" << endl;
		cout << "Exiting Program...\n\n"<<endl;
		exit(2);}

	char paramNames[10];
	year  = new int[numTimes];
	month = new int[numTimes];
	day   = new int[numTimes];
	hour  = new int[numTimes];

	Rain  = new double[numTimes];
	
	for (int cnt = 0; cnt<numParams; cnt++) {
		readDataFile >> paramNames;
	}
	
	for (int count = 0;count<numTimes;count++) {
		for (int ct = 0;ct<numParams;ct++) {
			if (ct==0) {
				readDataFile >> year[count];}
			else if (ct==1) {
				readDataFile >> month[count];}
			else if (ct==2) {
				readDataFile >> day[count];}
			else if (ct==3) {
				readDataFile >> hour[count];
			}
			else if (ct==4) {
				readDataFile >> tempo;
				if (tempo < 0  || tempo > 200)
					Rain[count] = 9999.99;
				else
					Rain[count] = tempo;
			}
		}
	}

	readDataFile.close();
	
	rainGauges[num].setYear(year);
	rainGauges[num].setMonth(month);
	rainGauges[num].setDay(day);
	rainGauges[num].setHour(hour);

	robustNess(Rain, numTimes);
	
	rainGauges[num].setRain(Rain);

	delete [] Rain; delete [] hour;
	delete [] year; delete [] month; delete [] day;

}

/***************************************************************************
**
** tRainfall::robustNess() Function
**
** This function is used to check variables for NO_DATA flag = 9999.99
** It searches the double array forward and backwards, substituting the
** NO_DATA value with the previously read valid entry.
**
***************************************************************************/
void tRainfall::robustNess(double *variable, int size) 
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
	return;
}

/***************************************************************************
**
** tRainfall::assignStationToNode() Function
**
** Obtains the reference latitude and longitude for each station in the grid
** projection specified for the basin as well as the id code for each 
** tRainGauge station. Calls tResample for obtaining the Thiessen polygons
** assignments of each node to a station.
**
***************************************************************************/
void tRainfall::assignStationToNode() 
{
	double *stationLong, *stationLat;
	int *stationID;
	
	stationID   = new int[numStations];
	stationLong = new double[numStations];
	stationLat  = new double[numStations];
	
	arraySize = gridPtr->getNodeList()->getActiveSize();
	
	for (int ct=0;ct<numStations;ct++) {
		stationID[ct] = rainGauges[ct].getStation();
		stationLong[ct] = rainGauges[ct].getLong();
		stationLat[ct] = rainGauges[ct].getLat();
	}
	
	assignedRain = respPtr->doIt(stationID,stationLong,stationLat,numStations);
	
	delete [] stationID; 
	delete [] stationLong; 
	delete [] stationLat;
	return;
}

/***************************************************************************
**
** tRainfall::setToNode() Function
**
** Functions used to set the values of tCNode for each time step and for
** each voronoi cell. Also spatially-averaged raingauge input if desired.
**
***************************************************************************/
void tRainfall::setToNode() 
{
	int ct = 0;
	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	
	cNode = nodeIter.FirstP();
	
	while(nodeIter.IsActive()) { 
		cNode->setRain(gaugeRain[ct]);
		cNode = nodeIter.NextP();
		ct++;
	}
	
	// If Raindistribution = 1, weighted average of rain gauge data
	// Reassignment to tCNode after areal weighting
	
	if (optMAP == 1) {
		int id = 0;
		double *curGauge;
		
		double sumRain = 0.0;
		double sumArea = 0.0;
		double maxRain = 250.0;  //Maximum valid rainfall (mm/hr)
		
		arraySize = gridPtr->getNodeList()->getActiveSize();
		curGauge = new double[arraySize];
		
		// Compute Weighted Mean Gauge Rainfall in Basin
		cNode = nodeIter.FirstP();
		while( nodeIter.IsActive() ) {
			curGauge[id] = cNode->getRain();
			if (curGauge[id] < 0.0 || curGauge[id] > maxRain*rainDt)
				curGauge[id] = 0.0;
			sumRain = sumRain + cNode->getVArea()*curGauge[id];   
			sumArea = sumArea + cNode->getVArea();
			cNode = nodeIter.NextP();
			id++; 
		}
		
		// Assign Weighted Mean Rainfall Values
		cNode = nodeIter.FirstP();
		while( nodeIter.IsActive() ) { 
			cNode->setRain( (sumRain/sumArea) / rainDt );  
			cNode = nodeIter.NextP();
		}
		delete [] curGauge;  
	}
	return;
}

/***************************************************************************
**
** tRainfall::setfState() Function
**
** Functions used to set the values of fState
**
***************************************************************************/
void tRainfall::setfState(int state) 
{ 
	fState = state;
}

/***************************************************************************
**
** tRainfall::writeRestart() Function
** 
** Called from tSimulator during simulation loop
** 
***************************************************************************/
void tRainfall::writeRestart(fstream & rStr) const
{
  BinaryWrite(rStr, searchRain);
  BinaryWrite(rStr, rainfallType);
  BinaryWrite(rStr, numStations);
  BinaryWrite(rStr, arraySize);
  BinaryWrite(rStr, hourlyTimeStep);
  BinaryWrite(rStr, numRains);
  BinaryWrite(rStr, optForecast);
  BinaryWrite(rStr, fState);
  BinaryWrite(rStr, optMAP);
  BinaryWrite(rStr, rainDt);

  if (rainfallType == 3) {
    for (int i = 0; i < 4; i++)
      BinaryWrite(rStr, currentTime[i]);
    for (int i = 0; i < arraySize; i++) {
      BinaryWrite(rStr, gaugeRain[i]);
      BinaryWrite(rStr, latitude[i]);
      BinaryWrite(rStr, longitude[i]);
    } 
    for (int i = 0; i < numStations; i++)
      rainGauges[i].writeRestart(rStr);
  }

  BinaryWrite(rStr, precLapseRate);
  BinaryWrite(rStr, aveMAP);
  BinaryWrite(rStr, cumMAP);
  BinaryWrite(rStr, climate);

  tStorm::writeRestart(rStr);
}

/***************************************************************************
**
** tRainfall::readRestart() Function
**
***************************************************************************/
void tRainfall::readRestart(fstream & rStr)
{
  BinaryRead(rStr, searchRain);
  BinaryRead(rStr, rainfallType);
  BinaryRead(rStr, numStations);
  BinaryRead(rStr, arraySize);
  BinaryRead(rStr, hourlyTimeStep);
  BinaryRead(rStr, numRains);
  BinaryRead(rStr, optForecast);
  BinaryRead(rStr, fState);
  BinaryRead(rStr, optMAP);
  BinaryRead(rStr, rainDt);

  if (rainfallType == 3) {
    for (int i = 0; i < 4; i++)
      BinaryRead(rStr, currentTime[i]);
    for (int i = 0; i < arraySize; i++) {
      BinaryRead(rStr, gaugeRain[i]);
      BinaryRead(rStr, latitude[i]);
      BinaryRead(rStr, longitude[i]);
    }
    for (int i = 0; i < numStations; i++)
      rainGauges[i].readRestart(rStr);
  }

  BinaryRead(rStr, precLapseRate);
  BinaryRead(rStr, aveMAP);
  BinaryRead(rStr, cumMAP);
  BinaryRead(rStr, climate);

  tStorm::readRestart(rStr);
}

//=========================================================================
//
//
//                     End of tRainfall.cpp
//
//
//=========================================================================
