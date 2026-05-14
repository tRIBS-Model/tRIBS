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
	optMAP = 0;

	// If data are used as the model forcing
	if (rainfallType == 1)
		Cout<<"Rainfall Source: \t\tGridded Rainfall [mm/hr]\n";
	
	if (rainfallType == 1) { 
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
	else if (rainfallType == 2) {
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
		Cout<<"\t\t(1) Gridded Rainfall"<<endl;
		Cout<<"\t\t(2) Rain Gauge Station Rainfall"<<endl;
		Cout << "Exiting Program...\n\n"<<endl;
		exit(1);
	}
	
}

// Destructor
tRainfall::~tRainfall() 
{
	gridPtr = nullptr;
	respPtr = nullptr;
	simCtrl = nullptr;
	
	if (rainfallType == 1) {
		if ( infile.is_open() )
			infile.close();
	}
	if (rainfallType == 2) {
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
    if ( infile.is_open() )
		infile.close();
	
	if (t->minuteRn || t->dtRain < 1)
		snprintf(mrainfileIn,sizeof(mrainfileIn), "%s%02d%02d%04d%02d%02d.%s", inputname,
				t->monthRn, t->dayRn, t->yearRn, t->hourRn, t->minuteRn, extension);
	else {
		snprintf(mrainfileIn,sizeof(mrainfileIn),"%s%02d%02d%04d%02d.%s", inputname,
				t->monthRn, t->dayRn, t->yearRn, t->hourRn, extension);
	}
	infile.open(mrainfileIn);

	// Check if file opened
    if ( !(infile.is_open()) )
		return 0;
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

	curRain = respPtr->doIt(mrainfileIn, 1);

	while( nodeIter.IsActive() ) {
		if (curRain[id] < 0.0 || curRain[id] > maxRain*t->getRainDT())
			curRain[id] = 0.0;
		sumRain = sumRain + cn->getVArea()*curRain[id];
		sumArea = sumArea + cn->getVArea();
		cn = nodeIter.NextP();
		id++;
	}

	id = 0;
	cn = nodeIter.FirstP();
	while( nodeIter.IsActive() ) {
		if (optMAP == 1)
			curRain[id]=sumRain/sumArea;
		if (rainfallType == 1)
			cn->setRain( curRain[id]/t->getRainDT());
		cn = nodeIter.NextP();
		id++;
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
	NewRainData(hourlyTimeStep, t);


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
void tRainfall::NewRainData(int time, tRunTimer *timer) 
{  
    // Begin block for validation of rainfall data input 
    // Get the timestamp from the first rain gauge's data for the current index
    int file_yr = rainGauges[0].getYear(time);
    int file_mo = rainGauges[0].getMonth(time);
    int file_dy = rainGauges[0].getDay(time);
    int file_hr = rainGauges[0].getHour(time);

    // Get the current simulation time from the timer object
    int sim_yr = timer->year;
    int sim_mo = timer->month;
    int sim_dy = timer->day;
    int sim_hr = timer->hour;

    // Compare the two timestamps
    if (file_yr != sim_yr || file_mo != sim_mo || file_dy != sim_dy || file_hr != sim_hr) {
        std::cerr << "\n\nFATAL ERROR in rain gauge data" << std::endl;
        std::cerr << "Timestamp mismatch detected during simulation." << std::endl;
        std::cerr << "Simulation expected data for: " << sim_yr << "/" << sim_mo << "/" << sim_dy 
                  << " " << sim_hr << ":00" << std::endl;
        std::cerr << "But found data for:           " << file_yr << "/" << file_mo << "/" << file_dy 
                  << " " << file_hr << ":00" << std::endl;
        std::cerr << "This indicates a gap, duplicate, or incorrect data in the rain gauge file(s)." << std::endl;
        std::cerr << "Exiting Program...\n\n" << std::endl;
        exit(1);
    }
    // End block for validation of rainfall data input
	
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
** Format for RainGaugeStation File (CSV):
**
** Header (required, 5 columns):
** ID,DataFile,Northing,Easting,Elevation
**
** Body (one row per station):
** ID        (int)     Station identifier
** DataFile  (string)  Path to the station MDF data file
** Northing  (double)  Station northing coordinate (basin projection)
** Easting   (double)  Station easting coordinate (basin projection)
** Elevation (double)  Station elevation in meters
**
***************************************************************************/
void tRainfall::readGaugeStat(char *stationfile)
{
	Cout<<"\nReading Rain Gauge Station File '"
	    << stationfile <<"'..."<<endl<<flush;

	ifstream readFile(stationfile);
	if (!readFile) {
		cout << "File "<<stationfile<<" not found." << endl;
		cout<<"Exiting Program...\n\n"<<endl;
		exit(1);
	}

	// Validate CSV header: ID,DataFile,Northing,Easting,Elevation
	std::string headerLine;
	std::getline(readFile, headerLine);
	{
		std::istringstream hss(headerLine);
		std::string token;
		int ncols = 0;
		while (std::getline(hss, token, ',')) ncols++;
		if (ncols != 5) {
			cerr << "\nError: Rain gauge station file '" << stationfile
			     << "' header must have 5 columns "
			        "(ID,DataFile,Northing,Easting,Elevation)." << endl;
			exit(1);
		}
	}

	// Read all data rows
	std::vector<std::string> rows;
	std::string line;
	while (std::getline(readFile, line)) {
		if (!line.empty()) rows.push_back(line);
	}
	readFile.close();

	numStations = static_cast<int>(rows.size());
	rainGauges  = new tRainGauge[numStations];
	assert(rainGauges != 0);

	for (int count = 0; count < numStations; count++) {
		std::istringstream ss(rows[count]);
		std::string token;
		char fileName[kName];

		std::getline(ss, token, ',');
		rainGauges[count].setStation(std::stoi(token));

		std::getline(ss, token, ',');
		strncpy(fileName, token.c_str(), kName - 1);
		fileName[kName - 1] = '\0';
		rainGauges[count].setFileName(fileName);

		std::getline(ss, token, ',');
		rainGauges[count].setLat(std::stod(token));

		std::getline(ss, token, ',');
		rainGauges[count].setLong(std::stod(token));

		std::getline(ss, token, ',');
		rainGauges[count].setElev(std::stod(token));
	}
}

/***************************************************************************
**
** tRainfall::readGaugeData() Function
**
**
** Reads and assigns data values to tRainGauge objects.
** File Format (CSV):
**
** Header (required, 5 columns):
** Year,Month,Day,Hour,Rain_mm/hr
**
** Body Lines (one row per timestep):
**      Year  (int)    4-digit year
**      Month (int)    Month (1-12)
**      Day   (int)    Day of month
**      Hour  (int)    Hour of day
**      Rain  (double) Rainfall rate in mm/hr; values <0 or >200 flagged as NO_DATA (9999.99)
**
***************************************************************************/
void tRainfall::readGaugeData(int num)
{
	char fileName[kMaxNameSize];
	snprintf(fileName, sizeof(fileName), "%s", rainGauges[num].getFileName());

	cout<<"\nReading RainGauge Data File '"<<fileName<<"'..."<<endl<<flush;

	ifstream readDataFile(fileName);
	if (!readDataFile) {
		cout << "\nFile " <<fileName<<" not found!" << endl;
		cout << "Exiting Program...\n\n"<<endl;
		exit(2);
	}

	// Validate CSV header: Year,Month,Day,Hour,Rain_mm/hr
	std::string headerLine;
	std::getline(readDataFile, headerLine);
	{
		std::istringstream hss(headerLine);
		std::string token;
		int ncols = 0;
		while (std::getline(hss, token, ',')) ncols++;
		if (ncols != 5) {
			cerr << "\nError: Rain gauge data file '" << fileName
			     << "' header must have 5 columns "
			        "(Year,Month,Day,Hour,Rain_mm/hr)." << endl;
			exit(1);
		}
	}

	std::vector<int>    Year, Month, Day, Hour;
	std::vector<double> Rain;

	std::string line;
	while (std::getline(readDataFile, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string token;

		std::getline(ss, token, ',');
		Year.push_back(std::stoi(token));
		std::getline(ss, token, ',');
		Month.push_back(std::stoi(token));
		std::getline(ss, token, ',');
		Day.push_back(std::stoi(token));
		std::getline(ss, token, ',');
		Hour.push_back(std::stoi(token));
		std::getline(ss, token, ',');
		double tempo = std::stod(token);
		Rain.push_back((tempo < 0 || tempo > 200) ? 9999.99 : tempo);
	}
	readDataFile.close();

	robustNess(Rain.data(), static_cast<int>(Rain.size()));

	rainGauges[num].setYear(Year);
	rainGauges[num].setMonth(Month);
	rainGauges[num].setDay(Day);
	rainGauges[num].setHour(Hour);
	rainGauges[num].setRain(Rain);
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
  BinaryWrite(rStr, optMAP);
  BinaryWrite(rStr, rainDt);

  if (rainfallType == 2) {
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
  BinaryRead(rStr, optMAP);
  BinaryRead(rStr, rainDt);

  if (rainfallType == 2) {
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
}

//=========================================================================
//
//
//                     End of tRainfall.cpp
//
//
//=========================================================================
