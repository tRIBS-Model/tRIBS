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
**  tHydroMetConvert.cpp:  Function file for tHydroMetConvert Class
**
***************************************************************************/

#include "src/tHydro/tHydroMetConvert.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//                  Section 1: tHydroMetConvert Constructor/Destructor
//
//
//=========================================================================

tHydroMetConvert::tHydroMetConvert(){}

tHydroMetConvert::tHydroMetConvert(tInputFile &infile)
{
	convertData = infile.ReadItem(convertData, "CONVERTDATA");
	
	if(convertData == 1 || convertData == 3){
		infile.ReadItem(mdiFile, "HYDROMETCONVERT");
		infile.ReadItem(mdfFilebase, "HYDROMETBASENAME");
		infile.ReadItem(sdfFile, "HYDROMETSTATIONS");     
	}
	else if(convertData == 2){
		infile.ReadItem(mdiFile, "GAUGECONVERT");
		infile.ReadItem(mdfFilebase, "GAUGEBASENAME");
		infile.ReadItem(sdfFile, "GAUGESTATIONS");
	}
	
	if(convertData != 3)
		initialize();
}

tHydroMetConvert::~tHydroMetConvert()
{
	if(convertData != 3){
		Cout<<"tHydroMetConvert Object has been destroyed..."<<endl;
		for(int ci=0;ci<numFiles;ci++)
			delete fNameArray[ci];
		for(int co=0;co<numFiles;co++)
			delete lNameArray[co];
		for(int ct=0;ct<numStations;ct++)
			delete sNameArray[ct];    
		for(int ce=0;ce<numParameters;ce++)
			delete pNameArray[ce];
		for(int ck=0;ck<Ncount;ck++){
			delete param[ck]; 
		}
		for(int cb=0;cb<Ncount;cb++){
			delete stat[cb];
		}
		delete [] elev; delete [] latitude; delete [] longitude;
		delete [] pNameArray; delete [] sNameArray; delete [] fNameArray; 
		delete [] date; delete [] hour; delete [] value; delete [] dateArray;
		delete [] lookFor; delete [] numTimes; delete [] dNameArray;
		delete [] param; delete [] stat; delete [] lNameArray;
	}
}

//=========================================================================
//
//
//                  Section 1: tHydroMetConvert Functions
//
//
//=========================================================================

/***************************************************************************
**
** tHydroMetConvert::callConvertRFC() Function
**
**
***************************************************************************/
void tHydroMetConvert::callConvertRFC()
{
	int option = 0;
	
	if(convertData == 1){
		cout<<"\nConverting HydroMet Data from RFC Format to tRIBS Format"<<endl;
		cout<<"-----------------------------------------------------------"<<endl;
	}
	else if(convertData == 2){
		cout<<"\nConverting RainGauge Data from RFC Format to tRIBS Format"<<endl;
		cout<<"-----------------------------------------------------------"<<endl;
	}
	
	readMDI();
	
	for(int ct=0;ct<numFiles;ct++){
		readLocData(ct);
		readPointData(ct);   
	}
	
	for(int cx=0;cx<numStations;cx++){
		if(convertData == 1)
			writeMDF(cx);
		else if(convertData == 2)
			writeGaugeMDF(cx);
	}
	
	if(convertData == 1)
		writeSDF(option);
	else if(convertData == 2)
		writeGaugeSDF(option);
	
	return;
}

/***************************************************************************
**
** tHydroMetConvert::callConvertDMIP() Function
**
**
***************************************************************************/
void tHydroMetConvert::callConvertDMIP()
{
	cout<<"\nConverting HydroMet Data from DMIP Format to tRIBS Format"<<endl;
	cout<<"-----------------------------------------------------------"<<endl;
	
	readAndWriteDMIP();   
	return; 
}

/***************************************************************************
**
** initialize() Function
**
***************************************************************************/
void tHydroMetConvert::initialize()
{
	param = new char*[kCount];
	stat = new char*[kCount];
	date = new double[kCount];
	hour = new double[kCount];
	value = new double[kCount];
	lookFor = new int[kCount];
	dateArray = new double[kCount];
	
	elev = new double[sCount];
	latitude = new double[sCount];
	longitude = new double[sCount];
	numTimes = new int[sCount];
	
	Ncount = 0;
	for(int ct=0;ct<kCount;ct++){
		date[ct] = 0.0;
		hour[ct] = 0.0;
		value[ct] = 0.0;
		lookFor[ct] = 0;
		dateArray[ct] = 0.0;
	}
	for(int st=0;st<sCount;st++){
		elev[st] = 0.0;
		latitude[st] = 0.0;
		longitude[st] = 0.0;
		numTimes[st] = 0;
	}
	return;
}

/***************************************************************************
**
** tHydroMetConvert::readMDI() Function
**
** Reads the Meteorologic Data Input (MDI) file used to specify which
** Operational RFC Point Data Files are to be opened, which stations 
** and parameters are to be read. All are character inputs. Fixed
** number of parameters per station, fixed number of stations per file.
**
** MDI File Structure:
** 
** #Files #Stations #Parameters
** MERGE or SEPARATE OPTION KEYWORD
** Name of File 1
** Name of Location File 1
** Name of File 2
** Name of Location File 2
** ...
** Name of Station 1 (#)
** Name of Station 2
** ..
** Name of Parameter 1    Number Station # -> Used for Merge only
** Name of Parameter 2    Number Station # -> # represents data belonging to
** ..                                         this particular station 
**
***************************************************************************/
void tHydroMetConvert::readMDI()
{
	char mergeOption[12];
	
	cout<<"\nIn readMDI..."<<endl<<flush;
	
	ifstream readMDI(mdiFile);
	if(!readMDI){
		cout<< "\nFile "<<mdiFile<<" not found..." <<endl;
		cout << "Exiting Program..."<<endl;
		exit(1);
	}
	
	readMDI >> numFiles;
	readMDI >> numStations;
	readMDI >> numParameters;
	readMDI >> mergeOption;
	
	if(strcmp(mergeOption,"MERGE")==0)
		optMerge = 1;
	else if(strcmp(mergeOption,"SEPARATE")==0)
		optMerge = 0;
	else{
		cout<<"\nHydroMet File Merge Option requires either";
		cout<<" MERGE or SEPARATE keywords"<<endl;
	}
	
	fNameArray = new char*[numFiles];
	lNameArray = new char*[numFiles];
	sNameArray = new char*[numStations];
	pNameArray = new char*[numParameters];
	dNameArray = new int[numParameters];
	
	for(int ct=0;ct<numFiles;ct++){   
		fNameArray[ct] = new char[kName];
		lNameArray[ct] = new char[kName];
		readMDI >> fNameArray[ct];
		readMDI >> lNameArray[ct];
	}
	for(int cti=0;cti<numStations;cti++){
		sNameArray[cti] = new char[10];
		readMDI >> sNameArray[cti];
	}
	for(int ctn=0;ctn<numParameters;ctn++){
		pNameArray[ctn] = new char[10];   
		readMDI >> pNameArray[ctn];
		readMDI >> dNameArray[ctn];		
	}
	readMDI.close();
	return;
}

/***************************************************************************
**
** tHydroMetConvert::readLocData() Function
**
** Reads the Operational RFC Location Data File containing Station ID's,
** and data describing the station: elevation, latitude and longitude
** (both in decimal degree). This data is stored and used in creating
** tRIBS *.sdf files (meteorologic station data files).
**
***************************************************************************/
void tHydroMetConvert::readLocData(int ct)
{
	char token[] = "|";
	char *elevString, *latString, *lngString;
	char *lName = new char[10];
	char *temp = new char[50];
	char *locRead = new char[kName];
	
	cout<<"\nIn readLocData..."<<endl<<flush;
	ifstream readLoc(lNameArray[ct]);
	if(!readLoc){
		cout << "\nFile "<<lNameArray[ct]<<" not found..." <<endl;
		cout << "Exiting Program..."<<endl;
		exit(2);
	}
	
	while(readLoc.getline(locRead,kName)){
		lName = strtok(locRead,token);   
		for(int scnt=0;scnt<numStations;scnt++){
			if(strcmp(lName, sNameArray[scnt])==0){
				temp = strtok(NULL,token);
				elevString = strtok(NULL,token);
				elev[scnt] = atof(elevString);
				latString = strtok(NULL,token);
				latitude[scnt] = atof(latString)/10000.0;
				lngString = strtok(NULL,token);
				longitude[scnt] = atof(lngString)/10000.0;
			}
		}
	} 
	readLoc.close();
	delete [] lName;
	delete [] temp;
	delete [] locRead;
	return;
}

/***************************************************************************
**
** tHydroMetConvert::readPointData() Function
**
** Reads the Operational RFC Point Data Files -- Tables delimited by "|"
** for each station/date/variable combination. Compares entry to the
** station and parameters specified in the MDI file. Saves the day, hour
** and parameter value for each entry.
**
***************************************************************************/
void tHydroMetConvert::readPointData(int ct)
{
	char *lineRead = new char[kName];
	char *dateString, *hourString, *valueString;
	char *yearString, *monString, *dayString, *timeString;
	char token[] = "|";
	char token2[] = "-";
	char *sName = new char[10];
	char *pName = new char[10];
	char *temp = new char[50];
	
	cout<<"\nIn readPointData..."<<endl<<flush;
	ifstream readPoint(fNameArray[ct]);
	if(!readPoint){
		cout<< "\nFile "<<fNameArray[ct]<<" not found..." <<endl;
		cout << "Exiting Program..."<<endl;
		exit(2);
	}
	
	while(readPoint >> lineRead){
		sName = strtok(lineRead,token); 
		for(int scnt=0;scnt<numStations;scnt++){
			if(strcmp(sName, sNameArray[scnt])==0){
				pName = strtok(NULL,token);
				for(int pcnt=0;pcnt<numParameters;pcnt++){
					if(strcmp(pName,pNameArray[pcnt])==0){
						for(int ij=0;ij<3;ij++)
							temp = strtok(NULL,token);
						dateString = strtok(NULL,token);
						if(dateString[0] == '9'){  
							hourString = strtok(NULL,token);
							valueString = strtok(NULL,token);
							char *years = new char[8];
							char t, t2, t3, t4, t5, t6;
							t = dateString[0]; t2 = dateString[1];
							t3 = dateString[2]; t4 = dateString[3];
							t5 = dateString[4]; t6 = dateString[5];
							years[0] = '1'; years[1] = '9';
							years[2] = t; years[3] = t2;
							years[4] = t3; years[5] = t4;
							years[6] = t5; years[7] = t6;
							date[Ncount] = atof(years);
							hour[Ncount] = atof(hourString);
							delete[] years;
						}
						else{
							yearString = strtok(dateString,token2);
							monString = strtok(NULL,token2);
							dayString = strtok(NULL,token2);
							timeString = strcat(strcat(yearString,monString),dayString);
							date[Ncount] = atof(timeString);
							readPoint >> lineRead;
							hourString = strtok(lineRead,token);
							valueString = strtok(NULL,token);
							hour[Ncount] = atof(hourString) * 10000;
						}
						value[Ncount] = atof(valueString);
						param[Ncount] = pNameArray[pcnt];
						stat[Ncount] = sNameArray[scnt];
						Ncount++;
					}
				} 
			}
		}
	}
	
	readPoint.close();
	delete [] sName;
	delete [] pName; 
	delete [] temp; 
	delete [] lineRead;
	return;
}

/***************************************************************************
**
** tHydroMetConvert::writeSDF() Function
**
** Writes the *.sdf file from the Operational RFC Point and Location Data
**
** Minor Modifications required to *.sdf outside of tRIBS:
**     1. Replace second latitude and longitude with those computed
**        for the projection of interest in the basin (i.e UTM, Equal-Area)
**        Necessarily equal to the projection used in input grids.
**        Not doing so, will cause failure of tResample for Point Stations
**     2. Replace GMT with the actual difference between the time zone
**        and the Greenwich Mean Time (-hour). Not doing so, will cause
**        failure of tEvapoTrans Radiation Scheme.
**
***************************************************************************/
void tHydroMetConvert::writeSDF(int Option)
{
	int numStatParams = 10;
	char gmt[] = "GMT";
	char blat[] = "BasinLat";
	char blong[] = "BasinLong";
	int numMetParams = 11;
	char command[kMaxNameLength];
	
	cout<<"\nIn writeSDF..."<<endl<<flush;
	
#ifdef ALPHA_64
    sprintf(command, "rm %s", sdfFile);
    system(command);
#elif defined LINUX_32
    unlink(sdfFile);
#elif defined WIN
#else 
    snprintf(command,sizeof(command), "rm %s", sdfFile);
    system(command);
#endif
	
	if(Option == 0){
		ifstream sdfIfs(sdfFile);  
		if(!sdfIfs){
			ofstream sdfOfs(sdfFile);  
			if(!sdfOfs.good()){
				cout<< "File "<<sdfFile<<" could not be opened!!!" <<endl;
				cout<<"Exiting Program..."<<endl;
				exit(1);}
			sdfOfs << numStations;
			sdfOfs << " " << numStatParams<<endl;
			for(int ct=0;ct<numStations;ct++){
				sdfOfs << ct+1 << " ";
				sdfOfs << mdfFilebase<<"_"<<sNameArray[ct]<< ".mdf ";
				sdfOfs << latitude[ct] << " ";
				sdfOfs << blat << " ";
				sdfOfs << longitude[ct] << " ";
				sdfOfs << blong << " ";
				sdfOfs << gmt << " ";
				sdfOfs << numTimes[ct] << " ";
				sdfOfs << numMetParams << " ";
				sdfOfs << elev[ct] << "\n";
			}
			sdfOfs.close();
		}
		else{
			cout<<"File "<<sdfFile<<" exists!"<<endl<<flush;
			cout<<"Exiting Program..."<<endl;
			exit(2);
		}
		sdfIfs.close();
	}
	else if(Option == 1){
		ofstream sdfOfs(sdfFile);
		if(!sdfOfs.good()){
			cout<< "File "<<sdfFile<<" could not be opened!!!" <<endl;
			cout<<"Exiting Program..."<<endl;
			exit(1);}
		sdfOfs << "1";
		sdfOfs << " " << numStatParams<<endl;
		sdfOfs << "1 ";
		sdfOfs << mdfFilebase<<"_merge"<< ".mdf ";
		sdfOfs << latitude[0] << " ";
		sdfOfs << blat << " ";
		sdfOfs << longitude[0] << " ";
		sdfOfs << blong << " ";
		sdfOfs << gmt << " ";
		sdfOfs << numTimes[0] << " ";
		sdfOfs << numMetParams << " ";
		sdfOfs << elev[0] << "\n";
		sdfOfs.close();
	}
	return;
}

/***************************************************************************
**
** tHydroMetConvert::writeGaugeSDF() Function
**
** Writes the *.sdf file from the Operational RFC Point data for the Rain
** Gauge files.
**
***************************************************************************/
void tHydroMetConvert::writeGaugeSDF(int Option)
{
	int numStatParams = 6;
	int numMetParams = 5;
	char command[kMaxNameLength];
	
#ifdef ALPHA_64
    sprintf(command, "rm %s", sdfFile);
    system(command);
#elif defined LINUX_32
    unlink(sdfFile);
#elif defined WIN
#else 
    snprintf(command,sizeof(command), "rm %s", sdfFile);
    system(command);
#endif
	
	cout<<"\nIn writeGaugeSDF..."<<endl<<flush;
	
	if(Option == 0){
		ifstream sdfIfs(sdfFile);  
		if(!sdfIfs){
			ofstream sdfOfs(sdfFile);  
			if(!sdfOfs.good()){
				cout<< "File "<<sdfFile<<" could not be opened!!!" <<endl;
				cout<<"Exiting Program..."<<endl;
				exit(1);}
			sdfOfs << numStations;
			sdfOfs << " " << numStatParams<<endl;
			for(int ct=0;ct<numStations;ct++){
				sdfOfs << ct+1 << " ";
				sdfOfs << mdfFilebase<<"_"<<sNameArray[ct]<< ".mdf ";
				sdfOfs << latitude[ct] << " ";
				sdfOfs << longitude[ct] << " ";
				sdfOfs << numTimes[ct] << " ";
				sdfOfs << numMetParams << " ";
			}
			sdfOfs.close();
		}
		else{
			cout<<"File "<<sdfFile<<" exists!"<<endl<<flush;
			cout<<"Exiting Program..."<<endl;
			exit(2);
		}
		sdfIfs.close();
	}
	else if(Option == 1){
		ofstream sdfOfs(sdfFile);
		if(!sdfOfs.good()){
			cout<< "File "<<sdfFile<<" could not be opened!!!" <<endl;
			cout<<"Exiting Program..."<<endl;
			exit(1);}
		sdfOfs << "1";
		sdfOfs << " " << numStatParams<<endl;
		sdfOfs << "1 ";
		sdfOfs << mdfFilebase<<"_merge"<< ".mdf ";
		sdfOfs << latitude[0] << " ";  
		sdfOfs << longitude[0] << " "; 
		sdfOfs << numTimes[0] << " ";
		sdfOfs << numMetParams << " ";
		sdfOfs.close();
	}
	return;
}

/***************************************************************************
**
** tHydroMetConvert::writeMDF() Function
**
** Writes the MDF files (Meteorologic Data Files) read directly by tRIBS
** This function reads the arrays created by readPointData() and manipulates
** them such as to print out all the parameters for each measurement hour.
** The algorithm attempts to optimize the search, sort and assigment of
** the parameters to the outputfile. No Data assigned 9999.99.
** 
** Parameters and Conversions:
**
** 	"PA" in inHG    converted to mb
** 	"TD" in F       converted to C
** 	"XR" in %       (no conversion)
** 	"XC" in tenths  (no conversion)
** 	"US" in mph     converted to m/s
** 	"TA" in F 	converted to C
** 	"TS" in F       converted to C
**
***************************************************************************/
void tHydroMetConvert::writeMDF(int ct)
{
	char fullName[kName];
	char extension[] = ".mdf"; 
	double currentDate, currentHour;
	double year, month, day; 
	double hourArray[24];
	char *pDate = new char[20];
	char tmpChar[4];
	char *dayStr, *monthStr, *tmpStr, *yearStr;
	
	//Compose output filename
	strcpy(fullName, mdfFilebase);
	strcat(fullName,"_");
	strcat(fullName,sNameArray[ct]);
	strcat(fullName,extension);
	
	cout<<"\nIn writeMDF "<<"for Station "<<ct+1<<"...."<<endl<<flush;
	
	ofstream mdfOfs(fullName);
	
	mdfOfs << "Y M D H ";
	mdfOfs << "PA ";
	mdfOfs << "TD ";
	mdfOfs << "XC ";
	mdfOfs << "US ";
	mdfOfs << "TA ";
	mdfOfs << "TS ";
	mdfOfs << "NR \n";
	
	//Look For rows corresponding to station 
	int numLines=0;
	for(int lp=0; lp<Ncount; lp++){
		if(strcmp(stat[lp],sNameArray[ct])==0){
			lookFor[numLines]=lp; 
			numLines++;}}
	
	//Extract Unique Dates from date array
	int size = 0;
	int nsize = 0;
	int notEqual = 0;
	for(int mn=0;mn<numLines;mn++){
		notEqual = 0;
		if(mn==0){
			dateArray[size] = date[lookFor[mn]];
			size++;}
		else{    
			nsize=0;
			for(int scr=0;scr<size;scr++){
				if(date[lookFor[mn]]!=dateArray[scr])
					notEqual++;}
			if(notEqual==size){
				dateArray[size] = date[lookFor[mn]];
				nsize++;}
			size = size+nsize;}
	}
	
	//Sort the Unique Dates, Assign dateArray
	insertSort(dateArray,size-1);
	
	//Assign hourArray
	for(int hr=0;hr<=23;hr++)
		hourArray[hr] = hr*10000; 
	
	//Initialize counters
	currentDate = dateArray[0]; 
	currentHour = hourArray[0];
	
	//While Loop for each Unique Date
	
	int totalSize = size*24;
	numTimes[ct] = totalSize;
	int hourCounter = 0;
	int dayCounter = 0;
	
	while( totalSize > 0 ){
		double lineWrite[11];
		//<--- Compose Year, Month, Day  
		//<--- Y2K Fix!
		snprintf(pDate,sizeof(pDate),"%f",currentDate);//WR--09192023: 'sprintf' is deprecated: This function is provided for compatibility reasons only.
		for(int i = 0; i<4; i++)
			tmpChar[i] = pDate[i];
		yearStr = tmpChar;;
		year = atol(yearStr);
		tmpStr = NULL; yearStr = NULL;
		
		char tmp[] = "";
		tmpStr = strncat(tmp, &(pDate[4]),1);
		monthStr = strncat(tmpStr, &(pDate[5]),1);
		month = atol(monthStr);
		
		tmpStr = NULL; monthStr = NULL;
		char tmp2[] = "";
		tmpStr = strncat(tmp2,&(pDate[6]), 1);
		dayStr = strncat(tmpStr, &(pDate[7]), 1);
		day = atol(dayStr);
		tmpStr = NULL; dayStr = NULL;
		
		//Create LineWrite for each Unique Date
		lineWrite[0] = (double)year;
		lineWrite[1] = (double)month;
		lineWrite[2] = (double)day; 
		lineWrite[3] = (double)(currentHour/10000);
		
		//No data value assigned 9999.999
		for(int sk=4;sk<=10;sk++)
			lineWrite[sk] = 9999.99;
		
		//cout<<"\nCurrentDate = "<< currentDate<<endl;
		
		for(int sc=0;sc<numLines;sc++){
			//cout<<"DateLookFor = "<<date[lookFor[sc]]<<endl;
			if(date[lookFor[sc]]==currentDate && hour[lookFor[sc]]==currentHour){
				if(strcmp(param[lookFor[sc]],"PA")==0)
					lineWrite[4] = value[lookFor[sc]]*33.8639; 
				if(strcmp(param[lookFor[sc]],"TD")==0)
					lineWrite[5] = (5.0/9.0)*(value[lookFor[sc]]-32.0);
				if(strcmp(param[lookFor[sc]],"XC")==0)
					lineWrite[6] = value[lookFor[sc]];  
				if(strcmp(param[lookFor[sc]],"US")==0)
					lineWrite[7] = value[lookFor[sc]]*0.44704;
				if(strcmp(param[lookFor[sc]],"TA")==0)
					lineWrite[8] = (5.0/9.0)*(value[lookFor[sc]]-32.0);
				if(strcmp(param[lookFor[sc]],"TS")==0)
					lineWrite[9] = (5.0/9.0)*(value[lookFor[sc]]-32.0);     
			}
		}
		
		mdfOfs << lineWrite[0] << " ";
		mdfOfs << lineWrite[1] << " ";
		mdfOfs << lineWrite[2] << " ";
		mdfOfs << lineWrite[3] << " ";
		mdfOfs << lineWrite[4] << " ";
		mdfOfs << lineWrite[5] << " ";
		mdfOfs << lineWrite[6] << " ";
		mdfOfs << lineWrite[7] << " ";
		mdfOfs << lineWrite[8] << " ";
		mdfOfs << lineWrite[9] << " ";
		mdfOfs << lineWrite[10] << "\n";
		
		hourCounter++;
		if(hourCounter==24){
			dayCounter++;
			currentDate = dateArray[dayCounter]; 
			hourCounter = 0;}
		currentHour = hourArray[hourCounter];
		totalSize--;
	}
	
	mdfOfs.close();
	for(int sz=0;sz<numLines;sz++)
		lookFor[sz] = 0;
	for(int st=0;st<size;st++)
		dateArray[st] = 0;
	
	delete [] pDate;
	return; 
}

/***************************************************************************
**
** tHydroMetConvert::writeGaugeMDF() Function
**
** Writes the MDF files (Meteorologic Data Files) read directly by tRIBS
** for the Rain Gauge input. Similar to writeMDF()
** 
** Parameters and Conversions:
**
** 	"PP" in inches/hr    converted to mm/hr (factor of 25.4)
**
***************************************************************************/
void tHydroMetConvert::writeGaugeMDF(int ct)
{
	char fullName[kName];
	char extension[] = ".mdf"; 
	double currentDate, currentHour;
	double year, month, day; 
	double hourArray[24];
	char *pDate = new char[20];
	char tmpChar[4];
	char *dayStr, *monthStr, *tmpStr, *yearStr;
	
	//Compose output filename
	strcpy(fullName, mdfFilebase);
	strcat(fullName,"_");
	strcat(fullName,sNameArray[ct]);
	strcat(fullName,extension);
	
	cout<<"\nIn writeGaugeMDF "<<"for Station "<<ct+1<<"...."<<endl<<flush;
	
	ofstream mdfOfs(fullName);
	
	mdfOfs << "Y M D H R \n";   //for rain gauge mdf
	
	//Look For rows corresponding to station 
	int numLines=0;
	for(int lp=0; lp<Ncount; lp++){
		if(strcmp(stat[lp],sNameArray[ct])==0){
			lookFor[numLines]=lp; 
			numLines++;}}
	
	//Extract Unique Dates from date array
	int size = 0;
	int nsize = 0;
	int notEqual = 0;
	for(int mn=0;mn<numLines;mn++){
		notEqual = 0;
		if(mn==0){
			dateArray[size] = date[lookFor[mn]];
			size++;}
		else{    
			nsize=0;
			for(int scr=0;scr<size;scr++){
				if(date[lookFor[mn]]!=dateArray[scr])
					notEqual++;}
			if(notEqual==size){
				dateArray[size] = date[lookFor[mn]];
				nsize++;}
			size = size+nsize;}
	}
	
	//Sort the Unique Dates, Assign dateArray
	insertSort(dateArray,size-1);
	
	//Assign hourArray
	for(int hr=0;hr<=23;hr++)
		hourArray[hr] = hr*10000; 
	
	//Initialize counters
	currentDate = dateArray[0]; 
	currentHour = hourArray[0];
	
	//While Loop for each Unique Date
	
	int totalSize = size*24;
	numTimes[ct] = totalSize;
	int hourCounter = 0;
	int dayCounter = 0;
	
	while( totalSize > 0 ){
		double lineWrite[11];
		
		//<--- Compose Year, Month, Day
        snprintf(pDate,sizeof (pDate),"%f",currentDate);//WR--09192023: 'sprintf' is deprecated: This function is provided for compatibility reasons only.
		for(int i = 0; i<4; i++)
			tmpChar[i] = pDate[i];
		yearStr = tmpChar;;
		year = atol(yearStr);
		tmpStr = NULL; yearStr = NULL;
		char tmp[] = "";
		tmpStr = strncat(tmp, &(pDate[2]),1);
		monthStr = strncat(tmpStr, &(pDate[3]),1);
		month = atol(monthStr);
		tmpStr = NULL; monthStr = NULL;
		char tmp2[] = "";
		tmpStr = strncat(tmp2,&(pDate[4]), 1);
		dayStr = strncat(tmpStr, &(pDate[5]), 1);
		day = atol(dayStr);
		tmpStr = NULL; dayStr = NULL;
		
		//Create LineWrite for each Unique Date
		lineWrite[0] = (double)year;
		lineWrite[1] = (double)month;
		lineWrite[2] = (double)day; 
		lineWrite[3] = (double)(currentHour/10000);
		
		//No data value assigned 9999.999
		for(int sk=4;sk<=4;sk++)     // for rain gauge mdf
			lineWrite[sk] = 9999.99;
		
		for(int sc=0;sc<numLines;sc++){
			if(date[lookFor[sc]]==currentDate && hour[lookFor[sc]]==currentHour){
				if(strcmp(param[lookFor[sc]],"PP")==0) {   // for rain gauge mdf
                    if (value[lookFor[sc]] >= 0) {
                        lineWrite[4] = value[lookFor[sc]] * 25.4;  //conversion from in to mm
                    } else {
                        lineWrite[4] = 9999.99;        //convert null values
                    }
                }
			}
		}
		
		mdfOfs << lineWrite[0] << " ";
		mdfOfs << lineWrite[1] << " ";
		mdfOfs << lineWrite[2] << " ";
		mdfOfs << lineWrite[3] << " ";
		mdfOfs << lineWrite[4] << "\n";
		
		hourCounter++;
		if(hourCounter==24){
			dayCounter++;
			currentDate = dateArray[dayCounter]; 
			hourCounter = 0;}
		currentHour = hourArray[hourCounter];
		totalSize--;
	}
	
	mdfOfs.close();
	for(int sz=0;sz<numLines;sz++)
		lookFor[sz] = 0;
	for(int st=0;st<size;st++)
		dateArray[st] = 0;
	
	delete [] pDate;
	return;
}

/***************************************************************************
**
** tHydroMetConvert::callMerge() Function
**
***************************************************************************/
void tHydroMetConvert::callMerge()
{
	char fullName[kName];
	char outName[kName];
	char extension[] = ".mdf"; 
	int numParms = 11;
	int columnRead, notUsedint, rows;
	double notUseddouble;
	double *pa, *td, *ta, *ws, *sc, *ts, *nr;
	int *dyear, *dday, *dmonth, *dhour;
	char labels[2];
	
	
	if(optMerge==0)
		cout<<"\nIn Merge Files Option...MERGE OFF..."<<endl<<flush;
	
	if(optMerge==1){   
		cout<<"\nIn Merge Files Option...MERGE ON..."<<endl<<flush;
		strcpy(outName, mdfFilebase);
		strcat(outName, "_merge");
		strcat(outName, extension);
		ofstream writeMerge(outName);
		
		dyear = new int[numTimes[0]];
		dmonth = new int[numTimes[0]];
		dday = new int[numTimes[0]];
		dhour = new int[numTimes[0]];
		
		pa = new double[numTimes[0]];
		td = new double[numTimes[0]];
		ta = new double[numTimes[0]];
		ws = new double[numTimes[0]];
		sc = new double[numTimes[0]];
		ts = new double[numTimes[0]];
		nr = new double[numTimes[0]];
		
		for(int sr=0;sr<numTimes[0];sr++){
			pa[sr] = 9999.99;
			td[sr] = 9999.99; 
			ta[sr] = 9999.99; 
			ws[sr] = 9999.99;
			sc[sr] = 9999.99;
			ts[sr] = 9999.99;
			nr[sr] = 9999.99;}
		
		for(int ct=0;ct<numStations;ct++){
			strcpy(fullName, mdfFilebase);
			strcat(fullName,"_");
			strcat(fullName,sNameArray[ct]);
			strcat(fullName,extension);
			for(int dt=0;dt<numParameters;dt++){
				if(dNameArray[dt]==(ct+1)){
					ifstream readData(fullName);
					if(!readData){
						cout << "File " <<fullName<<" not found!!!" << endl;
						cout<<"Exiting Program..."<<endl;
						exit(1);}
					for(int na=0;na<numParms;na++){
						readData >> labels;
						if(strcmp(labels,pNameArray[dt])==0)
							columnRead = na;}
					rows=0; 
					for(int lp=0;lp<numTimes[ct];lp++){
						if(ct==0){
							readData >> dyear[rows];
							readData >> dmonth[rows];
							readData >> dday[rows]; 
							readData >> dhour[rows];}
						else{
							for(int i=0;i<4;i++)
								readData >> notUsedint;}
						if(columnRead==4)
							readData >> pa[rows];
						else
							readData >> notUseddouble;
						if(columnRead==5)
							readData >> td[rows];
						else 
							readData >> notUseddouble;
						if(columnRead==6)
							readData >> sc[rows];
						else 
							readData >> notUseddouble;
						if(columnRead==7)
							readData >> ws[rows];
						else 
							readData >> notUseddouble;
						if(columnRead==8)
							readData >> ta[rows];
						else 
							readData >> notUseddouble;
						if(columnRead==9)
							readData >> ts[rows];
						else 
							readData >> notUseddouble;
						if(columnRead==10)
							readData >> nr[rows];
						else 
							readData >> notUseddouble;
						rows++;
					}
					readData.close();
				}
			}
		}
		writeMerge << "Y M D H ";
		writeMerge << "PA ";
		writeMerge << "TD ";
		writeMerge << "XC ";
		writeMerge << "US ";
		writeMerge << "TA ";
		writeMerge << "TS ";
		writeMerge << "NR \n";
		
		for(int rw=0;rw<rows;rw++){
			writeMerge << dyear[rw] << " ";
			writeMerge << dmonth[rw]<< " ";
			writeMerge << dday[rw]<< " ";
			writeMerge << dhour[rw]<< " ";
			writeMerge << pa[rw]<< " ";
			writeMerge << td[rw]<< " ";
			writeMerge << sc[rw]<< " ";
			writeMerge << ws[rw]<< " ";
			writeMerge << ta[rw]<< " ";
			writeMerge << ts[rw]<< " ";
			writeMerge << nr[rw]<< "\n";}  
		writeMerge.close();
		
		char sdfFileOLD[kName]; 
		char command[kName];
		strcpy(sdfFileOLD,sdfFile);
		strcat(sdfFileOLD,".old");
		rename(sdfFile, sdfFileOLD);
		
		int option = 1;
		writeSDF(option);
		
		delete [] pa; delete [] ta; delete [] ts; delete [] td; delete [] ws;
		delete [] sc; delete [] nr;
		delete [] dyear; delete [] dmonth; delete [] dday; delete [] dhour;
	}
	return;
}

/***************************************************************************
**
** tHydroMetConvert::callGaugeMerge() Function
**
***************************************************************************/
void tHydroMetConvert::callGaugeMerge()
{
	char fullName[kName];
	char outName[kName];
	char extension[] = ".mdf"; 
	int numParms = 5;  //for rain gauges
	int columnRead, notUsedint, rows;
	double notUseddouble;
	double *pp;
	int *dyear, *dday, *dmonth, *dhour;
	char labels[2];
	
	if(optMerge==0)
		cout<<"\nIn Merge Files Option...MERGE OFF..."<<endl<<flush;
	
	if(optMerge==1){   
		cout<<"\nIn Merge Files Option...MERGE ON.."<<endl<<flush;
		strcpy(outName, mdfFilebase);
		strcat(outName, "_merge");
		strcat(outName, extension);
		ofstream writeMerge(outName);
		
		dyear = new int[numTimes[0]];
		dmonth = new int[numTimes[0]];
		dday = new int[numTimes[0]];
		dhour = new int[numTimes[0]];
		
		pp = new double[numTimes[0]];  //rain gauge
		
		for(int sr=0;sr<numTimes[0];sr++)
			pp[sr] = 9999.99; 
		
		for(int ct=0;ct<numStations;ct++){
			strcpy(fullName, mdfFilebase);
			strcat(fullName,"_");
			strcat(fullName,sNameArray[ct]);
			strcat(fullName,extension);
			for(int dt=0;dt<numParameters;dt++){
				if(dNameArray[dt]==(ct+1)){
					ifstream readData(fullName);
					if(!readData){
						cout << "File " <<fullName<<" not found!!!" << endl;
						cout<<"Exiting Program..."<<endl;
						exit(1);}
					for(int na=0;na<numParms;na++){
						readData >> labels;
						if(strcmp(labels,pNameArray[dt])==0)
							columnRead = na;}
					rows=0; 
					for(int lp=0;lp<numTimes[ct];lp++){
						if(ct==0){
							readData >> dyear[rows];
							readData >> dmonth[rows];
							readData >> dday[rows]; 
							readData >> dhour[rows];}
						else{
							for(int i=0;i<4;i++)
								readData >> notUsedint;}
						if(columnRead==4)
							readData >> pp[rows];
						else
							readData >> notUseddouble;
						rows++;
					}
					readData.close();
				}
			}
		}
		writeMerge << "Y M D H R \n";
		
		for(int rw=0;rw<rows;rw++){
			writeMerge << dyear[rw] << " ";
			writeMerge << dmonth[rw]<< " ";
			writeMerge << dday[rw]<< " ";
			writeMerge << dhour[rw]<< " ";
			writeMerge << pp[rw]<< "\n";
		}  
		writeMerge.close();
		
		char sdfFileOLD[kName]; 
		char command[kName];
		strcpy(sdfFileOLD,sdfFile);
		strcat(sdfFileOLD,".old");
		rename(sdfFile, sdfFileOLD);
		
		int option = 1;
		writeSDF(option);
		
		delete[] dyear; delete[] dmonth; delete[] dday; delete[] dhour;
		delete[] pp;
	}
	return;
}

/***************************************************************************
**
** tHydroMetConvert::readAndWriteDMIP() Function
**
** Reads the DMIP Observed HydroMet Data Files for a basin used in
** the HYDROMETCONVERT keyword (ie. 07197000.dat for Baron Fork)
**
** Different variables than for RFC data:
**  DMIP = "Y M D H PA VP XC US TA IS IL"
**  RFC  = "Y M D H PA TD XC US TA TS NR"
** 
** Assign VP in same fashion as TD or RH. Ignore IS, IL radiation
**
***************************************************************************/
void tHydroMetConvert::readAndWriteDMIP()
{
	int nFields = 9;
	int nParams = 11;
	double rain, temp, swrad, lwrad;
	double press, vpress, wind, dens;
	int cyear, cmonth, cday, chour;
	double cover = 9999.99;
	char cdate[12];
	char metName[kName];
	char gaugeName[kName];
	char lineRead[kName];
	char extension[] = ".mdf";
	int numLines = 0;
	char *temp1 = new char[4];
	char *temp2 = new char[2];
	char *temp3 = new char[2];
	char *temp4 = new char[2];
	
	cout<<"\nIn readAndWriteDMIP..."<<endl<<flush;
	
	ifstream readDMIPfile(mdiFile);
	
	strcpy(metName, mdfFilebase);     //MDF file name
	strcat(metName, extension);
	
	strcpy(gaugeName, mdfFilebase); //MDF gauge file name
	strcat(gaugeName,"_gauge");
	strcat(gaugeName,extension);
	
	ofstream writeDMIPfile(metName);
	
	ofstream writeDMIPgauge(gaugeName);
	
	writeDMIPfile << "Y M D H PA VP XC US TA IS IL\n";
	writeDMIPgauge << "Y M D H R\n";
	
	if(!readDMIPfile){
		cout<< "\nFile "<<mdiFile<<" not found..." <<endl;
		cout << "Exiting Program..."<<endl;
		exit(1);
	}
	
	//Find out length of file
	while(readDMIPfile.getline(lineRead,kName)){
		numLines++;
	}
	
	readDMIPfile.close();
	
	readDMIPfile.open(mdiFile);
	
	for(int nl = 0;nl<numLines;nl++){
		for(int i = 0; i<nFields;i++){    //Read from file
			if(i==0)
				readDMIPfile >> cdate;
			if(i==1)
				readDMIPfile >> rain;
			if(i==2)
				readDMIPfile >> temp;
			if(i==3)
				readDMIPfile >> swrad;
			if(i==4)
				readDMIPfile >> lwrad;
			if(i==5)
				readDMIPfile >> dens;
			if(i==6)
				readDMIPfile >> press;
			if(i==7)
				readDMIPfile >> vpress;
			if(i==8)
				readDMIPfile >> wind;
		}
		
		//Parse date string
		temp1[0] = cdate[0];
		temp1[1] = cdate[1];
		temp1[2] = cdate[2]; 
		temp1[3] = cdate[3];
		cyear = atoi(temp1);
		
		temp2[0] = cdate[4];
		temp2[1] = cdate[5];
		cmonth = atoi(temp2);
		
		temp3[0] = cdate[6];
		temp3[1] = cdate[7];
		cday = atoi(temp3); 
		
		temp4[0] = cdate[8];
		temp4[1] = cdate[9];
		chour = atoi(temp4);
		
		for(int j = 0; j<nParams;j++){       //Write to file
			if(j==0){                           //Unit conversions
				writeDMIPfile << cyear<<" ";
				writeDMIPgauge << cyear<<" ";
			}
			if(j==1){
				writeDMIPfile << cmonth<<" ";
				writeDMIPgauge << cmonth<<" ";
			}
			if(j==2){
				writeDMIPfile << cday<<" ";
				writeDMIPgauge << cday<<" ";
			}
			if(j==3){
				writeDMIPfile << chour-1<<" ";
				writeDMIPgauge << chour-1<<" ";
			}
			if(j==4){
				writeDMIPfile << (press*10.0)<<" ";
				writeDMIPgauge << rain <<" ";
			}
			if(j==5)
				writeDMIPfile << (vpress*10.0)<<" ";
			if(j==6)
				writeDMIPfile << cover<<" ";
			if(j==7)
				writeDMIPfile << wind<<" ";
			if(j==8)
				writeDMIPfile << temp<<" ";
			if(j==9)
				writeDMIPfile << swrad<<" ";
			if(j==10)
				writeDMIPfile << lwrad;
		}
		writeDMIPfile << "\n";
		writeDMIPgauge << "\n";
	}
	readDMIPfile.close();
	writeDMIPfile.close();
	writeDMIPgauge.close();
	
	delete[] temp1;
	delete[] temp2;
	delete[] temp3;
	delete[] temp4;
	
	return;
}

/***************************************************************************
**
** tHydroMetConvert::insertSort() Function
**
***************************************************************************/
void tHydroMetConvert::insertSort(double a[], int N)
{
	int i,j;
	double v;
	
	for(i=2;i<=N;i++){
		v=a[i];
		j = i;
		while(a[j-1]>v){
			a[j] = a[j-1];
			j--;
		}
		a[j] = v;
	}  
	return;
}


//=========================================================================
//
//
//                    End of tHydroMetConvert.cpp
//
//
//=========================================================================
