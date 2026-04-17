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
**  tPreProcess.cpp:  Function file for tPreProcess Class
**
***************************************************************************/

#include "src/Headers/globalIO.h"
#include "src/tSimulator/tPreProcess.h"

//=========================================================================
//
//
//                  Section 1: tPreProcess Constructor/Destructor
//
//
//=========================================================================

tPreProcess::tPreProcess()
{

}

tPreProcess::tPreProcess(tInputFile &infile) {
	CheckInputFile(infile);
}

tPreProcess::~tPreProcess() 
{
	Cout<<"tPreProcess Object has been destroyed..."<<endl;
}

//=========================================================================
//
//
//                  Section 2: tPreProcess Functions
//
//
//=========================================================================


/***************************************************************************
**
** tPreProcess::CheckInputFile() Function
**
** Function calls the input parameters in the *.in file and checks that
** each is there and has valid type of input (double, int, string). It
** calls functions in tInputFile for reading. No additional checks made
** to ensure pathnames are correct or data is valid. 
**
** NOTE: Upon adding new keywords to an *.in file, this function must be
** updated to reflect the changes. 
**
***************************************************************************/
void tPreProcess::CheckInputFile(tInputFile &infile) 
{
	double tempVariable = 0.0;
	int optmesh, optrain, optrock, optmet;

	int optres;// JECR 2015
	int optsoil;// JorgeGiuseppe 2015

	// SKY2008Snow from AJR2007
	int optradshelt; //, optwindshelt;

	int optfrcst, optgw;
   int optpar, optgraph, optrest, optv;
	char tempString[kName];

	// SKYnGM2008LU: added by AJR 2007
	int optlu;

	// SKYnGM2008LU
	int optluinterp;

	// Commented out several items for compatibility of existing data sets -VIVA //WR reverted 08282023

    // BEGIN Move tControl Arguments to .in file WR 08282023
    IterReadItem(infile, tempVariable,"OPTGROUNDWATER"); //   Cout<<"\t-G    Run groundwater model"<<endl;
    IterReadItem(infile, tempVariable,"OPTSPATIAL"); //  Cout<<"\t-R    Write intermediate states (spatial output)"<<endl;
    // END
	IterReadItem(infile, tempString,"STARTDATE");     //Run and time parameters
	IterReadItem(infile, tempVariable,"RUNTIME");
	IterReadItem(infile, tempVariable,"TIMESTEP");
	IterReadItem(infile, tempVariable,"GWSTEP");
	IterReadItem(infile, tempVariable,"METSTEP");
	IterReadItem(infile, tempVariable,"RAININTRVL");
	IterReadItem(infile, tempVariable,"OPINTRVL");
	IterReadItem(infile, tempVariable,"SPOPINTRVL");
	IterReadItem(infile, tempVariable,"INTSTORMMAX");
	//IterReadItem(infile, tempVariable,"RAINSEARCH");
	IterReadItem(infile, tempVariable, "TLINKE");
	IterReadItem(infile, tempVariable,"BASEFLOW");       //Flow parameters
	IterReadItem(infile, tempVariable,"KINEMVELCOEF");
	IterReadItem(infile, tempVariable,"FLOWEXP");
	IterReadItem(infile, tempVariable,"CHANNELROUGHNESS");
	
	tempVariable = IterReadItem(infile,tempVariable,"CHANNELWIDTHCOEFF");
	if (tempVariable <= 0) {
		IterReadItem(infile,tempVariable,"CHANNELWIDTH");
	}
	else {
		IterReadItem(infile, tempVariable,"CHANNELWIDTHEXPNT");
		IterReadItem(infile, tempString,  "CHANNELWIDTHFILE");
		IterReadItem(infile, tempVariable,"WIDTHINTERPOLATION");
	}
	
	IterReadItem(infile, tempVariable,"OPTEVAPOTRANS");   //Options
	IterReadItem(infile, tempVariable,"OPTINTERCEPT");
	IterReadItem(infile, tempVariable,"GFLUXOPTION");
	//IterReadItem(infile, tempVariable,"OPTRUNON");

	// SKY2008Snow from AJR2007
	double optsnow = IterReadItem(infile, tempVariable,"OPTSNOW");
	if (optsnow == 1) {
		// Check if the keyword exists before trying to read it
		if (infile.IsItemIn("SNOWFILENAME")) {
			IterReadItem(infile, tempString, "SNOWFILENAME");
			CheckFileExists(tempString, "SNOWFILENAME");
		} else {
			// If it's not there, we don't crash. 
			// We just print a warning that defaults will be used.
			cout << "\n-----------------------------------------------------------------" << endl;
			cout << " WARNING: 'SNOWFILENAME' not specified in the .in file." << endl;
			cout << " tRIBS will use the default hardcoded snow parameters." << endl;
			cout << "-----------------------------------------------------------------\n" << endl;
		}
	}
	IterReadItem(infile, tempVariable,"OPTRADSHELT");


	optmesh=optrain=optrock=optmet=optfrcst=optgw=0; //Int
  	optpar=optgraph=optrest=0;

	optres=0; // JECR 2015
	optsoil=0; // JorgeGiuseppe 2015
	
	optmesh  = IterReadItem(infile, optmesh, "OPTMESHINPUT");
	optrain  = IterReadItem(infile, optrain, "RAINSOURCE");
	optmet   = IterReadItem(infile, optmet,  "METDATAOPTION");
	
	// SKYnGM2008LU: added by AJR 2007
	optlu = IterReadItem(infile,optlu, "OPTLANDUSE");

	// SKYnGM2008LU
	if (optlu == 1) {
		optluinterp = IterReadItem(infile,optluinterp, "OPTLUINTERP");
	}

	optrock  = IterReadItem(infile, optrock ,"OPTBEDROCK");
	optfrcst = IterReadItem(infile, optfrcst,"FORECASTMODE");
	//optgw    = IterReadItem(infile, optgw,   "OPTGWFILE");
	
	optres = IterReadItem(infile, tempVariable,"OPTRESERVOIR"); // JECR 2015
	optsoil = IterReadItem(infile, tempVariable,"OPTSOILTYPE"); // JorgeGiuseppe 2015	

	if (optmesh == 1) {
		IterReadItem(infile, tempString,  "INPUTDATAFILE");
	}
	else if (optmesh == 2) {
		IterReadItem   (infile, tempString,"POINTFILENAME");
		CheckFileExists(tempString,"POINTFILENAME");
	}
	else if (optmesh == 3) {
		IterReadItem(infile, tempString,"ARCINFOFILENAME");
	}
	
	/****************** Start of modifications by JECR 2015 *********************/
	if (optres == 1) {	
		IterReadItem   (infile, tempString,"RESPOLYGONID");	//Reservoir polygon ID
		CheckFileExists(tempString,"RESPOLYGONID");

		IterReadItem   (infile, tempString,"RESDATA");    //Reservoir parameters
		CheckFileExists(tempString,"RESDATA");
	}

	if (optsoil == 1) {	 //JorgeGiuseppe2015
		IterReadItem   (infile, tempString,"SCGRID");    //File with soil grid paths
		CheckFileExists(tempString,"SCGRID");	
	}
	
	/******************** End of modifications by JECR 2015 *********************/ 

	IterReadItem   (infile, tempString,"SOILTABLENAME");    //Watershed grids
	CheckFileExists(tempString,"SOILTABLENAME"); 
	
	IterReadItem   (infile, tempString,"SOILMAPNAME");
	CheckFileExists(tempString,"SOILMAPNAME");
	
	IterReadItem   (infile, tempString,"LANDTABLENAME");
	CheckFileExists(tempString,"LANDTABLENAME");
	
	IterReadItem   (infile, tempString,"LANDMAPNAME");
	CheckFileExists(tempString,"LANDMAPNAME");
	
	if (!optgw) {
		IterReadItem   (infile, tempString,"GWATERFILE");   //Groundwater: input as grid
		CheckFileExists(tempString,"GWATERFILE");
	}

	// SKY2008Snow from AJR2007
	optradshelt  = IterReadItem(infile, optradshelt ,"OPTRADSHELT");
	if (optradshelt){
		IterReadItem(infile,tempString,"DEMFILE");
		CheckFileExists(infile,tempString,"DEMFILE");
	}

	if (optrain == 1) {                //Rainfall data
		IterReadItem(infile, tempString,"RAINFILE");
		IterReadItem(infile, tempString,"RAINEXTENSION");
		IterReadItem(infile, tempVariable,"RAINDISTRIBUTION");
	}
	else if (optrain == 2) {
		IterReadItem   (infile, tempString,"GAUGESTATIONS");
		CheckFileExists(tempString,"GAUGESTATIONS");
		IterReadItem   (infile, tempVariable,"RAINDISTRIBUTION");
	}
	
	if (optrock == 0) {                                  //Bedrock data
		IterReadItem(infile,tempVariable,"DEPTHTOBEDROCK");
	}
	else if (optrock == 1) {
		IterReadItem   (infile, tempString,"BEDROCKFILE");
		CheckFileExists(tempString,"BEDROCKFILE");
	}
	
	if (optmet == 1) {                                          //Meteorological data
		IterReadItem   (infile, tempString,"HYDROMETSTATIONS");
		CheckFileExists(tempString,"HYDROMETSTATIONS");
	}
	else if (optmet == 2) {
		IterReadItem   (infile, tempString,"HYDROMETGRID");
		CheckFileExists(tempString,"HYDROMETGRID");
	}

	// SKYnGM2008LU: added by AJR 2007
	if (optlu == 1) {
		IterReadItem( infile, tempString,"LUGRID");
		CheckFileExists(tempString, "LUGRID");
	}

	// SKY2008Snow from AJR2007
	IterReadItem(infile,tempString,"TEMPLAPSE"); //Lapse rates
	IterReadItem(infile,tempString,"PRECLAPSE");
	IterReadItem(infile,tempVariable,"HILLALBOPT");



	IterReadItem(infile, tempString,"OUTFILENAME");        //Output
#ifndef PARALLEL_TRIBS
	CheckPathNameCorrect(tempString, "OUTFILENAME");
#endif
	
	IterReadItem(infile, tempString,"OUTHYDROFILENAME");
#ifndef PARALLEL_TRIBS
	CheckPathNameCorrect(tempString, "OUTHYDROFILENAME");
#endif
	
	IterReadItem(infile, tempString,"NODEOUTPUTLIST");
	IterReadItem(infile, tempString,"HYDRONODELIST");
	IterReadItem(infile, tempString,"OUTLETNODELIST");

	if (optfrcst != 0 ) {                //Forecasting
		IterReadItem(infile, tempVariable,"FORECASTTIME");
		IterReadItem(infile, tempVariable,"FORECASTLEADTIME");
		IterReadItem(infile, tempVariable,"FORECASTLENGTH");
	}
	else if (optfrcst == 1)
		IterReadItem(infile, tempString,  "FORECASTFILE");
	else if (optfrcst == 3)
		IterReadItem(infile, tempVariable,"CLIMATOLOGY");
	
   // Restart options
   optrest = IterReadItem(infile, optrest, "RESTARTMODE");
   if (optrest > 0) {
     IterReadItem(infile, tempString, "RESTARTDIR");
     IterReadItem(infile, tempString, "RESTARTFILE");
     if (optrest == 1 || optrest == 3) {
       IterReadItem(infile, tempVariable, "RESTARTINTRVL");
     }
   }

   // Parallel and graph file options
   optpar = IterReadItem(infile, optpar, "PARALLELMODE");
   if (optpar > 0) {
     optgraph = IterReadItem(infile, optgraph, "GRAPHOPTION");
     if (optgraph > 0) {
       IterReadItem(infile, tempString, "GRAPHFILE");
     }
   }

	Cout<<"\nInput File Keywords Checked..."<<endl<<flush;
	return;
}

/***************************************************************************
**
** tPreProcess::CheckFileExists() Function
**
** Function to check if referenced pathname of file in the string KEYWORDS
** exists in the indicated location. No checking performed of file validity
** or structure, just the presence. 
**
***************************************************************************/
void tPreProcess::CheckFileExists(char* filename, const char* keyword)
{
	ifstream readFile(filename);
	if (!readFile) {
		cerr<<"\nFile '"<<filename<<"' for parameter '"<<keyword<<"' not found. Exiting."<<endl;
		exit(1);
	}
	return;
}

/***************************************************************************
**
** tPreProcess::CheckPathNameCorrect() Function
**
** Function to check the validity if referenced pathname 
**
***************************************************************************/
void tPreProcess::CheckPathNameCorrect(char* filename, const char* keyword)
{
	char bname[kName];
	char zero[] = "0";
	snprintf(bname, sizeof(bname), "%s%s", filename, zero);

	ofstream testFile(bname);
	if (!testFile) {
		cerr<<"\nOutput path '"<<filename<<"' for parameter '"<<keyword<<"' is not writable. Exiting."<<endl;
		exit(1);
	}
	testFile.close();
	remove(bname);
	return;
}

/***************************************************************************
**
** tPreProcess::IterReadItem() Functions
**
** Function to check if a parameter 'itemCode' exists in the input file
** No checking performed of file validity or structure, just the presence
**
***************************************************************************/
double tPreProcess::IterReadItem(tInputFile &infile, double datType,
								 const char *itemCode)
{
	datType = infile.ReadItem(datType, itemCode);
	Cout<<"Parameter =   "<<itemCode<<"\t\t\t"<<datType<<endl;
	if (datType < -999000.) {
		cerr<<"\nRequired parameter '"<<itemCode<<"' is missing or invalid in the input file. Exiting."<<endl;
		exit(1);
	}
	return datType;
}

/***************************************************************************
**
** tPreProcess::IterReadItem() Functions
**
** Function to check if a parameter 'itemCode' exists in the input file
** No checking performed of file validity or structure, just the presence
**
***************************************************************************/
int tPreProcess::IterReadItem(tInputFile &infile, int datType,
							  const char *itemCode)
{
	datType = infile.ReadItem(datType, itemCode);
	Cout<<"Parameter =   "<<itemCode<<"\t\t\t"<<datType<<endl;
	if (datType == -9999) {
		cerr<<"\nRequired parameter '"<<itemCode<<"' is missing or invalid in the input file. Exiting."<<endl;
		exit(1);
	}
	return datType;
}

/***************************************************************************
**
** tPreProcess::IterReadItem() Functions
**
** Function to check if a parameter 'itemCode' exists in the input file
** No checking performed of file validity or structure, just the presence
**
***************************************************************************/
void tPreProcess::IterReadItem(tInputFile &infile, char * theString,
							   const char *itemCode)
{
	char errr[] = "-999";
	infile.ReadItem(theString, itemCode);
	Cout<<"Parameter =   "<<itemCode<<"\t\t\t"<<theString<<endl;
	if (!strcmp(theString, errr)) {
		cerr<<"\nRequired parameter '"<<itemCode<<"' is missing or invalid in the input file. Exiting."<<endl;
		exit(1);
	}
	return;
}

//=========================================================================
//
//
//                    End of tPreProcess.cpp
//
//
//=========================================================================
