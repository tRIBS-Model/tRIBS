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

tPreProcess::tPreProcess(SimulationControl *simCtrPtr, tInputFile &infile) {
	
	simCtrl = simCtrPtr;
	
	if (simCtrl->Check_label == 'Y') {
		CheckInputFile(infile);
	}
	
	convertData = infile.ReadItem(convertData, "CONVERTDATA");
	
	if (convertData == 1) {
		tHydroMetConvert hydroMetInput(infile);
		hydroMetInput.callConvertRFC();
		hydroMetInput.callMerge();
		cout<<"\n-----------------------------------------------"<<endl;
		cout<<"tRIBS HydroMet Pre-Processor completed"<<endl;
		cout<<"Use tRIBS with Convert Data = 0 for Model Runs"<<endl;
		cout<<"Exiting Program..."<<endl;
		cout<<"-----------------------------------------------"<<endl;
		exit(1);
	}
	else if (convertData == 2) {
		tHydroMetConvert rainGaugeInput(infile);
		rainGaugeInput.callConvertRFC();
		rainGaugeInput.callGaugeMerge();
		cout<<"\n-----------------------------------------------"<<endl;
		cout<<"tRIBS RainGauge Pre-Processor completed"<<endl;
		cout<<"Use tRIBS with Convert Data = 0 for Model Runs"<<endl;
		cout<<"Exiting Program..."<<endl;
		cout<<"-----------------------------------------------"<<endl;
		exit(2);
	}
	else if (convertData == 3) {
		tHydroMetConvert hydroMetInput(infile);
		hydroMetInput.callConvertDMIP();
		cout<<"\n-----------------------------------------------"<<endl;
		cout<<"tRIBS HydroMet Pre-Processor completed... "<<endl;
		cout<<"MDF Files Created for Met and Rain Gauge Data..."<<endl;
		cout<<"Use tRIBS with Convert Data = 0 for Model Runs"<<endl;
		cout<<"Must create SDF file based on Station Data for use in tRIBS"<<endl;
		cout<<"Exiting Program..."<<endl;
		cout<<"-----------------------------------------------"<<endl;
		exit(3);
	}
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
	int optmesh, optrain, optrock, optconv, optmet;

	int optres;// JECR 2015
	int optsoil;// JorgeGiuseppe 2015

	// SKY2008Snow from AJR2007
	int optradshelt; //, optwindshelt;

	int optfrcst, optstoch, optgw;
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
    IterReadItem(infile, tempVariable,"OPTINTERHYDRO"); // Cout<<"\t-H    Write intermediate hydrographs (.mrf)"<<endl;
    IterReadItem(infile, tempVariable,"OPTHEADER"); //  Cout<<"\t-M    Do NOT Write headers in pixel/hydrograph/voronoi output files"<<endl;
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
	IterReadItem(infile, tempVariable,"VELOCITYCOEF");
	IterReadItem(infile, tempVariable,"VELOCITYRATIO");
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
	IterReadItem(infile, tempVariable,"OPTSNOW");
	IterReadItem(infile, tempVariable,"MINSNTEMP");
	IterReadItem(infile, tempVariable,"SNLIQFRAC"); // Added by CJC 2020
	IterReadItem(infile, tempVariable,"OPTRADSHELT");


	optmesh=optrain=optrock=optconv=optmet=optfrcst=optstoch=optgw=0; //Int
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

	optconv  = IterReadItem(infile, optconv, "CONVERTDATA");

	optrock  = IterReadItem(infile, optrock ,"OPTBEDROCK");
	optfrcst = IterReadItem(infile, optfrcst,"FORECASTMODE");
	optstoch = IterReadItem(infile, optstoch,"STOCHASTICMODE");
	//optgw    = IterReadItem(infile, optgw,   "OPTGWFILE");
	
	optres = IterReadItem(infile, tempVariable,"OPTRESERVOIR"); // JECR 2015
	optsoil = IterReadItem(infile, tempVariable,"OPTSOILTYPE"); // JorgeGiuseppe 2015	

	if (optmesh == 1) {
		IterReadItem(infile, tempString,  "INPUTDATAFILE");
		IterReadItem(infile, tempVariable,"INPUTTIME");
	}
	else if (optmesh == 2) {
		IterReadItem   (infile, tempString,"POINTFILENAME");
		CheckFileExists(infile, tempString,"POINTFILENAME");
	}
	else if (optmesh == 3) {
		IterReadItem(infile, tempString,"ARCINFOFILENAME");
	}
	
	/****************** Start of modifications by JECR 2015 *********************/
	if (optres == 1) {	
		IterReadItem   (infile, tempString,"RESPOLYGONID");	//Reservoir polygon ID
		CheckFileExists(infile, tempString,"RESPOLYGONID");

		IterReadItem   (infile, tempString,"RESDATA");    //Reservoir parameters
		CheckFileExists(infile, tempString,"RESDATA");
	}

	if (optsoil == 1) {	 //JorgeGiuseppe2015
		IterReadItem   (infile, tempString,"SCGRID");    //File with soil grid paths
		CheckFileExists(infile, tempString,"SCGRID");	
	}
	
	/******************** End of modifications by JECR 2015 *********************/

	IterReadItem   (infile, tempString,"SOILTABLENAME");    //Watershed grids
	CheckFileExists(infile, tempString,"SOILTABLENAME"); 
	
	IterReadItem   (infile, tempString,"SOILMAPNAME");
	CheckFileExists(infile, tempString,"SOILMAPNAME");
	
	IterReadItem   (infile, tempString,"LANDTABLENAME");
	CheckFileExists(infile, tempString,"LANDTABLENAME");
	
	IterReadItem   (infile, tempString,"LANDMAPNAME");
	CheckFileExists(infile, tempString,"LANDMAPNAME");
	
	if (!optgw) {
		IterReadItem   (infile, tempString,"GWATERFILE");   //Groundwater: input as grid
		CheckFileExists(infile, tempString,"GWATERFILE");
	}

	// SKY2008Snow from AJR2007
	optradshelt  = IterReadItem(infile, optradshelt ,"OPTRADSHELT");
	if (optradshelt){
		IterReadItem(infile,tempString,"DEMFILE");
		CheckFileExists(infile,tempString,"DEMFILE");
	}

	if (optrain == 1 || optrain == 2) {                //Rainfall data
		IterReadItem(infile, tempString,"RAINFILE");
		IterReadItem(infile, tempString,"RAINEXTENSION");
		IterReadItem(infile, tempVariable,"RAINDISTRIBUTION");
	}
	else if (optrain == 3) {
		IterReadItem   (infile, tempString,"GAUGESTATIONS");
		CheckFileExists(infile, tempString,"GAUGESTATIONS");
		IterReadItem   (infile, tempVariable,"RAINDISTRIBUTION");
	}
	
	if (optrock == 0) {                                  //Bedrock data
		IterReadItem(infile,tempVariable,"DEPTHTOBEDROCK");
	}
	else if (optrock == 1) {
		IterReadItem   (infile, tempString,"BEDROCKFILE");
		CheckFileExists(infile, tempString,"BEDROCKFILE");
	}
	
	if (optconv == 1) {
		IterReadItem   (infile, tempString,"HYDROMETCONVERT");   //Hydromet data
		CheckFileExists(infile, tempString,"HYDROMETCONVERT");
		IterReadItem   (infile, tempString,"HYDROMETBASENAME");
	}
	else if (optconv == 2) {
		IterReadItem   (infile, tempString,"GAUGECONVERT");
		CheckFileExists(infile, tempString,"GAUGECONVERT");
		IterReadItem   (infile, tempString,"GAUGEBASENAME");
	}
	
	if (optmet == 1) {
		IterReadItem   (infile, tempString,"HYDROMETSTATIONS");
		CheckFileExists(infile, tempString,"HYDROMETSTATIONS");
	}
	else if (optmet == 2) {
		IterReadItem   (infile, tempString,"HYDROMETGRID");
		CheckFileExists(infile, tempString,"HYDROMETGRID");
	}

	// SKYnGM2008LU: added by AJR 2007
	if (optlu == 1) {
		IterReadItem( infile, tempString,"LUGRID");
		CheckFileExists(infile, tempString, "LUGRID");
	}

	// SKY2008Snow from AJR2007
	IterReadItem(infile,tempString,"TEMPLAPSE"); //Lapse rates
	IterReadItem(infile,tempString,"PRECLAPSE");
	IterReadItem(infile,tempVariable,"HILLALBOPT");



	IterReadItem(infile, tempString,"OUTFILENAME");        //Output
#ifndef PARALLEL_TRIBS
	CheckPathNameCorrect(infile, tempString, "OUTFILENAME");
#endif
	
	IterReadItem(infile, tempString,"OUTHYDROFILENAME");
#ifndef PARALLEL_TRIBS
	CheckPathNameCorrect(infile, tempString, "OUTHYDROFILENAME");
#endif
	
	IterReadItem(infile, tempString,"OUTHYDROEXTENSION");
	//IterReadItem(infile,tempString,"RIBSHYDOUTPUT");
	
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
	
	if (optstoch != 0 && optstoch != 6) {         //Stochastic Rainfall
		IterReadItem(infile, tempVariable,"PMEAN");
		IterReadItem(infile, tempVariable,"STDUR");
		IterReadItem(infile, tempVariable,"ISTDUR");
		if (optstoch != 1)
			IterReadItem(infile,tempVariable, "SEED");
		if (optstoch == 3 ||  optstoch == 4 || optstoch == 5) {
			IterReadItem(infile, tempVariable,"PERIOD");
			IterReadItem(infile, tempVariable,"MAXPMEAN");
			IterReadItem(infile, tempVariable,"MAXSTDURMN");
			IterReadItem(infile, tempVariable,"MAXISTDURMN");
		}
	}
	
	if (optstoch == 6) {
		IterReadItem   (infile, tempString,"WEATHERTABLENAME");
		CheckFileExists(infile, tempString,"WEATHERTABLENAME");
	}
	
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

   // Visualization options
   optv = IterReadItem(infile, optv, "OPTVIZ");
   if (optv > 0) {
       IterReadItem(infile, tempString,"OUTVIZFILENAME");
#ifndef PARALLEL_TRIBS
       CheckPathNameCorrect(infile, tempString, "OUTVIZFILENAME");
#endif
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
void tPreProcess::CheckFileExists(tInputFile &infile,
								  char* filename, const char* keyword) 
{
	char strg[kName];
	char tempString[kName];
	int InpStatus = 0;
	
	while ( !InpStatus ) {
		ifstream readFile(filename);
		if (!readFile) {
			cout<<"\nFile "<<filename<<" for parameter "<<keyword<<" not found..."<<endl;
			cout<<"\nCorrect the .in file and type 'y'"
				<<"\n\n>>";
			cin>>strg;
			infile.CloseOldAndOpenNew(infile.GetInFileName());
			IterReadItem(infile, tempString, keyword);
			strcpy(filename, tempString);
		}
		else
			InpStatus = 1;
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
void tPreProcess::CheckPathNameCorrect(tInputFile &infile, char* filename,
									   const char* keyword)
{
	char strg[kName];
	char bname[kName];
	char tempString[kName];
	char lsl[] = "ls -l |";
	char cat[] = "cat";
	char rm[]  = "rm";
	char larger[] = ">";
	char zero[] = "0";
	
	int InpStatus = 0;
	
	while ( !InpStatus ) {   
		// temporary file output 
		snprintf(bname, sizeof(bname), "%s%s", filename, zero);
		
		ofstream readFile(bname);
		
		if (!readFile) {
			cout<<"\nPathname '"<<filename<<"' for parameter '"<<keyword<<"' is incorrect"<<endl;
			cerr<<"\nCorrect the .in file and type 'y'"
				<<"\n\n>>";
			cin>>strg;
			infile.CloseOldAndOpenNew(infile.GetInFileName());
			IterReadItem(infile,tempString,keyword);
			strcpy(filename, tempString);
		} 
		else {
			InpStatus = 1;
#ifdef ALPHA_64
			remove(bname);
#elif defined LINUX_32
			unlink(bname);
#elif defined WIN
			remove(bname);
#else 
			remove(bname);
#endif
			
		}
	}
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
	char strg[kName];
	int InpStatus = 0;
	
	while ( !InpStatus ) {
		datType = infile.ReadItem(datType, itemCode);
		Cout<<"Parameter =   "<<itemCode<<"\t\t\t"<<datType<<endl;
		
		if (datType < -999000.) {
			cerr<<"\nThe input parameter is either not specified"
			<<"\nor wrong. Correct the .in file and type 'y'"
			<<"\n\n>>";
			cin>>strg;
			infile.CloseOldAndOpenNew(infile.GetInFileName());
		}
		else
			InpStatus = 1;
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
	char strg[kName];
	int InpStatus = 0;
	
	while ( !InpStatus ) {
		datType = infile.ReadItem(datType, itemCode);
		Cout<<"Parameter =   "<<itemCode<<"\t\t\t"<<datType<<endl;
		
		if (datType == -9999) {
			cout<<"\nThe input parameter is either not specified"
			<<"\n\tor wrong. Correct the .in file and type 'y'"
			<<"\n\n>>";
			cin>>strg;
			infile.CloseOldAndOpenNew(infile.GetInFileName());
		}
		else
			InpStatus = 1;
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
	char strg[kName];
	char errr[] = "-999";
	int InpStatus = 0;
	
	while ( !InpStatus ) {
		infile.ReadItem(theString, itemCode);
		Cout<<"Parameter =   "<<itemCode<<"\t\t\t"<<theString<<endl;
		
		if (!strcmp(theString, errr)) {
			cout<<"\nThe input parameter is either not specified"
			<<"\nor wrong. Correct the .in file and type 'y'"
			<<"\n\n>>";
			cin>>strg;
			infile.CloseOldAndOpenNew(infile.GetInFileName());
		}
		else
			InpStatus = 1;
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
