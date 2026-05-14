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
**  tInvariant.cpp:   Functions for tInvariant classes (see tInvariant.h)
**
**
***************************************************************************/

#include "src/Headers/globalIO.h"
#include "src/tRasTin/tInvariant.h"

//=========================================================================
//
//
//                  Section 1: tInvariant GenericSoilData Functions
//
//
//=========================================================================

/***************************************************************************
**
**  GenericSoilData(tInputFile *)
**
**  Constructor and Destructor for GenericSoilData Class
**
***************************************************************************/
GenericSoilData::GenericSoilData(tMesh<tCNode> *mesh, 
								 tInputFile *input, tResample *resamp)
{  
	int id;
	double *sc;
    
    // Changes added by Giuseppe in August 2016 to allow reading soil grids
    int stdgrid_opt = 100; // Assign a random value
    stdgrid_opt = input->ReadItem(stdgrid_opt, "OPTSOILTYPE" );
  		
	// The last one calls tResample functions
	input->ReadItem(soilTable, "SOILTABLENAME"); // input table
	input->ReadItem(soilGrid,  "SOILMAPNAME");   // reads input file
	double *tmp = resamp->doIt(soilGrid, 2);     // resamples grid
	
	tCNode *cn;
	tMeshListIter <tCNode> niter ( mesh->getNodeList() );
	id = 0;
	for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
		cn->setSoilID((int)tmp[id]);  //sets soil ID to tCNode 
		id++;
	}
	
	const int kSoilColumns = 12; // ID + 11 soil properties
	ifstream Inp0(soilTable);
	if (!Inp0) {
		cout <<"\nFile "<<soilTable<<" not found!"<<endl;
		cout <<"Exiting Program...\n\n"<<endl;
		exit(2);
	}

	std::string header;
	std::getline(Inp0, header);
	int colCount = 1;
	for (char c : header)
		if (c == ',') ++colCount;
	if (colCount != kSoilColumns) {
		cout << "\nError: Soil table '" << soilTable << "' has " << colCount
		     << " columns, expected " << kSoilColumns << " (ID + 11 parameters)." << endl;
		cout << "Exiting Program...\n\n" << endl;
		exit(2);
	}

	std::vector<std::vector<double>> rows;
	std::string line;
	while (std::getline(Inp0, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string token;
		std::vector<double> row(kSoilColumns);
		for (int j = 0; j < kSoilColumns; j++) {
			std::getline(ss, token, ',');
			row[j] = std::stod(token);
		}
		rows.push_back(row);
	}
	Inp0.close();

	numClass = static_cast<int>(rows.size());
	assert(numClass > 0);

	int maxSoilID = 0;
	for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP())
		maxSoilID = std::max(maxSoilID, cn->getSoilID());
	if (maxSoilID > numClass) {
		cout << "\nError: Soil raster contains class ID " << maxSoilID
		     << " but soil table '" << soilTable << "' only has "
		     << numClass << " rows." << endl;
		cout << "Exiting Program...\n\n" << endl;
		exit(2);
	}

	SoilClass = new SoilType* [numClass];
	assert(SoilClass != nullptr);

	sc = new double [kSoilColumns];
	assert(sc != nullptr);

	for (int i = 0; i < numClass; i++) {
		for (int j = 0; j < kSoilColumns; j++)
			sc[j] = rows[i][j];
		SoilClass[i] = new SoilType(sc, kSoilColumns);
		assert(SoilClass[i] != nullptr);
	}
	delete [] sc;

    // Giuseppe 2016 - Assign the soil properties to the nodes via the set functions
    id = 0;
    for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
        setSoilPtr( cn->getSoilID() );
        
        cn->setKs(getSoilProp(1)); // Surface hydraulic conductivity
        cn->setThetaS(getSoilProp(2)); // Saturation moisture content
        cn->setThetaR(getSoilProp(3)); // Residual moisture content
        cn->setPoreSize(getSoilProp(4)); // Pore-size distribution index
        cn->setAirEBubPres(getSoilProp(5)); // Air entry bubbling pressure
        cn->setDecayF(getSoilProp(6)); // Decay parameter in the exp
        cn->setSatAnRatio(getSoilProp(7)); // Anisotropy ratio (saturated)
        cn->setUnsatAnRatio(getSoilProp(8)); // Anisotropy ratio (unsaturated)
        cn->setPorosity(getSoilProp(9)); // Porosity
        cn->setVolHeatCond(getSoilProp(10)); // Volumetric Heat Conductivity
        cn->setSoilHeatCap(getSoilProp(11)); // Soil Heat Capacity
        
        id++;
    }
    
    //}                                 // Comment out later Giuseppe 2016
    //else if (stdgrid_opt == 1)        // Comment out later Giuseppe 2016
    if (stdgrid_opt == 1)
    {
        input->ReadItem(scfile, "SCGRID");

        cout<<"\nReading Soil Cover Data Grid File: ";
        cout<< scfile<<"..."<<endl<<flush;

        ifstream readFile(scfile);
        if (!readFile) {
            cout << "\nFile "<<scfile<<" not found!" << endl;
            cout << "Exiting Program...\n\n"<<endl;
            exit(1);
        }

        std::string header;
        std::getline(readFile, header);
        int colCount = 1;
        for (char c : header)
            if (c == ',') ++colCount;
        if (colCount != 3) {
            cout << "\nError: Soil grid file '" << scfile << "' has " << colCount
                 << " columns, expected 3 (Variable,BasePath,FileExtension)." << endl;
            cout << "Exiting Program...\n\n" << endl;
            exit(1);
        }

        std::vector<std::vector<std::string>> rows;
        std::string line;
        while (std::getline(readFile, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::vector<std::string> row(3);
            for (int j = 0; j < 3; j++)
                std::getline(ss, row[j], ',');
            rows.push_back(row);
        }
        readFile.close();

        int numParameters = static_cast<int>(rows.size());
        SCgridBaseNames.resize(numParameters);
        SCgridParamNames.resize(numParameters);
        SCgridExtNames.resize(numParameters);
        SCgridName.resize(numParameters);

        for (int ct=0;ct<numParameters;ct++) {

            const std::string &paramName = rows[ct][0];
            const std::string &baseName  = rows[ct][1];
            const std::string &extName   = rows[ct][2];

            SCgridParamNames[ct] = make_unique<char[]>(paramName.length() + 1);
            strcpy(SCgridParamNames[ct].get(), paramName.c_str());

            if ( (strcmp(SCgridParamNames[ct].get(),"KS")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"TS")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"TR")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"PI")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"PB")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"FD")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"AR")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"UA")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"PO")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"VH")!=0) &&
                (strcmp(SCgridParamNames[ct].get(),"SH")!=0) ) {
                cout << "\nA soil cover parameter name in the SC gdf file is an unexpected one."<<endl;
                cout << "\nExpected variables: KS,TS,TR,PI,PB,FD,AR,UA,PO,VH or SH" << endl;
                cout << "\tCheck and re-run the program" << endl;
                cout << "\nExiting Program..."<<endl<<endl;
                exit(1);
            }

            SCgridBaseNames[ct] = make_unique<char[]>(baseName.length() + 1);
            strcpy(SCgridBaseNames[ct].get(), baseName.c_str());

            if (strcmp(SCgridBaseNames[ct].get(),"NO_DATA")==0) {
                Cout << "\nCannot use NO_DATA for SC Grids"<<endl;
                Cout << "\nExiting Program..."<<endl<<endl;
                exit(1);
            }

            SCgridExtNames[ct] = make_unique<char[]>(extName.length() + 1);
            strcpy(SCgridExtNames[ct].get(), extName.c_str());

            SCgridName[ct] = make_unique<char[]>(100);
            strcpy(SCgridName[ct].get(), baseName.c_str());
            strcat(SCgridName[ct].get(), ".");
            strcat(SCgridName[ct].get(), extName.c_str());
            
            if (strcmp(SCgridParamNames[ct].get(),"KS")==0)
            {
                // Resample Ks grid
                cout << "\nResampling Ks grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                
                //tCNode *cn; /WR Debug both are previously defined
                // tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setKs(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"TS")==0)
            {
                // Resample ThetaS grid
                cout << "\nResampling ThetaS grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                //tCNode *cn;
                //tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setThetaS(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"TR")==0)
            {
                // Resample ThetaR grid
                cout << "\nResampling ThetaR grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                //tCNode *cn;
                //tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setThetaR(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"PI")==0)
            {
                // Resample Pore Index grid
                cout << "\nResampling Pore Index grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022

                //tCNode *cn;
                //tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setPoreSize(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"PB")==0)
            {
                // Resample Air E. Bubbling P grid
                cout << "\nResampling Air Entry Bubbling Pressure grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                //tCNode *cn;
                //tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {

                    cn->setAirEBubPres(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"FD")==0)
            {
                // Resample Decay grid
                cout << "\nResampling Decay Exponent grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                //tCNode *cn;
                // tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setDecayF(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"AR")==0)
            {
                // Resample Decay grid
                cout << "\nResampling Saturated Anisotropy Ratio grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                //tCNode *cn;
                //tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setSatAnRatio(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"UA")==0)
            {
                // Resample Decay grid
                cout << "\nResampling Unsaturated Anisotropy Ratio grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                //tCNode *cn;
                //tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setUnsatAnRatio(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"PO")==0)
            {
                // Resample Porosity grid
                cout << "\nResampling Porosity grid..." << endl;;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                //tCNode *cn;
                //tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setPorosity(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"VH")==0)
            {
                // Resample Volumetric Heat grid
                cout << "\nResampling Volumetric Heat Conductivity grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                //tCNode *cn;
                //tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setVolHeatCond(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
            if (strcmp(SCgridParamNames[ct].get(),"SH")==0)
            {
                // Resample Soil Heat grid
                cout << "\nResampling Soil Heat Capacity grid..." << endl;
                tmp = resamp->doIt(SCgridName[ct].get(), 1);     // resamples grid // Bug fix to replace hardcoded indexing CJC 2022
                //tCNode *cn;
                //tMeshListIter <tCNode> niter ( mesh->getNodeList() );
                id = 0;
                for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
                    cn->setSoilHeatCap(tmp[id]);  //sets soil ID to tCNode
                    id++;
                }
            }
            
        }
        
    }
    
    
}

GenericSoilData::~GenericSoilData() {   
	for (int iprop=0;iprop < numClass;iprop++)
		delete SoilClass[iprop];
	delete [] SoilClass;
}

/***************************************************************************
**
**  GenericSoilData:: Get and Set Soil Ptr
**
**  SetSoilParameters: To set (reset) soil parameter values
**
***************************************************************************/
void GenericSoilData::SetSoilParameters(tMesh<tCNode> *mesh,
										tResample *resamp, tInputFile &infile, int option)
{
	int id;
	
	if ( option ) {     // Needs a new soil resampling 
		Cout<<"\nResampling Soils......"<<endl<<flush;
		infile.ReadItem(soilGrid,  "SOILMAPNAME"); //Reads input file
		double *tmp = resamp->doIt(soilGrid, 2);   //Resamples grid
		
		tCNode *cn;
		tMeshListIter <tCNode> niter ( mesh->getNodeList() );
		id = 0;
		for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP())  {
			cn->setSoilID((int)tmp[id]);  //Sets soil ID to tCNode 
			id++;
		}
	}
	
	infile.ReadItem(soilTable, "SOILTABLENAME"); // input table
	const int kSoilColumns = 12; // ID + 11 soil properties
	ifstream Inp0(soilTable);
	if (!Inp0) {
		cout <<"File "<<soilTable<<" not found!"<<endl;
		cout <<"Exiting Program...\n\n"<<endl;
		exit(2);
	}

	std::string header;
	std::getline(Inp0, header);
	int colCount = 1;
	for (char c : header)
		if (c == ',') ++colCount;
	if (colCount != kSoilColumns) {
		cout << "\nError: Soil table '" << soilTable << "' has " << colCount
		     << " columns, expected " << kSoilColumns << " (ID + 11 parameters)." << endl;
		cout << "Exiting Program...\n\n" << endl;
		exit(2);
	}

	std::vector<std::vector<double>> rows;
	std::string line;
	while (std::getline(Inp0, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string token;
		std::vector<double> row(kSoilColumns);
		for (int j = 0; j < kSoilColumns; j++) {
			std::getline(ss, token, ',');
			row[j] = std::stod(token);
		}
		rows.push_back(row);
	}
	Inp0.close();

	int rowCount = static_cast<int>(rows.size());

	if ( option ) {         // Needs a new soil resampling
		for (int iprop=0;iprop < numClass;iprop++)
			delete SoilClass[iprop];
		delete [] SoilClass;

		numClass = rowCount;
		assert(numClass > 0);

		SoilClass = new SoilType* [numClass];
		assert(SoilClass != 0);

		double *sc = new double [kSoilColumns];
		assert(sc != 0);

		for (int i=0; i < numClass; i++) {
			for (int j=0; j < kSoilColumns; j++)
				sc[j] = rows[i][j];
			SoilClass[i] = new SoilType(sc, kSoilColumns);
			assert(SoilClass[i] != 0);
		}
		delete [] sc;
	}
	else if ( !option ) {
		if (rowCount != numClass) {
			cout<<"\nWarning! Number of classes does not correspond to previous";
			cout <<"\nProceeding with latter number of classes"<<endl<<flush;
		}
		int readCount = std::min(numClass, rowCount);
		for (int i=0; i < readCount; i++) {
			for (int j=0; j < kSoilColumns; j++)
				SoilClass[i]->setProperty( j, rows[i][j] );
		}
	}
	}

void GenericSoilData::printSoilPars() 
{
	cerr<<numClass<<" "<<SoilClass[0]->numProps<<endl;
	for (int i=0; i < numClass; i++) {
		for (int j=0; j < SoilClass[i]->numProps; j++)
			cerr<<SoilClass[i]->sProperty[j]<<" ";
		cerr<<endl;
	}
	return;
}

void GenericSoilData::setSoilPtr(int soilID) 
{
	if (soilID < 0 || soilID > numClass) { 
		cout <<"\nError: In setSoilPtr: soilclass > numClass or < 0"<<endl;
		cout <<"SoilClass = "<<soilID<<endl;
		cout <<"Exiting Program...\n\n"<<endl;
		exit(10);
	}
	currClass = soilID - 1;   // '-1' In order to make the corresponding 
							  //      indices start from '0'
	return;
}

double GenericSoilData::getSoilProp(int propID) 
{
	if (propID < 0 || propID > SoilClass[currClass]->numProps) { 
		cout <<"\nError: In getSoilProp: propID > numProps or < 0"<<endl;
		cout <<"PropID = "<<propID<<endl;
		cout <<"Exiting Program...\n\n"<<endl;
		exit(10);
	}
	return(SoilClass[currClass]->getProperty(propID)); 
}

//=========================================================================
//
//
//                  Section 2: tInvariant SoilType Functions
//
//
//=========================================================================

/***************************************************************************
**
**  SoilType:: Constructor and Destructor
**
***************************************************************************/
SoilType::SoilType(double *a, int n) 
{
	numProps = n;
	assert(n > 0);
	sProperty = new double [n];
	assert(sProperty != 0);
	for (int i=0; i < numProps; i++)
		sProperty[i] = a[i];
}

SoilType::SoilType() 
{
	sProperty = nullptr;
	numProps  = -999;
}

SoilType::~SoilType() 
{
	if (sProperty) 
		delete [] sProperty;
	numProps = 0;
}

/***************************************************************************
**
**  SoilType:: get and set Property
**
***************************************************************************/

double SoilType::getProperty(int id) {
	return (sProperty[id]);
}

void SoilType::setProperty(int id, double param) {
	sProperty[id] = param;
}


//=========================================================================
//
//
//                  Section 3: tInvariant GenericLandData Functions
//
//
//=========================================================================

/***************************************************************************
**
**  GenericLandData(tInputFile *)
**
**  Constructor and Destructor for GenericLandData Class
**
***************************************************************************/
GenericLandData::GenericLandData(tMesh<tCNode> *mesh,
								 tInputFile *input, tResample *resamp)
{
	int id;
	double *lc;
	
	input->ReadItem(landTable, "LANDTABLENAME");   // input table
	input->ReadItem(landGrid,  "LANDMAPNAME");     // reads input file
	double *tmp = resamp->doIt(landGrid, 2);       // resamples grid
	
	tCNode *cn;
	tMeshListIter <tCNode> niter ( mesh->getNodeList() );
	id = 0;
	for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP())  {
		cn->setLandUse((int)tmp[id]);  //<- sets land use ID to tCNode 
		id++;
	}
	
	const int kLandColumns = 13; // ID + 12 land use properties
	ifstream Inp0(landTable);
	if (!Inp0) {
		cout <<"\nFile "<<landTable<<" not found!"<<endl;
		cout <<"Exiting Program...\n\n"<<endl;
		exit(2);
	}

	std::string header;
	std::getline(Inp0, header);
	int colCount = 1;
	for (char c : header)
		if (c == ',') ++colCount;
	if (colCount != kLandColumns) {
		cout << "\nError: Land use table '" << landTable << "' has " << colCount
		     << " columns, expected " << kLandColumns << " (ID + 12 parameters)." << endl;
		cout << "Exiting Program...\n\n" << endl;
		exit(2);
	}

	std::vector<std::vector<double>> rows;
	std::string line;
	while (std::getline(Inp0, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string token;
		std::vector<double> row(kLandColumns);
		for (int j = 0; j < kLandColumns; j++) {
			std::getline(ss, token, ',');
			row[j] = std::stod(token);
		}
		rows.push_back(row);
	}
	Inp0.close();

	numClass = static_cast<int>(rows.size());
	assert(numClass > 0);

	int maxLandID = 0;
	for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP())
		maxLandID = std::max(maxLandID, cn->getLandUse());
	if (maxLandID > numClass) {
		cout << "\nError: Land use raster contains class ID " << maxLandID
		     << " but land use table '" << landTable << "' only has "
		     << numClass << " rows." << endl;
		cout << "Exiting Program...\n\n" << endl;
		exit(2);
	}

	LandClass = new LandType* [numClass];
	assert(LandClass != 0);

	lc = new double [kLandColumns];
	assert(lc != 0);

	for (int i = 0; i < numClass; i++) {
		for (int j = 0; j < kLandColumns; j++)
			lc[j] = rows[i][j];
		LandClass[i] = new LandType(lc, kLandColumns);
		assert(LandClass[i] != 0);
	}
	delete [] lc;
}

GenericLandData::~GenericLandData() {   
	for (int iprop=0;iprop < numClass;iprop++)
		delete LandClass[iprop];
	delete [] LandClass;
}

/***************************************************************************
**
**  GenericLandData:: Get and Set Land Ptr
**
**  SetLandParameters: To set (reset) land parameter values
**
***************************************************************************/
void GenericLandData::SetLtypeParameters(tMesh<tCNode> *mesh,
										 tResample *resamp, tInputFile &infile, int option)
{
	int id;
	
	if ( option ) {              // Needs a new soil resampling 
		Cout<<"\nResampling Landuse......"<<endl<<flush;
		infile.ReadItem(landGrid,  "LANDMAPNAME");  // Reads input file
		double *tmp = resamp->doIt(landGrid, 2);    // Resamples grid
		
		tCNode *cn;
		tMeshListIter <tCNode> niter ( mesh->getNodeList() );
		id = 0;
		for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP())  {
			cn->setLandUse((int)tmp[id]);  // Sets land use ID to tCNode 
			id++;
		}
	}
	
	infile.ReadItem(landTable, "LANDTABLENAME"); // Input table
	const int kLandColumns = 13; // ID + 12 land use properties
	ifstream Inp0(landTable);
	if (!Inp0) {
		cout <<"File "<<landTable<<" not found!!!"<<endl;
		cout<<", tInvariant Error"<<endl;
		exit(2);
	}

	std::string header;
	std::getline(Inp0, header);
	int colCount = 1;
	for (char c : header)
		if (c == ',') ++colCount;
	if (colCount != kLandColumns) {
		cout << "\nError: Land use table '" << landTable << "' has " << colCount
		     << " columns, expected " << kLandColumns << " (ID + 12 parameters)." << endl;
		cout << "Exiting Program...\n\n" << endl;
		exit(2);
	}

	std::vector<std::vector<double>> rows;
	std::string line;
	while (std::getline(Inp0, line)) {
		if (line.empty()) continue;
		std::istringstream ss(line);
		std::string token;
		std::vector<double> row(kLandColumns);
		for (int j = 0; j < kLandColumns; j++) {
			std::getline(ss, token, ',');
			row[j] = std::stod(token);
		}
		rows.push_back(row);
	}
	Inp0.close();

	int rowCount = static_cast<int>(rows.size());

	if ( option ) { // Needs a new landuse resampling
		for (int iprop=0;iprop < numClass;iprop++)
			delete LandClass[iprop];
		delete [] LandClass;

		numClass = rowCount;
		assert(numClass > 0);
		LandClass = new LandType* [numClass];
		assert(LandClass != 0);

		double *lc = new double [kLandColumns];
		assert(lc != 0);

		for (int i=0; i < numClass; i++) {
			for (int j=0; j < kLandColumns; j++)
				lc[j] = rows[i][j];
			LandClass[i] = new LandType(lc, kLandColumns);
			assert(LandClass[i] != 0);
		}
		delete [] lc;
	}
	else if ( !option ) {
		if (rowCount != numClass) {
			cout<<"\nWarning! Number of classes does not correspond to previous";
			cout <<"\nProceeding with latter number of classes"<<endl<<flush;
		}
		int readCount = std::min(numClass, rowCount);
		for (int i=0; i < readCount; i++) {
			for (int j=0; j < kLandColumns; j++)
				LandClass[i]->setProperty( j, rows[i][j] );
		}
	}
	return;
}

void GenericLandData::setLandPtr(int landID) 
{
	if (landID < 0 || landID > numClass) { 
		cout <<"\nError: In setLandlPtr: landclass > numClass or < 0"<<endl;
		cout <<"LandClass = "<<landID<<endl;
		cout <<"Exiting Program...\n\n"<<endl;
		exit(10);
	}
	currClass = landID - 1;  // '-1' In order to make the corresponding 
							 //      indices start from '0'
	return;
}

double GenericLandData::getLandProp(int propID) {
	//if (propID < 0 || propID > LandClass[currClass]->numProps) { 
	if (propID < 0 || propID >= LandClass[currClass]->numProps) { // GMnSKY2008MLE
		cout <<"\nError: In getLandProp: propID > numProps or < 0"<<endl;
		cout <<"PropID = "<<propID<<endl;
		cout <<"Exiting Program...\n\n"<<endl;
		exit(10);
	}
	return(LandClass[currClass]->getProperty(propID)); 
}

//=========================================================================
//
//
//                  Section 4: tInvariant LandType Functions
//
//
//=========================================================================

/***************************************************************************
**
**  LandType:: Constructor and Destructor
**
***************************************************************************/
LandType::LandType(double *a, int n) 
{
	numProps  = n;
	assert(n > 0);
	lProperty = new double [n];
	assert(lProperty != nullptr);
	
	for (int i=0; i < numProps; i++)
		lProperty[i] = a[i];
}

LandType::LandType() 
{
	lProperty = nullptr;
	numProps  = -999;
}

LandType::~LandType() 
{
	if (lProperty) 
		delete [] lProperty;
	numProps = 0;
}

/***************************************************************************
**
**  LandType:: Set and Get Property
**
***************************************************************************/

double LandType::getProperty(int id) {
	if ((id < 0) || ( id >= numProps)) {
	cout << "LandType getProperty call warning: id is" << id <<endl;
	}
	return (lProperty[id]);
}

void LandType::setProperty(int id, double param) {
	lProperty[id] = param;
}

//=========================================================================
//
//
//                      End of tInvariant.cpp
//
//
//=========================================================================

