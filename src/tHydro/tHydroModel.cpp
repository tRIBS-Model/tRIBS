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
**  tHydroModel.cpp:   Function file for tHydroModel Class (see tHydroModel.h)
**
***************************************************************************/

#include "src/tHydro/tHydroModel.h"
#include "src/Headers/globalIO.h"

#ifdef PARALLEL_TRIBS
#include "src/tParallel/tParallel.h"
#include "src/tGraph/tGraph.h"
#endif

//=========================================================================
//
//
//                  Section 1: tHydroModel Constructors/Destructors
//
//
//=========================================================================

tHydroModel::tHydroModel(SimulationControl *simCtrPtr, tMesh<tCNode> *gptr,
						 tInputFile &infile, tResample *resamp,
						 tWaterBalance *balance, tRunTimer *t)
{
	// Initialization
	NwtOld = NwtNew = 0.;
	MuOld = MuNew = 0.;
	MiOld = MiNew = 0.;
	NfOld = NfNew = 0.;
	NtOld = NtNew = 0.;
	RuOld = RuNew = 0.;
	RiOld = RiNew = 0.;
	nodeList = nullptr;

	gridPtr = gptr;
	simCtrl = simCtrPtr;
	balPtr = balance;
	timer = t;

	// Soil and Land Use Objects
	Cout<<"\nResampling soil and land use grids..."<<endl;
	soilPtr = new GenericSoilData(gptr, &infile, resamp);
	landPtr = new GenericLandData(gptr, &infile, resamp);
	SetHydroMVariables( infile, resamp, 0);

}

/*************************************************************************
**
**  tHydroModel::SetHydroMVariables(tInputFile &, tResample *)
**
**  Sets up basic data members for tHydroModel object
**
*************************************************************************/
void tHydroModel::SetHydroMVariables(tInputFile &infile,
									 tResample *resamp, int keep)
{
	char nodeFile[kMaxNameSize]; //Base Name of the file w/ node IDs
	char baseName[kMaxNameSize];
	char currName[kMaxNameSize+20];

	IntStormMAX = infile.ReadItem(IntStormMAX, "INTSTORMMAX" );
	EToption    = infile.ReadItem(EToption, "OPTEVAPOTRANS");
	Ioption     = infile.ReadItem(Ioption, "OPTINTERCEPT");
	gFluxOption = infile.ReadItem(gFluxOption, "GFLUXOPTION");
	BRoption    = infile.ReadItem(BRoption, "OPTBEDROCK");
	percolationOption = infile.ReadItem(percolationOption, "OPTPERCOLATION");

	if (infile.IsItemIn( "OPTGWFILE" ))
		GWoption = infile.ReadItem(GWoption, "OPTGWFILE");
	else
		GWoption = 0; //Default option

	if (infile.IsItemIn( "OPTRUNON" ))
		RunOnoption = infile.ReadItem(RunOnoption, "OPTRUNON");
	else
		RunOnoption = 0; //Default option

	// SKY2008Snow from AJR2007
	SnOpt = infile.ReadItem(SnOpt, "OPTSNOW");

	infile.ReadItem(gwatfile, "GWATERFILE" );
	infile.ReadItem(nodeFile, "HYDRONODELIST");

	RdstrOption = 0;

	// Bedrock optional input in meters
	if (!BRoption) {
		DtoBedrock = infile.ReadItem(DtoBedrock, "DEPTHTOBEDROCK");
		DtoBedrock = DtoBedrock*1000.0; //convert to mm
	}
	else if (BRoption == 1) {
		infile.ReadItem(bedrockfile, "BEDROCKFILE");
	}
	else {
		cout<<"\nWarning: Only two options available for bedrock depth:"<<endl;
		cout<<"\t(0) Read in uniform value from DEPTHTOBEDROCK"<<endl;
		cout<<"\t(1) Read in ASCII grid values from BEDROCKFILE"<<endl;
		cout<<"Your option is "<<BRoption<<"\tAssumed UNIFORM case..."<<endl;
		DtoBedrock = infile.ReadItem(DtoBedrock, "DEPTHTOBEDROCK");
		DtoBedrock = DtoBedrock*1000.0; //convert to mm
	}
	
	// Read the soil depths. If not in the file, use the original defaults.
	if (infile.IsItemIn( "SURFACESOILDEPTH" )) {
		surfaceSoilDepth = infile.ReadItem(surfaceSoilDepth, "SURFACESOILDEPTH");
	} else {
		surfaceSoilDepth = 100.0; // Default to 100mm
	}
	// If a decision made to keep the state vars don't change anything,
	// keep vars from previous run. Otherwise, re-initialize everything
	if ( !keep )
		InitSet(resamp);
	else
		InitIntegralVars();

	if (nodeList != nullptr)
		delete [] nodeList;
	SetHydroNodes(nodeFile);

	Stok = 0.;
	TotRain = 0.;
	TotGWchange = 0.;
	TotMoist = 0.;
}


tHydroModel::~tHydroModel()
{
	gridPtr = nullptr;
	delete soilPtr;
	delete landPtr;
	if (nodeList != nullptr)
		delete [] nodeList;
	Cout <<"tHydroModel Object has been destroyed..."<<endl;
}

/*************************************************************************
**
**  tHydroModel::HydroNodesExist()
**
**  If there is an active hydronode list, returns TRUE, FALSE otherwise
**
*************************************************************************/
int tHydroModel::HydroNodesExist()
{
	if (nodeList)
		return 1;
	else
		return 0;
}

/*************************************************************************
**
**  tHydroModel::SetHydroNodes(char *nodeFile)
**
**  Reads in the nodes IDs for which output at every time step will be
**  produced through the use of the HYDRONODELIST keyword
**
*************************************************************************/
void tHydroModel::SetHydroNodes(char *nodeFile)
{
	ifstream Inp(nodeFile);
	if (!Inp) {
		Cout <<"\nHydroNodeList File "<<nodeFile<<" not found."<<endl;
		Cout<<"\tWarning: The file with hydro node IDs does not "
			<<"exist. No output will be written."<<endl;
		return;
	}
	int cnt = 0;
	Inp>>numNodes;
	nodeList = new int[numNodes];
	for (int i = 0; i < numNodes; i++) {
		Inp>>nodeList[i];
		if (nodeList[i] >= 0 &&
			nodeList[i] <= gridPtr->getNodeList()->getActiveSize())
			cnt++;
	}
	Inp.close();

	if (!cnt) {
		delete [] nodeList;
		nodeList = nullptr;
	}
	else
		Cout <<"\nHydroNodeList has been set up"<<endl;
}

//=========================================================================
//
//
//                  Section 2: tHydroModel InitSet Function
//
//
//=========================================================================

/*************************************************************************
**
**  tHydroModel::InitSet(tResample *resamp)
**
**  Initializes nodes at the beginning of simulation, reads GW file
**
*************************************************************************/
void tHydroModel::InitSet(tResample *resamp)
{
	int id = 0;
	char wish = 'Z';
	char keep = 'Z';
	char filein[80];
	char yesno[20];
	double *tmp;
	double *bedrock;
	double bedRock = 0.0;

	tCNode * cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );

	// Print out options used for Hydrologic Processes
	Cout<<"\nHydrologic Options: "<<endl;
	Cout<<"\nEvapotranspiration Option: \t"<<EToption<<endl;
	Cout<<"Interception Option: \t\t"<< Ioption<<endl;
	Cout<<"Ground Heat Flux Option: \t"<< gFluxOption<<endl;
	Cout<<"Bedrock Depth Option: \t\t" << BRoption<<endl;


	// Groundwater initial file option for g
	if (!GWoption) {
		// Resample ASCII grid
		Cout<<"\nResampling groundwater table grid..."<<endl;
		tmp = resamp->doIt(gwatfile, 1); // resamples input GW ASCII grid
	}
	else if (GWoption == 1) {



		ifstream source( gwatfile );
        if (!source.is_open()) {
            cerr << "Failed to open file: " << gwatfile << std::endl;
            exit(EXIT_FAILURE);
        }

		int Vcnt = (gridPtr->getNodeList()->getActiveSize());
		tmp = new double[Vcnt];
		assert(tmp != nullptr);
		for (int i=0; i < Vcnt; i++) {
			source>>ID>>NwtNew;
			tmp[i] = NwtNew;
		}


		source.close();
	}
	else {
		Cout<<"\nWarning: Only two options available for GW initial file:"<<endl;
		Cout<<"\t(0) Read in ASCII grid from GWATERFILE"<<endl;
		Cout<<"\t(1) Read in values from specified Voronoi File"<<endl;
		Cout<<"Your option is "<<GWoption<<"\tAssumed gridded case..."<<endl;
		Cout<<"\nResampling groundwater table grid...\n"<<endl;
		tmp = resamp->doIt(gwatfile, 1); // resamples input GW
	}


	//Assign Water Table to tCNode
	int id2 = 0;
	for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() ){
		NwtNew = tmp[id2];              //Initial GW in mm

		cn->setNwtOld(NwtNew);
		cn->setNwtNew(NwtNew);
		id2++;
	}

	// Resamples input BR ASCII grid
	if (BRoption == 1) {
		Cout<<"\nResampling depth to bedrock grid..."<<endl;
		bedrock = resamp->doIt(bedrockfile,1);
	}

	// Calculate Basin Area
	BasArea = 0.0;
	for (cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP())
		BasArea += cn->getVArea();

	// Loop through nodes to get Initial GW conditions
	id = 0;
	for (cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP())
	{
        // Giuseppe 2016 - Begin changes to allow reading soil properties from grids
        //soilPtr->setSoilPtr( cn->getSoilID() );
        //        Ksat    = soilPtr->getSoilProp(1);  // Surface hydraulic conductivity
        //        Ths     = soilPtr->getSoilProp(2);  // Saturation moisture content
        //        Thr     = soilPtr->getSoilProp(3);  // Residual moisture content
        //        PoreInd = soilPtr->getSoilProp(4);  // Pore-size distribution index
        //        Psib    = soilPtr->getSoilProp(5);  // Air entry bubbling pressure
        //        F       = soilPtr->getSoilProp(6);  // Decay parameter in the exp
        //        Ar      = soilPtr->getSoilProp(7);  // Anisotropy ratio (saturated)
        //        UAr     = soilPtr->getSoilProp(8);  // Anisotropy ratio (unsaturated)
        //        porosity = soilPtr->getSoilProp(9); // Porosity
        Ksat = cn->getKs();  // Surface hydraulic conductivity
        Ths = cn->getThetaS(); // Saturation moisture content
        Thr = cn->getThetaR(); // Residual moisture content
        PoreInd = cn->getPoreSize(); // Pore-size distribution index
        Psib = cn->getAirEBubPres(); // Air entry bubbling pressure

        if (Psib >= 0) {
            std::cerr << "Error: Psib is zero or positive!" << std::endl;
            exit(1);
        }

        F = cn->getDecayF(); // Decay parameter in the exp
        Ar = cn->getSatAnRatio(); // Anisotropy ratio (saturated)
        UAr = cn->getUnsatAnRatio(); // Anisotropy ratio (unsaturated)
        porosity = cn->getPorosity(); // Porosity
        // Giuseppe 2016 - End changes to allow reading soil properties from grids
		
        Eps = 3 + 2/PoreInd;

		NwtNew = cn->getNwtOld();              //Initial GW in mm

		if (BRoption == 1)
			bedRock = bedrock[id]*1000.;     //Convert meters into mm
		else if (!BRoption)
			bedRock = DtoBedrock;

		//Temporary Lines added by Ricardo Mantilla
		//cout<<"bedRockBefore: id: "<<id<<" depth: "<<bedRock<<endl;
		//bedRock=(cn->getZ()-minElev)*1000;
		//cout<<"bedRockAfter: id: "<<id<<" depth: "<<bedRock<<endl;

		if (NwtNew > bedRock)
			NwtNew = bedRock - 1.0;

		if (NwtNew < 0.0) {
			NwtNew = bedRock - 1.0;
			cout<<"\nWarning: Resampled WT depth < 0 : assigned to "
				<<"DEPTHTOBEDROCK. Check boundary polygons, "
				<<"NODE ID = "<<id<<endl<<flush;
		}

		if (NwtNew >= fabs(Psib))
			MiNew = get_Total_Moist(NwtNew);
		else if ( NwtNew >= 0.0 && NwtNew < fabs(Psib) ) {
			MiNew  = 0.;
			NwtNew = 0.;
		}
		else
			cout<<"\nInitialization: Warning: Nwt < 0: Nwt = "<<NwtNew<<"\n\n";

		// Set Initial Variables
		cn->setNwtOld(NwtNew);      //Unsaturated zone Members
		cn->setNwtNew(NwtNew);
		cn->setMuOld(MiNew);
		cn->setMiOld(MiNew);
		cn->setMuNew(MiNew);
		cn->setMiNew(MiNew);

		cn->setNfOld(0.0);
		cn->setNfNew(0.0);
		cn->setNtOld(0.0);
		cn->setNtNew(0.0);
		cn->setQpin(0.0);
		cn->setQpout(0.0);
		cn->setRiOld(0.0);
		cn->setRiNew(0.0);
		cn->setRuOld(0.0);
		cn->setRuNew(0.0);

		cn->setGwaterChng(0.0);
		cn->setSrf_Hr(0.0);
		cn->setsrf(0.0);
		cn->setsbsrf(0.0);
		cn->sethsrf(0.0);
		cn->setpsrf(0.0);
		cn->setesrf(0.0);

		cn->setInterceptLoss(0.0);       //Interception Members
		cn->setNetPrecipitation(0.0);
		cn->setCanStorage(0.0);
		cn->setPotEvap(0.0);
		cn->setActEvap(0.0);
		cn->setStormLength( 0, 0 );
		cn->setCumIntercept(0.0);
		cn->setEvapoTrans(0.0);          //Evapotranspiration Members
		cn->setEvapWetCanopy(0.0);
		cn->setEvapDryCanopy(0.0);
		cn->setEvapSoil(0.0);
		cn->setSoilMoisture(0.0);
		cn->setSoilMoistureSC(0.0);
		cn->setSoilMoistureUNSC(0.0);
		cn->setRootMoisture(0.0);
		cn->setRootMoistureSC(0.0);
		cn->setTransmiss(0.0);
		cn->setAirPressure(0.0);         //Meteorological Members
		cn->setDewTemp(0.0);
		cn->setRelHumid(0.0);
		cn->setSkyCover(0.0);
		cn->setWindSpeed(0.0);
		cn->setAirTemp(0.0);
		cn->setSurfTemp(0.0);
		cn->setSoilTemp(0.0);
		cn->setNetRad(0.0);
		cn->setGridET(0.0);
		cn->setShortRadIn(0.0);
		cn->setLongRadIn(0.0);
		cn->setLongRadOut(0.0);
		cn->setGFlux(0.0);
		cn->setHFlux(0.0);
		cn->setLFlux(0.0);

		cn->setCanopyStorVol(0.0); // SKYnGM2008LU: tWaterBalance Member

		// SKYnGM2008LU: Land Use Members (for both static and dynamic)
		landPtr->setLandPtr( cn->getLandUse() );
		P_LU    = landPtr->getLandProp(1);
		S_LU    = landPtr->getLandProp(2);
		K_LU    = landPtr->getLandProp(3);
		b2_LU   = landPtr->getLandProp(4);
		Al_LU   = landPtr->getLandProp(5);
		h_LU    = landPtr->getLandProp(6);
		Kt_LU   = landPtr->getLandProp(7);
		Rs_LU   = landPtr->getLandProp(8);
		V_LU    = landPtr->getLandProp(9);
		LAI_LU  = landPtr->getLandProp(10);
		// CJC2025 Stress Thresholds
		SE_LU   = landPtr->getLandProp(11);
		ST_LU   = landPtr->getLandProp(12);
		RZ_LU   = landPtr->getLandProp(13);

		// Handle Rootzone Depth Input
		if (RZ_LU >= 9999.99) { RZ_LU = 1000.0; } else { RZ_LU *= 1000.0; } // If user's rootzone depth is no data (9999.99) set to 1m default
		RZ_LU   = std::max(RZ_LU, 100.0);
		RZ_LU   = std::min(RZ_LU, bedRock);

		cn->setRootZoneDepth(RZ_LU);
		cn->setRootZoneDepthInPrevGrid(RZ_LU);
		cn->setRootZoneDepthInUntilGrid(RZ_LU);
		cn->setLandUseAlb(Al_LU);
		cn->setLandUseAlbInPrevGrid(Al_LU);
		cn->setLandUseAlbInUntilGrid(Al_LU);
		cn->setThroughFall(P_LU);
		cn->setThroughFallInPrevGrid(P_LU);
		cn->setThroughFallInUntilGrid(P_LU);
		cn->setVegHeight(h_LU);
		cn->setVegHeightInPrevGrid(h_LU);
		cn->setVegHeightInUntilGrid(h_LU);
		cn->setStomRes(Rs_LU);
		cn->setStomResInPrevGrid(Rs_LU);
		cn->setStomResInUntilGrid(Rs_LU);
		cn->setVegFraction(V_LU);
		cn->setVegFractionInPrevGrid(V_LU);
		cn->setVegFractionInUntilGrid(V_LU);
		cn->setCanFieldCap(S_LU);
		cn->setCanFieldCapInPrevGrid(S_LU);
		cn->setCanFieldCapInUntilGrid(S_LU);
		cn->setDrainCoeff(K_LU);
		cn->setDrainCoeffInPrevGrid(K_LU);
		cn->setDrainCoeffInUntilGrid(K_LU);
		cn->setDrainExpPar(b2_LU);
		cn->setDrainExpParInPrevGrid(b2_LU);
		cn->setDrainExpParInUntilGrid(b2_LU);
		cn->setOptTransmCoeff(Kt_LU);
		cn->setOptTransmCoeffInPrevGrid(Kt_LU);
		cn->setOptTransmCoeffInUntilGrid(Kt_LU);
		cn->setLeafAI(LAI_LU);
		cn->setLeafAIInPrevGrid(LAI_LU);
		cn->setLeafAIInUntilGrid(LAI_LU);
		// CJC2025 Stress Thresholds
		cn->setEvapThresh(SE_LU);
		cn->setEvapThreshInPrevGrid(SE_LU);
		cn->setEvapThreshInUntilGrid(SE_LU);
		cn->setTransThresh(ST_LU);
		cn->setTransThreshInPrevGrid(ST_LU);
		cn->setTransThreshInUntilGrid(ST_LU);

		// SKY2008Snow
		cn->setLandFact(0.0);

		cn->setIntStormVar(0.0);
		cn->setBedrockDepth(bedRock);
		cn->setBasinArea(BasArea);

		cn->setHlevel(0.0);              //Routing Members
		cn->setQstrm(0.0);
		cn->setFlowVelocity(0.0);

		cn->setSoilMoisture(ComputeSurfSoilMoist(surfaceSoilDepth));
		cn->setRootMoisture(ComputeSurfSoilMoist(RZ_LU));
		cn->setSoilMoistureSC(cn->getSoilMoisture()/Ths);
		cn->setRootMoistureSC(cn->getRootMoisture()/Ths);

        cn->setSoilCutoff(get_Upper_Moist(surfaceSoilDepth,bedRock)/surfaceSoilDepth); // volumetric soil moisture content of soil zone if water table reaches bedrock
        cn->setRootCutoff(get_Upper_Moist(RZ_LU,bedRock)/RZ_LU); // volumetric soil moisture content of root zone if water table reaches bedrock
        cn->setRootCutoffDepth(RZ_LU); // cache the rootzone depth the cutoff was computed for



		if (NwtNew)
			cn->setSoilMoistureUNSC(cn->getMuNew()/NwtNew/Ths);
		else
			cn->setSoilMoistureUNSC(1.0);

		// Set Initial Water Balance Storage
		cn->setUnSaturatedStorage(MiNew/1000*cn->getVArea()); // WR 12192023: convert to m3
		cn->setSaturatedStorage(Ths*(bedRock-NwtNew)*cn->getVArea()/1000.0);

		// Make sure that DtoBedrock assigned proper value
		DtoBedrock = bedRock;

		id++;
	}

	// Deallocate memory if necessary
	if (wish == 'y')
		delete [] tmp;

	// Initialize non-state variables (cumulative or integrated)
	InitIntegralVars();

	// Variables to compute relative input from each factor - Viva (4/25/2004)
	fSoi100=fTop100=fClm100=fGW100=dM100=dMRt=mTh100=mThRt=0.0;
}

/*************************************************************************
**
**  tHydroModel::InitIntegralVars()
**
**  Set integral variables to zero for new run
**
*************************************************************************/
void tHydroModel::InitIntegralVars()
{
	tCNode *cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );
	for ( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() ) {

		// SKY2008Snow from AJR2007
        // Giuseppe 2016 - Begin changes to allow reading soil properties from grids
        //soilPtr->setSoilPtr( cn->getSoilID() );
        //Ths = soilPtr->getSoilProp(2);
        Ths = cn->getThetaS(); // Saturation moisture content
        // Giuseppe 2016 - End changes to allow reading soil properties from grids

		cn->satOccur    = 0;
		cn->hsrfOccur   = 0.0;
		cn->sbsrfOccur  = 0.0;
		cn->psrfOccur   = 0.0;
		cn->satsrfOccur = 0.0;
		cn->RechDisch   = 0.0;

		cn->setAvSoilMoisture(floor(cn->getSoilMoistureSC()*1.0E+4));
		cn->addAvSoilMoisture(cn->getRootMoistureSC()*1.0E-1);

		cn->setAvEvapFract(0.0);
		cn->setAvET(0.0);

		// SKYnGM2008LU: Land Use Members (for both static and dynamic)
		landPtr->setLandPtr( cn->getLandUse() );
		P_LU    = landPtr->getLandProp(1);
		S_LU    = landPtr->getLandProp(2);
		K_LU    = landPtr->getLandProp(3);
		b2_LU   = landPtr->getLandProp(4);
		Al_LU   = landPtr->getLandProp(5);
		h_LU    = landPtr->getLandProp(6);
		Kt_LU   = landPtr->getLandProp(7);
		Rs_LU   = landPtr->getLandProp(8);
		V_LU    = landPtr->getLandProp(9);
		LAI_LU  = landPtr->getLandProp(10);
		// CJC2025 Stress Thresholds
		SE_LU  = landPtr->getLandProp(11);
		ST_LU  = landPtr->getLandProp(12);
		RZ_LU  = landPtr->getLandProp(13);

		// Handle Rootzone Depth Input
		if (RZ_LU >= 9999.99) { RZ_LU = 1000.0; } else { RZ_LU *= 1000.0; } // If user's rootzone depth is no data (9999.99) set to 1m default
		RZ_LU  = std::max(RZ_LU, 100.0);
		RZ_LU  = std::min(RZ_LU, cn->getBedrockDepth());

		cn->setAvRootZoneDepth(RZ_LU);
		cn->setAvThroughFall(P_LU);
		cn->setAvCanFieldCap(S_LU);
		cn->setAvDrainCoeff(K_LU);
		cn->setAvDrainExpPar(b2_LU);
		cn->setAvLandUseAlb(Al_LU);
		cn->setAvVegHeight(h_LU);
		cn->setAvOptTransmCoeff(Kt_LU);
		cn->setAvStomRes(Rs_LU);
		cn->setAvVegFraction(V_LU);
		cn->setAvLeafAI(LAI_LU);
		cn->setAvEvapThresh(SE_LU);
		cn->setAvTransThresh(ST_LU);

		CheckMoistureContent( cn );
	}
}

/*************************************************************************
**
**  tHydroModel::CheckMoistureContent(tCNode *cn)
**
**  Checks if initial moisture content corresponds to a correct value
**
*************************************************************************/
void tHydroModel::CheckMoistureContent(tCNode *cn)
{
    double Mdelt;
	// An additional check is made for elements
	// that behave in some strange fashion
	SetupNodeUSZ( cn );

	Mdelt = get_Total_Moist(NwtOld);
	if (fabs(MiOld-Mdelt) >= 1.0E-3) {
		if (simCtrl->Verbose_label == 'Y') {
            double dM;

			cout <<"Warning! Initial. moist. imb., ID = "<<ID
			<<"; MiOld = "<<MiOld<<"; MActual = "<<Mdelt<<endl;
			cout << "Node #" << cn->getID()<<" ("<<cn->getX()<<","
				<< cn->getY() << ")"<<endl<<flush;
			cout << "-------------------------" <<'\n';
			cout <<"Ksat = "<<Ksat
				<<"; Ths = "<<Ths
				<<"; Thr = "<<Thr
				<<";\nPoreInd = "<<PoreInd
				<<"; Eps = "<<Eps
				<<"; Psib = "<<Psib
				<<";\nF = "<<F
				<<"; Ar = "<<Ar
				<<"; UAr = "<<UAr<<";\n";
			cout << "NwtNew = " << cn->getNwtNew()<<";\n";
			cout << "MiNew  = " << cn->getMiNew() <<";\n";
			cout << "MuNew  = " << cn->getMuNew() <<";\n";
			cout << "RuNew  = " << cn->getRuNew() <<";\n";
			cout << "NfNew  = " << cn->getNfNew() <<";\n";
			cout << "NtNew = " << cn->getNtNew()  <<";\n";
			cout << endl << flush;

			// Re-initialize such element
			dM = MuOld - MiOld;
			MiNew = Mdelt;
			MuNew = MiNew + dM;
			cn->setMuNew(MuNew);
			cn->setMiNew(MiNew);
			cn->setMuOld(MuNew);
			cn->setMiOld(MiNew);
		}
	}
	}

/*************************************************************************
**
**  tHydroModel::SetupNodeUSZ(tCNode *cn)
**
**  Sets up node 'cn' at the beginning of USZ simulation
**
*************************************************************************/
void tHydroModel::SetupNodeUSZ(tCNode *cn)
{
	ID = cn->getID();

	soilPtr->setSoilPtr( cn->getSoilID() ); //set current s.type

	// Get soil hydraulic properties
    // Giuseppe 2016 - Begin changes to allow reading soil properties from grids
    //    Ksat    = soilPtr->getSoilProp(1);  // Surface hydraulic conductivity
    //    Ths     = soilPtr->getSoilProp(2);  // Saturation moisture content
    //    Thr     = soilPtr->getSoilProp(3);  // Residual moisture content
    //    PoreInd = soilPtr->getSoilProp(4);  // Pore-size distribution index
    //    Psib    = soilPtr->getSoilProp(5);  // Air entry bubbling pressure
    //    F       = soilPtr->getSoilProp(6);  // Decay parameter in the exp
    //    Ar      = soilPtr->getSoilProp(7);  // Anisotropy ratio (saturated)
    //    UAr     = soilPtr->getSoilProp(8);  // Anisotropy ratio (unsaturated)
    //    porosity = soilPtr->getSoilProp(9); // Porosity
    Ksat = cn->getKs();  // Surface hydraulic conductivity
    Ths = cn->getThetaS(); // Saturation moisture content
    Thr = cn->getThetaR(); // Residual moisture content
    PoreInd = cn->getPoreSize(); // Pore-size distribution index
    Psib = cn->getAirEBubPres(); // Air entry bubbling pressure
    F = cn->getDecayF(); // Decay parameter in the exp
    Ar = cn->getSatAnRatio(); // Anisotropy ratio (saturated)
    UAr = cn->getUnsatAnRatio(); // Anisotropy ratio (unsaturated)
    porosity = cn->getPorosity(); // Porosity
    // Giuseppe 2016 - End changes to allow reading soil properties from grids
	Eps = 3 + 2/PoreInd;
	
	// Get dynamic variables from prior time step
	NwtOld = NwtNew= cn->getNwtOld();  // nwt
	MuOld  = MuNew = cn->getMuOld();   // mu
	MiOld  = MiNew = cn->getMiOld();   // mti
	NtOld  = NtNew = cn->getNtOld();   // nt
	NfOld  = NfNew = cn->getNfOld();   // nf
	RuOld  = RuNew = cn->getRuOld();   // ru
	RiOld  = RiNew = cn->getRiOld();   // ri

	QpOut          = cn->getQpout();   // QpOut mm^3/hour
	QpIn           = cn->getQpin();    // QpIn  mm^3/hour

	// Calculated fluxes as depth
	QpOut = QpOut*1.0E-6/(cn->getVArea());  // mm/hour = on a horizontal plane
	QpIn  = QpIn *1.0E-6/(cn->getVArea());  // mm/hour = on a horizontal plane

	IntStormVar    = cn->getIntStormVar(); // interstorm variable

	// Re-set runoff which is accumulated over 'EtIStep' time interval
	if ( !(fmod((timer->getCurrentTime() - timer->getTimeStep()),
				timer->getEtIStep())) )
		cn->setSrf_Hr(0.0);
}

/*************************************************************************
**
**  tHydroModel::MaxSoilETRate(tCNode *, double)
**
**  Maximum soil-sourced ET rate [mm hr^-1, horizontal plane] the column
**  can deliver this step. tEvapoTrans calls this before committing ET so
**  that soil evaporation and transpiration are supply-limited at
**  formation: the column cannot dry below the hydrostatic profile with
**  the water table pinned at bedrock (Newton() clamps Nwt at DtoBedrock),
**  so any demand beyond this rate would be silently refused by the
**  unsaturated update while remaining reported as ET, creating phantom
**  water.
**
*************************************************************************/
double tHydroModel::MaxSoilETRate(tCNode *cn, double dt)
{
	double DtoB = cn->getBedrockDepth();

	// Brooks-Corey retention parameters (as SetupNodeUSZ sets them)
	Ths     = cn->getThetaS();
	Thr     = cn->getThetaR();
	PoreInd = cn->getPoreSize();
	Psib    = cn->getAirEBubPres();

	double NwtOldL = cn->getNwtOld();
	double MuOldL  = cn->getMuOld();

	// Extractable water above the water-table-at-bedrock floor [mm],
	// slope-normal
	double extractable = MuOldL + Ths*(DtoB - NwtOldL) - get_Total_Moist(DtoB);
	if (extractable < 0.0)
		extractable = 0.0;

	// Slope factor, matching UnSaturatedZone()'s Cos. R1 carries ET as
	// Ractual*Cos, so the horizontal ET rate the column can supply is
	// extractable / (dt*Cos).
	double edgeSlope = cn->getFlowEdg()->getSlope();
	double Cos = (edgeSlope > 0.0) ? 1.0/sqrt(1.0 + edgeSlope*edgeSlope) : 1.0;

	return extractable/(dt*Cos);
}

//=========================================================================
//
//
//                  Section 3: tHydroModel Unsaturated Code
//
//
//=========================================================================

/*************************************************************************
**
**  tHydroModel::UnsaturatedZone( double )
**
**  The unsaturated zone dynamics algorithm: explicit scheme
**
*************************************************************************/
void tHydroModel::UnSaturatedZone(double dt)
{

	if (simCtrl->Verbose_label == 'Y')
		cout<<"->Unsaturated zone simulation... "<<endl<<flush;

	tCNode * cn;
	tCNode * dnode;
	tCNode * cdest;
	tEdge  * ce;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );
	tMeshListIter<tEdge>  edgIter( gridPtr->getEdgeList() );

	double Kunsat;
	double Ractual;
	double xxsrf;
	double Mperch, Mdelt, Mdva, AA, BB;
	double qn;
	double Nstar, NwtNext, NfNext;
	double ThSurf;
	double EvapSoi, EvapVeg;
	double airTemp, Ta_hi, Ta_lo, alphat;
	double mbR1used, totSurf;

   Ractual = xxsrf = 0.0;        //SMM - added 08132008
   Mperch = Mdelt = Mdva = AA = BB = 0.0; //SMM - added 08132008
   qn = Nstar = NwtNext = NfNext = 0.0;   //SMM - added 08132008
   ThSurf = 0;

	// SKY2008Snow from AJR2007
	double routeWE, snWE; //RINEHART 2007 @ NMT

	double RADAR = 0.1;   //Threshold for RADAR data, mm/hour
	int Pixel_State;

	enum {Storm_Evol, WTStaysAtSurf, WTDropsFromSurf, WTGetsToSurf,
        Storm_Unsat_Evol, Perched_Evol, Perched_SurfSat,
        StormToInterTransition, ExactInitial, IntStormBelow};

	// Start of Simulation Loop for Unsaturated Zone
	//---------------------------------------------

	int id = 0;
	cn = nodIter.FirstP();
	while ( nodIter.IsActive() ) {

#ifdef PARALLEL_TRIBS
  // Receive Qpin from incoming outlet node(s) on another processor
  // if this is a stream head node
  if (tGraph::isUpstreamNode(cn)) {
    id = cn->getReach();
    tGraph::receiveQpin(id, cn);

    // If runon option set, receive flux node values
    if (RunOnoption)
        tGraph::receiveRunFlux(cn);
  }
#endif

		// Setup basic variables
		SetupNodeUSZ( cn );

		// Initialize Runoff
		srf=hsrf=esrf=psrf=satsrf=sbsrf=0.0;

		// Get geometry properties
		ce = cn->getFlowEdg();

		// Slope factors. Originally: alpha = atan(slope); Cos = cos(alpha);
		// Sin = sin(alpha). Rewritten with the identities
		// cos(atan(s)) = 1/sqrt(1+s^2) and sin(atan(s)) = s/sqrt(1+s^2)
		// to avoid three expensive trig calls per node
		double edgeSlope = ce->getSlope();
		if (edgeSlope > 0.0) {
			Cos = 1.0/sqrt(1.0 + edgeSlope*edgeSlope);
			Sin = edgeSlope*Cos;
		}
		else {
			Cos = 1.0;
			Sin = 0.0;
		}
		
		// Get Bedrock depth for computing Nwt
		DtoBedrock = cn->getBedrockDepth(); // Added by CJC2020
		
		// Get Actual Rainfall after ET and I
		EvapSoi = cn->getEvapSoil();
		EvapVeg = cn->getEvapDryCanopy();

		// No DEW is allowed in soil
		if (EvapSoi < 0.0)
			EvapSoi = 0.0;
		if (EvapVeg < 0.0)
			EvapVeg = 0.0;

		// Get Actual Rainfall after ET and I
		if (Ioption == 0 && EToption == 0)
			Ractual = cn->getRain();

		else if (Ioption == 0 && EToption != 0)
			Ractual = cn->getRain() - EvapSoi - EvapVeg;

		else if (Ioption != 0 && EToption == 0) {
			cout<<"\nError: Cannot Have Rutter Model (on) and ET (off)"<<endl;
			Ractual = cn->getRain();
		}
		else if (Ioption != 0 && EToption != 0) {
            Ractual = cn->getNetPrecipitation() - EvapSoi - EvapVeg;
        }

		// Runon calculations (if on), runon value [mm hr^-1]
		qrunon = 0.0;
		if (RunOnoption)
			qrunon = GetCellRunon( cn, dt );

		// Calculate accumulated rain [m^3]
		TotRain += (Ractual*dt*cn->getVArea()/1000.);

		// SKY2008Snow from AJR2007
		if (SnOpt) {
			snWE = cn->getLiqWE() + cn->getIceWE();
			routeWE = cn->getLiqRouted();
			if ((snWE > 1e-3) || (routeWE > 0.)) {
                //WR 12182023 removed EvapVeg from route water following above approach and because transpiration can still occur w/snow
				Ractual = 10.0*routeWE-EvapVeg; //have to convert to mm // Changed from R to Ractual CJC2020
			}
		}

        // Adjust saturated hydraulic conductivity based on air
//        // temperature to account for frozen soil effects CJC2020
//        Ta_hi = 8;
//        Ta_lo = 4;
//        alphat = 0.004;
//        airTemp = cn->getAirTemp();
//        if (airTemp <= Ta_lo) {
//            Ksat = Ksat*alphat;
//        }
//        else if ((airTemp < Ta_hi) && (airTemp > Ta_lo)) {
//            Ksat = Ksat*alphat + (airTemp - Ta_lo)*((Ksat - Ksat*alphat)/(Ta_hi - Ta_lo));
//        }
		
		// Total rate of in/outfluxes [mm hr^-1]
		R1 = (Ractual + (QpIn-QpOut))*Cos + qrunon;  // Moved this line to be below routeWE  calculation CJC2020
		R = Ractual; // Moved this line to be below routeWE  calculation CJC2020


		// Step1: Compute K_unsaturated and rainfall
		//---------------------------------------------
		if (Ksat != 0.0) {
            ThSurf = get_InitMoist_depthZ(0.);
        }

		if (RuOld == 0.0) {
            Kunsat = Ksat * pow(((ThSurf - Thr) / (Ths - Thr)), Eps); //Initialized profile
        }
		else {
            Kunsat = RuOld;  // There is wetted wedge
        }


		// Step2: Determine the node state
		//---------------------------------------------

		// Splitting into different situations
		Pixel_State = -1000;
		if ( (NwtOld==0.) && (R1>=0.) ) { // WT initially at surface & stays there
            Pixel_State = WTStaysAtSurf;
        }

		if ( (NwtOld==0.) && (R1<0.) ) {   // WT initially at surface & drops
            Pixel_State = WTDropsFromSurf;
        }

		if ((NwtOld>0.)&&((R1*dt)>=(NwtOld*Ths-MuOld))&&(NfOld==0.||NfOld==NwtOld)) {
            Pixel_State = WTGetsToSurf;       // WT initially at some depth, reaches surface
        }

		if (Pixel_State==-1000) {        // None of the above
			MuNew = MuOld + dt*R1;  	   // Calculate MuNew from forcing

			if (MuNew < MiOld) {           // Apply the pertinent interstorm equation
                Pixel_State = IntStormBelow;
            }
			if (fabs(MuNew - MiOld) < 1.0E-6) {
                Pixel_State = ExactInitial; // Exactly Initialized State
            }
			if (MuNew > MiOld) {
				if (R1 > 0.0) {

					// Check if moisture redistribution is needed
					if (R > RADAR ) {

						// Must destroy the edge and redistribute moisture
						if (IntStormVar >= IntStormMAX && NfOld > 0.0) {
							if (MuNew >= NwtOld*Ths) {
								NwtNew = 0.0;
								MiNew = 0.0;
								sbsrf = (MuNew - NwtOld*Ths)/dt;
							}
							else {
								// There is enough space above Nt
								Mdelt = Ths*NwtOld - MuOld;
								NwtNew = Newton(Mdelt, NwtOld);
								MiNew = get_Total_Moist(NwtNew);
							}
							NwtOld = NwtNew;
							MiOld = MuOld = MuNew = MiNew;
							NtOld = NtNew = 0.0;
							NfOld = NfNew = 0.0;
							RiOld = RiNew = 0.0;
							RuOld = RuNew = 0.0;
							IntStormVar   = 0.0;
						}
					}

					// Proceed with old/new values
					if (NtOld == NfOld) {
                        Pixel_State = Storm_Unsat_Evol;
                    }

					if ((NtOld < NfOld) && (NtOld != 0.0)) {
                        Pixel_State = Perched_Evol;
                    }
					if ((NtOld == 0.0) && (NfOld > 0.0)) {
						ThRiNf = get_InitMoist_depthZ(NfOld);
						ThReNf = Ths;
						set_Suction_Term(NfOld);
						qn = Ksat*F*NfOld/expm1(F*NfOld)*(Cos + G/NfOld); //WR debug 01112024: Converted all instances of exp(x)-1 to expm1(x) to address loss in precision

						xxsrf = qn - RiOld*Cos; // Net infiltration rate
						if (R1 >= xxsrf) {
                            Pixel_State = Perched_SurfSat;
                        }
						else {
                            Pixel_State = StormToInterTransition;
                        }
					}

					// --------------------------------------------------
					// The following handles two aspects:
					//       a) Rainfall is continuing and Nf has reached
					//          water table -> redistribute moisture
					//       b) Interstorm  subsurface lateral fluxes
					//          should not cause the formation of new Nf
					//       -> NtOld = NwtOld, moisture is redistributed
					if ((NfOld == NwtOld)||(R<=RADAR && IntStormVar>0 && NfOld==0.0))
						Pixel_State = Storm_Evol;
      }
				if (R1 <= 0.0)
					Pixel_State = StormToInterTransition;
    }
			if (Ksat == 0.0)
				Pixel_State = Perched_SurfSat;
  }

		// Print old values of state variables
		PrintOldVars( cn, ce, Ractual, Pixel_State );

		// Step 3: Go to the cases
		//---------------------------------------------

		switch (Pixel_State) {


			//----------------
			case WTStaysAtSurf:

				// Saturated column: all net influx leaves as runoff, so
				// psrf + sbsrf must equal R1 (>= 0 by the state entry
				// condition). R carries ET with its sign, a negative R
				// (ET on a laterally-supplied saturated node) is debited
				// against the inflow, not dropped; the pair is clamped
				// so both components stay non-negative. CJC2026
				if (QpIn > QpOut)
					psrf = (QpIn-QpOut)*Cos;
				else
					psrf = 0.0;

				sbsrf = R1 - psrf;
				if (sbsrf < 0.0) {
					psrf += sbsrf;
					sbsrf = 0.0;
				}

				MiNew=0.0;
				MuNew=0.0;
				RiNew=0.0;
				RuNew=0.0;
				NfNew=0.0;
				NtNew=0.0;
				NwtNew=0.0;
				QpOut=0.0;

				// If ExactInitial during interstorm period
				// use 0.1 mm/hour is a threshold value for RADAR

                if (R <= RADAR) {
					IntStormVar += dt;
				}
                else {
                    if (IntStormVar > 0.0)
                        IntStormVar = 0.0;
                }
					break;


				//----------------
			case WTGetsToSurf:

				// Partitioning into sbsrf & psrf
				// Column fills to the surface: runoff is the net influx
				// beyond the remaining deficit, psrf + sbsrf =
				// R1 - (NwtOld*Ths - MuOld)/dt (>= 0 by the state entry
				// condition). Lateral inflow fills the deficit first and
				// its excess exits as psrf; the rain-side R carries ET
				// with its sign so a negative R is debited against the
				// lateral excess, not dropped; the pair is clamped so
				// both components stay non-negative. CJC2026
				totSurf = R1 + MuOld/dt - NwtOld*Ths/dt;
				if (QpIn > QpOut)
					psrf = (QpIn - QpOut)*Cos + MuOld/dt - NwtOld*Ths/dt;
				else
						psrf = 0.0;
				if (psrf < 0.0)
					psrf = 0.0;
				sbsrf = totSurf - psrf;
				if (sbsrf < 0.0) {
					psrf += sbsrf;
					sbsrf = 0.0;
				}

				MuNew=0.0;
				RiNew=0.0;
				RuNew=0.0;
				NfNew=0.0;
				NtNew=0.0;
				NwtNew=0.0;
				MiNew=0.0;
				QpOut=0.0;

				if (IntStormVar > 0.0)
					IntStormVar = 0.0;

                break;


				//----------------
			case WTDropsFromSurf:

				R1 = (R + (QpIn - QpOut))*Cos;
				if (R1 < 0.0)
					NwtNew = Newton(fabs(R1*dt), 0.0);

                MiNew = get_Total_Moist(NwtNew); //WR Clang-Tidy: Misleading indentation: statement is indented too deeply. Maybe should be in scope of if statement.

				MuNew = MiNew;
				RiNew = 0.0;
				RuNew = 0.0;
				QpOut=0.0;

				IntStormVar += dt;

                if(IntStormVar >= IntStormMAX) {
					NfNew = 0.0;
					NtNew = 0.0;
				}else {
                    NfNew=NwtNew;
                    NtNew=NwtNew;
                }
                break;


				//----------------
			case ExactInitial:

				// Moisture conditions correspond to initialized state
				MiNew=MiOld;
				MuNew=MiNew;
				NwtNew=NwtOld;
				RiNew=0.0;
				RuNew=0.0;
				QpOut=0.0;

				if (R <= 0.0)
					IntStormVar += dt;

                if (NfOld == NwtOld) {
                    if (IntStormVar >= IntStormMAX) {
							NfNew = 0.0;
							NtNew = 0.0;
						}
                    else {
							NfNew=NwtNew;
							NtNew=NwtNew;
						}
					}
                else {
                    NfNew = 0.0;
                    NtNew = 0.0;
                }

                break;


				//----------------
			case IntStormBelow:

				R1 = (R + (QpIn - QpOut))*Cos;


                Mdelt = Ths*NwtOld - (MuOld + R1*dt); // Pixel moisture deficit
				NwtNew = Newton(Mdelt, NwtOld);


				RiNew = 0.0;
				RuNew = 0.0;
				MiNew = get_Total_Moist(NwtNew);
				MuNew = MiNew;

				IntStormVar += dt;  // <- NOVA

				// If preceding state was StormToInterTransition
				// assign fronts to the surface
                if ((NfOld < NwtOld) && (NfOld > 0.0)) {
					NfNew = NtNew = 0.0;
				}
                // If preceding state was {IntStormBelow, ExactInitial, WTDrops}
                else {
                    if (IntStormVar < IntStormMAX) {
                        if (NfOld == 0.0)
                            NfNew = NtNew = 0.0;
                        else
                            NfNew = NtNew = NwtNew;
						}
                    else
                        NfNew = NtNew = 0.0;
					}
					QpOut = 0.0;
				break;


				//----------------
			case StormToInterTransition:

				if (R <= RADAR) {
					IntStormVar += dt;
				}
				else
					IntStormVar = 0.0;

				R1 = (R + (QpIn - QpOut))*Cos;
				MuNew  = MuOld + R1*dt;
				NwtNew = NwtOld;
				MiNew  = MiOld;
				RiNew  = RiOld;

				// The case of deep interstorm period:
				// T >= Tmax (w/o rain on surface)
				if (IntStormVar >= IntStormMAX && RdstrOption) {

					// Destroy edge and redistribute moisture
					if (MuNew >= NwtOld*Ths) {
						NwtNew = 0.0;
						sbsrf = (MuNew - NwtOld*Ths)/dt;
						MiNew = 0.0;
					}
					// There is enough space above Nt
					else  {
						Mdelt = Ths*NwtOld - MuNew;
						NwtNew = Newton(Mdelt, NwtOld);
						MiNew = get_Total_Moist(NwtNew);
					}
					MuNew = MiNew;
					NtNew = 0.0;
					NfNew = 0.0;
					RiNew = 0.0;
					RuNew = 0.0;
				}
                // Rainfall hiatus or interstorm beginning:
                // T < Tmax (w/o rain on surface), wedge is moving down

                else {
                    Mdelt = get_Lower_Moist(NfOld, NwtNew);
                    Mdelt = MuNew - Mdelt;          //Amount of water in the wetted edge
                    AA = Eps/F*(Ths-Thr)*(1.0 - exp(-F*NfOld/Eps)) + Thr*NfOld;

                    // Correct variables in imbalance (artifacts of Perched_Evol case)
                    if (Mdelt > AA && NtOld == NfOld) {
                        Mdelt = AA - 1.0E-4;

                        if (simCtrl->Verbose_label == 'Y') {
                            cout<<"\nWarning: IMBALANCE: Correction by "<<Mdelt - AA<<" mm"<<endl;
                            cout<<"ID = "<<ID<<endl;
							}

						}


						// CASE 1 : UNSATURATED EDGE
						if (Mdelt <= AA) {

							if (R1 == 0.0) {
								RuNew = RuOld;
								Nstar = (log(Ksat/RuNew))/F;
							}
							else {   // [ R1 < 0 ]
								RuNew = get_RechargeRate(Mdelt, NfOld);
								Nstar = (log(Ksat/RuNew))/F;
							}

							// MOVE FRONTS: Wedge gets very close to the MiOld profile (0.05%)
							if ((100.*(MuNew-MiOld)/MiOld) <= 0.05) {
								NfNew = NtNew = 0.0;   //Pixel is "ready" for another rainfall

								if ((MuNew - MiOld) > 1.0E-4) {
									Mdelt = Ths*NwtOld - MuNew;
									NwtNew = Newton(Mdelt, NwtOld);
								}
								MiNew = get_Total_Moist(NwtNew);
								MuNew = MiNew;
								RiNew = RuNew = 0.0;
							}

							// Wetted wedge is still SIGNIFICANT (RuNew>RiNew >> 0.1%)
							else {
								ThRiNf = get_InitMoist_depthZ(NfOld);
								ThReNf = get_EdgeMoist_depthZ(NfOld);

								// Gravity cannot be dominant any longer --
								// The idea is to just redistribute the moisture
								// in the profile by raising the water table
								if (ThReNf <= ThRiNf) {
									Mdelt = Ths*NwtOld - (MuOld + R1*dt);
									NwtNew = Newton(Mdelt, NwtOld);
									NtNew = NfNew = 0.0;
									RuNew = RiNew = 0.0;
									MiNew = get_Total_Moist(NwtNew);
									MuNew = MiNew;
									QpOut = 0.0;
								}
								else {
									set_Suction_Term(NfOld);
									qn = RuNew*Cos + Ksat*exp(-F*NfOld)*G/NfOld;
									NfNew = NfOld + dt*(qn - RiNew*Cos)/(ThReNf-ThRiNf);
									NtNew = NfNew;


									if (NfNew < (NwtNew+Psib)) {
										Mdelt = get_Lower_Moist(NfNew, NwtNew);
										Mdelt = MuNew - Mdelt;
										RuNew = get_RechargeRate(Mdelt, NfNew);
										AA = Eps/F*(Ths-Thr)*(1.0 - exp(-F*NfNew/Eps)) + Thr*NfNew;

										// If the wedge perches, re-define variables
										if (Mdelt > AA) {
											NtNew = (log(Ksat/RuNew))/F;
											RuNew = Ksat*exp(-F*NtNew);
										}
									}
									// Wetting front hits the capillary fringe
									else if (NfNew >= (NwtNew+Psib))  {
										if ((MuNew - MiOld) > 1.0E-4) {
											Mdelt = Ths*NwtOld - MuNew;
											NwtNew = Newton(Mdelt, NwtOld);
										}
										MiNew = get_Total_Moist(NwtNew);
										MuNew = MiNew;
										RiNew = RuNew = 0.0;

										cdest = (tCNode*)ce->getDestinationPtrNC();
										NwtNext = cdest->getNwtOld();
										NfNext  = cdest->getNfOld();

										if ((NfNext == NwtNext) && (IntStormVar < IntStormMAX))
											NfNew = NtNew = NwtNew;
										else
											NfNew = NtNew = 0.0;
										RuNew = 0.0;
									}
								}
							}
						}

						// CASE 2 :  PERCHED & SURFACE SATURATED WEDG

						// a.) Let's move the wetting front
						else {
							ThRiNf = get_InitMoist_depthZ(NfOld);
							ThReNf = Ths;
							set_Suction_Term(NfOld);
							qn = Ksat*F*(NfOld-NtOld)/(exp(F*NfOld)-exp(F*NtOld));
							qn *= (Cos + G/NfOld);

							NfNew = NfOld + dt*(qn-RiNew*Cos)/(Ths-ThRiNf);

							// Check possible situations with Nf & Nwt
							if (fabs(NfNew - (NwtNew+Psib)) <= 1.0E-3) {
								NwtNew = NtOld-Psib;
								if (NwtNew > DtoBedrock)
									NwtNew = DtoBedrock;
								cdest = (tCNode*)ce->getDestinationPtrNC();
								NwtNext = cdest->getNwtOld();
								NfNext  = cdest->getNfOld();

								if ((NfNext == NwtNext) && (IntStormVar < IntStormMAX))
									NfNew = NwtNew;
								else
									NfNew = 0.0;
							}
							else if (NfNew > (NwtNew+Psib)) {
								NwtNew = NtOld-Psib;
								if (NwtNew > DtoBedrock)
									NwtNew = DtoBedrock;
								cdest = (tCNode*)ce->getDestinationPtrNC();
								NwtNext = cdest->getNwtOld();
								NfNext  = cdest->getNfOld();

								if ((NfNext == NwtNext) && (IntStormVar < IntStormMAX))
									NfNew = NwtNew;
								else
									NfNew = 0;
							}

							// b.) Let's move the top front

							// The wetting front has reached the water table
							if ((NfNew == 0.0) || (NfNew == NwtNew)) {
								if (MuNew >= NwtOld*Ths) {
									NwtNew = 0.0;
									NtNew = NfNew = 0.0;
									sbsrf = (MuNew - NwtOld*Ths)/dt;
									MuNew = MiNew = 0.0;
									RuNew = RiNew = 0.0;
								}
								else  {    //Enough space above Nt
									Mdelt = Ths*NwtOld - MuNew;
									NwtNew = Newton(Mdelt, NwtOld);
									RiNew = RuNew = 0.0;
									MiNew = get_Total_Moist(NwtNew);
									MuNew = MiNew;
									if (NfNew == NwtNew)
										NfNew = NtNew = NwtNew;
									else
										NfNew = NtNew = 0.0;
								}
							}

							//  The wetting front is still above the water table
							else if ((NfNew > 0.0) && (NfNew != NwtNew)) {
								Mdelt = Eps/F*(Ths-Thr)*(1.0 - exp(-F*NtOld/Eps)) + Thr*NtOld;
								Mdva  = Eps/F*(Ths-Thr)*(1.0 - exp(-F*NfNew/Eps)) + Thr*NfNew;

								// Max of water that can be lost without
								// transition to unsaturated state
								AA = (Mdelt+(NfNew-NtOld)*Ths - Mdva);

								// Next state is unsaturated
								if (R1 < (qn - AA/dt)) {
									Mperch = get_Lower_Moist(NfNew, NwtNew);
									Mperch = MuNew - Mperch; // water left in the wetted wedge
									RuNew = get_RechargeRate(Mperch, NfNew);
									NtNew = NfNew;

									if (Mdelt >= Mdva)
										MuNew -= Mdelt - Mdva + 1.0E-4;
								}
								// Next state is either surfuce or perched saturated
								else {
									BB = -dt*(qn - R1)/(Ths-Thr) - Eps/F*exp(-F*NtOld/Eps);
									AA = -exp(F*(BB-NtOld)/Eps);
									// AA below the Lambert W branch point (-1/e) has no
									// real W0 value; clamp to the branch-point limit
									// W0(-1/e) = -1 instead of driving the Halley
									// iteration into a non-convergent region.
									// Original: NtNew = NtOld + (Eps/F*LambertW(AA) - BB);
									if (AA < (-1/exp(1.0))) {
										NtNew = NtOld + (Eps/F*(-1.0) - BB);
										static int warnLambert1 = 0;
										if (++warnLambert1 <= 10)
											cout<<"\nWarning: Value for LAMBERT f-n is too small\n\n";
										else if (warnLambert1 == 11)
											cout<<"\nWarning: Value for LAMBERT f-n is too small (further warnings suppressed)\n\n";
									}
									else
										NtNew = NtOld + (Eps/F*LambertW(AA) - BB);
									RuNew = Ksat*exp(-F*NtNew);
								}
							}
	}

						QpOut = 0.0;
						if ((NtNew == NfNew) && (NfNew != NwtNew) && (NfNew > 10.))
							QpOut = get_UnSat_LateralFlow(NfNew, RuNew, ce->getVEdgLen()*1000.);

						if ((NfNew-NtNew) > 1.0)
							QpOut = get_Sat_LateralFlow(NtNew, NfNew, RuNew, ce->getVEdgLen()*1000.);

						QpOut *= Sin;

      } // Matches *else* - the case of T < Tmax
					break;


				//---------------
			case Storm_Evol:

				R1 = (R + (QpIn-QpOut))*Cos;
				Mdelt = Ths*NwtOld - (MuOld + R1*dt);
				NwtNew = Newton(Mdelt, NwtOld);
				RiNew = RuNew = 0.0;
				MiNew = get_Total_Moist(NwtNew);
				MuNew = MiNew;
				QpOut = 0.0;

				if (R <= RADAR) {
					IntStormVar += dt;

					if (NfOld == NwtOld) {
						if (IntStormVar >= IntStormMAX)
							NfNew = NtNew = 0.0;
						else
							NfNew=NtNew=NwtNew;
					}
					else
						NfNew = NtNew = 0.0;
				}

                else {
                    if (IntStormVar > 0.0)
                        IntStormVar = 0.0;
                    NfNew=NtNew=NwtNew;
                }
					break;


				//--------------------
			case Storm_Unsat_Evol:

				if (R <= RADAR)
					IntStormVar += dt;
				else
					IntStormVar = 0.0;

				R1 = (R + (QpIn-QpOut))*Cos;

				NwtNew= NwtOld;
				MuNew = MuOld + R1*dt;
				MiNew = MiOld;
				RiNew = RiOld = 0.0;

				// CASE 1 : NO FRONTS YET / START FROM INITIALIZED STATE

				if (NfOld == 0.0 && NtOld == 0.0) {

					ThRiNf = get_InitMoist_depthZ(0.0);
					// Only gravitational flow
					if (R < Kunsat && R1/Cos < Kunsat) {
						ThReNf = pow((R1/Ksat),(1/Eps))*(Ths-Thr)+Thr;
						qn = R1;
						NfNew = dt*qn/(ThReNf-ThRiNf);
					}
					else {
						if (R1 < Ksat)  {
							ThReNf = Thr + (Ths-Thr)*pow((R1/Ksat),(1/Eps));
							set_Suction_Term(0.0);
							if (fabs(ThReNf-ThRiNf) < 1.0E-4)
								qn = R1;
							else
								qn = Ksat*exp(-F*R1/(ThReNf-ThRiNf))*G*(ThReNf-ThRiNf)/R1 + R1;
							NfNew = dt*qn/(Ths-ThRiNf);
						}
						else  {   // [ R1 >= Ksat ]
							ThReNf = Ths;
							set_Suction_Term(0.0);
							qn = Ksat*G*(ThReNf-ThRiNf)/Ksat + R1;
							NfNew  = dt*qn/(ThReNf-ThRiNf);
						}
					}

					// Artificial forcing to a larger value
					if (NfNew > 0.0 && NfNew < 32000.) {
						if (NfNew < 1.0)
							NfNew = 1.0;
						NtNew  = NfNew;
					}
					else
						NtNew = NfNew = NwtOld + 100.;

					// Calculation of Requivalent:
					if (NfNew < NwtOld) {
						Mdelt = get_Lower_Moist(NfNew, NwtNew);
						Mdelt = MuNew - Mdelt;
						AA = Eps/F*(Ths-Thr)*(1.0 - exp(-F*NfNew/Eps)) + Thr*NfNew;

						if (Mdelt >= AA) {
							// CASE 1
							if (fabs(Mdelt - AA) <= 1.0E-6) {
								NtNew = NfNew;
								NfNew += 1.0E-5;
								RuNew = Ksat*exp(-F*NtNew);
							}
							// CASE 2
							else {
								if (Mdelt >= NfNew*Ths) {
									ThRiNf = get_InitMoist_depthZ(NfNew);
									NfNew += (Mdelt - NfNew*Ths)/(Ths - ThRiNf);
									if (NfNew >= (NwtNew + Psib)) {
										NfNew = NtNew = NwtNew = 0.0;
										MuNew = MiNew = 0.0;
										RuNew = RiNew = 0.0;
									}
									else {
										NtNew = 0.0;
										RuNew = Ksat*F*NfNew/expm1(F*NfNew);
										AA = get_Lower_Moist(NfNew, NwtNew);
										BB = NfNew*Ths;
										MuNew = AA + BB;
									}
								}
								// Top front is NOT at the surface
								else if (Mdelt < NfNew*Ths) {
									RuNew = get_RechargeRate(Mdelt, NfNew);
									NtNew = (log(Ksat/RuNew))/F;
									NfNew += 1.0E-5;
								}
							}
						}
						// CASE 2
						else
							RuNew = get_RechargeRate(Mdelt, NfNew);


						if ((RiNew-RuNew) < 1.0E-2 && (RiNew - RuNew) > 0.0)
							RuNew = RiNew;
						else if ((RiNew-RuNew) > 1.0E-2) {
							static int warnRuRi = 0;
							if (++warnRuRi <= 10)
								cout<<"\nWarning: UNSAT: RuNew < RiNew: id = "<<cn->getID()<<"\n";
							else if (warnRuRi == 11)
								cout<<"\nWarning: UNSAT: RuNew < RiNew (further warnings suppressed)\n";
						}

						if (RuNew > Ksat) {
							static int warnRuKsat = 0;
							if (++warnRuKsat <= 10)
								cout<<"\nWarning: UNSAT: RuNew > Ksat: id = "
								<<cn->getID()<<"\n\tRuNew = "<<RuNew<<"\tKsat = "<<Ksat<<"\n";
							else if (warnRuKsat == 11)
								cout<<"\nWarning: UNSAT: RuNew > Ksat (further warnings suppressed)\n";
						}
                    }
				}

                // CASE 2 : FRONTS BELOW THE SURFACE
                //          A) Situation of surface perching

                else if (NfOld > 0.0 && NtOld > 0.0) {
                    Mdelt = get_Lower_Moist(NfOld, NwtNew);
                    Mdelt = MuNew - Mdelt;

                    // The following case can be caused by too high dt
                    if (Mdelt >= NfOld*Ths) {
                        ThRiNf = get_InitMoist_depthZ(NfOld);
                        ThReNf = Ths;
                        set_Suction_Term(NfOld);
                        qn = Ksat*F*NfOld/expm1(F*NfOld);

                        // Use SurfSatModel functions
                        if ((qn*(Cos+G/NfOld)-RiOld*Cos) <= R1) {
                            qn *= (Cos + G/NfOld);
                            if (R1 >= (qn - RiOld*Cos)) {
                                hsrf += (R1 - (qn - RiOld*Cos));
                                R1 -= hsrf;
                            }
                            else {
                                static int warnStateDef1 = 0;
                                if (++warnStateDef1 <= 10)
                                    cout<<"\nWarning: Wrong state def.: R < Keqviv-RiOld*Cos: id = "
                                    <<cn->getID()<<"\n\n";
                                else if (warnStateDef1 == 11)
                                    cout<<"\nWarning: Wrong state def.: R < Keqviv-RiOld*Cos (further warnings suppressed)\n\n";
                            }

                            NtNew = 0.0;
                            MuNew = MuOld + R1*dt;
                            NfNew = NfOld + dt*(qn-RiNew*Cos)/(Ths-ThRiNf);
                            RuNew = Ksat*F*NfNew/expm1(F*NfNew);

                            // Check to see whether WT has to be modified
                            if ((fabs(NfNew-(NwtNew+Psib))<=1.0E-3) || (NfNew>(NwtNew+Psib))) {
                                NwtNew = MiNew = 0.0;
                                RuNew = RiNew = 0.0;
                                Mperch = MuNew - NwtOld*Ths;
                                if (Mperch > 0.0)
                                    sbsrf += Mperch/dt;
                                MuNew = 0.0;
                                NtNew = NfNew = 0.0;
                            }
                        }
                        // Keep the unsaturated edge
                        else {
                            ThReNf = get_EdgeMoist_depthZ(NfOld);
                            set_Suction_Term(NfOld);
                            qn = RuOld*Cos + Ksat*exp(-F*NfOld)*G/NfOld;
                            NfNew = NfOld + dt*(qn-RiNew*Cos)/(Ths-ThRiNf);
                            NtNew = NfNew;
                            if (NfNew < (NwtNew+Psib)) {
                                Mdelt = get_Lower_Moist(NfNew, NwtNew);
                                Mdelt = MuNew - Mdelt;
                                RuNew = get_RechargeRate(Mdelt, NfNew);
                            }
                            else {
                                if (MuNew > Ths*NwtNew) {
                                    sbsrf += (MuNew  - Ths*NwtNew);
                                    NfNew = NtNew = NwtNew = 0.0;
                                    MuNew = MiNew = 0.0;
                                    RuNew = RiNew = 0.0;
                                }
                            }
                        }
                    }
                    // CASE 2 : FRONTS BELOW THE SURFACE
                    //          B) unsaturated wedge evolution
                    else {
                        RuNew = get_RechargeRate(Mdelt, NfOld);
                        Nstar = (log(Ksat/RuNew))/F;

                        // The following situation occurs when a new slug of water
                        // added to storage above NfOld exceeds the equivalent
                        // moisture amount corresponding to NfOld depth but still
                        // less than the maximum moisture capacity: NfOld*Ths
                        // * The method used here to redefine Ntop and Ru is _approximate_
                        // * More accurate way would be using LambertW function

                        if (Nstar <= NfOld) {
                            NtNew = Nstar;
                            NfNew = NfOld+1.0E-5;
                        }
                        // We also need to move the wetting front deeper but
                        // it is inconvenient because this requires re-writing
                        // Perched_Sat case. With another code implementation,
                        // we would call "Perched_Evol" function
                        else {  // NfOld < Nstar - edge is still to be saturated
                            ThRiNf = get_InitMoist_depthZ(NfOld);
                            ThReNf = get_EdgeMoist_depthZ(NfOld);
                            if (ThReNf <= ThRiNf) { // not too much water has been added

                                // Water table is too close to the surface
                                // The idea is to just redistribute the moisture
                                // in the profile raising the local water table
                                Mdelt = Ths*NwtOld - (MuOld + R1*dt);
                                NwtNew = Newton(Mdelt, NwtOld);
                                RuNew = RiNew = 0.0;
                                MiNew = get_Total_Moist(NwtNew);
                                MuNew = MiNew;
                                QpOut = 0.0;

                                if (NwtOld < 2.0*fabs(Psib)) {
                                    NtNew = NfNew = NwtNew;
                                }
                                else {
                                    NtNew = NfNew = 0.0;
									}
                            }

                            else {
                                set_Suction_Term(NfOld);
                                qn = RuNew*Cos + Ksat*exp(-F*NfOld)*G/NfOld;
                                NfNew  = NfOld + dt*(qn - RiNew*Cos)/(ThReNf-ThRiNf);
                                NtNew  = NfNew;
                                if (NfNew < (NwtNew+Psib)) {
                                    Mdelt = get_Lower_Moist(NfNew, NwtNew);
                                    Mdelt = MuNew - Mdelt;
                                    RuNew = get_RechargeRate(Mdelt, NfNew);
                                }
                            }
                        } // matches NfOld < Nstar
                    }
                }
					// Redistribute moisture along the profile:
					// get to initialized state, raise water table
					if ((NfNew > (NwtNew+Psib)) && (NfNew != NwtNew)) {

						Mdelt = Ths*NwtNew - MuNew;
						NwtNew = Newton(Mdelt, NwtOld);
						RiNew = RuNew = 0.0;
						MuNew = MiNew = get_Total_Moist(NwtNew);
						NfNew = NtNew = NwtNew;
					}

					// Get the subsurface flux downstream
					QpOut = 0.0;
				if ((NtNew == NfNew) && (NfNew != NwtNew) && (NfNew > 10.0))
					QpOut = get_UnSat_LateralFlow(NfNew, RuNew, ce->getVEdgLen()*1000.0);

                if ((NfNew-NtNew) > 1.0)
                    QpOut = get_Sat_LateralFlow(NtNew, NfNew, RuNew, ce->getVEdgLen()*1000.0);

                QpOut *= Sin;
				break;


				//---------------
			case Perched_Evol:

				if (R <= RADAR)
					IntStormVar += dt;
				else
					IntStormVar = 0.0;

				R1 = (R + (QpIn-QpOut))*Cos;

				NwtNew = NwtOld;
				MuNew = MuOld + R1*dt;
				MiNew = MiOld;
				RiNew = RiOld;

				// The problem here is that the new coming moisture is not
				// taken into account up to the point where we define NtNew
				// We define the flux qn based on the OLD position of Nf and Nt
				// Given a sufficiently small time step this should not be an issue

				ThRiNf = get_InitMoist_depthZ(NfOld);
				ThReNf = Ths;
				set_Suction_Term(NfOld);

				qn = Ksat*F*(NfOld-NtOld)/(exp(F*NfOld)-exp(F*NtOld));
				qn *= (Cos + G/NfOld);

				// Wetting front evolution
				NfNew = NfOld + dt*(qn-RiNew*Cos)/(Ths-ThRiNf);

				// Check for possible situations with Nf and Nwt
				if (fabs(NfNew - (NwtNew+Psib)) <= 1.0E-3) {
					NwtNew = NtOld-Psib;
					if (NwtNew > DtoBedrock)
						NwtNew = DtoBedrock;
					NfNew = NwtNew;
				}
                else if (NfNew > (NwtNew+Psib)) {
                    NwtNew = NtOld-Psib;
                    if (NwtNew > DtoBedrock)
                        NwtNew = DtoBedrock;
                    NfNew = NwtNew;
                    R1 += qn - ((NwtOld+Psib-NfOld)*(Ths-ThRiNf)/dt + RiNew*Cos);
                }

                // Wetting front has reached the water table
                if (NfNew == NwtNew) {

                    if (MuNew >= NwtOld*Ths) {
                        NwtNew = NtNew = NfNew = 0.0;
                        sbsrf = (MuNew - NwtOld*Ths)/dt;
                        MuNew = MiNew = 0.0;
                        RuNew = RiNew = 0.0;
                    }
                    else  {
                        Mdelt = Ths*NwtOld - MuNew;
                        NwtNew = Newton(Mdelt, NwtOld);
                        RiNew = RuNew = 0.0;
                        MiNew = get_Total_Moist(NwtNew);
                        MuNew = MiNew;
                        NtNew = NfNew = NwtNew;
                    }
                }
                // If NfNew is still above water table

                else if (NfNew != NwtNew) {
                    Mdelt = Eps/F*(Ths-Thr)*(1.0 - exp(-F*NtOld/Eps)) + Thr*NtOld;
                    Mdva  = Eps/F*(Ths-Thr)*(1.0 - exp(-F*NfNew/Eps)) + Thr*NfNew;
                    // Max amount of water that can be lost without
                    // transforming the state to unsaturated state
                    AA = (Mdelt+(NfNew-NtOld)*Ths - Mdva);

                    // Next state is unsaturated
                    if ((R1+RiNew*Cos) < (qn - AA/dt)) {
                        Mperch = get_Lower_Moist(NfNew, NwtNew);
                        Mperch = MuNew - Mperch;
                        RuNew = get_RechargeRate(Mperch, NfNew);
                        NtNew = NfNew;
                    }
                    //  Next state is either Surf_Sat or Perch_Sat
                    else {
                        xxsrf = (R1 - qn)*dt + Mdelt;
                        // The influx is high enough to fill space above Ntop
                        if (xxsrf >= NtOld*Ths) {
                            NtNew = 0.0;
                            hsrf += (xxsrf - NtOld*Ths)/dt;
                            MuNew -= hsrf*dt;
                            RuNew = Ksat*F*NfNew/expm1(F*NfNew);
                        }
                        else {
                            BB = -dt*(qn - R1)/(Ths-Thr) - Eps/F*exp(-F*NtOld/Eps);
                            AA = -exp(F*(BB-NtOld)/Eps);
                            // AA below the Lambert W branch point (-1/e) has no
                            // real W0 value; clamp to the branch-point limit
                            // W0(-1/e) = -1 instead of driving the Halley
                            // iteration into a non-convergent region.
                            // Original: NtNew = NtOld + (Eps/F*LambertW(AA) - BB);
                            if (AA < (-1.0/exp(1.0))) {
                                NtNew = NtOld + (Eps/F*(-1.0) - BB);
                                static int warnLambert2 = 0;
                                if (++warnLambert2 <= 10)
                                    cout<<"\n\n\t\tWarning: Value for LAMBERT f-n is too small\n\n";
                                else if (warnLambert2 == 11)
                                    cout<<"\n\n\t\tWarning: Value for LAMBERT f-n is too small (further warnings suppressed)\n\n";
                            }
                            else
                                NtNew = NtOld + (Eps/F*LambertW(AA) - BB);
                            RuNew = Ksat*exp(-F*NtNew);
                        }

                        // Only for the cases when Nf has not reached WT

                        if (NwtNew != 0.0) {
                            AA = get_Lower_Moist(NfNew, NwtNew);
                            BB=(NfNew-NtNew)*Ths + Eps/F*(Ths-Thr)*(1.0-exp(-F*NtNew/Eps))+Thr*NtNew;

                            if (MuNew > (AA + BB)) {
                                ThRiNf = get_InitMoist_depthZ(NfNew);
                                // To redistribute imbalance
                                NfNew += (MuNew - (AA+BB))/(Ths - ThRiNf);

                                if (NfNew >= (NwtNew+Psib)) {
                                    // Front reached the capillary fringe:
                                    // collapse the wedge into the water table
                                    // conservatively, solve for the WT whose
                                    // hydrostatic profile holds exactly MuNew,
                                    // exporting any above-saturation excess.
                                    // The original code resets here (MuNew =
                                    // get_Total_Moist(NtNew-Psib) for NtNew>0,
                                    // MuNew = 0 with full saturation for
                                    // NtNew==0) discarded the difference and
                                    // destroyed wedge water. CJC2026
                                    if (MuNew >= NwtOld*Ths) {
                                        sbsrf += (MuNew - NwtOld*Ths)/dt;
                                        NwtNew = NfNew = NtNew = 0.0;
                                        MuNew = MiNew = 0.0;
                                        RiNew = RuNew = 0.0;
                                    }
                                    else {
                                        Mdelt = Ths*NwtOld - MuNew;
                                        NwtNew = Newton(Mdelt, NwtOld);
                                        MiNew = get_Total_Moist(NwtNew);
                                        MuNew = MiNew;
                                        NfNew = NtNew = NwtNew;
                                        RiNew = RuNew = 0.0;
                                    }
                                }
                                else {
                                    // The wedge excess was absorbed by
                                    // advancing NfNew above; MuNew already
                                    // holds it. The original code, MuNew = AA+BB
                                    // reset destroyed that excess. CJC2026
                                    if (NtNew == 0.0)
                                        RuNew = Ksat*F*NfNew/expm1(F*NfNew);
									}
								}
							}
						}
					}

					QpOut = 0.0;
				if ((NtNew == NfNew) && (NfNew != NwtNew) && (NfNew > 10.0))
					QpOut = get_UnSat_LateralFlow(NfNew, RuNew, ce->getVEdgLen()*1000.0);

                if ((NfNew-NtNew) > 1.0)
                    QpOut = get_Sat_LateralFlow(NtNew, NfNew, RuNew, ce->getVEdgLen()*1000.0);

                QpOut *= Sin;

                break;


				//-------------------
			case Perched_SurfSat:

				if (R <= RADAR )
					IntStormVar += dt;
				else
					IntStormVar = 0.0;

				if (Ksat == 0.0) {
					MuNew=MiNew=0.0;
					RiNew=RuNew=0.0;
					NfNew=NtNew=0.0;
					NwtNew=NwtOld;
					QpOut=0.0;
					hsrf+=R*Cos;
				}
					else  {
						NwtNew = NwtOld;

						// Pearched surface saturation has already developed
						if ((NtOld == 0.0) && (NfOld > 0.0) && (NfOld != NwtOld)) {
							ThRiNf = get_InitMoist_depthZ(NfOld);
							ThReNf = Ths;
							set_Suction_Term(NfOld);
							qn = Ksat*F*NfOld/expm1(F*NfOld);
							qn *= (Cos + G/NfOld);

							// The newcoming moisture is not taken into account
							// define the flux based on old state variables

							R1 = (R + (QpIn-QpOut))*Cos;

							// Infiltration excess runoff is generated
							if (R1 >= qn) {
								hsrf += (R1 - qn);
								R1 -= hsrf;
							}

							else {
								static int warnStateDef2 = 0;
								if (++warnStateDef2 <= 10)
									cout<<"\nWarning: WRONG state def.: R < Keqviv-RiOld*Cos: id = "
									<<cn->getID()<<"\n\n";
								else if (warnStateDef2 == 11)
									cout<<"\nWarning: WRONG state def.: R < Keqviv-RiOld*Cos (further warnings suppressed)\n\n";
                            }

							NtNew = 0.0;
							RiNew = RiOld;
							MiNew = MiOld;
							MuNew = MuOld + R1*dt;
							NfNew = NfOld + dt*(qn-RiNew*Cos)/(Ths-ThRiNf);
						}

						if (NfNew >= 1.0) {
							RuNew = Ksat*F*NfNew/expm1(F*NfNew);
							if (RuNew > Ksat)
								RuNew = Ksat;
						}
						else
							RuNew = Ksat;

						// Check to see if WT is needed to be modified
						if ((fabs(NfNew-(NwtNew+Psib)) <= 1.0E-3) || (NfNew > (NwtNew+Psib))) {
							Mperch = MuNew - NwtOld*Ths;
							if (Mperch >= 0.0) {
								sbsrf += Mperch/dt;
								NfNew = NtNew = NwtNew = 0.0;
								MiNew = MuNew = 0.0;
								RiNew = RuNew = 0.0;
							}
							else {
								NwtNew = Newton(fabs(Mperch), 0.0);
								MuNew = MiNew = get_Total_Moist(NwtNew);
								RiNew = RuNew = 0.0;
								NfNew=NtNew=NwtNew;
							}
						}

						// Only for the case when Nf has not reached water table
						if (NwtNew != 0.0 && NfNew >= 1.0 && NfNew != NwtNew) {
							AA = get_Lower_Moist(NfNew, NwtNew);
							BB = NfNew*Ths;
							if (MuNew >= (AA+BB)) {
								ThRiNf = get_InitMoist_depthZ(NfNew);
								NfNew += (MuNew - (AA+BB))/(Ths-ThRiNf);
								if (NfNew >= (NwtNew+Psib)) {
									if (MuNew >= (NwtOld*Ths))
										sbsrf += (MuNew - NwtOld*Ths);
									NwtNew = 0.0;
									RiNew = RuNew = 0.0;
									MiNew = MuNew = 0.0;
									NfNew = NtNew = 0.0;
								}
								else {
									MuNew = AA+BB;
									RuNew = Ksat*F*NfNew/expm1(F*NfNew);
								}
							}
						}
						QpOut = 0.0;
						if ((NfNew-NtNew) > 10.0) {
							QpOut = get_Sat_LateralFlow(NtNew, NfNew, RuNew, ce->getVEdgLen()*1000.0);
							QpOut *= Sin;
						}
					}
					break;

            default:
                cerr<<"No Case Identified in tHydroModel::UnSaturatedZone"<<endl; //WR debug: Clang warning
  }

		// End of case handling for unsat+perched timestep

		// Step 4: Print Statements
		//---------------------------------------------

		if (simCtrl->Verbose_label == 'Y') {
			if (NtNew<0 || NtOld<0) { cout <<"\nWarning: Top Front < 0\n\n";
				cout<<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld
					<<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout<<"NtNew = "<< NtNew << "\tR1 = " << R1 <<"; id = "<<cn->getID()<<endl;}

			if (NfNew<0 || NfOld<0) { cout<<"\nWarning: Wetting Front < 0\n\n";
				cout<<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld
					<<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout<<"NfNew = " <<NfNew<<"\tR1 = "<<R1<<"; id = "<<cn->getID()<<endl;}

			if (NtNew>NwtNew) { cout <<"\nWarning: Top Front > WT depth\n\n";
				cout<<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld
					<<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout<<"NtNew = " << NtNew <<", NwtNew = " << NwtNew
					<<"\tR1 = " << R1 <<"; id = "<<cn->getID()<<endl; }

			if (NfNew>NwtNew) { cout<<"\nWarning: Top Front > WT depth\n\n";
				cout<<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld
					<<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout<<"NfNew = " << NfNew <<", NwtNew = " << NwtNew
					<<"\tR1 = " << R1 <<"; id = "<<cn->getID()<<endl; }

			if (NtNew>NfNew) { cout<<"\nWarning: Top Front > Wet front\n\n";
				cout<<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld
					<<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout<<"NtNew = "<<NtNew<<", NfNew = "<<NfNew
					<<"\tR1 = "<<R1 <<"; id = "<<cn->getID()<<endl; }

			if (NwtNew<0 || NwtOld<0) { cout<<"\nWarning: Water Table < 0\n\n";
				cout<<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld
					<<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout<<"NwtNew = "<<NwtNew<<"\tR1 = "<<R1<<"; id = "<<cn->getID()<<endl; }

			if (hsrf<0.0)  cout<<"\nWarning: RUNOFF - hsrf  < 0: id = "<<cn->getID()<<"\n";
			if (sbsrf<0.0) cout<<"\nWarning: RUNOFF - sbsrf < 0: id = "<<cn->getID()<<"\n";
			if (psrf<0.0)  cout<<"\nWarning:  RUNOFF - psrf  < 0: id = "<<cn->getID()<<"\n";

			if (hsrf>999999.) {
				cout<<"\nWarning! RUNOFF comp.- hsrf  "
				<<"> 999999: id = "<<cn->getID()<<"\n\t\t-> Assigned to zero\n\n";
				hsrf = 0.0;
			}
			if (sbsrf>999999.) {
				cout<<"\nWarning! RUNOFF comp.- sbsrf  "
				<<"> 999999: id = "<<cn->getID()<<"\n\t\t-> Assigned to zero\n\n";
				sbsrf = 0.0;
			}
			if (psrf>999999.) {
				cout<<"\nWarning! RUNOFF comp.- psrf  "
				<<"> 999999: id = "<<cn->getID()<<"\n\t\t-> Assigned to zero\n\n";
				psrf = 0.0;
			}

			if (MuNew < MiNew)
				cout<<"\nWarning: MuNew < MiNew -- ENDDD --: id = "<<cn->getID()<<"\n\n";

			if ((RuNew < RiNew) && (RuNew != 0.0))
				cout<<"\nWarning: RuNew < RiNew: id = "<<cn->getID()<<"\n\n";

			if (RuNew < 0.0)
				cout<<"\nWarning: RuNew < 0: id = "<<cn->getID()<<"\n\n";

			if (RiNew > 0.0) {
				cout<<"\nWarning: NON-Equilibrium rate - RiNew>0: id =";
				cout<<cn->getID()<<"\n\n";}

			if ((NfNew < NfOld) && (NfNew != 0.0) &&
				(NfNew != NwtNew) && (NfOld != NwtOld))
				cout<<"\nWarning: NfNew < NfOld: id = "<<cn->getID()<<"\n\n";
		}


		// Adding runoff contributions from Unsat Zone dynamics
		esrf = psrf;
		esrf = esrf/Cos;

		srf = hsrf + sbsrf + psrf;

		// 'Cos' factor in order to transform
		// to the rate for a horizontal plane
		srf /= Cos;
		hsrf  /= Cos;
		sbsrf /= Cos;
		psrf  /= Cos;

		// Update Variables
		double dM1, dM2, ThSurf0, ThRoot0, QpOut0;

		ThSurf0 = cn->getSoilMoistureSC();
		ThRoot0 = cn->getRootMoistureSC();
		QpOut0  = cn->getQpout();

		cn->setNwtNew(NwtNew);
		cn->setMuNew(MuNew);
		cn->setMiNew(MiNew);
		cn->setNfNew(NfNew);
		cn->setNtNew(NtNew);
		cn->setRiNew(RiNew);
		cn->setRuNew(RuNew);
		cn->setQpout(QpOut);
		cn->setIntStormVar(IntStormVar);

		cn->addSrf_Hr(srf*dt); // WR debug 02062024-- this scaled by dt so appears that it is only a total (mm), but its re-set every hour it the total runoff in an hour as written to .pixel file
		cn->setsrf(srf*dt);
		cn->addCumSrf(srf*dt);
		cn->setsbsrf(sbsrf*dt);
		cn->sethsrf(hsrf*dt);
		cn->setpsrf(psrf*dt);
		cn->setesrf(esrf*dt);

		// Soil moisture in the top 10 cm
		ThSurf = ComputeSurfSoilMoist(surfaceSoilDepth);
		cn->setSoilMoisture( ThSurf );
		cn->setSoilMoistureSC( ThSurf/Ths );

		// Soil moisture in the unsaturated zone
		if (NwtNew)
			cn->setSoilMoistureUNSC( MuNew/NwtNew/Ths );
		else
			cn->setSoilMoistureUNSC( 1.0 );

		// Need to divide by the total # of time steps elapsed
		AA = (double)timer->getElapsedSteps(timer->getCurrentTime());
		// The integer part - surface SM
		BB = floor(cn->getAvSoilMoisture())*1.0E-4;
		// The decimal part - root SM
		Mdelt = (cn->getAvSoilMoisture() - floor(cn->getAvSoilMoisture()))*1.0E+1;
		cn->setAvSoilMoisture(0.0);
		cn->setAvSoilMoisture(floor((BB*AA + ThSurf/Ths)/(AA+1)*1.0E+4));

		// Estimate average root soil moisture
		{
			double rzDepth = std::min(cn->getRootZoneDepth(), DtoBedrock);
			ThSurf = ComputeSurfSoilMoist(rzDepth);
			cn->setRootMoisture( ThSurf );
			cn->setRootMoistureSC( ThSurf/Ths );
			// RootCutoff = get_Upper_Moist(rzDepth, DtoBedrock)/rzDepth depends
			// only on the rootzone depth, bedrock depth and the static soil
			// properties, so recompute it only when the rootzone depth changes
			// (e.g. dynamic land-use grids) instead of every time step: the
			// pow() calls inside get_Upper_Moist showed up in runtime profiles
			if (rzDepth != cn->getRootCutoffDepth()) {
				cn->setRootCutoff( get_Upper_Moist(rzDepth, DtoBedrock) / rzDepth );
				cn->setRootCutoffDepth(rzDepth);
			}
		}
		cn->addAvSoilMoisture((Mdelt*AA + cn->getRootMoisture()/Ths)/(AA+1.0)*1.0E-1);

		// Set other variables
		cn->setRecharge((NwtNew-NwtOld)*Ths/(Cos*dt));
		cn->setUnSatFlowOut(QpOut*1.0E-6/(cn->getVArea()));
		cn->setUnSatFlowIn(QpIn*1.0E-6/(cn->getVArea()));

		//Saturation frequencies
		if (hsrf>0.0)
			cn->hsrfOccur=cn->hsrfOccur+floor(hsrf*1.0E+3)+1.0E-6;
		if (sbsrf>0.0)
			cn->sbsrfOccur=cn->sbsrfOccur+floor(sbsrf*1.0E+3)+1.0E-6;
		if (psrf>0.0)
			cn->psrfOccur=cn->psrfOccur+floor(psrf*1.0E+3)+1.0E-6;

		// Saturation occurrence counter. The original test was simply
		//     if (ComputeSurfSoilMoist(1.0)/Ths > 0.999) satOccur++;
		// i.e. the average moisture in the top 1 mm exceeds 99.9% of
		// saturation. ComputeSurfSoilMoist is expensive (pow/exp/log), so in
		// the most common state like initial hydrostatic profile (no wetting
		// front) with a deep enough water table it is skipped: moisture
		// only increases with depth, so the top 1 mm cannot be near
		// saturation when NwtNew >= NwtSatThresh (see ComputeNwtSatThreshold
		// for the derivation). The counter values are unchanged.
		if (cn->getNwtSatThresh() < 0.0)
			cn->setNwtSatThresh(ComputeNwtSatThreshold());
		if (!((NfNew == 0.0 || NfNew == NwtNew) && NwtNew >= cn->getNwtSatThresh())) {
			if (ComputeSurfSoilMoist(1.0)/Ths > 0.999)
				cn->satOccur = cn->satOccur + 1;
		}

		// The term with 'sbsrf' is not exactly right...
		cn->RechDisch=cn->RechDisch+((NwtOld-NwtNew)+sbsrf*dt*Cos/Ths)*1.0E-3;

		// Set QpIn to downslope cell
		dnode = (tCNode *)ce->getDestinationPtrNC();
		dnode->addQpin(QpOut);

#ifdef PARALLEL_TRIBS
  // Send Qpin from reach outlet node to
  // reach head node, if on another processor
  int idin = cn->getReach();
  int idout = dnode->getReach();
  if (idin != idout && tGraph::hasDownstream(idin)) {
    tGraph::sendQpin(idin, dnode, QpOut);

      // If runon option is set, send srf to node above
      // outlet to head
      if (RunOnoption)
          tGraph::sendRunFlux(cn);
  }

#endif

		// The following units are [m^3]
		Stok += srf*dt*cn->getVArea()/1000.0;
        TotGWchange += (NwtNew-NwtOld)*Ths*cn->getVArea()/(Cos*1000.0);
		TotMoist += (MuNew-MuOld)*cn->getVArea()/(Cos*1000.0);

		// Code related to soil moisture study by Valerio Noto
		double AreaF, fs, ft, fc;
        fs = ft = fc = 0;

		QpOut0 *= ((1.0E-6)/(cn->getVArea()));

		// GW rise
		if (NwtNew < NwtOld) {
			fs = fabs(MuNew+Ths*(NwtOld-NwtNew)-MuOld)/dt-(QpIn-QpOut0)*Cos+(EvapSoi+EvapVeg);
			ft = fabs(QpIn-QpOut0)*Cos;
			fc = fabs(EvapSoi+EvapVeg);
		}
		// GW drop
		else if (NwtNew > NwtOld) {
			if (RuOld > 0)
				fs = RuOld;
			else
				fs = cn->getNetPrecipitation();
			ft = fabs(QpIn-QpOut0)*Cos;
			fc = fabs(EvapSoi+EvapVeg);
		}
		// GW is zero
		else if (NwtOld == 0.0 && NwtNew == 0.0) {
			ft = (QpIn-QpOut0)*Cos;
			fc = fabs(EvapSoi+EvapVeg);
			if (ft-fc > 0.0)
				fs = 0.0;
			else
				fs = fabs(ft+fc);
			ft = fabs(ft);
		}
		// GW does not change
		else if (fabs(NwtNew-NwtOld) < 1.0E-9) {
			fs = RuNew;
			ft = fabs(QpIn-QpOut0)*Cos;
			fc = fabs(EvapSoi+EvapVeg);
		}

		// Relative contribution to the basin-scale
		AreaF = (cn->getVArea())/BasArea;
		fs *= AreaF;
		ft *= AreaF;
		fc *= AreaF;

		fSoi100 += fs;
		fTop100 += ft;
		fClm100 += fc;
		dM100   += ((cn->getSoilMoistureSC()) - ThSurf0)*AreaF*Ths*100.;
		dMRt    += ((cn->getRootMoistureSC()) - ThRoot0)*AreaF*Ths*1000.;
		mTh100  += ((cn->getSoilMoistureSC()) + ThSurf0)/2.0*AreaF;
		mThRt   += ((cn->getRootMoistureSC()) + ThRoot0)/2.0*AreaF;

		/*
			if (Pixel_State == Storm_Evol)
  	     cout << "\tCASE: *** Storm_Evol ***\n";
         else if (Pixel_State == WTStaysAtSurf)
  	     cout << "\tCASE: *** WTStaysAtSurf ***\n";
         else if (Pixel_State == WTDropsFromSurf)
  	     cout << "\tCASE: *** WTDropsFromSurf ***\n";
         else if (Pixel_State == WTGetsToSurf)
  	     cout << "\tCASE: *** WTGetsToSurf ***\n";
         else if (Pixel_State == Storm_Unsat_Evol)
  	     cout << "\tCASE: *** Storm_Unsat_Evol ***\n";
         else if (Pixel_State == Perched_Evol)
  	     cout << "\tCASE: *** Perched_Evol ***\n";
         else if (Pixel_State == Perched_SurfSat)
  	     cout << "\tCASE: *** Perched_SurfSat ***\n";
         else if (Pixel_State == StormToInterTransition)
  	     cout << "\tCASE: *** StormToInterTransition ***\n";
         else if (Pixel_State == ExactInitial)
  	     cout << "\tCASE: *** ExactInitial ***\n";
         else if (Pixel_State == IntStormBelow)
  	     cout << "\tCASE: *** IntStormBelow ***\n";

		 cout<<"\tQpIn = "<<QpIn<<"  QpOut0 = "<<QpOut0<<endl;
		 cout<<"\tdm100 = "<<((cn->getSoilMoistureSC()) - ThSurf0)*Ths*100.
		 <<"  fSoi100 = "<<fs
		 <<"  fTop100 = "<<ft
		 <<"  fClm100 = "<<fc
		 <<"  NfOld = "<<NfOld<<"  NwtOld = "<<NwtOld<<"  NwtNew = "<<NwtNew
		 <<endl<<flush;
		 */

		PrintNewVars( cn, Ractual );

		R = R1 = -999.0;
		id++;

		cn = nodIter.NextP();

  } // End of long while node loop

#ifdef PARALLEL_TRIBS
  // Send Qstrm, NwtOld, and NfOld to upstream reach outlet nodes
  tGraph::sendOverlap();
  // Receive as Qstrm, NwtOld, and NfOld from downstream reach head nodes
  tGraph::receiveOverlap();
#endif

	if (simCtrl->Verbose_label == 'Y') {
		// Output water balance information every 100 hours or at the end
		AA = timer->getCurrentTime();
		if (nodeList ||
			(AA/100.-floor(AA/100.)) < (timer->getTimeStep()-timer->getTimeStep()/2.)/100. ||
			timer->RemainingTime() == 0.0) {

			cout<<"\t\tRUNOFF = "<<Stok<<" M^3"<<endl<<flush;
			cout<<"\t\tTOTRAIN = "<<TotRain<<" M^3"<<endl<<flush;
			cout<<"\t\tMOISTCHN = "<<TotMoist<<" M^3"<<endl<<flush;
			cout<<"\t\tGWCHANGE = "<<TotGWchange<<" M^3"<<endl<<flush;
		}
	}

	// In order to properly take care of the actual generated runoff in the element
	if (RunOnoption) {
		cn = nodIter.FirstP();
		while ( nodIter.IsActive() ) {
			SetCellRunon( cn, cn->getSrf()/dt, cn->getRunOn()/dt, dt, 1 );
			cn = nodIter.NextP();
		}
	}
}

//=========================================================================
//
//
//                  Section 4: tHydroModel Get/Set Functions
//
//
//=========================================================================

/*************************************************************************
**
**  tHydroModel::PowPsibOverNwt( double )
**
**  Returns pow((-Psib/Nwt), PoreInd), the Brooks-Corey term shared by
**  get_Total_Moist and get_Upper_Moist. Those functions are called
**  several times per node per time step with the same water table depth
**  Nwt, so the last result is cached and reused: pow() is expensive and
**  dominated the runtime profile. The cache keys on all three inputs
**  (Nwt, Psib, PoreInd), so a stale value can never be returned. The
**  result is always identical to calling pow() directly.
**
*************************************************************************/
double tHydroModel::PowPsibOverNwt(double Nwt)
{
	if (Nwt != cachedPowNwt || Psib != cachedPowPsib ||
		PoreInd != cachedPowPoreInd) {
		cachedPowNwt = Nwt;
		cachedPowPsib = Psib;
		cachedPowPoreInd = PoreInd;
		cachedPowVal = pow((-Psib/Nwt), PoreInd);
	}
	return cachedPowVal;
}

/*************************************************************************
**
**  tHydroModel::get_Total_Moist( double )
**
**  Allows to get the total moisture value in the initial profile
**
*************************************************************************/
double tHydroModel::get_Total_Moist(double Nwt)
{
	double dM = 0.0;
	if (Nwt >= fabs(Psib)) {
		if (PoreInd >(1.0-1.0E-6) && PoreInd <(1.0+1.0E-6))
			dM = Thr*(Nwt+Psib)-fabs(Psib)*(Ths-Thr)*log((-Psib)/Nwt)-Ths*Psib;
		else
			dM = Thr*(Nwt+Psib)-Ths*Psib-(Ths-Thr)/(PoreInd-1.0)*
				(Psib + Nwt*PowPsibOverNwt(Nwt)); // = pow((-Psib/Nwt),PoreInd), cached
	}
	else if (Nwt == 0.0)
		dM = 0.0;
	else {
		if (simCtrl->Verbose_label == 'Y') {
			cout<<"\nWarning: Improper use of 'get_Total_Moist'!"<<endl;
			cout<<"\nNwt = "<<Nwt<<";\tNODE ID = "<<ID<<endl<<endl;
		}
		// Instead of changing this in the code
		NwtNew = 0.0;
	}

	return dM;
}

/*************************************************************************
**
**  tHydroModel::get_Upper_Moist( double, double )
**
**  Allows to get total moisture in the initial profile above the point
**
**  Arguments: Nf  - depth above which the total moisture content
**                   is to be evaluated
**             Nwt - water table depth of the initial moisture profile
**
*************************************************************************/
double tHydroModel::get_Upper_Moist(double Nf, double Nwt)
{
	double dM;
	if (Nf <= (Nwt+Psib)) {
		if (PoreInd >(1.0-1.0E-6) && PoreInd <(1.0+1.0E-6))
			dM = Thr*Nf + fabs(Psib)*(Ths-Thr)*log(Nwt/(Nwt-Nf));
		else
			dM = Thr*Nf - (Ths-Thr)/(PoreInd-1.0)*
				((Nf-Nwt)*pow((Psib/(Nf-Nwt)),PoreInd) + Nwt*PowPsibOverNwt(Nwt));
	}
	else {
		if (PoreInd >(1.0-1.0E-6) && PoreInd <(1.0+1.0E-6))
			dM = get_Total_Moist(Nwt) + Ths*Psib + Ths*(Nf-(Nwt+Psib));
		else
			dM = Thr*(Nwt+Psib) - (Ths-Thr)/(PoreInd-1.0)*
				(Psib + Nwt*PowPsibOverNwt(Nwt)) + Ths*(Nf-(Nwt+Psib));
	}
	return dM;
}

/*************************************************************************
**
**  tHydroModel::get_Lower_Moist( double, double )
**
**  Allows to get total moisture in the initial profile below point Nf
**
*************************************************************************/
double tHydroModel::get_Lower_Moist(double Nf, double Nwt) const
{
	double dM;
	if (Nf <= (Nwt+Psib)) {
		if (PoreInd >(1.0-1.0E-6) && PoreInd <(1.0+1.0E-6))
			dM = Thr*(Nwt+Psib-Nf)-fabs(Psib)*(Ths-Thr)*log((-Psib)/(Nwt-Nf))-Ths*Psib;
		else
			dM = Thr*(Nwt+Psib-Nf) - Ths*Psib - (Ths-Thr)/(PoreInd-1.0)*
				(Psib + (Nwt-Nf)*pow((Psib/(Nf-Nwt)),PoreInd));
	}
	else
		dM = Ths*(Nwt-Nf);
	return dM;
}

/*************************************************************************
**
**  tHydroModel::get_InitMoist_depthZ( double )
**
**  Allows to get a moisture value at depth Nf
**
*************************************************************************/
double tHydroModel::get_InitMoist_depthZ(double Nf) const
{
	if (Nf >= (NwtOld+Psib))
		return Ths;
	else
		return (Thr+(Ths-Thr)*pow((Psib/(Nf-NwtOld)),PoreInd));
}

/*************************************************************************
**
**  tHydroModel::get_EdgeMoist_depthZ( double )
**
**  Allows to get a moisture value in the penetrating unsaturated edge
**
*************************************************************************/
double tHydroModel::get_EdgeMoist_depthZ(double Nf) const
{
	return(Thr + exp(F*Nf/Eps)*(Ths-Thr)*pow((RuNew/Ksat),(1.0/Eps)));
}

/*************************************************************************
**
**  tHydroModel::set_Suction_Term( double )
**
**  Sets value for the suction term in the modified Green-Ampt model
**
**  Note: ThRiNf & ThReNf must be defined before the function run
**
**  Computes, with PoreIndF = PoreInd*exp(-F*Nf/2) being the pore-size
**  distribution index adjusted for the conductivity decay with depth:
**
**     Se(theta) = ((theta - Thr)/(Ths - Thr)) ^ (3 + 1/PoreIndF)
**     G = -Psib*(Se0 - SeIn)/(3*PoreIndF + 1)        [mm]
**
**  where SeIn uses theta = ThRiNf, Se0 uses theta = ThReNf, and Se0 = 1
**  in the saturated phase (ThReNf == Ths). The exp() term and the Se
**  exponent are evaluated once and reused.
**
*************************************************************************/
void tHydroModel::set_Suction_Term(double Nf)
{
	double PoreIndF = PoreInd*exp(-F*Nf/2.);
	double SeExp    = 3. + 1./PoreIndF;

	SeIn = pow(((ThRiNf-Thr)/(Ths-Thr)), SeExp);

	if (ThReNf == Ths) { // <--- Saturated phase
		G = -Psib*(1. - SeIn)/(3.*PoreIndF + 1.);  // UNITS: mm
	}
	else {
		Se0  = pow(((ThReNf-Thr)/(Ths-Thr)), SeExp);
		G = -Psib*(Se0 - SeIn)/(3.*PoreIndF + 1.);  // UNITS: mm
	}
}

/*************************************************************************
**
**  tHydroModel::get_RechargeRate( double,double )
**
*************************************************************************/
double tHydroModel::get_RechargeRate(double dM, double Nf) const
{
	double BB;
	BB = (dM - Thr*Nf)/(Ths-Thr);
	BB = BB*F/(Eps*expm1(F*Nf/Eps));
	return(Ksat*pow(BB, Eps));
}

/*************************************************************************
**
**  tHydroModel::get_UnSat_LateralFlow(double, double, double)
**
*************************************************************************/
double tHydroModel::get_UnSat_LateralFlow(double Nf, double Ru,
										  double veLength) const
{
	return(Nf*Ru*(UAr-1.0)*veLength); // UNITS: mm^3/hour
}

/*************************************************************************
**
**  tHydroModel::get_Sat_LateralFlow(double, double, double, double )
**
*************************************************************************/
double tHydroModel::get_Sat_LateralFlow(double Nt, double Nf,
										double Ru, double veLength) const
{
	double BB;
	if (Nt == 0.0)
		BB = (Ksat*UAr/F*(1.-exp(-F*Nf))-Ksat*F*Nf*Nf/expm1(F*Nf))*veLength;
	else {
		BB = Nt*Ru*(UAr-1.0) + Ksat*UAr/F*(exp(-F*Nt)-exp(-F*Nf));
		BB = (BB - Ksat*F*(Nf-Nt)*(Nf-Nt)/(exp(F*Nf)-exp(F*Nt)))*veLength;
	}
	return BB;  // UNITS: mm^3/hour
}

/*************************************************************************
**
**  tHydroModel::getTransmissivityFinD( double )
**
*************************************************************************/
double tHydroModel::getTransmissivityFinD(double Nwt) const
{
	return( Ar * Ksat * (exp(-F*Nwt) - exp(-F*DtoBedrock))/F );
}

/*************************************************************************
**
**  tHydroModel::getTransmissivityInfD( double )
**
*************************************************************************/
double tHydroModel::getTransmissivityInfD(double Nwt) const
{
	return( Ar * Ksat * exp(-F*Nwt)/F );
}

//=========================================================================
//
//
//                  Section 5: tHydroModel Reset Functions
//
//
//=========================================================================

/*************************************************************************
**
**  tHydroModel::Reset( )
**
*************************************************************************/
void tHydroModel::Reset()
{
	tCNode * cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList());
	cn = nodIter.FirstP();
	while ( nodIter.IsActive() ) {
		cn->setNwtOld(cn->getNwtNew());
		cn->setMuOld(cn->getMuNew());
		cn->setMiOld(cn->getMiNew());
		cn->setNfOld(cn->getNfNew());
		cn->setNtOld(cn->getNtNew());
		cn->setRuOld(cn->getRuNew());
		cn->setRiOld(cn->getRiNew());
		cn->setQpin(0.0);
		cn->setsrf(0.0);
		cn->setsbsrf(0.0);
		cn->sethsrf(0.0);
		cn->setpsrf(0.0);
		cn->setsatsrf(0.0);
		cn->setGwaterChng( 0.0 );
		cn = nodIter.NextP();
	}
}

/*************************************************************************
**
**  tHydroModel::ResetGW( )
**
*************************************************************************/
void tHydroModel::ResetGW()
{
	tCNode * cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList());
	cn = nodIter.FirstP();
	while ( nodIter.IsActive() ) {
		cn->setNwtOld(cn->getNwtNew());
		cn->setMuOld(cn->getMuNew());
		cn->setMiOld(cn->getMiNew());
		cn->setNfOld(cn->getNfNew());
		cn->setNtOld(cn->getNtNew());
		cn->setRuOld(cn->getRuNew());
		cn->setRiOld(cn->getRiNew());
		cn->setesrf(0.0);
		cn->setsatsrf(0.0);
		cn->setGwaterChng( 0.0 );
		cn = nodIter.NextP();
	}
	}

//=========================================================================
//
//
//                  Section 6: tHydroModel ComputeFluxes Functions
//
//
//=========================================================================


/***************************************************************************
**
**  tHydroModel::ComputeFluxesNodes1D()
**
**  This function simply takes the steepest direction in GW topography
**  for the current node. The flux is computed based on the same assumptions
**  as before. The difference though is that the direction of flow is SINGLE
**  and it is determined looping through the *nodes*, not the *edges*.
**
***************************************************************************/
void tHydroModel::ComputeFluxesNodes1D()
{
	tCNode * cn;
	tCNode * cnorg;
	tCNode * cndest;
	tEdge  *ce, *firstedg, *nbredg;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );

	int    ctr;
	double slp;        			    //Steepest slope found so far
	double Cos1, Cos2;
	double Width, Transmissivity, WTSlope;
	double thisWTElevation, nextWTElevation;  //Absolute elevation of WT, m abs.

	for ( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() ) {
		WTSlope = -999.0;
		cnorg = cn;
		// Originally Cos1 = cos(atan(slope)); identical to 1/sqrt(1+s^2)
		// but avoids two expensive trig calls per node
		double slpOrg = (cnorg->getFlowEdg())->getSlope();
		Cos1 = 1.0/sqrt(1.0 + slpOrg*slpOrg);

        // Giuseppe 2016 - Begin changes to allow reading soil properties from grids
		//soilPtr->setSoilPtr( cnorg->getSoilID() );
        //Psib = soilPtr->getSoilProp(5);
        Psib = cnorg->getAirEBubPres(); // Giuseppe 2016 - Changes for soil grids
		// Giuseppe 2016 - End changes to allow reading soil properties from grids

		// Depending on whether the water table is at the surface or not,
		// define the gradient of the GW head
		if (cnorg->getNwtOld() == 0.0)
			thisWTElevation = (cnorg->getZ()) - ((-Psib)/(Cos1*1000.0));
		else
			thisWTElevation = (cnorg->getZ()) - (cnorg->getNwtOld()/(Cos1*1000.0));

		firstedg = cnorg->getFlowEdg();
		nbredg = firstedg;
		ce = firstedg;

		ctr = 0;
		while ( ce != firstedg || ctr == 0)
		{
			cndest = (tCNode *)ce->getDestinationPtrNC();
			if ( (cnorg->getBoundaryFlag()  != kOpenBoundary) &&
				 (cndest->getBoundaryFlag() != kOpenBoundary) &&
				 (cnorg->getBoundaryFlag()  != kClosedBoundary) &&
				 (cndest->getBoundaryFlag() != kClosedBoundary) ) {

				// Originally Cos2 = cos(atan(slope)), see Cos1 above
				double slpDest = (cndest->getFlowEdg())->getSlope();
				Cos2 = 1.0/sqrt(1.0 + slpDest*slpDest);

                // Giuseppe 2016 - Begin changes to allow reading soil properties from grids
				//                soilPtr->setSoilPtr( cndest->getSoilID() );
                //                Psib = soilPtr->getSoilProp(5);
                Psib = cndest->getAirEBubPres(); // Giuseppe 2016 - Changes for soil grids
				// Giuseppe 2016 - End changes to allow reading soil properties from grids

				if (cndest->getNwtOld() == 0.0)
					nextWTElevation = cndest->getZ() - ((-Psib)/(Cos2*1000.0));
				else
					nextWTElevation = cndest->getZ() - (cndest->getNwtOld()/(Cos2*1000.0));

				// Compute only positive fluxes
            DtoBedrock = cnorg->getBedrockDepth(); //SMM - 09232008
				if (thisWTElevation > nextWTElevation &&
					cnorg->getNwtOld()  <= DtoBedrock &&
					cndest->getNwtOld() <= DtoBedrock) {

					slp = (thisWTElevation - nextWTElevation)/ce->getLength();
					if (WTSlope < slp) {
						WTSlope = slp;
						nbredg = ce;
					}
				}
			}
			ce = ce->getCCWEdg();
			ctr++;
			if ( ctr > 50 ) {
				cout << "\nMesh error in Node "<<cn->getID()<<endl;
			}
		} //End of looping through edges


		if ( WTSlope > 0.0 ) {
			cndest = (tCNode *)nbredg->getDestinationPtrNC();

            // Get soil hydraulic properties
			// Giuseppe 2016 - Begin changes to allow reading soil properties from grids
            //            soilPtr->setSoilPtr( cnorg->getSoilID() );
            //            Ksat    = soilPtr->getSoilProp(1);  // Surface hydraulic conductivity
            //            Ths     = soilPtr->getSoilProp(2);  // Saturation moisture content
            //            Thr     = soilPtr->getSoilProp(3);  // Residual moisture content
            //            PoreInd = soilPtr->getSoilProp(4);  // Pore-size distribution index
            //            Psib    = soilPtr->getSoilProp(5);  // Air entry bubbling pressure
            //            F       = soilPtr->getSoilProp(6);  // Decay parameter in the exp
            //            Ar      = soilPtr->getSoilProp(7);  // Anisotropy ratio (saturated)
            //            UAr     = soilPtr->getSoilProp(8);  // Anisotropy ratio (unsaturated)
            //            porosity = soilPtr->getSoilProp(9); // Porosity
            Ksat = cnorg->getKs();  // Surface hydraulic conductivity
            Ths = cnorg->getThetaS(); // Saturation moisture content
            Thr = cnorg->getThetaR(); // Residual moisture content
            PoreInd = cnorg->getPoreSize(); // Pore-size distribution index
            Psib = cnorg->getAirEBubPres(); // Air entry bubbling pressure
            F = cnorg->getDecayF(); // Decay parameter in the exp
            Ar = cnorg->getSatAnRatio(); // Anisotropy ratio (saturated)
            UAr = cnorg->getUnsatAnRatio(); // Anisotropy ratio (unsaturated)
            porosity = cnorg->getPorosity(); // Porosity
			// Giuseppe 2016 - End changes to allow reading soil properties from grids
            Eps = 3 + 2/PoreInd;

			DtoBedrock = cnorg->getBedrockDepth();  //Local variable

			// Transmissivity (depth averaged quantity) (MM^2/HOUR)
			Transmissivity = getTransmissivityFinD( cnorg->getNwtOld() );

			// Width in the direction of flow (MM)
			Width = nbredg->getVEdgLen()*1000.0;

			// Unconfined aquifer HGL (MM^3/HOUR)
			QOut = Transmissivity * Width * WTSlope;

			// Store outgoing flux from origin
			cnorg->addGwaterChng(QOut);

			// Record incoming flux to destination
			cndest->addGwaterChng(-QOut);

			cnorg->setTransmiss(Transmissivity);
		}
		else {
			// Giuseppe 2016 - Begin changes to allow reading soil properties from grids
            //            soilPtr->setSoilPtr( cnorg->getSoilID() );
            //            Ksat    = soilPtr->getSoilProp(1);
            //            F       = soilPtr->getSoilProp(6);
            //            Ar      = soilPtr->getSoilProp(7);
            Ksat = cnorg->getKs();  // Surface hydraulic conductivity
            F = cnorg->getDecayF(); // Decay parameter in the exp
            Ar = cnorg->getSatAnRatio(); // Anisotropy ratio (saturated)
			// Giuseppe 2016 - End changes to allow reading soil properties from grids

			Transmissivity = Ar * Ksat * exp(-F*cnorg->getNwtOld())/F;
			cnorg->setTransmiss(Transmissivity);
			DtoBedrock = cnorg->getBedrockDepth();
		}
	}
	}

/***************************************************************************
**
**  tHydroModel::ComputeFluxesEdgesND()
**
**  This function uses multiple direction approach, i.e. the routine
**  partitions GW flow in all downstream neighbors of a given Voronoi
**  cell.  Width of the flow is determined based on the Voronoi side length
**
***************************************************************************/
void tHydroModel::ComputeFluxesEdgesND()
{
	tCNode * cnorg;
	tCNode * cndest;
	tEdge  * ce;
	tMeshListIter<tEdge>  edgIter( gridPtr->getEdgeList() );

	double Cos1, Cos2;
	double Width, Transmissivity, WTSlope;
	double thisWTElevation, nextWTElevation; //Absolute elevation of WT, m abs.
	double deficit; 			   //Average deficit between two edges.

	for ( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() ) {
		// Destination and Origin Nodes
		cnorg  = (tCNode *)ce->getOriginPtrNC();
		cndest = (tCNode *)ce->getDestinationPtrNC();

		// Excluding calculation of flux to the outlet point
		if ( (cnorg->getBoundaryFlag() != kOpenBoundary) &&
			 (cndest->getBoundaryFlag() != kOpenBoundary)
			 &&  (cnorg->getBoundaryFlag() != kClosedBoundary) &&
			 (cndest->getBoundaryFlag() != kClosedBoundary) ) {

			// Slopes for subsurface fl. Originally cos(atan(slope));
			// identical to 1/sqrt(1+s^2) but avoids four expensive trig
			// calls per edge
			double slpOrg  = (cnorg->getFlowEdg())->getSlope();
			double slpDest = (cndest->getFlowEdg())->getSlope();
			Cos1 = 1.0/sqrt(1.0 + slpOrg*slpOrg);
			Cos2 = 1.0/sqrt(1.0 + slpDest*slpDest);

            // Giuseppe 2016 - Begin changes to allow reading soil properties from grids
			//            soilPtr->setSoilPtr( cnorg->getSoilID() );
            //            Psib = soilPtr->getSoilProp(5);
            Psib = cnorg->getAirEBubPres(); // Air entry bubbling pressure
			// Giuseppe 2016 - End changes to allow reading soil properties from grids

			// Depending on whether the water table is at the surface or not,
			// define the gradient of the GW head
			if (cnorg->getNwtOld() == 0.0) {
                thisWTElevation = (cnorg->getZ()) - ((-Psib) / (Cos1 * 1000.0));
            }
			else {
                thisWTElevation = (cnorg->getZ()) - (cnorg->getNwtOld() / (Cos1 * 1000.0));
            }

			// Giuseppe 2016 - Begin changes to allow reading soil properties from grids
            //            soilPtr->setSoilPtr( cndest->getSoilID() );
            //            Psib = soilPtr->getSoilProp(5);
            Psib = cndest->getAirEBubPres(); // Air entry bubbling pressure
			// Giuseppe 2016 - End changes to allow reading soil properties from grids

			if (cndest->getNwtOld() == 0.0) {
                nextWTElevation = (cndest->getZ()) - ((-Psib) / (Cos2 * 1000.0));
            }
			else {
                nextWTElevation = (cndest->getZ()) - (cndest->getNwtOld() / (Cos2 * 1000.0));
            }

			// Compute only positive fluxes
         DtoBedrock = cnorg->getBedrockDepth(); //SMM - 09232008
			if (thisWTElevation > nextWTElevation &&
				cnorg->getNwtOld() <= DtoBedrock &&
				cndest->getNwtOld() <= DtoBedrock) {

				// Giuseppe 2016 - Begin changes to allow reading soil properties from grids
                //soilPtr->setSoilPtr( cnorg->getSoilID() );
                // Get soil hydraulic properties
                //                Ksat    = soilPtr->getSoilProp(1);  // Surface hydraulic conductivity
                //                Ths     = soilPtr->getSoilProp(2);  // Saturation moisture content
                //                Thr     = soilPtr->getSoilProp(3);  // Residual moisture content
                //                PoreInd = soilPtr->getSoilProp(4);  // Pore-size distribution index
                //                Psib    = soilPtr->getSoilProp(5);  // Air entry bubbling pressure
                //                F       = soilPtr->getSoilProp(6);  // Decay parameter in the exp
                //                Ar      = soilPtr->getSoilProp(7);  // Anisotropy ratio (saturated)
                //                UAr     = soilPtr->getSoilProp(8);  // Anisotropy ratio (unsaturated)
                //                porosity = soilPtr->getSoilProp(9); // Porosity
                Ksat = cnorg->getKs();  // Surface hydraulic conductivity
                Ths = cnorg->getThetaS(); // Saturation moisture content
                Thr = cnorg->getThetaR(); // Residual moisture content
                PoreInd = cnorg->getPoreSize(); // Pore-size distribution index
                Psib = cnorg->getAirEBubPres(); // Air entry bubbling pressure
                F = cnorg->getDecayF(); // Decay parameter in the exp
                Ar = cnorg->getSatAnRatio(); // Anisotropy ratio (saturated)
                UAr = cnorg->getUnsatAnRatio(); // Anisotropy ratio (unsaturated)
                porosity = cnorg->getPorosity(); // Porosity
				// Giuseppe 2016 - End changes to allow reading soil properties from grids
                Eps = 3 + 2/PoreInd;

				DtoBedrock = cnorg->getBedrockDepth();  //Local variable

				deficit  = cnorg->getNwtOld();
				deficit += cndest->getNwtOld();
				deficit  = deficit/2;           // Average Water Table depth

				WTSlope = (thisWTElevation - nextWTElevation)/ce->getLength();

				// Transmissivity (depth averaged quantity) (MM^2/HOUR)
				Transmissivity = getTransmissivityFinD( cnorg->getNwtOld() );

				// Width in the direction of flow (MM)
				Width = ce->getVEdgLen()*1000.0;

				// Constrain the GW gradient
				if (WTSlope > 1.0) {
                    WTSlope = 1.0;
                }

                // Unconfined aquifer HGL (MM^3/HOUR)
				QOut = Transmissivity * Width * WTSlope;

				// Compute Voronoi polygon shape factor and constrain the model dynamics
				if (Width) {
					deficit = cndest->getVArea()/(Width*Width*10.0E-6);
					if (deficit <= 0.1)
						QOut *= deficit;
				}

				if (cnorg->getID() == -1 || cndest->getID() == -1) {
					if (simCtrl->Verbose_label == 'Y') {
						cout<<"ORIGIN Node ID = "<<cnorg->getID()<<"  ("<<cnorg->getX()<<","
						<<cnorg->getY()<<")"<<endl<<flush;
						cout<<"DESTIN Node ID = "<<cndest->getID()<<"  ("<<cndest->getX()<<","
							<<cndest->getY()<<")"<<endl<<flush;
						cout<<"\tthisWTElevation = "<<thisWTElevation
							<<" m\tnextWTElevation = "<<nextWTElevation<<" m"<<endl<<flush;
						cout<<"\tFluxin Origin BEFORE: "<<cnorg->getGwaterChng()*1.0E-9
							<<" m^3/hr"<<endl<<flush;
						cout<<"\tFluxin Destin BEFORE: "<<cndest->getGwaterChng()*1.0E-9
							<<" m^3/hr"<<endl<<flush;
						cout<<"\tWidth = "<<Width<<" mm\tWTSlope = "<<WTSlope<<"\tTransmissivity = "
							<<Transmissivity*1.0E-6<<" m^2/hr\tQOut = "
							<<QOut*1.0E-9<<" m^3/hr"<<endl<<flush;
					}
				}

            // No need to add unless not 0.0 SMM - 09232008
            if (QOut > 0.0 || QOut < 0.0) {

				    // Record outgoing flux from origin
				    cnorg->addGwaterChng(QOut);

				    // Record incoming flux to destination
				    cndest->addGwaterChng(-QOut);

            }

				cnorg->setTransmiss(Transmissivity);
			}
			else {
				// Giuseppe 2016 - Begin changes to allow reading soil properties from grids
                //                soilPtr->setSoilPtr( cnorg->getSoilID() );
                //                Ksat    = soilPtr->getSoilProp(1);
                //                F       = soilPtr->getSoilProp(6);
                //                Ar      = soilPtr->getSoilProp(7);
                Ksat = cnorg->getKs();  // Surface hydraulic conductivity
                F = cnorg->getDecayF(); // Decay parameter in the exp
                Ar = cnorg->getSatAnRatio(); // Anisotropy ratio (saturated)
				// Giuseppe 2016 - End changes to allow reading soil properties from grids

				Transmissivity = Ar * Ksat * exp(-F*cnorg->getNwtOld())/F;
				cnorg->setTransmiss(Transmissivity);
			}
		}

		/*
			// <><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
			// Artifical case to let the planar HILLSLOPE drain through
			// the lower boundary that is composed of the stream nodes

			// Note: For the following, one needs to change the iteration loop to:
			// for ( ce=edgIter.FirstP(); !( edgIter.AtEnd() ); ce=edgIter.NextP() ) {

			else if ( (cnorg->getBoundaryFlag() == kStream) &&
					  (cndest->getBoundaryFlag() == kClosedBoundary) ) {
				int cflag;
				double al, d1, x1, y1, x2, y2;
				tArray<double> xy(2), xy1(2);
				tCNode *cnn;
				tEdge  *firstedg, *curedg;

				d1  = ce->getLength();
				xy  = cnorg->get2DCoords();
				xy1 = cndest->get2DCoords();
				x1  = xy1[0]-xy[0];
				y1  = xy1[1]-xy[1];
				x2  = 0.0;
				y2  = d1;

				// Compute the angle (radians) between the
				// reference vector (North) and flow edge
				al = acos((x1*x2 + y1*y2)/(d1*d1));
				if (xy1[0] < xy[0])
					al = 8*atan(1)-al;

				// If the aspect angle is around pi, accept this direction:
				if (al > (4*atan(1)-0.1) && al < (4*atan(1)+0.1)) {

					cflag = 0;
					firstedg = curedg = cnorg->getFlowEdg();
					while (!cflag) {
						cnn = (tCNode*)curedg->getDestinationPtrNC();
						d1 = curedg->getLength();
						xy  = cnorg->get2DCoords();
						xy1 = cnn->get2DCoords();
						x1  = xy1[0]-xy[0];
						y1  = xy1[1]-xy[1];
						x2  = 0.0;
						y2  = d1;
						alpha = acos((x1*x2 + y1*y2)/(d1*d1));
						if (xy1[0] < xy[0])
							alpha = 8*atan(1)-alpha;
						curedg = curedg->getCCWEdg();

						// If the current spoke goes North (upstream), accept it
						if (fabs(alpha) < 1E-1 || curedg == firstedg)
							cflag = 1;
					}

					// Approximate the hydraulic head in node 'cnorg' with the
					// hydraulic head in the node 'cnn': the one North of  'cnorg'
					alpha = atan( (cnn->getFlowEdg())->getSlope() );
					Cos1 = cos(alpha);
					Cos2 = 1;

					soilPtr->setSoilPtr( cnn->getSoilID() );
					Psib = soilPtr->getSoilProp(5);
					if (cnn->getNwtOld() == 0.0)
						thisWTElevation = cnn->getZ() - ((-Psib)/(Cos1*1000));
					else
						thisWTElevation = cnn->getZ() - (cnn->getNwtOld()/(Cos1*1000));

					soilPtr->setSoilPtr( cnorg->getSoilID() );
					Psib = soilPtr->getSoilProp(5);
					if (cnorg->getNwtOld() == 0.0)
						nextWTElevation = cnorg->getZ() - ((-Psib)/(Cos2*1000));
					else
						nextWTElevation = cnorg->getZ() - (cnorg->getNwtOld()/(Cos2*1000));

					// Approximate the gradient and make sure it is positive
					// i.e. the flux discharges through the node 'cnorg'
					WTSlope = ( thisWTElevation - nextWTElevation )/d1;

					if (WTSlope < 0)
						WTSlope = cnn->getFlowEdg()->getSlope();

					soilPtr->setSoilPtr( cnorg->getSoilID() );
					Ksat    = soilPtr->getSoilProp(1);
					F       = soilPtr->getSoilProp(6);
					Ar      = soilPtr->getSoilProp(7);
					DtoBedrock = cnorg->getBedrockDepth();

					// ..... Transmissivity (depth averaged quantity) (MM^2/HOUR) .....
					Transmissivity = getTransmissivityFinD( cnorg->getNwtOld() );

					Width = ce->getVEdgLen()*1000;           // (MM)
					Width = 25*1000;
					QOut = Transmissivity * Width * WTSlope; // (MM^3/HOUR)


					if (cnorg->getID() == 365) {
						cout<<"\t### WTSLOPE ### = "<<WTSlope<<endl;
						cout<<"\tFluxin Origin BEFORE: "
							<<cnorg->getGwaterChng()*1.0E-9<<" m^3/hr"<<endl<<flush;
						cout<<"\tWidth = "<<Width
							<<" mm\tWTSlope = "    <<WTSlope
							<<"\tTransmissivity = "<<Transmissivity*1.0E-6
							<<" m^2/hr\tQOut = "   <<QOut*1.0E-9<<" m^3/hr"<<endl<<flush;
					}

					cnorg->addGwaterChng(QOut);
					cnorg->setTransmiss(Transmissivity);

				}
			}
		// <><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
		*/
			}

	}

//=========================================================================
//
//
//                  Section 7: tHydroModel Saturated Zone
//
//
//=========================================================================

/*************************************************************************
**
**  tHydroModel::SetupNodeSZ(tCNode *cn)
**
**  Sets up node 'cn' at the beginning of SZ simulation
**
*************************************************************************/
void tHydroModel::SetupNodeSZ(tCNode *cn)
{
	ID = cn->getID();
	// Giuseppe 2016 - Begin changes to allow reading soil properties from grids
    //soilPtr->setSoilPtr( cn->getSoilID() );
    
    // Get soil hydraulic properties
    //    Ksat    = soilPtr->getSoilProp(1);  // Surface hydraulic conductivity
    //    Ths     = soilPtr->getSoilProp(2);  // Saturation moisture content
    //    Thr     = soilPtr->getSoilProp(3);  // Residual moisture content
    //    PoreInd = soilPtr->getSoilProp(4);  // Pore-size distribution index
    //    Psib    = soilPtr->getSoilProp(5);  // Air entry bubbling pressure
    //    F       = soilPtr->getSoilProp(6);  // Decay parameter in the exp
    //    Ar      = soilPtr->getSoilProp(7);  // Anisotropy ratio (saturated)
    //    UAr     = soilPtr->getSoilProp(8);  // Anisotropy ratio (unsaturated)
    Ksat = cn->getKs();  // Surface hydraulic conductivity
    Ths = cn->getThetaS(); // Saturation moisture content
    Thr = cn->getThetaR(); // Residual moisture content
    PoreInd = cn->getPoreSize(); // Pore-size distribution index
    Psib = cn->getAirEBubPres(); // Air entry bubbling pressure
    F = cn->getDecayF(); // Decay parameter in the exp
    Ar = cn->getSatAnRatio(); // Anisotropy ratio (saturated)
    UAr = cn->getUnsatAnRatio(); // Anisotropy ratio (unsaturated)
    porosity = cn->getPorosity(); // Porosity
	// Giuseppe 2016 - End changes to allow reading soil properties from grids
	Eps = 3 + 2/PoreInd;

	// Get dynamic variables
	NwtOld = NwtNew= cn->getNwtOld();
	MuOld  = MuNew = cn->getMuOld();
	MiOld  = MiNew = cn->getMiOld();
	NtOld  = NtNew = cn->getNtOld();
	NfOld  = NfNew = cn->getNfOld();
	RuOld  = RuNew = cn->getRuOld();
	RiOld  = RiNew = cn->getRiOld();
}

/*************************************************************************
**
**  tHydroModel::SaturatedZone(double dtGW)
**
*************************************************************************/
void tHydroModel::SaturatedZone(double dtGW)
{

#ifdef PARALLEL_TRIBS
  // Exchange Nwt (send localFlux NwtNew)
  tGraph::sendNwt();
  tGraph::receiveNwt();
#endif

	// Calculate the fluxes between cells
	ComputeFluxesEdgesND();

#ifdef PARALLEL_TRIBS
  // Exchange groundwater (send remoteFlux getGwaterChng)
  tGraph::sendGroundWater();
  tGraph::receiveGroundWater();
#endif

	// Calculate changes in water table level
	tCNode * cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );

	int cnt = 0;
	int Couple_State;
	double Area, Nstar;
	double Mdelt, dM1, dM2, dM3, ThSurf, AA, BB;
	double ThSurf0, ThRoot0, fsum, AreaF;

	double gwdm = 0.0;
	double mth100 = 0.0;
	double mthrt = 0.0;
	double AreaGW = 0.0;
	int gwcnt = 0;
	int id = 0;

	for ( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() ) {

		// Setup the node
		SetupNodeSZ( cn );

		// Initialize runoff
		srf = satsrf = 0.0;

		// Get geometry. Originally Cos = cos(atan(slope)); identical to
		// 1/sqrt(1+s^2) but avoids two expensive trig calls per node
		double edgeSlope = cn->getFlowEdg()->getSlope();
		Cos = (edgeSlope > 0.0) ? 1.0/sqrt(1.0 + edgeSlope*edgeSlope) : 1.0;

		Area = cn->getVArea();   // M^2;
		
		// Get Bedrock depth for computing Nwt
		DtoBedrock = cn->getBedrockDepth(); // added by CJC2020

		// Calculate approximate water table depth (Cos to get actual area)
		NwtNew = NwtOld + dtGW*(cn->getGwaterChng()*1.0E-6)*Cos/(Area*Ths);
		if (NwtNew < 0.0) {
			satsrf = fabs(NwtNew*Ths)/dtGW;
			NwtNew=0.0;
		}

		// Potential States
		enum {GW_Exfiltrate, GW_IntStorm_Like, GW_Initial, GW_Positive_Bal};
		Couple_State = -1000;

		if (NwtNew == NwtOld) {
            Couple_State=GW_Initial;
        }
		else if (NwtNew-NwtOld > 1.0E-3) {
            Couple_State=GW_IntStorm_Like;
        }
		else if (NwtOld-NwtNew > 1.0E-3) {
            Couple_State=GW_Positive_Bal;
        }
		if (MuOld > NwtNew*Ths) {
            Couple_State=GW_Exfiltrate;
        }
		if (Couple_State==-1000) {
            Couple_State=GW_Initial;
        }
		
		// State Switch Statements
		//---------------------------------------------


		switch (Couple_State) {
			//----------------
			case GW_Exfiltrate:
				satsrf += (MuOld - NwtNew*Ths)/dtGW;
				NwtNew=0;
				MuNew=MiNew=0.0;
				RiNew=RuNew=0.0;
				NfNew=NtNew=0.0;
				break;


				//----------------
			case GW_Initial:
				NwtNew=NwtOld;
				MuNew=MuOld;
				MiNew=MiOld;
				NfNew=NfOld;
				NtNew=NtOld;
				RuNew=RuOld;
				RiNew=RiOld;
				break;


				//----------------
			case GW_IntStorm_Like:

				// Redefine Nwt according to the new moisture deficit in the element
				// The amount of water that is extracted
				Mdelt = (NwtNew - NwtOld)*Ths;

				if (NwtOld > 0.0)
					NwtNew = Newton((Ths*NwtOld - (MiOld-Mdelt)), NwtNew);
				else
					NwtNew = Newton((Ths*NwtOld - (MiOld-Mdelt)), 0.0);

				MiNew = get_Total_Moist(NwtNew);
				RiNew = 0.0;

				// Initialized state at the beginning
				if ((NfOld == 0.0) || (NfOld == NwtOld)) {
					RuNew = 0.0;
					MuNew = MiNew;
					if (NfOld == 0.0)
						NfNew=NtNew=0.0;
					else
						NfNew=NtNew=NwtNew;
				}
					// There was a wetted wedge of moisture
					else {
						dM1 = get_Upper_Moist(NfOld, NwtNew);
						dM2 = MuOld - MiOld;
						Mdelt = dM1 + dM2;
						dM3 = Eps/F*(Ths-Thr)*(1.0 - exp(-F*NfOld/Eps)) + Thr*NfOld;
						MuNew = MiNew + dM2;

						// It is an unsaturated wedge
						if (Mdelt < dM3) {
							NtNew = NfNew = NfOld;
							// Wetting front in the same position
							// but the wedge has become "thinner"
							RuNew = get_RechargeRate(Mdelt, NfOld);
						}
						// Perched Saturated wedge or Surface Perched case
						else {
							// Mdelt is the amount of water we
							// need to subtract from the wedge
							Mdelt = get_Upper_Moist(NfOld, NwtOld);
							Mdelt -= dM1;

							// Perched Saturated wedge
							if ((NtOld > 0.0) && (NtOld < NfOld)) {
								// Note: it makes more sense to substract moisture
								// from below the wetting front rather from the top:
								// effect on surface runoff should be smaller
								NfNew = Newton(Mdelt, NwtNew, NfOld, 0);

								if (NfNew <= NtOld) {
									if ((NtOld-NfNew) >= 10.0) {
										cout<<"GW_Intstorm_Like: Warning: ";
										cout<<"Difference higher then the thresh.!!"<<endl;
									}
									else
										NfNew = NtOld+0.001;
								}
								NtNew = NtOld;
								RuNew = RuOld;
							}
							// Surface Perched case
							else if (NtOld == 0.0 && NfOld > 0.0) {
								NfNew = Newton(Mdelt, NwtNew, NfOld, 0);
								NtNew = NtOld;
								RuNew = Ksat*F*NfNew/expm1(F*NfNew);
							}
							if (NtNew > NfOld) {
								cout<<"\nWarning: Incorrect case definition:";
								cout<<" NtNew > NfOld: id = "<<cn->getID()<<"\n";
							}
						}
					}
					break;


				//----------------
			case GW_Positive_Bal:

				// Water table has not reached the wetting front yet
				if (((NwtNew+Psib) > NfOld) || (NfOld == NwtOld)) {

					MiNew = get_Total_Moist(NwtNew);
					dM1 = get_Lower_Moist(NwtNew, NwtOld);
					dM2 = MiOld - dM1;
					// Moisture imbalance: we must account for it
					Mdelt = dM1 - (MiNew-dM2);

					if (Mdelt < 0.0) {
						cout<<"\n\t\tGround Water Model: Warning! Mdelt < 0\n";
						cout<<"\t\tNwtOld = "<<NwtOld<<"  NwtNew = "<<NwtNew<<endl;
						cout<<"\t\tMiOld = " <<MiOld<<"  MiNew = " <<MiNew<<endl;
						cout<<"\t\tMdelt = " <<Mdelt<<endl;
						cout<<"\t\tid = "<<cn->getID()<<"\n";
					}

					// Redefine Nwt
					NwtNew = Newton((Ths*NwtNew - (MiNew+Mdelt)), NwtNew);
					MiNew = get_Total_Moist(NwtNew);

					// The element has been in an initialized state
					if ((NfOld == 0.0) || (NfOld == NwtOld)) {
						MuNew = MiNew;
						RuNew = RiNew = 0.0;

						if (NfOld == NwtOld)
							NfNew = NtNew = NwtNew;
						else
							NfNew = NtNew = 0.0;
					}

					else if ((NfOld > 0.0) && (NfOld != NwtOld)) {

						// After relocation of the water table, it could reach the
						// wetting front. If Nwt has not reached:
						if ((NwtNew+Psib) > NfOld) {
							Mdelt = get_Upper_Moist(NfOld, NwtNew);
							Mdelt += (MuOld - MiOld); //Amount of SM above Nf

							// Element was in the unsaturated state
							if (NtOld == NfOld) {

								// The re-adjusted moisture profile lead to perching
								// of moisture in the top layer: adjust everything else
								if (Mdelt >= NfOld*Ths) {
									// Total influx from the saturated zone
									Mdelt = -dtGW*(cn->getGwaterChng()*1.0E-6)*Cos/Area;
									NwtNew = Newton((Ths*NwtOld - (MuOld+Mdelt)), (NtOld-Psib));
									MuNew = MiNew = get_Total_Moist(NwtNew);
									RiNew = RuNew = 0.0;
									NfNew=NtNew=NwtNew;
								}
								else {
									RuNew = get_RechargeRate(Mdelt, NfOld);
									Nstar = (log(Ksat/RuNew))/F;
									// If the wedge becomes perched
									if (Nstar <= NfOld) {
										NtNew = Nstar;
										NfNew = NfOld+1.0E-5;
									}
									// It is still unsaturated
									else {
										NfNew = NtNew = NfOld;
									}
									MuNew = MiNew + (MuOld - MiOld);
									RiNew = 0.0;
								}
							}

							// Element was in perched or surface saturated state
							else if (NtOld < NfOld) {

								if (Mdelt >= NfOld*Ths) {
									Mdelt -= NfOld*Ths;
									dM3 = NwtOld;
									NwtOld = NwtNew;

									// Redistribute imbalance
									dM1 = Ths*(NwtNew+Psib-NfNew);
									//dM2 = get_Z1Z2_Moist(NfNew, NwtNew+Psib, NwtNew);
									// SKY2008Snow from AJR2007
									dM2 = get_Z1Z2_Moist(NfNew, NwtNew+Psib, NwtNew,0);

									if ((dM1-dM2-Mdelt) > 1.0E-6)
										NfNew = Newton(Mdelt, NwtNew, NfOld, 1);
									else
										NfNew = NwtNew+Psib+1.0E-6;
									NtNew = NtOld;
								}
								else {
									dM1 = get_Upper_Moist(NfOld, NwtNew);
									dM2 = get_Upper_Moist(NfOld, NwtOld);
									Mdelt = dM1 - dM2;
									if (Mdelt < 0.0)
										cout<<"\nWarning: GW Model: Mdelt< 0:id = "<<cn->getID()<<"\n";
									dM3 = NwtOld;
									NwtOld = NwtNew;
									// Redistribute imbalance
									dM1 = Ths*(NwtNew+Psib-NfNew);
									//dM2 = get_Z1Z2_Moist(NfNew, NwtNew+Psib, NwtNew);
									// SKY2008Snow from AJR2007
									dM2 = get_Z1Z2_Moist(NfNew, NwtNew+Psib, NwtNew,0);

									if ((dM1-dM2-Mdelt) > 1.0E-6)
										NfNew = Newton(Mdelt, NwtNew, NfOld, 1);
									else
										NfNew = NwtNew+Psib+1.0E-6;
									NtNew = NtOld;
								}

								// There is groundwater - unsaturated zone interaction
								if (NfNew >= (NwtNew+Psib)) {
									if (NtOld > 0.0) {
										Mdelt = -dtGW*(cn->getGwaterChng()*1.0E-6)*Cos/Area;
										NwtOld = dM3;
										NwtNew = Newton((Ths*NwtOld - (MuOld+Mdelt)),(NtOld-Psib));
										MuNew = MiNew = get_Total_Moist(NwtNew);
										RiNew = RuNew = 0.0;
										NfNew=NtNew=NwtNew;
									}
									else if (NtNew == 0.0) {
										NwtNew = NfNew = NtNew = 0.0;
										MuNew = MiNew = 0.0;
										RiNew = RuNew = 0.0;
									}
								}
								// Wetting front has not reached water table yet
								else {
									MuNew = MiNew + (MuOld - MiOld);
									RiNew = 0.0;
									if (NtNew == 0.0)
										RuNew = Ksat*F*NfNew/expm1(F*NfNew);
									else
										RuNew = RuOld;
								}
							}
						}

						// Defined Nwt has reached NfOld
						else {
							Mdelt = -dtGW*(cn->getGwaterChng()*1.0E-6)*Cos/Area;
							NwtNew = Newton((Ths*NwtOld - (MuOld+Mdelt)), NwtOld);
							MuNew = MiNew = get_Total_Moist(NwtNew);
							RiNew=RuNew=0.0;
							NfNew=NtNew=NwtNew;
						}
					}
				} // Matches (NwtNew+Psib) > NfOld || (NfOld == NwtOld)


				// New water table has risen and reached NfOld
				else if ( ((NwtNew+Psib) <= NfOld) && (NfOld != NwtOld) ) {
					MiNew = get_Total_Moist(NwtNew);
					dM1 = get_Lower_Moist(NwtNew, NwtOld);
					dM2 = MiOld - dM1;
					Mdelt = dM1 - (MiNew-dM2); //Moisture imbalance

					if (Mdelt < 0) {
						cout<<"\n\t\tGround Water Model: Warning! Mdelt < 0\n";
						cout<<"\t\tNwtOld = "<<NwtOld<<"  NwtNew = "<<NwtNew<<endl;
						cout<<"\t\tMiOld = " <<MiOld<<"  MiNew = " <<MiNew<<endl;
						cout<<"\t\tMdelt = " <<Mdelt<<endl;
						cout<<"\t\tid = "<<cn->getID()<<"\n";
					}

					// This only works for (dM1 > (MiNew-dM2))
					if (NtOld == NfOld)
						NwtNew = Newton((Ths*NwtOld - (MuOld+Ths*(NwtOld-NwtNew))), NwtNew);
					else if (NtOld < NfOld)
						NwtNew = Newton((Ths*NwtOld - (MuOld+Ths*(NwtOld-NwtNew))), (NtOld-Psib));

					MuNew = MiNew = get_Total_Moist(NwtNew);
					RuNew = RiNew = 0.0;
					NtNew = NwtNew;
					NfNew = NwtNew;
				}
				break;
            default:
                cerr<<"No case identified in tHydroModel::SaturatedZone"<<endl;
                break;

		}   //End of Switch Statements


		// Print out Statements
		if (simCtrl->Verbose_label == 'Y') {
			if (NtNew<0.0 || NtOld<0.0) {
				cout <<"\nWarning: Top Front < 0\n\n";
				cout <<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld;
				cout <<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout <<"NtNew = "<<NtNew<<"\nid = "<<cn->getID();
				cout <<"\tQIn = "<<QIn<<"  QOut = "<<QOut<<endl;
			}

			if (NfNew<0.0 || NfOld<0.0) {
				cout << "\nWarning: Wetting Front < 0\n\n\n";
				cout <<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld;
				cout <<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout <<"NfNew = "<<NfNew<<"\nid = "<<cn->getID();
				cout <<"\tQIn = "<<QIn<<"  QOut = "<<QOut<<endl;
			}

			if (NtNew>NwtNew)  {
				cout << "\nWarning: Top Front > WT depth\n\n";
				cout <<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld;
				cout <<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout <<"NtNew = "<<NtNew<<", NwtNew = "<<NwtNew<<"\nid = "<<cn->getID();
				cout <<"\tQIn = "<<QIn<<"  QOut ="<<QOut<<endl;
			}

			if (NfNew>NwtNew)  {
				cout << "\nWarning: Top Front > WT depth\n\n";
				cout <<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld;
				cout <<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout <<"NfNew = "<<NfNew<<", NwtNew = "<<NwtNew<<"\nid = "<<cn->getID();
				cout << "\tQIn = "<<QIn<<"  QOut = "<<QOut<<endl;
			}

			if (NtNew>NfNew) {
				cout << "\nWarning: Top Front > Wet front\n\n";
				cout <<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld;
				cout <<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout <<"NtNew = "<<NtNew<<", NfNew = "<<NfNew<<"\nid = "<<cn->getID();
				cout<<"\tQIn = "<<QIn<<"  QOut = "<<QOut<<endl;
			}

			if (NwtNew<0.0 || NwtOld<0.0) {
				cout << "\nWarning: Water table < 0\n\n";
				cout <<"NwtOld: "<<NwtOld<<", MiOld: "<<MiOld<<", MuOld: "<<MuOld;
				cout <<", NtOld= "<<NtOld<<", NfOld= "<<NfOld<<endl;
				cout << "NwtNew = " << NwtNew <<"\nid = "<<cn->getID();
				cout << "\tQIn = " << QIn << "  QOut = " << QOut << endl;
			}

			if (satsrf<0.0) {
				cout <<"\nWarning: RUNOFF component (satsrf) < 0: id = ";
				cout <<cn->getID()<<"\n\n";
			}

			if (satsrf>999999.0) {
				cout<<"\nWarning: RUNOFF comp.- satsrf  "
				<<"> 999999: id = "<<cn->getID()<<"\n\t\t-> Assigned to zero\n\n";
				satsrf = 0.0;
			}

			if (MuNew<MiNew) {
				cout <<"\nWarning: Total moisture content (MuNew) < MiNew: id = ";
				cout <<cn->getID()<<"\n\n";
			}
		}

		// GW contribution to surface runoff
		esrf = psrf+satsrf;
		esrf = esrf/Cos;
		satsrf /= Cos;
		srf = satsrf;

		// Print the variables of interest
		PrintNewGWVars( cn, Couple_State );

		// Update Variables
		ThSurf0 = cn->getSoilMoistureSC();
		ThRoot0 = cn->getRootMoistureSC();

		cn->setNwtNew(NwtNew);
		cn->setMuNew(MuNew);
		cn->setMiNew(MiNew);
		cn->setNfNew(NfNew);
		cn->setNtNew(NtNew);
		cn->setRuNew(RuNew);
		cn->setRiNew(RiNew);

		cn->addSrf_Hr(srf*dtGW);// WR debug 02062024, this is not in hours?
		cn->setsrf(cn->getSrf()+srf*dtGW); // this is a total depth in mm#
		cn->addCumSrf(cn->getSrf()); // Added line compared to old tRIBS version, not documented CJC 2022 //WRdebug 02062024: Above line was called, then called again in addCumSrf doubling the
		cn->setsatsrf(satsrf*dtGW);
		cn->setesrf(esrf*dtGW);

		if (satsrf>0.0) {
            cn->satsrfOccur = cn->satsrfOccur + floor(satsrf * 1.0E+3) + 1.0E-6;
        }

		// Saturation occurrence counter: same gated test as in
		// UnSaturatedZone(), see the comment there. Original test:
		//     if (ComputeSurfSoilMoist(1.0)/Ths > 0.999) satOccur++;
		if (cn->getNwtSatThresh() < 0.0)
			cn->setNwtSatThresh(ComputeNwtSatThreshold());
		if (!((NfNew == 0.0 || NfNew == NwtNew) && NwtNew >= cn->getNwtSatThresh())) {
			if (ComputeSurfSoilMoist(1.0)/Ths > 0.999)
				cn->satOccur = cn->satOccur + 1;
		}

		cn->RechDisch=cn->RechDisch+((NwtOld-NwtNew)+satsrf*dtGW*Cos/Ths)*1.0E-3;

		// Soil moisture in the top 10 cm
		ThSurf = ComputeSurfSoilMoist(surfaceSoilDepth);
		cn->setSoilMoisture( ThSurf );
		cn->setSoilMoistureSC( ThSurf/Ths );

		// Soil moisture in the unsaturated zone
		if (NwtNew) {
            cn->setSoilMoistureUNSC(MuNew / NwtNew / Ths);
        }
		else {
            cn->setSoilMoistureUNSC(1.0);
        }

		// Need to divide by the total # of time steps elapsed
		AA = (double)timer->getElapsedSteps(timer->getCurrentTime());
		// The integer part --surface SM--
		BB = floor(cn->getAvSoilMoisture())*1.0E-4;
		// The decimal part --root SM--
		Mdelt = (cn->getAvSoilMoisture() - floor(cn->getAvSoilMoisture()))*1.0E+1;
		cn->setAvSoilMoisture(0.0);
		cn->setAvSoilMoisture(floor((BB*AA + ThSurf/Ths)/(AA+1)*1.0E+4));

		// Estimate average root soil moisture
		{
			double rzDepth = std::min(cn->getRootZoneDepth(), DtoBedrock);
			ThSurf = ComputeSurfSoilMoist(rzDepth);
			cn->setRootMoisture( ThSurf );
			cn->setRootMoistureSC( ThSurf/Ths );
			// Recompute RootCutoff only when the rootzone depth changes
			if (rzDepth != cn->getRootCutoffDepth()) {
				cn->setRootCutoff( get_Upper_Moist(rzDepth, DtoBedrock) / rzDepth );
				cn->setRootCutoffDepth(rzDepth);
			}
		}
		cn->addAvSoilMoisture((Mdelt*AA + cn->getRootMoisture()/Ths)/(AA+1.0)*1.0E-1);

		cn->setRecharge((NwtNew-NwtOld)*Ths/(Cos*dtGW));
		// The following are in [M^3]
		Stok += srf*dtGW*cn->getVArea()/1000.0;
		TotMoist += (MuNew - MuOld)*cn->getVArea()/(Cos*1000.0);
        TotGWchange += (NwtNew-NwtOld)*Ths*cn->getVArea()/(Cos*1000.0);

		// Estimate and output relative factors
		if (fabs((cn->getSoilMoistureSC()) - ThSurf0) > 1.0E-6) {
			AA = 1.0;
			fGW100 += fabs(dtGW*(cn->getGwaterChng()*1.0E-6)/Area)*AA*Area;
			gwdm   += ((cn->getSoilMoistureSC()) - ThSurf0)*Ths*100.0*Area;
			AreaGW += Area;
			gwcnt++;
			/*
				cout<<"\t\tfGW100 = "<<dtGW*(cn->getGwaterChng()*1E-6)/Area
			 <<"\tdm100 = "<<((cn->getSoilMoistureSC()) - ThSurf0)*Ths*100.
			 <<"\tNwtOld = "<<NwtOld
			 <<endl<<flush; */
		}

		AreaF = Area/BasArea;
		dMRt += ((cn->getRootMoistureSC()) - ThRoot0)*Ths*1000.0*AreaF;

		// Mean soil moisture within the GW time interval over the domain
		mth100 += ((cn->getSoilMoistureSC()) + ThSurf0)/2.0*AreaF;
		mthrt  += ((cn->getRootMoistureSC()) + ThRoot0)/2.0*AreaF;
	
		id++;

		// if (cn->getNwtOld() > DtoBedrock) cnt++;
	}
	// cout<<"In total "<<cnt<<" cells have WT > BEDROCK\n"<<endl;

	if (gwcnt > 0) {
		fGW100 /= AreaGW;
		gwdm   /= AreaGW;
	}

	dM100   += gwdm;
	mTh100  /= (dtGW/(timer->getTimeStep()));
	mThRt   /= (dtGW/(timer->getTimeStep()));
	mTh100   = (mTh100+mth100)/2.0;
	mThRt    = (mThRt+mthrt)/2.0;
	fSoi100 /= (dtGW/(timer->getTimeStep()));
	fTop100 /= (dtGW/(timer->getTimeStep()));
	fClm100 /= (dtGW/(timer->getTimeStep()));

	// Could also be:
	//      /= (gridPtr->getNodeList()->getActiveSize());

	fsum = (fabs(fSoi100)+fabs(fTop100)+fabs(fClm100)+fabs(fGW100));
	if (!fsum) fsum = 1.0E-6;

	//fctout<<setprecision(10)<<timer->getCurrentTime()<<"  "
	//    <<dM100<<"  "<<dMRt<<"  "
	//	  <<mTh100<<"  "<<mThRt<<"  "
	//	  <<fsum<<"  "
	//	  <<fSoi100/fsum<<"  "<<fTop100/fsum<<"  "
	//	  <<fClm100/fsum<<"  "<<fGW100/fsum<<endl;

	fSoi100=fTop100=fClm100=fGW100=dM100=dMRt=mTh100=mThRt=0.0;

#ifndef PARALLEL_TRIBS
    if (simCtrl->Verbose_label == 'Y') { //WR 08292023 change conditions from nodelist to Verbose
		cout<<"\t\tRUNOFF = "<<Stok<<" M^3"<<endl<<flush;
		cout<<"\t\tMOISTCHN = "<<TotMoist<<" M^3"<<endl<<flush;
		cout<<"\t\tGWCHANGE = "<<TotGWchange<<" M^3"<<endl<<flush;
	}
#endif

	} // End of SaturatedZone Routine...


//=========================================================================
//
//
//                  Section 8: tHydroModel Surface Soil Moisture
//
//
//=========================================================================

/*************************************************************************
**
**  tHydroModel::ComputeSurfSoilMoist(double Z)
**
**  The function computes average soil moisture value for the top
**  soil layer of thickness Z (mm). Can be used only at the beginning or
**  at the end of 'UnsaturatedZone' or 'SaturatedZone' functions when
**  the state variables are set.
**
*************************************************************************/
double tHydroModel::ComputeSurfSoilMoist(double Z)
{
	double dM = 0;
    double Nstar = 0;

	if (NwtNew == 0.0)
		dM = Z*Ths;

	else if ((NwtNew > 0.0) && (NfNew==0.0 || NfNew == NwtNew)) //Init. prof.
		dM = get_Upper_Moist(Z, NwtNew);

	else if (NfNew > 0.0 && NtNew==NfNew && NfNew < NwtNew) { //Unsat edge
		if (fabs(RuNew) > 1.0E-6) {
			Nstar = log(Ksat/RuNew)/F;
			if (NfNew >= Z)
				dM = Eps/F*(Ths-Thr)*(exp(F*(Z-Nstar)/Eps)-exp(-F*Nstar/Eps))+Thr*Z;
			else {
				dM  = get_Upper_Moist(Z, NwtNew);
				dM += (MuNew - MiNew);
			}
		}
		else
			dM = get_Upper_Moist(Z, NwtNew);
	}
	else if (NfNew > 0.0 && NtNew<NfNew && NtNew > 0.0) { //Perched_sat edge
		Nstar = NtNew;
		if (NtNew >= Z)
			dM = Eps/F*(Ths-Thr)*(exp(F*(Z-Nstar)/Eps)-exp(-F*Nstar/Eps))+Thr*Z;
		else {
			if (NfNew >= Z) { //Z is between the top front and the wet front
				dM  = Eps/F*(Ths-Thr)*(1.0-exp(-F*NtNew/Eps))+Thr*NtNew;
				dM += ((Z-NtNew)*Ths);
			}
			else if (NfNew < Z) {
				dM  = get_Upper_Moist(Z, NwtNew);
				dM += (MuNew - MiNew);
			}
		}
	}
	else if (NfNew > 0.0 && NtNew == 0.0) { //Surf_sat edge
		if (NfNew >= Z)
			dM = Z*Ths;
		else if (NfNew < Z) {
			dM  = get_Upper_Moist(Z, NwtNew);
			dM += (MuNew - MiNew);
		}
	}
	return (dM/Z); // returns average soil moisture over depth Z
}

/*************************************************************************
**
**  tHydroModel::ComputeNwtSatThreshold()
**
**  Returns the water table depth (mm) above which the top 1 mm of soil
**  can possibly be saturated, used to skip the expensive satOccur test
**  ComputeSurfSoilMoist(1.0)/Ths > 0.999 when it must be false.
**
**  Derivation, valid for the initial hydrostatic profile (no wetting
**  front), where moisture only increases with depth toward the water
**  table, so the top 1 mm average cannot exceed the moisture at 1 mm:
**
**     theta(1mm) = Thr + (Ths-Thr)*(-Psib/(Nwt-1))^PoreInd  >  0.999*Ths
**
**  rearranges to Nwt < NwtSatThresh with
**
**     NwtSatThresh = 1 + (-Psib)*((Ths-Thr)/(0.999*Ths-Thr))^(1/PoreInd)
**
**  Depends only on static soil properties: computed once per node and
**  cached in tCNode.
**
*************************************************************************/
double tHydroModel::ComputeNwtSatThreshold() const
{
	double denom = 0.999*Ths - Thr;
	if (denom <= 0.0)
		return 1.0E+9; // degenerate soil parameters: never skip the full test
	return 1.0 + (-Psib)*pow((Ths-Thr)/denom, 1.0/PoreInd);
}

/*************************************************************************
**
** tHydroModel::get_Z1Z2_Moist(double Z1, double Z2, double Nwt)
**
** Allows to get a moisture content in the initial profile between depths
** Z1 and Z2 with water table depth equal to Nwt
**
*************************************************************************/
//double tHydroModel::get_Z1Z2_Moist(double Z1, double Z2, double Nwt)
// SKY2008Snow from AJR2007
double tHydroModel::get_Z1Z2_Moist(double Z1, double Z2, double Nwt, int mark)
{
	double dM;

	// SKY2008Snow from AJR2007
	mark++;

	if (Z1 >= Z2) {
		if (simCtrl->Verbose_label == 'Y') {
			cout<<"\nWARNING: Wrong use of the function get_Z1Z2_Moist()!"
			<<" \nArgument Z2 must be larger than Z1!"<<endl<<flush;
		}
	}

	if (Z2 <= (Nwt+Psib+1.0E-3))
		if (PoreInd >(1.0-1.0E-6) && PoreInd <(1.0+1.0E-6))
			dM = Thr*(Z2-Z1)-fabs(Psib)*(Ths-Thr)*log((Nwt-Z2)/(Nwt-Z1));
		else
			dM = Thr*(Z2-Z1)-(Ths-Thr)/(PoreInd-1.0)*
				((Z2-Nwt)*pow((-Psib/(Nwt-Z2)),PoreInd) -
				 (Z1-Nwt)*pow((-Psib/(Nwt-Z1)),PoreInd));

	else if (Z2 > (Nwt+Psib+1.0E-3) && Z1 < (Nwt+Psib+1.0E-3)) {
		//dM  = get_Z1Z2_Moist(Z1, Nwt+Psib, Nwt);
		// SKY2008Snow from AJR2007
		dM  = get_Z1Z2_Moist(Z1, Nwt+Psib, Nwt, mark);
		dM += (Z2-(Nwt+Psib))*Ths;
	}

	else
		dM = (Z2-Z1)*Ths;
	return dM;
}

/***************************************************************************
**
**  tHydroModel::GetCellRunon() and SetCellRunon()
**
**  Simple runon scheme that allows to account for the surface runon
**  onto the element. 'Cos' of the element is assumed to be accounted for
**  in both 'runoff' and 'runon'.
**
***************************************************************************/
double tHydroModel::GetCellRunon( tCNode *cn, double DTUS )
{
	double runon;
	tCNode *cnn;
	tEdge  *firstedg;
	tEdge  *curedg;

	runon = 0.0;

	// Originally Cos = cos(atan(slope)) = 1/sqrt(1+s^2), avoids trig calls
	double edgeSlope = cn->getFlowEdg()->getSlope();
	Cos = (edgeSlope > 0.0) ? 1.0/sqrt(1.0 + edgeSlope*edgeSlope) : 1.0;

	firstedg = cn->getFlowEdg();
	curedg = firstedg->getCCWEdg();
	while (curedg != firstedg) {
		cnn = (tCNode*)curedg->getDestinationPtrNC();
		if (cnn->getBoundaryFlag() != kClosedBoundary &&
			cnn->getBoundaryFlag() != kOpenBoundary) {

			// Find upslope element contributing to the current element
			if (cnn->getFlowEdg()->getDestinationPtrNC() == (tNode*)cn &&
				cnn->getSrf() > 0.0) {

				// Account for area and slope differences [mm hour^-1]
				// 'Cos' of 'cnn' has been accounted for in setSrf()
				runon += ((cnn->getSrf())/DTUS)*
				(cnn->getVArea()) / ((cn->getVArea())/Cos);

				// Account for area and slope differences [mm hour^-1]
				// 'Cos' of 'cnn' has been accounted for in setSrf()
				runon += ((cnn->getSrf())/DTUS)*(cnn->getVArea()) / ((cn->getVArea())/Cos);
			}
		}
		curedg = curedg->getCCWEdg();
	}
	if (runon < 1.0E-6)
		runon = 0.0;
	return runon;
}

void tHydroModel::SetCellRunon( tCNode *cn, double runoff,
								double runon, double DTUS, int eOpt )
{
	double del_srf, Cos1, rf;
	tCNode *cnn;
	tEdge  *firstedg;
	tEdge  *curedg;

	// Set runon to the cell [mm]
	cn->setRunOn(runon*DTUS);

	if (!runon) {
		cn->setRunOn(0.0);
		return;
	}

	// If runoff is higher than runon - element contributes to runoff
	if (runoff >= runon) {
		if (eOpt)
			cn->setsrf((runoff - runon)*DTUS);
		else
			rf = (runoff - runon);
		return;
	}

	// Since runon is higher, set runoff to zero
	if (eOpt)
		cn->setsrf(0.0);
	else
		rf = 0.0;

	// Difference in runon and runoff [mm hr^-1]
	del_srf = runon - runoff;

	firstedg = cn->getFlowEdg();
	curedg = firstedg->getCCWEdg();
	while (curedg != firstedg) {
		cnn = (tCNode*)curedg->getDestinationPtrNC();
		if (cnn->getBoundaryFlag() != kClosedBoundary &&
			cnn->getBoundaryFlag() != kOpenBoundary) {

			// Find upslope element contributing to the current element
			if (cnn->getFlowEdg()->getDestinationPtrNC() == (tNode*)cn &&
				cnn->getSrf() > 0.0) {
				// Originally Cos1 = cos(atan(slope)) = 1/sqrt(1+s^2)
				double slpNbr = cnn->getFlowEdg()->getSlope();
				Cos1 = (slpNbr > 0.0) ? 1.0/sqrt(1.0 + slpNbr*slpNbr) : 1.0;

				// Account for area and slope differences [mm hour^-1]
				runoff = ((cnn->getSrf())/DTUS);

				if (!eOpt) {
					Stok  -= cnn->getSrf()*
					cnn->getVArea()/(1000.);

					Stok  += (runoff*(1.0 - del_srf/(runon))*DTUS)*
						cnn->getVArea()/(1000.);
				}
				// If upslope cell produces runoff - reduce it
				// proportionally to its contribution to 'runon'
				if (eOpt)
					cnn->setsrf(runoff*(1.0 - del_srf/(runon))*DTUS);
			}
		}
		curedg = curedg->getCCWEdg();
	}
	}

//=========================================================================
//
//
//                  Section 9: tHydroModel Lambert Function
//
//
//=========================================================================

/***************************************************************************
**
**  tHydroModel::LambertW( double )
**
**  Compute the Lambert W function of z.  This function satisfies
**  W(z).*exp(W(z)) = z, and can thus be used to express solutions of
**  transcendental equations involving exponentials or logarithms. n must
**  be integer, and specifies the branch of W to be computed; W(z) is a
**  shorthand for W(0,z), the principal branch. Branches 0 and -1 are the
**  only ones that can take on non-complex values. If either n or z are
**  non-scalar, the function is mapped to each element; both may be
**  non-scalar provided their dimensions agree. This implementation should
**  return values within 2.5*eps of its counterpart in Maple V, release 3
**  or later.  Please report any discrepancies to the author, Nici Schraudolph
**  <nic@idsia.ch>. For further details, see: Corless, Gonnet, Hare, Jeffrey,
**  and Knuth (1996), "On the Lambert W Function", Advances in Computational
**  Mathematics 5(4):329-359.
**
***************************************************************************/
double tHydroModel::LambertW(double z)
{
	fcomplex w, v, tt, ttt, t, p;
	double xx, xy;
	double c, f;
	int n;

	if (z == 0.0)     //  LambertW function in 0 is 0
		return 0.0;

	xx = 2.0*exp(1.0)*z + 2.0;

	if (xx < 0.0) {
		xx = fabs(xx);
		xx = sqrt(xx);
		setComplex(&w, -1.0, xx);
	}
	else {
		xx = sqrt(xx) - 1.0;
		setComplex(&w, xx, 0);
	}

	xx = fabs(z);
	xx = log(xx);
	if (z < 0.0)
		setComplex(&v, xx, PI);
	else
		setComplex(&v, xx, 0);

	if (z > 1.0) {
		xy = log(xx);
		setComplex(&tt, xy, 0);
		setCompNumber(&v, (Csub(v,tt)));
	}
	else if (z == 1.0)
		setComplex(&tt, 0, 0);
	else if ( z < 1.0 && z > 0.0) {
		xy = log(fabs(xx));
		setComplex(&tt, xy, PI);
		setCompNumber(&v, (Csub(v,tt)));
	}
	else if ( z < 0.0 && z > -1.0) {
		xx = atan(v.i/ v.r);         	//negative value
		xy = log(fabs(v.r / cos(xx)));
		xx += PI;                     	//positive imaginary part
		setComplex(&tt, xy, xx);
		setCompNumber(&v, (Csub(v,tt)));
	}

	// Choose strategy for initial guess
	c = fabs(z + 1.0/exp(1.0));
	if (c > 1.45) {
		c = 1.0;
		setCompNumber(&w, v);
	}
	else
		c = 0.0;

	// Halley iteration
	n = 0;
	setComplex(&t, 1000, 1000);

	while (n < 15 && (fabs(t.r) > xy || fabs(t.i) > xx)) {
		xx = exp(Cr(&w));
		setComplex(&p, (xx*cos(Ci(&w))), (xx*sin(Ci(&w))));

		setCompNumber(&t, (Cmul(w, p)));
		subtReal(&t, z);

		if (w.i == 0 && w.r == -1)
			f = 0.0;
		else
			f = 1.0;

		setCompNumber(&tt, Cmul(p, CaddFF(&w, f)));
		setCompNumber(&ttt,Csub(tt,Cdiv(Cmul(t,RCmul(0.5,CaddFF(&w, 2.0))),CaddFF(&w, f))));
		setCompNumber(&t, Cdiv(RCmul(f,t), ttt));
		setCompNumber(&w, Csub(w, t));

		xy = (2.48*LAMBEPS)*(1.0+fabs(w.r));
		xx = (2.48*LAMBEPS)*(1.0+fabs(w.i));
		n++;
	}
	return w.r;
}

//=========================================================================
//
//
//                  Section 9: tHydroModel Newton and Polyn Function
//
//
//=========================================================================


/***************************************************************************
**
**  tHydroModel::Newton (double,double)
**
**  Returns a value of water table depth based on the solution of water
**  balance equation.  Arguments: dM - moisture deficit in the unsaturated
**  zone calculated from previous state, forcing, and lateral fluxes.
**  Nwt is old value of water table depth: in the case Nwt > 0, it is used
**  as an initial guess.
**
***************************************************************************/
double tHydroModel::Newton(double dM, double Nwt)
{

	int ITMAX = 30;
	double DMIN   = 1.0E-5; //Minimum slope of the derivative function
	double DX_TOL = 1.0E-5; //Tolerance for dx
	int i;
	double C1, C2, C3;
	double fvalue, fdvalue, x(0), dx;
	double xinit, xup, Nwt_estim;

	C1 = Ths-Thr;
	C2 = C1*pow((-Psib),PoreInd)/(PoreInd-1.0);
	C3 = C1*Psib/(PoreInd-1.0) + Psib*C1 - dM;
	
	if (Nwt == 0.0) {
		// This is an extremum point (minimum), needs to be shifted
		xinit = dM*(PoreInd-1.0)/((Ths-Thr)*PoreInd) - Psib;
		xinit += 50.0;
		if (xinit <= fabs(Psib))    // If less than allowable
			xinit = fabs(Psib) + 50.0;
	}
	// A better initial estimate will be using NwtOld
	else
		xinit = Nwt;

	// Depending on what pore-size distribution index
	// use either a more complicated or a simpler solver
	if (PoreInd < 5.0) {
		if ((Ths*Nwt - get_Total_Moist(Nwt)) > dM)  // water table rises
			xup = ceil(Nwt);
		else                                        // water table drops
			xup = DtoBedrock + 1.0;
		x = rtsafe_mod(C1, C2, C3, dM/(Ths-Thr)+fabs(Psib), xup, DX_TOL, xinit);
	}

	// If the flux too high, then the value of water table can numerically
	// go below 'DtoBedrock', just let it go down to a such depth
	if (PoreInd >= 5.0 || fabs(x-xinit) <= 1.0E-7  ) {
		for (x = xinit, i=0; i<ITMAX; x-=dx, i++) {
			polyn(x, fvalue, fdvalue, C1, C2, C3);

			if (fabs(fdvalue) < DMIN) {
				if (simCtrl->Verbose_label == 'Y') {
					cout <<"\nWarning! NEWTON 1: derivative fell below minimum!"
					<<" ID = "<<ID<<endl<<endl;
					// return(x); // commented out by CJC2020
				}
				Nwt_estim = x; // added by CJC2020
				break; // added by CJC2020
			}

			dx = fvalue/fdvalue;

			if (fabs(dx) < DX_TOL) {
				Nwt_estim = x; // added by CJC2020
				break; // added by CJC2020
				// return(x); // commented out by CJC2020
			}
		}

		if (simCtrl->Verbose_label == 'Y') {
			cout<<"\nWarning: NEWTON 1: NO convergence in "<<ITMAX<<" steps\n";
			cout<<"\nERROR = "<<fvalue
				<<";  Nwt_initial = "<<Nwt
				<<";  dM = "<<dM
				<<";  Nwt_estim = "<<x
				<<";  ID = "<<ID
				<<"\nInitial value is kept..."
				<<endl<<endl<<flush;
		}
		
		if ((fabs(fdvalue) < DMIN) || (fabs(dx) < DX_TOL))  { // added by CJC2020
			x = Nwt_estim;
		}

		else { // added by CJC2020
			x = Nwt;
		}
	}
	if (x > DtoBedrock)
		x = DtoBedrock;
	return (x);
}

/***************************************************************************
**
**  tHydroModel::Newton (double, double, double, int)
**
**  Returns value of depth in the soil moisture profile (above Nf)
**  which allows extraction from the wetted profile 'dM' amount of
**  soil moisture (subtraction/addition from above/below the wetting front).
**  The solution is based on water balance in the profile. The function
**  is called when water table drops/rises and there is a need to slightly
**  re-adjust position of the wetting front to satisfy moisture balance
**  Arguments: dM - amount of moisture required to "subtract"
**  Nwt NEW value of the water table depth, Nf - old position of the
**  wetting front, used as initial guess, and flag: 0/1 - WT drop/rise case
**
***************************************************************************/
double tHydroModel::Newton(double dM, double Nwt, double Nf, int flag) const
{
	int ITMAX = 30;
	double DMIN   = 1.0E-5; // MINIMUM SLOPE OF THE DERIVATIVE F_N
	double DX_TOL = 1.0E-5; // TOLERANCE FOR dx
	int i;
	double C1, C2, C3;
	double fvalue, fdvalue, x, dx;
	double xinit;

	// dM, Nwt and Nf have to be passed to the function
	if (!flag)
		C1 = Ths-Thr;
	else
		C1 = -(Ths-Thr);
	C2 = C1*pow((-Psib),PoreInd)/(PoreInd-1);
	C3 = C2/pow((Nwt-Nf),(PoreInd-1)) - C1*Nf + dM;

	xinit = Nf; //Initial guess for the function

	for (x = xinit, i=0; i<ITMAX; x-=dx, i++) {
		polyn(x, Nwt, fvalue, fdvalue, C1, C2, C3);
		if (fabs(fdvalue) < DMIN) {
			if (simCtrl->Verbose_label == 'Y') {
				cout<<"\nWarning NEWTON: Derivative fell below minimum." <<endl;
				cout<<"ID = "<<ID<<"; fdvalue = "<<fdvalue<<endl;
				return(x);
			}
		}
		dx = fvalue/fdvalue;
		if (fabs(dx) < DX_TOL) {
			return(x);
		}
	}

	if (simCtrl->Verbose_label == 'Y') {
		cout<<"\nWarning: NEWTON 2: NO convergence in "<<ITMAX<<" steps\n";
		cout<<"ERROR = "<<fvalue <<";  flag  = "<<flag <<";  Nwt_initial = "<<Nwt
			<<";  Nf_initial = "<<Nf <<";  dM = "<<dM <<";  ID = "<<ID
			<<endl<<flush;
	}

	x = Nf;

	return (x);
}

/*****************************************************************************
**
**  tHydroModel::rtsafe_mod
**
**  Finds a root of the polynomial which lies in the range [x1 and x2]
**  starting from initial guess - xguess. Accuracy of estimation is xacc.
**  c1, c2, c3 are the polynomial coefficients
**
*****************************************************************************/
#define MAXITER 30
double tHydroModel::rtsafe_mod(double c1, double c2, double c3,
                               double x1, double x2, double xacc, double xguess)
{
	int j;
	double df(0.0), dx, dxold, f, fh(0.0), fl(0.0); // Giuseppe -- GMnSKY2008MLE
	double temp, xh, xl, rts;

	// 'fl' & 'fh' below are the evaluation function values
	// corresponding to arguments 'x1' and 'x2'

	polyn(x1,fl,df,c1,c2,c3);
	polyn(x2,fh,df,c1,c2,c3);

	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		if (simCtrl->Verbose_label == 'Y') {
			cout<<"`\nWarning: Root must be bracketed by negative and positive f_n values!"<<endl;
			cout<<"fl = "<<fl<<"; fh = "<<fh<<"; xguess = "<<xguess<<"; ID = "<<ID<<endl;
		}
		return xguess;
	}
	if (fl == 0.0) return x1;

	if (fh == 0.0) return x2;

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

	polyn(rts,f,df,c1,c2,c3);

	for (j=1; j <= MAXITER; j++) {  // Loop over allowed iterations
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) // Bisect if Newton is out of range
			|| (fabs(2.0*f) > fabs(dxold*df))) {    // or not decreasing fast enough
			dxold = dx;
			dx = 0.5*(xh-xl);
			rts = xl+dx;
			if (xl == rts) return rts; //Change in root is
		}                            // negligible, take it
		else {
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts; // Convergence criterion

		polyn(rts,f,df,c1,c2,c3);
		if (f < 0.0) // Maintain the bracket on the root
			xl=rts;
		else
			xh=rts;
	}

	if (simCtrl->Verbose_label == 'Y') {
		cout<<"\nWarning in NEWTON 3: NO convergence in "<<MAXITER<<" steps\n";
		cout<<"\nERROR = "<<f
			<<";  Nwt_initial = "<<xguess
			<<";  Nwt_estim = "<<rts
			<<";  ID = "<<ID
			<<"\nInitial value is kept..."
			<<endl<<endl<<flush;
	}
	return xguess;
}
#undef MAXITER

/***************************************************************************
**
**  tHydroModel::polyn
**
**  Both return value of the function in non-closed form (fv) as well as
**  the value of its derivative (dv). Used in Newton iterative procedure
**  to find root of an equation.
**
**  Evaluates the polynomial and its derivative:
**     fv = C1*x^PoreInd + C3*x^(PoreInd-1) + C2
**     dv = C1*PoreInd*x^(PoreInd-1) + C3*(PoreInd-1)*x^(PoreInd-2)
**
**  PERFORMANCE NOTE: pow() is expensive and this function runs inside the
**  Newton/rtsafe iteration loops (up to ~30 calls per node per time step).
**  The three powers above differ only by whole numbers, so instead of five
**  pow() calls we compute x^(PoreInd-2) once and multiply up by x:
**     x^(PoreInd-1) = x^(PoreInd-2) * x
**     x^PoreInd     = x^(PoreInd-1) * x
**  The math is essentially identical to the equations above.
**
***************************************************************************/
void tHydroModel::polyn(double x, double& fv, double& dv,
                        double C1, double C2, double C3) const
{
	double xPm2 = pow(x,(PoreInd-2.0)); // x^(PoreInd-2)
	double xPm1 = xPm2*x;               // x^(PoreInd-1)
	double xP   = xPm1*x;               // x^PoreInd
	fv = C1*xP + C3*xPm1 + C2;
	dv = C1*PoreInd*xPm1 + C3*(PoreInd-1.0)*xPm2;
}

/***************************************************************************
**
**  Second form, written in terms of the depth difference (Nwt - x):
**     fv = C1*x*(Nwt-x)^(PoreInd-1) + C3*(Nwt-x)^(PoreInd-1) - C2
**     dv = C1*(Nwt-x)^(PoreInd-1) - C1*(PoreInd-1)*x*(Nwt-x)^(PoreInd-2)
**          - C3*(PoreInd-1)*(Nwt-x)^(PoreInd-2)
**
**  Same change as above: (Nwt-x)^(PoreInd-2) is computed with a
**  single pow() call and (Nwt-x)^(PoreInd-1) is obtained by multiplying
**  by (Nwt-x), replacing five pow() calls with one.
**
***************************************************************************/
void tHydroModel::polyn(double x, double Nwt, double& fv, double& dv,
						double C1, double C2, double C3) const
{
	double y    = Nwt - x;
	double yPm2 = pow(y,(PoreInd-2.0)); // (Nwt-x)^(PoreInd-2)
	double yPm1 = yPm2*y;               // (Nwt-x)^(PoreInd-1)
	fv = C1*x*yPm1 + C3*yPm1 - C2;
	dv = C1*yPm1 - C1*(PoreInd-1.0)*x*yPm2 - C3*(PoreInd-1.0)*yPm2;
}


//=========================================================================
//
//
//                       Section 10: Print Variables
//
//
//=========================================================================


/*************************************************************************
**
**  tHydroModel::PrintOldVars(tCNode *cn)
**
**  Prints OLD dynamic variables for the current node if it is in the
**  list of nodes for which the output has to be made
**
*************************************************************************/
void tHydroModel::PrintOldVars(tCNode *cn, tEdge *ce, double Ractual,
                               int Pixel_State)
{
	if (simCtrl->Verbose_label == 'Y') {
		if ( nodeList ) {
			enum {Storm_Evol, WTStaysAtSurf, WTDropsFromSurf, WTGetsToSurf,
				Storm_Unsat_Evol, Perched_Evol, Perched_SurfSat,
				StormToInterTransition, ExactInitial, IntStormBelow};

			for (int ii=0; ii < numNodes; ii++) {
				if (cn->getID() == nodeList[ii]) {
					cout <<"------------------------->> BEFORE:\nNode #"
					<<cn->getID()<<" ("<<cn->getX()<<","
					<<cn->getY()<< ")"<<endl<<flush;
					if (Pixel_State == Storm_Evol)
						cout << "CASE: Storm_Evol \n";
					else if (Pixel_State == WTStaysAtSurf)
						cout << "CASE: WTStaysAtSurf \n";
					else if (Pixel_State == WTDropsFromSurf)
						cout << "CASE: WTDropsFromSurf \n";
					else if (Pixel_State == WTGetsToSurf)
						cout << "CASE: WTGetsToSurf \n";
					else if (Pixel_State == Storm_Unsat_Evol)
						cout << "CASE: Storm_Unsat_Evol\n";
					else if (Pixel_State == Perched_Evol)
						cout << "CASE: Perched_Evol \n";
					else if (Pixel_State == Perched_SurfSat)
						cout << "CASE: Perched_SurfSat \n";
					else if (Pixel_State == StormToInterTransition)
						cout << "CASE: StormToInterTransition \n";
					else if (Pixel_State == ExactInitial)
						cout << "CASE: ExactInitial \n";
					else if (Pixel_State == IntStormBelow)
						cout << "CASE: *** IntStormBelow \n";
					cout << "-------------------------" <<'\n';
					cout << "tBOUND.: " <<cn->getBoundaryFlag()<<'\n';
					cout << "Area  = "  <<cn->getVArea()<<" m^2"<<'\n';
					cout << "DESTIN ID = "<<ce->getDestinationPtrNC()->getID()<<"; Z = "
						<< ce->getDestinationPtrNC()->getZ()<<endl<<flush;
					cout << "Edge lgth: " <<ce->getLength()<<" m"<<endl<<flush;
					cout << "Flw width = "<<ce->getVEdgLen()<<" m"<<'\n';
					cout << "Cos   = "  <<Cos<<'\n';
					cout << "Sin   = "  <<Sin<<'\n';
					cout << "Rain = "  <<cn->getRain()
						<< "; NetRain = "<<Ractual<<"; R1 = "<<R1<<endl;
					cout << "NwtOld = " <<cn->getNwtOld()<< '\n';
					cout << "MiOld  = " <<cn->getMiOld()<< '\n';
					cout << "MuOld  = " <<cn->getMuOld()<< '\n';
					cout << "RuOld  = " <<cn->getRuOld()<< '\n';
					cout << "NfOld  = " <<cn->getNfOld()<< '\n';
					cout << "NtOld  = " <<cn->getNtOld()<< '\n';
					cout << "Qpin   = " <<cn->getQpin()*1.0E-6/(cn->getVArea())<<" mm/h"<< '\n';
					cout << "Qpout  = " <<cn->getQpout()*1.0E-6/(cn->getVArea())<<" mm/h"<< '\n';
					cout << "SurfTemp = "<<cn->getSurfTemp()<<" oC"<<endl;
					cout << "SoilTemp = "<<cn->getSoilTemp()<<" oC"<<endl;
					cout << "AirTemp = "<<cn->getAirTemp()<<" oC"<<endl;
					cout << "IntStormVar = "<<cn->getIntStormVar()<<" hrs"<<endl;
					cout<<endl<<flush;
				}
			}
		}
	}
	}

/*************************************************************************
**
**  tHydroModel::PrintNewVars(tCNode *cn)
**
**  Prints NEW dynamic variables for the current node if it is in the
**  list of nodes for which the output has to be made
**
*************************************************************************/
void tHydroModel::PrintNewVars(tCNode *cn, double Ractual)
{
	if (simCtrl->Verbose_label == 'Y') {
		if ( nodeList ) {
			for (int ii=0; ii < numNodes; ii++) {
				if (cn->getID() == nodeList[ii]) {
					cout << "------------------------->> AFTER:\nNode #" << cn->getID()
					<<" ("<<cn->getX()<<","<< cn->getY() << ")"<<endl<<flush;
					cout << "-------------------------" <<'\n';
					cout <<"Ksat = "<<Ksat
						<<"; Ths = "<<Ths
						<<"; Thr = "<<Thr
						<<";\nPoreInd = "<<PoreInd
						<<"; Eps = "<<Eps
						<<"; Psib = "<<Psib
						<<";\nF = "<<F
						<<"; Ar = "<<Ar
						<<"; UAr = "<<UAr<<";\n";
					cout << "Rain = "<<Ractual<<"; R1 = "<<R1<<endl;
					cout << "NwtNew = " <<cn->getNwtNew()<< '\n';
					cout << "MiNew  = " <<cn->getMiNew()<< '\n';
					cout << "MuNew  = " <<cn->getMuNew()<< '\n';
					cout << "RuNew  = " <<cn->getRuNew()<< '\n';
					cout << "NfNew  = " <<cn->getNfNew()<< '\n';
					cout << "NtNew  = " <<cn->getNtNew()<< '\n';
					cout << "Qpin   = " <<cn->getQpin()*1.0E-6/(cn->getVArea())<<" mm/h"<<'\n';
					cout << "Qpout  = " <<cn->getQpout()*1.0E-6/(cn->getVArea())<<" mm/h"<<'\n';
					cout << "Srf    = " <<cn->getSrf()<<" mm"<<endl<<flush;
					cout << "Srf_Hr = " <<cn->getSrf_Hr()<<" mm"<<endl<<flush;
					cout << "SurfTemp = "<<cn->getSurfTemp()<<" oC;\n";
					cout << "IntStormVar = "<<cn->getIntStormVar()<<" hrs"<<endl;
					cout<<endl<<flush;
				}
			}
		}
	}
	}

/*************************************************************************
**
**  tHydroModel::PrintNewGWVars(tCNode *cn, int Pixel_State)
**
**  Prints new dynamic variables for the current node if it is in the
**  list of nodes for which the output has to be made: after GW loop
**
*************************************************************************/
void tHydroModel::PrintNewGWVars(tCNode *cn, int Pixel_State)
{
	enum {GW_Exfiltrate, GW_IntStorm_Like, GW_Initial, GW_Positive_Bal};

	if (simCtrl->Verbose_label == 'Y') {
		if ( nodeList ) {
			for (int ii=0; ii < numNodes; ii++) {
				if (cn->getID() == nodeList[ii]) {
					cout<<"------------------------->> AFTER GW MODEL:\nNode #"
					<<cn->getID()<<endl;
					if (Pixel_State == GW_Exfiltrate)
						cout << "CASE: GW_Exfiltrate  \n";
					else if (Pixel_State == GW_IntStorm_Like)
						cout << "CASE: GW_IntStorm_Like \n";
					else if (Pixel_State == GW_Initial)
						cout << "CASE: GW_Initial \n";
					else if (Pixel_State == GW_Positive_Bal)
						cout << "CASE: GW_Positive_Bal \n";

					cout<<"\tNwtOld = "<<cn->getNwtOld()
						<<"; MiOld = "<<cn->getMiOld()
						<<"; MuOld = "<<cn->getMuOld()
						<<"; RuOld = "<<cn->getRuOld()
						<<"; NtOld = "<<cn->getNtOld()
						<<"; NfOld = "<<cn->getNfOld()<<";"<<endl
						<<"Transmissivity = "<<cn->getTransmiss()*1.0E-6<<" m^2/hr\n"
						<<"FLUX = "<<cn->getGwaterChng()*1.0E-9<<" m^3/hr\n"
						<<" NwtNew = "<<NwtNew<<";\n"
						<<" MiNew = "<<MiNew<<";\n"
						<<" MuNew = "<<MuNew<<";\n"
						<<" RuNew = "<<RuNew<<";\n"
						<<" NtNew = "<<NtNew<<";\n"
						<<" NfNew = "<<NfNew<<";\n"
						<<" Srf   = "<<srf<<" mm;\n"
						<<" Srf_Hour = "<<cn->getSrf_Hr()<<" mm;\n"
						<<endl<<flush;
				}
			}
		}
	}
	}

//=========================================================================
//
//
//                       End of tHydroModel.cpp
//
//
//=========================================================================


