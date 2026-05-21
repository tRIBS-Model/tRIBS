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
**  tOutput.cpp: Functions for output objects for classes tOutput and 
**               tCOutput (see tOutput.h)
**
*************************************************************************/

#include "src/tInOut/tOutput.h"
#include "src/Headers/globalIO.h"
#include "src/Headers/Inclusions.h"

#ifdef PARALLEL_TRIBS
#include "src/tGraph/tGraph.h"
#include "src/tParallel/tParallel.h"
#endif

//=========================================================================
//
//
//                  Section 1: tOutput Constructors/Destructors
//
//
//=========================================================================

/*************************************************************************
**
**  tOutput::Constructor (1)
**
**  The constructor takes two arguments, a pointer to the grid mesh and
**  a reference to an open input file. It reads the base name for the
**  output files from the input file, and opens and initializes these.
**
*************************************************************************/
template< class tSubNode >
tOutput<tSubNode>::tOutput(SimulationControl *simCtrPtr,
						   tMesh<tSubNode> * gridPtr, tInputFile &infile,
						   tResample *resamp)
{
	assert(gridPtr > 0);
	g = gridPtr;  
	respPtr = resamp;
	simCtrl = simCtrPtr;

	
	Cout<<"\nOutput Files:"<<endl<<endl;
	infile.ReadItem(baseName, "OUTFILENAME" );    //basename of output
	infile.ReadItem(nodeFile, "NODEOUTPUTLIST");  //pathname of node output
	
#ifdef PARALLEL_TRIBS
  // Nodes, edges, triangles, and Z files are only written
  // by the Master node
  if (tParallel::isMaster()) {
#endif

	char nodesext[10] = ".nodes";
	char edgesext[10] = ".edges";
	char trisext[10] = ".tri";
	char zofsext[10] = ".z";
	
	CreateAndOpenFileSingle( &nodeofs, nodesext );
	CreateAndOpenFileSingle( &edgofs,  edgesext );
	CreateAndOpenFileSingle( &triofs,  trisext );
	CreateAndOpenFileSingle( &zofs,    zofsext );

	nodeofs.setf( ios::fixed, ios::floatfield);
	zofs.setf   ( ios::fixed, ios::floatfield);

#ifdef PARALLEL_TRIBS
  }
#endif

	ReadNodeOutputList();
	CreateAndOpenPixel();
	dynvars = nullptr;
}

/*************************************************************************
**
**  tOutput::Constructor (2)
**
**  The constructor takes two arguments, a pointer to the grid mesh and
**  a reference to an open input file. It reads the base name for the
**  output files from the input file, and opens and initializes these.
**
*************************************************************************/
template< class tSubNode >
tOutput<tSubNode>::tOutput( SimulationControl *simCtrPtr,
							tMesh<tSubNode> * gridPtr, 
                            tInputFile &infile, tResample *resamp,
                            tRunTimer * timptr )
{
	assert(gridPtr != nullptr); //Updated to new c++ standards -WR
	g = gridPtr;        
	timer = timptr;    
	respPtr = resamp;
	simCtrl = simCtrPtr;

	
	Cout<<"\nOutput Files:"<<endl<<endl;
	infile.ReadItem( baseName, "OUTFILENAME" );          
	infile.ReadItem( nodeFile, "NODEOUTPUTLIST");  
	
#ifdef PARALLEL_TRIBS
  // Nodes, edges, triangles, and Z files are only written
  // by the Master node
  if (tParallel::isMaster()) {
#endif

	char nodesext[10] = ".nodes";
	char edgesext[10] = ".edges";
	char trisext[10] = ".tri";
	char zofsext[10] = ".z";
	
	CreateAndOpenFileSingle( &nodeofs, nodesext);
	CreateAndOpenFileSingle( &edgofs,  edgesext);
	CreateAndOpenFileSingle( &triofs,  trisext);
	CreateAndOpenFileSingle( &zofs,    zofsext);
	
	nodeofs.setf( ios::fixed, ios::floatfield);
	zofs.setf   ( ios::fixed, ios::floatfield);
	
#ifdef PARALLEL_TRIBS
  }
#endif

	ReadNodeOutputList(); 
	CreateAndOpenPixel();
	dynvars = nullptr;
}

template< class tSubNode >
tOutput<tSubNode>::~tOutput()
{
	g     = NULL;
	timer = NULL;
	if (nodeList)
		delete [] nodeList;
	if (uzel)
		delete [] uzel;
	if (pixinfo)
		delete [] pixinfo;
	if (dynvars)
		delete [] dynvars;
	Cout<<"tOutput Object has been destroyed..."<<endl<<flush;
}

//=========================================================================
//
//
//                  Section 2: tOutput Open and Write Functions
//
//
//=========================================================================

/*************************************************************************
**
**  tOutput::CreateAndOpenFile
**
**  Opens the output file stream pointed by theOFStream, giving it the
**  name <baseName><extension>, and checks to make sure that the ofstream
**  is valid.
**
**  Input:  theOFStream -- ptr to an ofstream object
**          extension -- file name extension (e.g., ".nodes")
**  Output: theOFStream is initialized to create an open output file
**  Assumes: extension is a null-terminated string, and the length of
**           baseName plus extension doesn't exceed kMaxNameSize+6
**           (ie, the extension is expected to be <= 6 characters)
**
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::CreateAndOpenFile( ofstream *theOFStream,
                                           char *extension )
{
	char fullName[kMaxNameSize+6];
	
	strcpy( fullName, baseName );
	strcat( fullName, extension );
	
#ifdef PARALLEL_TRIBS
// Add processor extension if running in parallel
  char procex[10];
  snprintf( procex,sizeof(procex), ".%-d", tParallel::getMyProc()); //WR--09192023: warning: 'sprintf' is deprecated: This function is provided for compatibility reasons only.  Due to security concerns inherent in the design of sprintf(3), it is highly recommended that you use snprintf(3) instead.
  strcat(fullName, procex);
#endif

	theOFStream->open( fullName );
	
	if ( !theOFStream->good() )
		cerr << "File "<<fullName<<" not created." << endl;
	
/*SMM
	Cout<<"Creating Output File: \t '"<<fullName<<"' "<<endl;
*/
	return;
}

/*************************************************************************
**
**  tOutput::CreateAndOpenFileSingle
**
**  Like CreateAndOpenFile but never appends a processor suffix.
**  Used in parallel mode so the master can write one merged output file.
**
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::CreateAndOpenFileSingle( ofstream *theOFStream,
                                                 char *extension )
{
	char fullName[kMaxNameSize+6];
	strcpy( fullName, baseName );
	strcat( fullName, extension );
	theOFStream->open( fullName );
	if ( !theOFStream->good() )
		cerr << "File " << fullName << " not created." << endl;
}

/*************************************************************************
**
**  tOutput::ReadNodeOutputList()
**
**  Opens and Reads the node list from a *.nol file whose structure is:
**
**  Number of Nodes
**  Node1 Node2 Node3 Node4 Node5 ...
**
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::ReadNodeOutputList() {
	
	ifstream readNOL(nodeFile);
	if (!readNOL) {
		Cout << "\nFile "<<nodeFile<<" not found..."<<endl;
		Cout<<"\tAttention: The specified file with node IDs does not "
			<<"exist.\n\t\t   No node output will be written."<<endl<<endl;
		numNodes = 0;
        nodeList = nullptr;
        uzel = nullptr;
        pixinfo = nullptr;
		return;
	}
	
	readNOL >> numNodes;
	nodeList = new int[numNodes];
	uzel = new tSubNode*[numNodes];
	pixinfo = new ofstream[numNodes];

#ifdef PARALLEL_TRIBS
  // Initialize to NULL
  for (int i = 0; i < numNodes; i++)
    uzel[i] = NULL;
#endif

	for (int i = 0; i < numNodes; i++) {
		readNOL >> nodeList[i]; 
	}
	
	readNOL.close();
	return;
}

/*************************************************************************
**
**  tOutput::CreateAndOpenPixel()
**
**  Write the header for the *.pixel output file.
** 
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::CreateAndOpenPixel()
{
	if ( nodeList ) {
		char pixelext[10] = ".pixel";
		char nodeNum[10], pixelnode[100];
		
      //SMM - Set interior nodes, added 08132008
      SetInteriorNode();

		for (int i = 0; i < numNodes; i++) {
#ifdef PARALLEL_TRIBS
         // Check if node is on this processor
         if ( (uzel[i] != NULL) && (nodeList[i] >= 0) ) {
#else
			if (nodeList[i] >= 0) {
#endif
				snprintf(nodeNum, sizeof(nodeNum), "%d",nodeList[i]);
				strcpy(pixelnode, nodeNum);
				strcat(pixelnode, pixelext);		 
				
				CreateAndOpenFile( &pixinfo[i], pixelnode );
				
				// Write Header
				pixinfo[i]<<"NodeID,"  //1
				<<"Time_hr,"
				<<"Nwt_mm,"
				<<"Nf_mm,"
				<<"Nt_mm,"             //5
				<<"Mu_mm,"
				<<"Mi_mm,"
				<<"QpOut_mm_h,"
				<<"QpIn_mm_h,"
				<<"Trnsm_m2_h,"        //10
				<<"GWflx_m3_h,"
				<<"Srf_Hour_mm,"
				<<"Rain_mm_h,"
				<<"SoilMoist_[],"
				<<"RootMoist_[],"      //15
				<<"AirT_oC,"
				<<"DewT_oC,"
				<<"SurfT_oC,"
				<<"SoilT_oC,"
				<<"Press_Pa,"          //20
				<<"RelHum_[],"
				<<"SkyCov_[],"
				<<"Wind_m_s,"
				<<"NetRad_W_m2,"
				<<"ShrtRadIn_W_m2,"    //25
				<<"ShortRadInSlope_W_m2,"    // JB2025 @ ASU
				<<"ShrtRadIn_dir_W_m2,"
				<<"ShrtRadIn_dif_W_m2,"
				<<"ShortAbsbVeg_W_m2,"
				<<"ShortAbsbSoi_W_m2," //30
				<<"LngRadIn_W_m2,"
				<<"LngRadOut_W_m2A,"
				<<"PotEvp_mm_h,"
				<<"ActEvp_mm_h,"
				<<"EvpTtrs_mm_h,"      //35
				<<"EvpWetCan_mm_h,"
				<<"EvpDryCan_mm_h,"
				<<"EvpSoil_mm_h,"
				<<"Gflux_W_m2,"
				<<"HFlux_W_m2,"        //40
				<<"Lflux_W_m2,"
				<<"NetPrecip_mm_hr,"
				<<"LiqWE_cm,"
				<<"IceWE_cm,"
				<<"SnWE_cm,"           //45
				<<"SnSub_cm,"
				<<"SnEvap_cm,"
				<<"U_kJ_m2,"
				<<"RouteWE_cm,"
				<<"SnTemp_C,"          //50
				<<"SurfAge_h,"
				<<"SnDepth_cm,"
				<<"SnDensity_kg_m3,"
				<<"DU_kJ_m2,"
				<<"snLHF_kJ_m2,"       //55
				<<"snSHF_kJ_m2,"
				<<"snGHF_kJ_m2,"
				<<"snPHF_kJ_m2,"
				<<"snRLout_kJ_m2,"
				<<"snRLin_kJ_m2,"      //60
				<<"snRSin_kJ_m2,"
				<<"Uerror_kJ_m2,"
				<<"IntSWEq_cm,"
				<<"IntSub_cm,"
				<<"IntSnUnload_cm,"     //65
				<<"CanStorage_mm,"
				<<"CumIntercept_mm,"
				<<"Interception_mm,"
				<<"Recharge_mm/hr,"
				<<"RunOn_mm,"          //70
				<<"Qstrm_m3_s,"
				<<"Hlevel_m,"
				<<"ThroughFall_[],"
				<<"CanFieldCap_mm,"
				<<"DrainCoeff_mm_hr,"  //75
				<<"DrainExpPar_1_mm,"
				<<"LandUseAlb_[],"
				<<"VegHeight_m,"
				<<"OptTransmCoeff_[],"
				<<"StomRes_s_m,"       //80
				<<"VegFraction[],"
				<<"LeafAI_[],"
				<<"RootZoneDepth_m"    //83
				<<"\n";
				
				pixinfo[i].setf( ios::right, ios::adjustfield );
				pixinfo[i].setf( ios::fixed, ios::floatfield);  
			}
		}
	}
	return;
}

/*************************************************************************
**
**  tOutput::CreateAndOpenDynVar()
**
**  Opens a number of files for selected dynamic variables 
** 
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::CreateAndOpenDynVar()
{
	if (g->getNodeList()->getActiveSize() < 0) {
		
		// Allocate memory for ofstream objects
		//dynvars = new ofstream[19];
		dynvars = new ofstream[36]; // SKY2008Snow, AJR2008
		
		char Nwt[10] = "_Nwt";
		char Nf[10]  = "_Nf";
		char Nt[10]  = "_Nt";
		char Mu[10]  = "_Mu";
		char Mi[10]  = "_Mi";
		char QOut[10] = "_QOut";
		char QIn[10]  = "_QIn";
		char Ru[10]    = "_Ru";
		char Trsm[10]  = "_Trsm";
		char GWflx[10] = "_GWflx";
		char Srf[10]   = "_Srf";
		char SMS[10]   = "_SMS";
		char SMRT[10]  = "_SMRT";
		char ET[10]    = "_ET";
		char Edc[10]  = "_Edc";
		char Ewc[10]  = "_Ewc";
		char Eso[10]  = "_Eso";
		char NetP[10] = "_NetP";
		char Qbot[10] = "_Qbot";
	
		// SKY2008Snow from AJR2007
		char SnWE[10] = "_SnWE";	    //added by AJR 2007 @ NMT
		char SnTempC[10] = "_snTempC";  //added by AJR 2007 @ NMT
		char IceWE[10] = "_IcWE";	    //added by AJR 2007 @ NMT
		char LiqWE[10] = "_LiWE";	    //added by AJR 2007 @ NMT
		char DU[10] = "_DU";	    //added by AJR 2007 @ NMT
		char Upack[10] = "_UTOT";	    //added by AJR 2007 @ NMT
		char snLHF[10] = "_sLHF";	    //added by AJR 2007 @ NMT
		char snSHF[10] = "_sSHF";	    //added by AJR 2007 @ NMT
		char snGHF[10] = "_sGHF";	    //added by AJR 2007 @ NMT
		char snPHF[10] = "_sPHF";	    //added by AJR 2007 @ NMT
		char snRLo[10] = "_sRLo";	    //added by AJR 2007 @ NMT
		char snRLi[10] = "_sRLi";	    //added by AJR 2007 @ NMT
		char snRSi[10] = "_sRSi";	    //added by AJR 2007 @ NMT
		char Uerr[10] = "_Uerr";	    //added by AJR 2007 @ NMT
		char IntSn[10] = "_IntSn";	    //added by AJR 2007 @ NMT
		char IntSub[10] = "_IntSub";    //added by AJR 2007 @ NMT
		char IntUnl[10] = "_IntUnl";    //added by AJR 2007 @ NMT


		//for (int i = 0; i < 19; i++) {
		for (int i = 0; i < 36; i++) { // SKY2008Snow from AJR2007
			if (i == 0)
				CreateAndOpenFile( &dynvars[i], Nwt);
			else if (i == 1)
				CreateAndOpenFile( &dynvars[i], Nf);
			else if (i == 2)
				CreateAndOpenFile( &dynvars[i], Nt);
			else if (i == 3)
				CreateAndOpenFile( &dynvars[i], Mu);
			else if (i == 4)
				CreateAndOpenFile( &dynvars[i], Mi);
			else if (i == 5)
				CreateAndOpenFile( &dynvars[i], QOut);
			else if (i == 6)
				CreateAndOpenFile( &dynvars[i], QIn);
			else if (i == 7)
				CreateAndOpenFile( &dynvars[i], Ru);
			else if (i == 8)
				CreateAndOpenFile( &dynvars[i], Trsm);
			else if (i == 9)
				CreateAndOpenFile( &dynvars[i], GWflx);
			else if (i == 10)
				CreateAndOpenFile( &dynvars[i], Srf);
			else if (i == 11)
				CreateAndOpenFile( &dynvars[i], SMS);
			else if (i == 12)
				CreateAndOpenFile( &dynvars[i], SMRT);
			else if (i == 13)
				CreateAndOpenFile( &dynvars[i], ET);
			else if (i == 14)
				CreateAndOpenFile( &dynvars[i], Edc);
			else if (i == 15)
				CreateAndOpenFile( &dynvars[i], Ewc);
			else if (i == 16)
				CreateAndOpenFile( &dynvars[i], Eso);
			else if (i == 17)
				CreateAndOpenFile( &dynvars[i], NetP);
			else if (i == 18)

			// SKY2008Snow from AJR2007
				CreateAndOpenFile( &dynvars[i], SnWE);	    //added by AJR 2007 @ NMT
			else if (i == 19)
				CreateAndOpenFile( &dynvars[i], SnTempC);   //added by AJR 2007 @ NMT
			else if (i == 20)
				CreateAndOpenFile( &dynvars[i], IceWE);	    //added by AJR 2007 @ NMT
			else if (i == 21)
				CreateAndOpenFile( &dynvars[i], LiqWE);	    //added by AJR 2007 @ NMT
			else if (i == 22)
				CreateAndOpenFile( &dynvars[i], DU);	    //added by AJR 2007 @ NMT
			else if (i == 23 )
				CreateAndOpenFile( &dynvars[i], Upack);	    //added by AJR 2007 @ NMT
			else if (i == 24)
				CreateAndOpenFile( &dynvars[i], snLHF);	    //added by AJR 2007 @ NMT
			else if (i == 25)
				CreateAndOpenFile( &dynvars[i], snSHF);	    //added by AJR 2007 @ NMT
			else if (i == 26)
				CreateAndOpenFile( &dynvars[i], snGHF);	    //added by AJR 2007 @ NMT
			else if (i == 27)
				CreateAndOpenFile( &dynvars[i], snPHF);	    //added by AJR 2007 @ NMT
			else if ( i == 28)
				CreateAndOpenFile( &dynvars[i], snRLo);	    //added by AJR 2007 @ NMT
			else if ( i == 29 )
				CreateAndOpenFile( &dynvars[i], snRLi);	    //added by AJR 2007 @ NMT
			else if ( i == 30 )
				CreateAndOpenFile( &dynvars[i], snRSi);	    //added by AJR 2007 @ NMT
			else if (i == 31)
				CreateAndOpenFile( &dynvars[i], Uerr);	    //added by AJR 2007 @ NMT
			else if (i == 32)
				CreateAndOpenFile( &dynvars[i], IntSn);	    //added by AJR 2007 @ NMT
			else if (i == 33)
				CreateAndOpenFile( &dynvars[i], IntSub);    //added by AJR 2007 @ NMT
			else if (i == 34)
				CreateAndOpenFile( &dynvars[i], IntUnl);    //added by AJR 2007 @ NMT
			else if (i == 35)

				CreateAndOpenFile( &dynvars[i], Qbot);
			dynvars[i].setf( ios::right, ios::adjustfield );
			dynvars[i].setf( ios::fixed, ios::floatfield);
		}
	}
	return;
}

/*************************************************************************
**
**  tOutput::WriteOutput
**
**  This function writes information about the mesh to four files called
**  name.nodes, name.edges, name.tri, and name.z, where "name" is a
**  name that the user has specified in the input file and which is
**  stored in the data member baseName.
**
**  Input: time -- time of the current output time-slice
**  Output: the node, edge, and triangle ID numbers are modified so that
**          they are numbered according to their position on the list
**  Assumes: the four file ofstreams have been opened by the constructor
**           and are valid
**
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::WriteOutput( double time )
{
	tNode * cn;
	tEdge * ce;
	tTriangle * ct;

	tMeshListIter<tSubNode> niter( g->getNodeList() );
	tMeshListIter<tEdge>    eiter( g->getEdgeList() );
	tListIter<tTriangle>    titer( g->getTriList() );
	
	int nnodes = g->getNodeList()->getSize();
	int nedges = g->getEdgeList()->getSize();
	int ntri   = g->getTriList()->getSize();
	
#ifdef PARALLEL_TRIBS
	// If running parallel, sum sizes on all processors
	int nGlobalActiveNodes = g->getNodeList()->getGlobalActiveSize();
	int nGlobalActiveEdges = g->getEdgeList()->getGlobalActiveSize();
	int nActiveNodes = g->getNodeList()->getActiveSize();
	int nActiveEdges = g->getEdgeList()->getActiveSize();

	cout<<"Proc " << tParallel::getMyProc()
	    <<": tOutput Characteristics:"<<endl;
  
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of nodes: \t\t"<<nnodes<<endl;
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of edges: \t\t"<<nedges<<endl;
  
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of active nodes: \t"<<nActiveNodes<<endl<<flush;
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of active edges: \t"<<nActiveEdges<<endl<<flush;
  
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of global active nodes: \t"<<nGlobalActiveNodes<<endl<<flush;
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of global active edges: \t"<<nGlobalActiveEdges<<endl<<flush;

#else

	int nActiveNodes = g->getNodeList()->getActiveSize();
	int nActiveEdges = g->getEdgeList()->getActiveSize();

	cout<<"\ntOutput Characteristics:"<<endl<<endl;
	cout<<"Number of nodes: \t\t"<<nnodes<<endl;
	cout<<"Number of edges: \t\t"<<nedges<<endl;
	cout<<"Number of triangles: \t\t"<<ntri<<endl;
	cout<<"Number of active nodes: \t"
		<<nActiveNodes<<endl<<flush;
	cout<<"Number of active edges: \t"
		<<nActiveEdges<<endl<<flush;

#endif

	if (nnodes > 0 && nedges > 0 && ntri > 0) {
#ifdef PARALLEL_TRIBS
  // Master node ONLY writes these files
  if (tParallel::isMaster()) {
#endif

	// Formating floating point numbers in output
	nodeofs.setf(ios::fixed, ios::floatfield);
	zofs.setf   (ios::fixed, ios::floatfield);
	edgofs.setf (ios::fixed, ios::floatfield);
	triofs.setf (ios::fixed, ios::floatfield);
	
	// Write 'node' file and 'z' file 
	nodeofs<<" "<<time<<"\n"<<nnodes<<"\n";
	zofs   <<" "<<time<<"\n"<<nnodes<<"\n";
	
	for ( cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() ) {
		nodeofs<<cn->getX()<<" "<<cn->getY()<<" "
		<<cn->getEdg()->getID()<<" "<<cn->getBoundaryFlag()<<"\n";
		zofs   <<cn->getZ()<<"\n";
	}
	
	// Write 'edge' file 
	edgofs<<" "<<time<<"\n"<<nedges<<"\n";
	for ( ce=eiter.FirstP(); !(eiter.AtEnd()); ce=eiter.NextP() )
		edgofs <<ce->getOriginPtrNC()->getID()<<" "
			<<ce->getDestinationPtrNC()->getID()<<" "
			<<ce->getCCWEdg()->getID()<<"\n";
	
	// Write 'triangle' file 
	int i;
	triofs<<" "<<time<<"\n"<<ntri<<"\n";
	for ( ct=titer.FirstP(); !(titer.AtEnd()); ct=titer.NextP() )  {
		for ( i=0; i<=2; i++ )
			triofs<<ct->pPtr(i)->getID()<<" ";
		for ( i=0; i<=2; i++ )  {
			if ( ct->tPtr(i) ) 
				triofs<<ct->tPtr(i)->getID()<<" ";
			else 
				triofs<<"-1 ";
		}
		triofs << ct->ePtr(0)->getID() << " " 
			<< ct->ePtr(1)->getID() << " " 
			<< ct->ePtr(2)->getID() << "\n";
	}
	
#ifdef PARALLEL_TRIBS
  }
#endif
	}

	// Set interior nodes
	SetInteriorNode();
	
	return;
}

/*************************************************************************
**
**  tOutput::SetInteriorNode()
**
**  Initializes pointers to basin interior outlets 
** 
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::SetInteriorNode()
{
	tSubNode * cnn;
	tMeshListIter<tSubNode> niter( g->getNodeList() );
	
	if (nodeList) {  // <--- If the node list in NOT empty
		for (int i=0; i < numNodes; i++) {
			if (nodeList[i] >= 0) {
#ifdef PARALLEL_TRIBS
           // The processor with the node of interest creates the file
           for( cnn=niter.FirstP(); niter.IsActive(); cnn=niter.NextP() ) {
#else
				for ( cnn=niter.FirstP(); !(niter.AtEnd()); cnn=niter.NextP() ) {
#endif
					if (cnn->getID() == nodeList[i]) {
						uzel[i] = cnn;   // <=== Defining node of interest 
						cout<<"\nNode of Interest ID: \t"<<nodeList[i]
							<<" has been set up..."<<endl<<flush;
					}
				}
			}
		}
	}
	return;
}

//=========================================================================
//
//
//                  Section 3: tOutput Void Write Functions
//
//
//=========================================================================
template< class tSubNode >
void tOutput<tSubNode>::WriteNodeData( double ) 
{
	cout<<"tOutput:WriteNodeData VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::WriteNodeData( double , tResample *)
{
	cout<<"tOutput:WriteNodeData VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::WriteDynamicVars( double )
{
	cout<<"tOutput:WriteDynamicVars VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::WritePixelInfo( double )
{
	cout<<"tOutput:WritePixelInfo VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::end_simulation()
{
	for (int i = 0; i < numNodes; i++) {
#ifdef PARALLEL_TRIBS
        // Check if node is on this processor
        if ((uzel[i] != NULL) && (nodeList[i] >= 0))
#endif
            pixinfo[i].close();

    }
}

//=========================================================================
//
//
//                  Section 4: tCOutput Derived Class (tRIBS Node)
//
//
//=========================================================================

/*************************************************************************
**
**  tCOutput::Constructors and Destructor
**
**  Constructor for the derived tCNode class. In addition to the 
**  inhereted constructor functions, it also opens a _voi file for
**  writing the voronoi mesh points as a ArcInfo generate file
**
*************************************************************************/
template< class tSubNode >
tCOutput<tSubNode>::tCOutput(SimulationControl *simCtrPtr, tMesh<tSubNode> *g, 
							 tInputFile &infile, tResample *resamp, 
							 tRunTimer *timptr) 
: tOutput<tSubNode>(simCtrPtr, g, infile, resamp, timptr)
{
	simCtrl = simCtrPtr;
	char nodeFileO[kMaxNameSize];

	// Set up interior outlets and open .qout files
	infile.ReadItem( outletName, "OUTHYDROFILENAME");
	infile.ReadItem( nodeFileO, "OUTLETNODELIST" );
	ReadOutletNodeList(nodeFileO);

#ifdef PARALLEL_TRIBS
  // Change the order when running in parallel,
  // so outlets are set per processor
  SetInteriorOutlet();
  CreateAndOpenOutlet();
#else
	CreateAndOpenOutlet();
	SetInteriorOutlet();
#endif

	{
		std::set<std::string> selection;
		char dynVarPath[kMaxNameSize];
		dynVarPath[0] = '\0';
		if (infile.IsItemIn("DYNVARFILE")) {
			infile.ReadItem(dynVarPath, "DYNVARFILE");
			ReadDynVarFile(dynVarPath, selection);
		}
		BuildDynVarTable(selection);
	}

	WriteNodeData( 0, resamp );
}

template< class tSubNode >
tCOutput<tSubNode>::tCOutput(SimulationControl *simCtrPtr, tMesh<tSubNode> *g,
							 tInputFile &infile, tResample *resamp )
: tOutput<tSubNode>(simCtrPtr, g, infile, resamp, this->timptr)
{
	BuildDynVarTable({});
}

template< class tSubNode >
tCOutput<tSubNode>::~tCOutput() 
{

	// GMnSKY2008MLE to fix memory leaks
	if (numOutlets > 0) {
		for (int j=0; j < numOutlets; j++) 
#ifdef PARALLEL_TRIBS
    if ( (Outlets[j] != NULL) && (OutletList[j] > 0) )
#endif
			outletinfo[j].close(); 
		delete [] OutletList; 
		delete [] Outlets; 
		delete [] outletinfo; 
	}       

    Cout<<"tCOutput Object has been destroyed..."<<endl<<flush;
}

/*************************************************************************
**
**  tCOutput::WriteNodeData()
**
**  Writes the voronoi vertices to a file *_voi in a format compatible
**  with ArcInfo generate files for input and transformation into
**  a polygon coverage. The output format should be readable by ArcInfo 
**  & Matlab
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteNodeData( double time )
{
	tSubNode *cn;
	tEdge *firstedg;
	tEdge *curedg;
	tMeshListIter<tSubNode> ni( this->g->getNodeList() );
	tArray<double> xy;
	char voiext[10] = "_voi";

	if (time == 0) {
#ifdef PARALLEL_TRIBS
		ostringstream buf;
		buf.setf(ios::fixed, ios::floatfield);
		cn = ni.FirstP();
		while (ni.IsActive()) {
			buf << cn->getID() << ',' << cn->getX() << ',' << cn->getY() << "\n";
			firstedg = cn->getFlowEdg();
			xy = firstedg->getRVtx();
			buf << xy[0] << ',' << xy[1] << "\n";
			curedg = firstedg->getCCWEdg();
			while (curedg != firstedg) {
				xy = curedg->getRVtx();
				buf << xy[0] << ',' << xy[1] << "\n";
				curedg = curedg->getCCWEdg();
			}
			buf << "END\n";
			cn = ni.NextP();
		}
		{
			string localStr = buf.str();
			int localLen = (int)localStr.size();
			int nprocs = tParallel::getNumProcs();
			vector<int> lengths(nprocs), offsets(nprocs);
			MPI_Gather(&localLen, 1, MPI_INT,
			           lengths.data(), 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
			if (tParallel::isMaster()) {
				int total = 0;
				for (int p = 0; p < nprocs; p++) { offsets[p] = total; total += lengths[p]; }
				vector<char> combined(total);
				MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
				            combined.data(), lengths.data(), offsets.data(), MPI_CHAR,
				            MASTER_PROC, MPI_COMM_WORLD);
				string combinedStr(combined.data(), total);
				vector<string> blocks;
				const string delim = "END\n";
				size_t pos = 0;
				while (pos < (size_t)total) {
					size_t end = combinedStr.find(delim, pos);
					if (end == string::npos) break;
					blocks.push_back(combinedStr.substr(pos, end - pos + delim.size()));
					pos = end + delim.size();
				}
				sort(blocks.begin(), blocks.end(), [](const string &a, const string &b) {
					return stoi(a) < stoi(b);
				});
				this->CreateAndOpenFileSingle(&vorofs, voiext);
				vorofs.setf(ios::fixed, ios::floatfield);
				for (const auto &blk : blocks) vorofs << blk;
				vorofs << "END\n";
				vorofs.close();
			} else {
				MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
				            nullptr, nullptr, nullptr, MPI_CHAR,
				            MASTER_PROC, MPI_COMM_WORLD);
			}
		}
#else
		this->CreateAndOpenFileSingle(&vorofs, voiext);
		vorofs.setf(ios::fixed, ios::floatfield);
		cn = ni.FirstP();
		while (ni.IsActive()) {
			vorofs << cn->getID() << ',' << cn->getX() << ',' << cn->getY() << "\n";
			firstedg = cn->getFlowEdg();
			xy = firstedg->getRVtx();
			vorofs << xy[0] << ',' << xy[1] << "\n";
			curedg = firstedg->getCCWEdg();
			while (curedg != firstedg) {
				xy = curedg->getRVtx();
				vorofs << xy[0] << ',' << xy[1] << "\n";
				curedg = curedg->getCCWEdg();
			}
			vorofs << "END\n";
			cn = ni.NextP();
		}
		vorofs << "END\n";
		vorofs.close();
#endif
	}
	return;
}

/*************************************************************************
**
**  tCOutput::WriteNodeData( double time, tResample *tresamp )
**
**  - Writes a file containing voronoi vertices: reads arrays
**    from the tResample object
**  
**  - Writes a file containing drainage areas for stream network
** 
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteNodeData( double time, tResample *tresamp )
{
	tSubNode *cn;
	tMeshListIter<tSubNode> ni( this->g->getNodeList() );
	int k, i = 0;
	char voiext[10]   = "_voi";
	char areaext[10]  = "_area";
	char widthext[10] = "_width";

	if (time == 0) {
#ifdef PARALLEL_TRIBS
		ostringstream voiBuf, areaBuf, widthBuf;
		voiBuf.setf(ios::fixed, ios::floatfield);
		areaBuf.setf(ios::fixed, ios::floatfield);
		widthBuf.setf(ios::fixed, ios::floatfield);

		if (tParallel::isMaster()) {
			areaBuf << "ID\tX\tY\tCArea\n";
			widthBuf << "ID\tX\tY\tWidth\tEdgL\tSlp\n";
		}

		cn = ni.FirstP();
		while (ni.IsActive()) {
			voiBuf << cn->getID() << ',' << cn->getX() << ',' << cn->getY() << "\n";
			for (k = 0; k < tresamp->nPoints[i]; k++)
				voiBuf << tresamp->vXs[i][k] << "," << tresamp->vYs[i][k] << "\n";
			voiBuf << "END\n";

			if (cn->getBoundaryFlag() == 3) {
				areaBuf << cn->getID()
				        << "\t" << cn->getX() << "\t" << cn->getY()
				        << "\t" << cn->getContrArea() << "\n";
				widthBuf << cn->getID()
				         << "\t" << cn->getX() << "\t" << cn->getY()
				         << "\t" << cn->getChannelWidth()
				         << "\t" << cn->getFlowEdg()->getLength()
				         << "\t" << cn->getFlowEdg()->getSlope() << "\n";
			}
			i++;
			cn = ni.NextP();
		}
		if (tGraph::hasLastReach()) {
			if (cn->getBoundaryFlag() == 2) {
				areaBuf << cn->getID()
				        << "\t" << cn->getX() << "\t" << cn->getY()
				        << "\t" << cn->getContrArea() << "\n";
				widthBuf << cn->getID()
				         << "\t" << cn->getX() << "\t" << cn->getY()
				         << "\t" << cn->getChannelWidth() << "\t0.0\t0.0\n";
			}
		}

		// Gather and write _voi: sort multi-line polygon blocks by node ID
		{
			string localStr = voiBuf.str();
			int localLen = (int)localStr.size();
			int nprocs = tParallel::getNumProcs();
			vector<int> lengths(nprocs), offsets(nprocs);
			MPI_Gather(&localLen, 1, MPI_INT,
			           lengths.data(), 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
			if (tParallel::isMaster()) {
				int total = 0;
				for (int p = 0; p < nprocs; p++) { offsets[p] = total; total += lengths[p]; }
				vector<char> combined(total);
				MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
				            combined.data(), lengths.data(), offsets.data(), MPI_CHAR,
				            MASTER_PROC, MPI_COMM_WORLD);
				string combinedStr(combined.data(), total);
				vector<string> blocks;
				const string delim = "END\n";
				size_t pos = 0;
				while (pos < (size_t)total) {
					size_t end = combinedStr.find(delim, pos);
					if (end == string::npos) break;
					blocks.push_back(combinedStr.substr(pos, end - pos + delim.size()));
					pos = end + delim.size();
				}
				sort(blocks.begin(), blocks.end(), [](const string &a, const string &b) {
					return stoi(a) < stoi(b);
				});
				this->CreateAndOpenFileSingle(&vorofs, voiext);
				vorofs.setf(ios::fixed, ios::floatfield);
				for (const auto &blk : blocks) vorofs << blk;
				vorofs << "END\n";
				vorofs.close();
			} else {
				MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
				            nullptr, nullptr, nullptr, MPI_CHAR,
				            MASTER_PROC, MPI_COMM_WORLD);
			}
		}

		// Gather and write _area: header from master, rows sorted by node ID
		{
			string localStr = areaBuf.str();
			int localLen = (int)localStr.size();
			int nprocs = tParallel::getNumProcs();
			vector<int> lengths(nprocs), offsets(nprocs);
			MPI_Gather(&localLen, 1, MPI_INT,
			           lengths.data(), 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
			if (tParallel::isMaster()) {
				int total = 0;
				for (int p = 0; p < nprocs; p++) { offsets[p] = total; total += lengths[p]; }
				vector<char> combined(total);
				MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
				            combined.data(), lengths.data(), offsets.data(), MPI_CHAR,
				            MASTER_PROC, MPI_COMM_WORLD);
				string combinedStr(combined.data(), total);
				size_t headerEnd = combinedStr.find('\n');
				string header = combinedStr.substr(0, headerEnd + 1);
				vector<string> rows;
				size_t pos = headerEnd + 1;
				while (pos < (size_t)total) {
					size_t end = combinedStr.find('\n', pos);
					if (end == string::npos) break;
					rows.push_back(combinedStr.substr(pos, end - pos + 1));
					pos = end + 1;
				}
				sort(rows.begin(), rows.end(), [](const string &a, const string &b) {
					return stoi(a) < stoi(b);
				});
				this->CreateAndOpenFileSingle(&drareaofs, areaext);
				drareaofs.setf(ios::fixed, ios::floatfield);
				drareaofs << header;
				for (const auto &row : rows) drareaofs << row;
				drareaofs.close();
			} else {
				MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
				            nullptr, nullptr, nullptr, MPI_CHAR,
				            MASTER_PROC, MPI_COMM_WORLD);
			}
		}

		// Gather and write _width: header from master, rows sorted by node ID
		{
			string localStr = widthBuf.str();
			int localLen = (int)localStr.size();
			int nprocs = tParallel::getNumProcs();
			vector<int> lengths(nprocs), offsets(nprocs);
			MPI_Gather(&localLen, 1, MPI_INT,
			           lengths.data(), 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
			if (tParallel::isMaster()) {
				int total = 0;
				for (int p = 0; p < nprocs; p++) { offsets[p] = total; total += lengths[p]; }
				vector<char> combined(total);
				MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
				            combined.data(), lengths.data(), offsets.data(), MPI_CHAR,
				            MASTER_PROC, MPI_COMM_WORLD);
				string combinedStr(combined.data(), total);
				size_t headerEnd = combinedStr.find('\n');
				string header = combinedStr.substr(0, headerEnd + 1);
				vector<string> rows;
				size_t pos = headerEnd + 1;
				while (pos < (size_t)total) {
					size_t end = combinedStr.find('\n', pos);
					if (end == string::npos) break;
					rows.push_back(combinedStr.substr(pos, end - pos + 1));
					pos = end + 1;
				}
				sort(rows.begin(), rows.end(), [](const string &a, const string &b) {
					return stoi(a) < stoi(b);
				});
				this->CreateAndOpenFileSingle(&widthsofs, widthext);
				widthsofs.setf(ios::fixed, ios::floatfield);
				widthsofs << header;
				for (const auto &row : rows) widthsofs << row;
				widthsofs.close();
			} else {
				MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
				            nullptr, nullptr, nullptr, MPI_CHAR,
				            MASTER_PROC, MPI_COMM_WORLD);
			}
		}
#else
		this->CreateAndOpenFileSingle(&vorofs, voiext);
		this->CreateAndOpenFileSingle(&drareaofs, areaext);
		this->CreateAndOpenFileSingle(&widthsofs, widthext);
		vorofs.setf(ios::fixed, ios::floatfield);
		drareaofs.setf(ios::fixed, ios::floatfield);
		widthsofs.setf(ios::fixed, ios::floatfield);
		drareaofs << "ID\tX\tY\tCArea\n";
		widthsofs << "ID\tX\tY\tWidth\tEdgL\tSlp\n";

		cn = ni.FirstP();
		while (ni.IsActive()) {
			vorofs << cn->getID() << ',' << cn->getX() << ',' << cn->getY() << "\n";
			for (k = 0; k < tresamp->nPoints[i]; k++)
				vorofs << tresamp->vXs[i][k] << "," << tresamp->vYs[i][k] << "\n";
			vorofs << "END\n";

			if (cn->getBoundaryFlag() == 3) {
				drareaofs << cn->getID()
				          << "\t" << cn->getX() << "\t" << cn->getY()
				          << "\t" << cn->getContrArea() << "\n";
				widthsofs << cn->getID()
				          << "\t" << cn->getX() << "\t" << cn->getY()
				          << "\t" << cn->getChannelWidth()
				          << "\t" << cn->getFlowEdg()->getLength()
				          << "\t" << cn->getFlowEdg()->getSlope() << "\n";
			}
			i++;
			cn = ni.NextP();
		}

		if (cn->getBoundaryFlag() == 2) {
			drareaofs << cn->getID()
			          << "\t" << cn->getX() << "\t" << cn->getY()
			          << "\t" << cn->getContrArea() << "\n";
			widthsofs << cn->getID()
			          << "\t" << cn->getX() << "\t" << cn->getY()
			          << "\t" << cn->getChannelWidth() << "\t0.0\t0.0\n";
		}
		vorofs << "END\n";
		vorofs.close();
		drareaofs.close();
		widthsofs.close();
#endif
	}
	return;
}

/*************************************************************************
**
**  tCOutput::WritePixelInfo()
**
**  Writes the dynamic variables of node of interest to a *.pixel file
**  The output format should be readable by ArcInfo & Matlab 
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WritePixelInfo( double time )
{ 
	if (this->numNodes > 0) {
		int hour, minute;
		char extension[20];
		
		hour   = (int)floor(time);
		minute = (int)((time-hour)*100);
        snprintf(extension,sizeof(extension),"%04d.%02d", hour, minute);
		
		// Writing to a file dynamic variables of node of interest  
		// The output format should be readable by ArcInfo & Matlab 
		for (int i = 0; i < this->numNodes; i++) {

#ifdef PARALLEL_TRIBS
  // Doesn't need to be less than active size
      if ( (this->uzel[i] != NULL) && (this->nodeList[i] >= 0) ) {
#else
			if ( this->uzel[i] && this->nodeList[i] < this->g->getNodeList()->getActiveSize()) {
#endif
				// CJC2025: Correct utputs by dividing by cos_slope
				// This code only runs for valid, local nodes.
				tEdge *flowEdge = this->uzel[i]->getFlowEdg();

				// Check if the edge itself is valid.
				double cos_slope = 1.0; // Default to 1.0 (no slope correction)
				if (flowEdge) {
					double slope_rad = atan(flowEdge->getSlope());
					cos_slope = cos(slope_rad);
					// Check to prevent division by zero, just in case.
					if (cos_slope < 1E-9) cos_slope = 1.E-9; 
				}
				
				this->pixinfo[i]<<this->nodeList[i]
				<<","<<extension
				/* 3 */   << "," << (this->uzel[i]->getNwtNew() / cos_slope) << ","

				<<setprecision(4)
				<< this->uzel[i]->getNfNew() / cos_slope << ","
				/* 5 */   << this->uzel[i]->getNtNew() / cos_slope << ","

				<< this->uzel[i]->getMuNew() / cos_slope << ","
				<< this->uzel[i]->getMiNew() / cos_slope << ","
				<< this->uzel[i]->getQpout()*1.E-6/this->uzel[i]->getVArea() << ","
				<< this->uzel[i]->getQpin() *1.E-6/this->uzel[i]->getVArea() << ","
				/* 10 */  << this->uzel[i]->getTransmiss()*1.E-6 << ","

				<< this->uzel[i]->getGwaterChng()*1.E-9 << ","
				<< this->uzel[i]->getSrf_Hr() << ","
				<< this->uzel[i]->getRain() << ","
				<< this->uzel[i]->getSoilMoistureSC() << ","
				/* 15 */  << this->uzel[i]->getRootMoistureSC() << ","

				<< this->uzel[i]->getAirTemp() << ","
				<< this->uzel[i]->getDewTemp() << ","
				<< this->uzel[i]->getSurfTemp() << ","
				<< this->uzel[i]->getSoilTemp() << ","
				/* 20 */  << this->uzel[i]->getAirPressure() << ","

				<< this->uzel[i]->getRelHumid() << ","
				<< this->uzel[i]->getSkyCover() << ","
				<< this->uzel[i]->getWindSpeed() << ","
				<< this->uzel[i]->getNetRad() << ","
				/* 25 */  << this->uzel[i]->getShortRadIn() << ","

				<< this->uzel[i]->getShortRadSlope() << "," // JB2025 @ ASU
				<< this->uzel[i]->getShortRadIn_dir() << ","
				<< this->uzel[i]->getShortRadIn_dif() << ","
				<< this->uzel[i]->getShortAbsbVeg() << ","
				/* 30 */  << this->uzel[i]->getShortAbsbSoi() << ","

				<< this->uzel[i]->getLongRadIn() << ","
				<< this->uzel[i]->getLongRadOut() << ","
				<< this->uzel[i]->getPotEvap() << ","
				<< this->uzel[i]->getActEvap() << ","
				/* 35 */  << this->uzel[i]->getEvapoTrans() << ","

				<< this->uzel[i]->getEvapWetCanopy() << ","
				<< this->uzel[i]->getEvapDryCanopy() << ","
				<< this->uzel[i]->getEvapSoil() << ","
				<< this->uzel[i]->getGFlux() << ","
				/* 40 */  << this->uzel[i]->getHFlux() << ","

				<< this->uzel[i]->getLFlux() << ","
				<< this->uzel[i]->getNetPrecipitation() << ","
				// SKY2008Snow from AJR2007
				<< this->uzel[i]->getLiqWE() << ","        //added by AJR 2007 @ NMT
				<< this->uzel[i]->getIceWE() << ","        //added by AJR 2007 @ NMT
				/* 45 */  << (this->uzel[i]->getLiqWE()+this->uzel[i]->getIceWE()) << "," //added by AJR 2007 @ NMT

				<< this->uzel[i]->getSnSub() << ","        // added by CJC2020
				<< this->uzel[i]->getSnEvap() << ","       // added by CJC2020
				<< this->uzel[i]->getUnode() << ","        //added by AJR 2007 @ NMT
				<< this->uzel[i]->getLiqRouted() << ","    //added by AJR 2007 @ NMT
				/* 50 */  << this->uzel[i]->getSnTempC() << ","    //added by AJR 2007 @ NMT

				<< this->uzel[i]->getCrustAge() << ","     //added by AJR 2007 @ NMT
				<< this->uzel[i]->getSnDepth() << ","      // added by CJC2025
				<< this->uzel[i]->getRhoSn() << ","        // added by CJC2025
				<< this->uzel[i]->getDU() << ","           //added by AJR 2007 @ NMT
				/* 55 */  << this->uzel[i]->getSnLHF() << ","      //added by AJR 2007 @ NMT

				<< this->uzel[i]->getSnSHF() << ","        //added by AJR 2007 @ NMT
				<< this->uzel[i]->getSnGHF() << ","        //added by AJR 2007 @ NMT
				<< this->uzel[i]->getSnPHF() << ","        //added by AJR 2007 @ NMT
				<< this->uzel[i]->getSnRLout() << ","      //added by AJR 2007 @ NMT
				/* 60 */  << this->uzel[i]->getSnRLin() << ","     //added by AJR 2007 @ NMT

				<< this->uzel[i]->getSnRSin() << ","       //added by AJR 2007 @ NMT
				<< this->uzel[i]->getUerror() << ","       //added by AJR 2007 @ NMT
				<< this->uzel[i]->getIntSWE() << ","       //added by AJR 2007 @ NMT
				<< this->uzel[i]->getIntSub() << ","       //added by AJR 2007 @ NMT
				/* 65 */  << this->uzel[i]->getIntSnUnload() << "," //added by AJR 2007 @ NMT

				<< this->uzel[i]->getCanStorage() << ","
				<< this->uzel[i]->getCumIntercept() << ","
				<< this->uzel[i]->getInterceptLoss() << ","
				/* 69 */  << this->uzel[i]->getRecharge() << ","
				<< this->uzel[i]->getRunOn() << ",";

				if (this->uzel[i]->getBoundaryFlag() == kStream)
					this->pixinfo[i]<<this->uzel[i]->getQstrm() << ","
					<< this->uzel[i]->getHlevel() << ",";
				else
					this->pixinfo[i]<<"0.0,0.0,";

				// SKYnGM2008LU
				this->pixinfo[i]<<setprecision(4)
				<< this->uzel[i]->getThroughFall() << ","
				<< this->uzel[i]->getCanFieldCap() << ","
				/* 75 */ << this->uzel[i]->getDrainCoeff() << ","
				<< this->uzel[i]->getDrainExpPar() << ","
				<< this->uzel[i]->getLandUseAlb() << ","
				<< this->uzel[i]->getVegHeight() << ","
				<< this->uzel[i]->getOptTransmCoeff() << ","
				/* 80 */ << this->uzel[i]->getStomRes() << ","
				<< this->uzel[i]->getVegFraction() << ","
				<< this->uzel[i]->getLeafAI() << ","
				/* 83 */ << this->uzel[i]->getRootZoneDepth() <<
				"\n" << flush; 

			}
		}
	}
}

/*************************************************************************
**
**  tCOutput::ReadDynVarFile()
**
**  Parses a single-row CSV file of column names into 'selection'.
**  On any error the set is left empty and the caller uses all columns.
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::ReadDynVarFile(const char *path,
                                        std::set<std::string> &selection)
{
	ifstream f(path);
	if (!f) {
		cerr << "\nWarning: DYNVARFILE '" << path
		     << "' not found. Using default dynamic output.\n";
		return;
	}
	string line;
	if (!getline(f, line) || line.empty()) {
		cerr << "\nWarning: DYNVARFILE '" << path
		     << "' is empty. Using default dynamic output.\n";
		return;
	}
	stringstream ss(line);
	string token;
	while (getline(ss, token, ',')) {
		size_t s = token.find_first_not_of(" \t\r\n");
		size_t e = token.find_last_not_of(" \t\r\n");
		if (s != string::npos)
			selection.insert(token.substr(s, e - s + 1));
	}
	Cout << "\n\tDynamic output filtered to " << selection.size()
	     << " variable(s) from: " << path << endl;
}

/*************************************************************************
**
**  tCOutput::BuildDynVarTable()
**
**  Populates activeDynCols with the full column set filtered by 'selection'.
**  An empty selection means all columns are included (default behavior).
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::BuildDynVarTable(const std::set<std::string> &selection)
{
	std::vector<DynVarCol> all = {
		{"Nwt", 5, [](tSubNode *cn) {
			double cs = cos(atan(cn->getFlowEdg()->getSlope()));
			return cn->getNwtNew() / (cs < 1e-9 ? 1e-9 : cs);
		}},
		{"Mu", 5, [](tSubNode *cn) {
			double cs = cos(atan(cn->getFlowEdg()->getSlope()));
			return cn->getMuNew() / (cs < 1e-9 ? 1e-9 : cs);
		}},
		{"Mi", 5, [](tSubNode *cn) {
			double cs = cos(atan(cn->getFlowEdg()->getSlope()));
			return cn->getMiNew() / (cs < 1e-9 ? 1e-9 : cs);
		}},
		{"Nf", 5, [](tSubNode *cn) {
			double cs = cos(atan(cn->getFlowEdg()->getSlope()));
			return cn->getNfNew() / (cs < 1e-9 ? 1e-9 : cs);
		}},
		{"Nt", 5, [](tSubNode *cn) {
			double cs = cos(atan(cn->getFlowEdg()->getSlope()));
			return cn->getNtNew() / (cs < 1e-9 ? 1e-9 : cs);
		}},
		{"Qpout",        5, [](tSubNode *cn) { return cn->getQpout() * 1.e-6 / cn->getVArea(); }},
		{"Qpin",         5, [](tSubNode *cn) { return cn->getQpin()  * 1.e-6 / cn->getVArea(); }},
		{"Srf",          4, [](tSubNode *cn) { return cn->getSrf_Hr(); }},
		{"Rain",         3, [](tSubNode *cn) { return cn->getRain(); }},
		{"ST",           3, [](tSubNode *cn) { return cn->getSnTempC(); }},
		{"IWE",          5, [](tSubNode *cn) { return cn->getIceWE(); }},
		{"LWE",          5, [](tSubNode *cn) { return cn->getLiqWE(); }},
		{"SnSub",        7, [](tSubNode *cn) { return cn->getSnSub(); }},
		{"SnEvap",       7, [](tSubNode *cn) { return cn->getSnEvap(); }},
		{"SnMelt",       7, [](tSubNode *cn) { return cn->getLiqRouted(); }},
		{"SnDepth",      5, [](tSubNode *cn) { return cn->getSnDepth(); }},
		{"Upack",        5, [](tSubNode *cn) { return cn->getUnode(); }},
		{"sLHF",         5, [](tSubNode *cn) { return cn->getSnLHF(); }},
		{"sSHF",         5, [](tSubNode *cn) { return cn->getSnSHF(); }},
		{"sGHF",         5, [](tSubNode *cn) { return cn->getSnGHF(); }},
		{"sPHF",         5, [](tSubNode *cn) { return cn->getSnPHF(); }},
		{"sRLo",         5, [](tSubNode *cn) { return cn->getSnRLout(); }},
		{"sRLi",         5, [](tSubNode *cn) { return cn->getSnRLin(); }},
		{"sRSi",         5, [](tSubNode *cn) { return cn->getSnRSin(); }},
		{"Uerr",         5, [](tSubNode *cn) { return cn->getUerror(); }},
		{"IntSWE",       5, [](tSubNode *cn) { return cn->getIntSWE(); }},
		{"IntSub",       5, [](tSubNode *cn) { return cn->getIntSub(); }},
		{"IntUnl",       5, [](tSubNode *cn) { return cn->getIntSnUnload(); }},
		{"SoilMoist",    3, [](tSubNode *cn) { return cn->getSoilMoistureSC(); }},
		{"RootMoist",    3, [](tSubNode *cn) { return cn->getRootMoistureSC(); }},
		{"CanStorage",   3, [](tSubNode *cn) { return cn->getCanStorage(); }},
		{"ActEvp",       3, [](tSubNode *cn) { return cn->getActEvap(); }},
		{"EvpSoil",      5, [](tSubNode *cn) { return cn->getEvapSoil(); }},
		{"ET",           5, [](tSubNode *cn) { return cn->getEvapoTrans(); }},
		{"GFlux",        3, [](tSubNode *cn) { return cn->getGFlux(); }},
		{"HFlux",        3, [](tSubNode *cn) { return cn->getHFlux(); }},
		{"LFlux",        3, [](tSubNode *cn) { return cn->getLFlux(); }},
		{"Qstrm",        3, [](tSubNode *cn) { return cn->getQstrm(); }},
		{"Hlev",         3, [](tSubNode *cn) { return cn->getHlevel(); }},
		{"FlwVlc",       3, [](tSubNode *cn) { return cn->getFlowVelocity(); }},
		{"ThroughFall",  5, [](tSubNode *cn) { return cn->getThroughFall(); }},
		{"CanFieldCap",  5, [](tSubNode *cn) { return cn->getCanFieldCap(); }},
		{"DrainCoeff",   5, [](tSubNode *cn) { return cn->getDrainCoeff(); }},
		{"DrainExpPar",  5, [](tSubNode *cn) { return cn->getDrainExpPar(); }},
		{"LandUseAlb",   5, [](tSubNode *cn) { return cn->getLandUseAlb(); }},
		{"VegHeight",    5, [](tSubNode *cn) { return cn->getVegHeight(); }},
		{"OptTransmCoeff",5,[](tSubNode *cn) { return cn->getOptTransmCoeff(); }},
		{"StomRes",      5, [](tSubNode *cn) { return cn->getStomRes(); }},
		{"VegFraction",  5, [](tSubNode *cn) { return cn->getVegFraction(); }},
		{"LeafAI",       5, [](tSubNode *cn) { return cn->getLeafAI(); }},
	};

	if (selection.empty()) {
		activeDynCols = std::move(all);
		return;
	}

	for (const auto &col : all) {
		if (selection.count(col.name))
			activeDynCols.push_back(col);
	}

	if (activeDynCols.empty()) {
		cerr << "\nWarning: No valid variable names matched in DYNVARFILE. "
		     << "Using default dynamic output.\n";
		activeDynCols = std::move(all);
	}
}

/*************************************************************************
**
**  tCOutput::WriteDynamicVars()
**
**  Writes a file containing dynamic variables for all the nodes
**  The option is turned on with the '-R' option in the command line
**  
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteDynamicVars( double time )
{
	tSubNode *cn;
	tMeshListIter<tSubNode> ni( this->g->getNodeList() );
	
#ifdef PARALLEL_TRIBS
   int nActiveNodes = this->g->getNodeList()->getGlobalActiveSize();
#else
   int nActiveNodes = this->g->getNodeList()->getActiveSize();
#endif
	int nnodes       = this->g->getNodeList()->getSize();
	
	int hour, minute;
	char extension[20];
	
	hour   = (int)floor(time);
	minute = (int)floor((time-hour)*60);

	// Write Header
	cout<<"\n\tHOUR = "<<hour<<"\tMINUTE = "<<minute<<"\n";

    snprintf(extension,sizeof(extension),".%04d_%02dd", hour, minute);

	ostringstream buf;
	ostream *pout;

#ifdef PARALLEL_TRIBS
	if (tParallel::isMaster()) {
		buf << "ID";
		for (const auto &col : activeDynCols)
			buf << ',' << col.name;
		if (time == 0)
			buf << ',' << "SoilID" << ',' << "LUseID";
		buf << '\n';
	}
	pout = &buf;
#else
	this->CreateAndOpenFile( &arcofs, extension );  //Opens file for writing

    // Write header
	arcofs << "ID";
	for (const auto &col : activeDynCols)
		arcofs << ',' << col.name;
	if (time == 0)
		arcofs << ',' << "SoilID" << ',' << "LUseID" << endl << flush;
	else
		arcofs << '\n';
	pout = &arcofs;
#endif

	cn = ni.FirstP();
	while (ni.IsActive()) {
		*pout << cn->getID();
		for (const auto &col : activeDynCols)
			*pout << ',' << setprecision(col.prec) << col.get(cn);
		if (time == 0)
			*pout << ',' << setprecision(0) << cn->getSoilID()
			      << ',' << setprecision(0) << cn->getLandUse();
		*pout << '\n';
		cn = ni.NextP();
	}

#ifdef PARALLEL_TRIBS
	{
		string localStr = buf.str();
		int localLen = (int)localStr.size();
		int nprocs = tParallel::getNumProcs();
		vector<int> lengths(nprocs), offsets(nprocs);
		MPI_Gather(&localLen, 1, MPI_INT,
		           lengths.data(), 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
		if (tParallel::isMaster()) {
			int total = 0;
			for (int i = 0; i < nprocs; i++) { offsets[i] = total; total += lengths[i]; }
			vector<char> combined(total);
			MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
			            combined.data(), lengths.data(), offsets.data(), MPI_CHAR,
			            MASTER_PROC, MPI_COMM_WORLD);
			string combinedStr(combined.data(), total);
			size_t headerEnd = combinedStr.find('\n');
			string header = combinedStr.substr(0, headerEnd + 1);
			vector<string> rows;
			size_t pos = headerEnd + 1;
			while (pos < (size_t)total) {
				size_t end = combinedStr.find('\n', pos);
				if (end == string::npos) break;
				rows.push_back(combinedStr.substr(pos, end - pos + 1));
				pos = end + 1;
			}
			sort(rows.begin(), rows.end(), [](const string &a, const string &b) {
				return stoi(a) < stoi(b);
			});
			this->CreateAndOpenFileSingle(&arcofs, extension);
			arcofs << header;
			for (const auto &row : rows) arcofs << row;
			arcofs.close();
		} else {
			MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
			            nullptr, nullptr, nullptr, MPI_CHAR,
			            MASTER_PROC, MPI_COMM_WORLD);
		}
	}
#else
	arcofs.close();
#endif


	
	// Call another output function only at the beginning
	// and end of simulation
	if ( time == 0.0 || this->timer->IsFinished() )
		WriteIntegrVars( time );
	
	return;
}

/*************************************************************************
**
**  WriteIntegrVars( double time )
**
**  Writes a file containing variables that contain integrated variables
**  The option is turned on with OPTSPATIAL = 1
**  
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteIntegrVars( double time )
{
	int hour, minute;
	int Occur, prec;
	char extension[20];
	double avRate, tmp1, tmp2;
	tSubNode *cn;
	tMeshListIter<tSubNode> ni( this->g->getNodeList() );
	
	hour   = (int)floor(time);
	minute = (int)floor((time-hour)*60);
	
	snprintf(extension, sizeof(extension), ".%04d_%02di", hour, minute);

	ostringstream buf;
	ostream *pout;

#ifdef PARALLEL_TRIBS
	if (tParallel::isMaster()) {
		buf << "ID" << ','                                     // 1
		    << "BndCd" << ','                                  // 2
		    << "Z" << ','                                      // 3
		    << "VAr" << ','                                    // 4
		    << "CAr" << ','                                    // 5
		    << "Curv" << ','                                   // 6
		    << "EdgL" << ','                                   // 7
		    << "Slp" << ','                                    // 8
		    << "FWidth" << ','                                 // 9
		    << "Aspect" << ','                                 // 10
		    << "SV" << ','                                     // 11
		    << "LV" << ','                                     // 12
		    << "AvSM" << ','                                   // 13
		    << "AvRtM" << ','                                  // 14
		    << "HOccr" << ','                                  // 15
		    << "HRt" << ','                                    // 16
		    << "SbOccr" << ','                                 // 17
		    << "SbRt" << ','                                   // 18
		    << "POccr" << ','                                  // 19
		    << "PRt" << ','                                    // 20
		    << "SatOccr" << ','                                // 21
		    << "SatRt" << ','                                  // 22
		    << "SoiSatOccr" << ','                             // 23
		    << "RchDsch" << ','                                // 24
		    << "AvET" << ','                                   // 25
		    << "EvpFrct" << ','                                // 26
		    << "cET" << ','                                    // 27
		    << "cEsoil" << ','                                 // 28
		    << "cLHF" << ','                                   // 29
		    << "cMelt" << ','                                  // 30
		    << "cSHF" << ','                                   // 31
		    << "cPHF" << ','                                   // 32
		    << "cRLIn" << ','                                  // 33
		    << "cRLo" << ','                                   // 34
		    << "cRSIn" << ','                                  // 35
		    << "cGHF" << ','                                   // 36
		    << "cUErr" << ','                                  // 37
		    << "cHrsSun" << ','                                // 38
		    << "cHrsSnow" << ','                               // 39
		    << "persTime" << ','                               // 40
		    << "peakWE" << ','                                 // 41
		    << "initTime" << ','                               // 42
		    << "peakTime" << ','                               // 43
		    << "cIntSub" << ','                                // 44
		    << "cSnSub" << ','                                 // 45
		    << "cSnEvap" << ','                                // 46
		    << "cIntUnl" << ','                                // 47
		    << "AvTF" << ','                                   // 48
		    << "AvCanFieldCap" << ','                          // 49
		    << "AvDrainCoeff" << ','                           // 50
		    << "AvDrainExpPar" << ','                          // 51
		    << "AvLUAlb" << ','                                // 52
		    << "AvVegHeight" << ','                            // 53
		    << "AvOTCoeff" << ','                              // 54
		    << "AvStomRes" << ','                              // 55
		    << "AvVegFract" << ','                             // 56
		    << "AvLeafAI" << ','                               // 57
		    << "AvEvapThresh" << ','                           // 58
		    << "AvTransThresh" << ','                          // 59
		    << "Bedrock_Depth_mm" << ','                       // 60
		    << "Ks" << ','                                     // 61
		    << "ThetaS" << ','                                 // 62
		    << "ThetaR" << ','                                 // 63
		    << "PoreSize" << ','                               // 64
		    << "AirEBubPress" << ','                           // 65
		    << "DecayF" << ','                                 // 66
		    << "SatAnRatio" << ','                             // 67
		    << "UnsatAnRatio" << ','                           // 68
		    << "Porosity" << ','                               // 69
		    << "VolHeatCond" << ','                            // 70
		    << "SoilHeatCap" << ','                            // 71
		    << "SoilID" << ','                                 // 72
		    << "LandUseID" << ','                              // 73
		    << "AvRootZoneDepth"                               // 74
		    << '\n';
	}
	pout = &buf;
#else
	this->CreateAndOpenFile(&intofs, extension);
	intofs << "ID" << ','                                     // 1
		   << "BndCd" << ','                                  // 2
		   << "Z" << ','                                      // 3
		   << "VAr" << ','                                    // 4
		   << "CAr" << ','                                    // 5
		   << "Curv" << ','                                   // 6
		   << "EdgL" << ','                                   // 7
		   << "Slp" << ','                                    // 8
		   << "FWidth" << ','                                 // 9
		   << "Aspect" << ','                                 // 10
		   << "SV" << ','                                     // 11
		   << "LV" << ','                                     // 12
		   << "AvSM" << ','                                   // 13
		   << "AvRtM" << ','                                  // 14
		   << "HOccr" << ','                                  // 15
		   << "HRt" << ','                                    // 16
		   << "SbOccr" << ','                                 // 17
		   << "SbRt" << ','                                   // 18
		   << "POccr" << ','                                  // 19
		   << "PRt" << ','                                    // 20
		   << "SatOccr" << ','                                // 21
		   << "SatRt" << ','                                  // 22
		   << "SoiSatOccr" << ','                             // 23
		   << "RchDsch" << ','                                // 24
		   << "AvET" << ','                                   // 25
		   << "EvpFrct" << ','                                // 26
		   << "cET" << ','                                    // 27
		   << "cEsoil" << ','                                 // 28
		   << "cLHF" << ','                                   // 29
		   << "cMelt" << ','                                  // 30
		   << "cSHF" << ','                                   // 31
		   << "cPHF" << ','                                   // 32
		   << "cRLIn" << ','                                  // 33
		   << "cRLo" << ','                                   // 34
		   << "cRSIn" << ','                                  // 35
		   << "cGHF" << ','                                   // 36
		   << "cUErr" << ','                                  // 37
		   << "cHrsSun" << ','                                // 38
		   << "cHrsSnow" << ','                               // 39
		   << "persTime" << ','                               // 40
		   << "peakWE" << ','                                 // 41
		   << "initTime" << ','                               // 42
		   << "peakTime" << ','                               // 43
		   << "cIntSub" << ','                                // 44
		   << "cSnSub" << ','                                 // 45
		   << "cSnEvap" << ','                                // 46
		   << "cIntUnl" << ','                                // 47
		   << "AvTF" << ','                                   // 48
		   << "AvCanFieldCap" << ','                          // 49
		   << "AvDrainCoeff" << ','                           // 50
		   << "AvDrainExpPar" << ','                          // 51
		   << "AvLUAlb" << ','                                // 52
		   << "AvVegHeight" << ','                            // 53
		   << "AvOTCoeff" << ','                              // 54
		   << "AvStomRes" << ','                              // 55
		   << "AvVegFract" << ','                             // 56
		   << "AvLeafAI" << ','                               // 57
		   << "AvEvapThresh" << ','                           // 58
		   << "AvTransThresh" << ','                          // 59
		   << "Bedrock_Depth_mm" << ','                       // 60
		   << "Ks" << ','                                     // 61
		   << "ThetaS" << ','                                 // 62
		   << "ThetaR" << ','                                 // 63
		   << "PoreSize" << ','                               // 64
		   << "AirEBubPress" << ','                           // 65
		   << "DecayF" << ','                                 // 66
		   << "SatAnRatio" << ','                             // 67
		   << "UnsatAnRatio" << ','                           // 68
		   << "Porosity" << ','                               // 69
		   << "VolHeatCond" << ','                            // 70
		   << "SoilHeatCap" << ','                            // 71
		   << "SoilID" << ','                                 // 72
		   << "LandUseID" << ','                              // 73
		   << "AvRootZoneDepth"                               // 74
		   << "\n";
	pout = &intofs;
#endif

	cn = ni.FirstP();
	while (ni.IsActive()) {

		*pout<<cn->getID()<<','                              //1
		     <<cn->getBoundaryFlag()<<','                    //2
		     <<setprecision(4)<<cn->getZ()<<','              //3
		     <<setprecision(7)<<cn->getVArea()<<','          //4
		     <<setprecision(7)<<cn->getContrArea()*1.E-6<<',' //5
		     <<setprecision(6)<<cn->getCurvature()<<','      //6
		     <<cn->getFlowEdg()->getLength()<<','            //7
		     <<cn->getFlowEdg()->getSlope()<<','             //8
		     <<cn->getFlowEdg()->getVEdgLen()<<','           //9
		     <<setprecision(4)<<cn->getAspect()<<','         //10
		     <<setprecision(7)<<cn->getSheltFact()<<','      //11
		     <<cn->getLandFact()<<','<<setprecision(4);      //12

		tmp1 = floor(cn->getAvSoilMoisture())*1.E-4;
		tmp2 = (cn->getAvSoilMoisture()-floor(cn->getAvSoilMoisture()))*1.E+1;
		*pout<<tmp1<<','<<                                   //13
		     tmp2<<',';                                      //14

		// -----------  separate runoff mechanism occurrence and rate
		Occur = (int)((cn->hsrfOccur-floor(cn->hsrfOccur))*1.E+6);
		if (Occur > 0) {
			avRate = floor(cn->hsrfOccur)/1000.0/Occur;
			prec = 3;
			if (avRate>100.0) prec++;
		} else {
			avRate = 0.0;
			prec = 0;
		}
		*pout<<setprecision(6)<<Occur<<                      //15
		     setprecision(prec)<<','<<avRate<<',';           //16

		// -----------
		Occur = (int)((cn->sbsrfOccur-floor(cn->sbsrfOccur))*1.E+6);
		if (Occur > 0) {
			avRate = floor(cn->sbsrfOccur)/1000.0/Occur;
			prec = 3;
			if (avRate>100.) prec++;
		} else {
			avRate = 0.0;
			prec = 0;
		}
		*pout<<setprecision(6)<<Occur<<                      //17
		     setprecision(prec)<<','<<avRate<<',';           //18

		// -----------
		Occur = (int)((cn->psrfOccur-floor(cn->psrfOccur))*1.E+6);
		if (Occur > 0) {
			avRate = floor(cn->psrfOccur)/1000./Occur;
			prec = 3;
			if (avRate>100.) prec++;
		} else {
			avRate = 0.0;
			prec = 0;
		}
		*pout<<setprecision(6)<<Occur<<                      //19
		     setprecision(prec)<<','<<avRate<<',';           //20

		// -----------
		Occur = (int)((cn->satsrfOccur-floor(cn->satsrfOccur))*1.E+6);
		if (Occur > 0) {
			avRate = floor(cn->satsrfOccur)/1000./Occur;
			prec = 3;
			if (avRate>100.) prec++;
		} else {
			avRate = 0.0;
			prec = 0;
		}
		*pout<<setprecision(6)<<Occur<<','                   //21
		     << setprecision(prec)<<avRate<<','              //22
		     << setprecision(6)<<cn->satOccur<<','           //23
		     << setprecision(3)<<cn->RechDisch<<','          //24
		     << setprecision(4)<<cn->getAvET()<<','          //25
		     << setprecision(4)<<cn->getAvEvapFract()<<','   //26
		     << setprecision(4)<<cn->getCumTotEvap() <<','   //27
		     << setprecision(4)<<cn->getCumBarEvap()<<','    //28
		     << setprecision(7)<<cn->getCumLHF()<<','        //29
		     << setprecision(7)<<cn->getCumMelt()<<','       //30
		     << setprecision(7)<<cn->getCumSHF()<<','        //31
		     << setprecision(7)<<cn->getCumPHF()<<','        //32
		     << setprecision(7)<<cn->getCumRLin()<<','       //33
		     << setprecision(7)<<cn->getCumRLout()<<','      //34
		     << setprecision(7)<<cn->getCumRSin()<<','       //35
		     << setprecision(7)<<cn->getCumGHF()<<','        //36
		     << setprecision(7)<<cn->getCumUerror()<<','     //37
		     << setprecision(7)<<cn->getCumHrsSun()<<','     //38
		     << setprecision(7)<<cn->getCumHrsSnow()<<','    //39
		     << setprecision(7)<<cn->getPersTimeMax()<<','   //40
		     << setprecision(7)<<cn->getPeakSWE()<<','       //41
		     << setprecision(7)<<cn->getInitPackTime()<<','  //42
		     << setprecision(7)<<cn->getPeakPackTime()<<','  //43
		     << setprecision(7)<<cn->getCumIntSub()<<','     //44
		     << setprecision(7)<<cn->getCumSnSub()<<','      //45
		     << setprecision(7)<<cn->getCumSnEvap()<<','     //46
		     << setprecision(7)<<cn->getCumIntUnl()<<','     //47
		     << setprecision(7)<<cn->getAvThroughFall()<<',' //48
		     << setprecision(7)<<cn->getAvCanFieldCap()<<',' //49
		     << setprecision(7)<<cn->getAvDrainCoeff()<<','  //50
		     << setprecision(7)<<cn->getAvDrainExpPar()<<',' //51
		     << setprecision(7)<<cn->getAvLandUseAlb()<<','  //52
		     << setprecision(7)<<cn->getAvVegHeight()<<','   //53
		     << setprecision(7)<<cn->getAvOptTransmCoeff()<<',' //54
		     << setprecision(7)<<cn->getAvStomRes()<<','     //55
		     << setprecision(7)<<cn->getAvVegFraction()<<',' //56
		     << setprecision(7)<<cn->getAvLeafAI()<<','      //57
		     << setprecision(7)<<cn->getAvEvapThresh()<<','  //58
		     << setprecision(7)<<cn->getAvTransThresh()<<',' //59
		     << setprecision(7)<<cn->getBedrockDepth()<<','  //60
		     << setprecision(7)<<cn->getKs()<<','            //61
		     << setprecision(7)<<cn->getThetaS()<<','        //62
		     << setprecision(7)<<cn->getThetaR()<<','        //63
		     << setprecision(7)<<cn->getPoreSize()<<','      //64
		     << setprecision(7)<<cn->getAirEBubPres()<<','   //65
		     << setprecision(7)<<cn->getDecayF()<<','        //66
		     << setprecision(7)<<cn->getSatAnRatio()<<','    //67
		     << setprecision(7)<<cn->getUnsatAnRatio()<<',' //68
		     << setprecision(7)<<cn->getPorosity()<<','      //69
		     << setprecision(7)<<cn->getVolHeatCond()<<','   //70
		     << setprecision(7)<<cn->getSoilHeatCap()<<','   //71
		     << setprecision(7)<<cn->getSoilID()<<','        //72
		     << setprecision(7)<<cn->getLandUse()<<','       //73
		     << setprecision(7)<<cn->getAvRootZoneDepth()    //74
		     << '\n';

		cn = ni.NextP();
	}

#ifdef PARALLEL_TRIBS
	{
		string localStr = buf.str();
		int localLen = (int)localStr.size();
		int nprocs = tParallel::getNumProcs();
		vector<int> lengths(nprocs), offsets(nprocs);
		MPI_Gather(&localLen, 1, MPI_INT,
		           lengths.data(), 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
		if (tParallel::isMaster()) {
			int total = 0;
			for (int i = 0; i < nprocs; i++) { offsets[i] = total; total += lengths[i]; }
			vector<char> combined(total);
			MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
			            combined.data(), lengths.data(), offsets.data(), MPI_CHAR,
			            MASTER_PROC, MPI_COMM_WORLD);
			string combinedStr(combined.data(), total);
			size_t headerEnd = combinedStr.find('\n');
			string header = combinedStr.substr(0, headerEnd + 1);
			vector<string> rows;
			size_t pos = headerEnd + 1;
			while (pos < (size_t)total) {
				size_t end = combinedStr.find('\n', pos);
				if (end == string::npos) break;
				rows.push_back(combinedStr.substr(pos, end - pos + 1));
				pos = end + 1;
			}
			sort(rows.begin(), rows.end(), [](const string &a, const string &b) {
				return stoi(a) < stoi(b);
			});
			this->CreateAndOpenFileSingle(&intofs, extension);
			intofs << header;
			for (const auto &row : rows) intofs << row;
			intofs.close();
		} else {
			MPI_Gatherv((char*)localStr.data(), localLen, MPI_CHAR,
			            nullptr, nullptr, nullptr, MPI_CHAR,
			            MASTER_PROC, MPI_COMM_WORLD);
		}
	}
#else
	intofs.close();
#endif
	return;
}

/*************************************************************************
**
**  WriteOutletInfo( double time )
**
**  Writes a file containing simulated streamflow and stage values
**  The output format should be readable by ArcInfo & Matlab 
**  
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteOutletInfo( double time )
{
	if (numOutlets > 0) {
		for (int i = 0; i < numOutlets; i++) {
#ifdef PARALLEL_TRIBS
			if ( (Outlets[i] != NULL) && (OutletList[i] > 0) ) {
#else
			if ( Outlets[i] && OutletList[i] < this->g->getNodeList()->getActiveSize()) {
#endif
				// CJC2025: Use fixed format with 4 decimal places for time
				outletinfo[i] << std::fixed << std::setprecision(4) << time << ","
				              << Outlets[i]->getQstrm() << ","
				              << Outlets[i]->getHlevel() << "\n";
			}
		}
	}
	return;
}

/*************************************************************************
**
**  tCOutput::ReadOutletList()
**
**  Opens and Reads the node list from a *.oul file whose structure is:
**
**  Number of Outlet Nodes
**  NodeID1 NodeID2 NodeID3 NodeID4 NodeID5 ...
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::ReadOutletNodeList(char *nodeFileO)
{
	ifstream readOUL(nodeFileO);
	if (!readOUL) {
		Cout<<"\n>>>Outlet Node File "<<nodeFileO<<" not found..."<<endl;
		Cout<<">>>No output for interior nodes will be written..."<<endl<<endl;
		numOutlets = 0;
		return;
	}
	
	readOUL>>numOutlets;
	OutletList  = new int[numOutlets];
	Outlets     = new tSubNode*[numOutlets];
	outletinfo  = new ofstream[numOutlets];
	for (int i = 0; i < numOutlets; i++)
		readOUL>>OutletList[i]; 
	
#ifdef PARALLEL_TRIBS
  // Initialize Outlets to NULL, used to determine local Outlets
  for (int i = 0; i < numOutlets; i++)
    Outlets[i] = NULL;
#endif

	readOUL.close();
	return;
}

/*************************************************************************
**
**  tOutput::CreateAndOpenOutlet()
**
**  Write the header for the *.iout output file.
** 
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::CreateAndOpenOutlet()
{
	if (numOutlets > 0)
	{
		if ( OutletList ) {
			char pixelext[10] = ".qout";
			char nodeNum[10], pixelnode[100];
			char fullName[kMaxNameSize+6];
		
			for (int i = 0; i < numOutlets; i++) {
#ifdef PARALLEL_TRIBS
       // Check if outlet node is on this processor
       if ( (Outlets[i] != NULL) && (OutletList[i] > 0) ) {
#else
				if (OutletList[i] >= 0) { //WR added = for single element case, where
#endif
					snprintf(nodeNum, sizeof(nodeNum),"_%d", OutletList[i]);
					strcpy(pixelnode, nodeNum);
					strcat(pixelnode, pixelext);
				
					strcpy( fullName, outletName );
					strcat( fullName, pixelnode );
				
					outletinfo[i].open( fullName );
				
					if ( !outletinfo[i].good() )
						cerr<<"File "<<fullName<<"can not be created.";
					
					// Write Header
					outletinfo[i] << "Time_hr,Qstrm_m3_s,Hlev_m\n";
				}
			}
		}
	}
	return;
}

/*************************************************************************
**
**  tOutput::SetInteriorOutlet()
**
**  Initializes pointers to basin interior outlets 
** 
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::SetInteriorOutlet()
{
	tSubNode * cnn;
	tMeshListIter<tSubNode> niter( this->g->getNodeList() );
	if ( OutletList ) {
		for (int i=0; i < numOutlets; i++) {
			if (OutletList[i] >= 0) { //WR added in = for single element case where node ID = 0
#ifdef PARALLEL_TRIBS
           // Each processor creates/writes outlet files for its points
           for (cnn=niter.FirstP(); niter.IsActive(); cnn=niter.NextP() ) {
#else
				for (cnn=niter.FirstP(); !(niter.AtEnd()); cnn=niter.NextP() ) {
#endif
					if (cnn->getID() == OutletList[i]) {
						Outlets[i] = cnn;   //Defining node of interest 
						Cout<<"\nOutlet of Interest ID: \t"<<OutletList[i]
							<<" has been set up..."<<endl<<flush;
					}
				}
			}
		}
	}
	return;
}

//=========================================================================
//
//
//                       End of tOutput.cpp
//
//
//=========================================================================
