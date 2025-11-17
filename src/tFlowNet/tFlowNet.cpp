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
**  tFlowNet.cpp: Functions for class tFlowNet (see tFlowNet.h) based on
**             CHILD and RIBS routines for Stream Net Routing
**
***************************************************************************/

#include "src/tFlowNet/tFlowNet.h"
#include "src/Headers/TemplDefinitions.h"
#include "src/Headers/globalIO.h"

#ifdef PARALLEL_TRIBS
#include "src/tParallel/tParallel.h"
#endif

//=========================================================================
//
//
//                  Section 1: tFlowNet Constructors/Destructors
//
//
//=========================================================================

tFlowNet::tFlowNet() 
{
	gridPtr = 0;
}

tFlowNet::tFlowNet(SimulationControl *simCtrPtr, tMesh<tCNode> *gridRef, 
				   tInputFile &infile, tRunTimer *timptr) 
{
	gridPtr = gridRef;
	timer   = timptr;
	simCtrl = simCtrPtr;
	BasArea = 0.0;
	
	SetFlowVariables( infile );
    
	// If the mesh was created by the MeshBuilder read FlowNet info from file
	int option = infile.ReadItem(option, "OPTMESHINPUT");
	if (option != 9) {

		Cout <<"\nCalculating slopes..."<< endl;
		CalcSlopes();
	
		Cout <<"\nInitializing flow directions..."<<endl;
		InitFlowDirs();
	
		Cout <<"\nComputing flow directions..."<<endl;
		FlowDirs();
	
		Cout <<"\nCorrecting sinks..."<<endl;
		FillLakes();
	
		Cout <<"\nSorting nodes by network order..."<<endl;
		SortNodesByNetOrder();
	
		Cout <<"\nSetting basin outlet and stream reaches..."<<endl;
		SetBasinOutlet();
		WeightedShortestPath();  
		SortStreamNodes();
		DrainAreaVoronoi();
		DeriveCurvature();
		DeriveStreamReaches( infile );
	
		Cout <<"\nInitialize velocities and time..."<<endl;
		setTravelVelocity(0.0);  
		initializeTravelTime();
	
		Cout <<"\nChecking drainage widths..."<<endl;
		CheckVDrainageWidths();
	
		Cout <<"\nSet reach numbers..."<<endl;
		SetReachInformation();
  }

  else {
		Cout <<"\nRead MeshBuilder flownet information..."<<endl;
		ReadFlowNetFromMeshBuilder();
	}
	
	Cout<<"\nHillslope velocity: \t\t"<<hillvel<<" m/sec"<<endl;
	Cout<<"Stream velocity: \t\t"<<streamvel<<" m/sec"<<endl;
	
	// Initialize Hydrograph class
	res = new tFlowResults( simCtrl, infile, timer, MaxTravel() );

	percolationOption = infile.ReadItem(percolationOption, "OPTPERCOLATION"); // ASM percolation option

}

tFlowNet::~tFlowNet() 
{
	gridPtr = NULL;
	timer   = NULL;
	delete res;
	Cout<<"tFlowNet Object has been destroyed..."<<endl<<flush;
}

/*****************************************************************************
**  
**  ReadFlowNetFromMeshBuilder()
**  
**  MeshBuilder already calculated and wrote out flow network
**
*****************************************************************************/

void tFlowNet::ReadFlowNetFromMeshBuilder()
{
	fstream flowStr("flow.meshb", ios::in | ios::binary);

	BinaryRead(flowStr, hillvel);
	BinaryRead(flowStr, streamvel);
	BinaryRead(flowStr, velratio);
	BinaryRead(flowStr, velcoef);
	BinaryRead(flowStr, flowexp);
	BinaryRead(flowStr, baseflow);
	BinaryRead(flowStr, flowout);
	BinaryRead(flowStr, maxttime);
	BinaryRead(flowStr, dist_hill_max);
	BinaryRead(flowStr, dist_stream_max);
	BinaryRead(flowStr, BasArea);

	// Get the actual node pointers from the save NodeTable lookup in mesh
	int id, count, size;
	tCNode* curnode;

	BinaryRead(flowStr, id);
	OutletNode = (tCNode*) gridPtr->getNodeFromID(id);

	BinaryRead(flowStr, size);
	for (int i = 0; i < size; i++) {
		BinaryRead(flowStr, id);
		curnode = (tCNode*) gridPtr->getNodeFromID(id);
		if (curnode != 0)
			HeadsLst.insertAtBack(curnode);
		else
			cout << "tFlowNet::ReadFlowNetFromMeshBuilder: Missing Head" << endl;
	}

	BinaryRead(flowStr, size);
	for (int i = 0; i < size; i++) {
		BinaryRead(flowStr, id);
		curnode = (tCNode*) gridPtr->getNodeFromID(id);
		if (curnode != 0)
			NodesLstH.insertAtBack(curnode);
		else
			cout << "tFlowNet::ReadFlowNetFromMeshBuilder: Missing Head" << endl;
	}

	BinaryRead(flowStr, size);
	for (int i = 0; i < size; i++) {
		BinaryRead(flowStr, id);
		curnode = (tCNode*) gridPtr->getNodeFromID(id);
		if (curnode != 0)
			NodesLstO.insertAtBack(curnode);
		else
			cout << "tFlowNet::ReadFlowNetFromMeshBuilder: Missing Outlet" << endl;
	}

	BinaryRead(flowStr, size);
	for (int i = 0; i < size; i++) {
		BinaryRead(flowStr, count);
		NNodes.insertAtBack(count);
	}
	flowStr.close();

	// Now we can delete the NodeTable
	gridPtr->deleteNodeTable();
}

//=========================================================================
//
//
//                  Section 2: tFlowNet Slopes/Direction Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  SetFlowVariables()
**  
**  Reads from InputFile and assigns basic varibales used in tFlowNet
**  As opposed to work conducted on mesh by other functions in constructor
**
*****************************************************************************/
void tFlowNet::SetFlowVariables(tInputFile &infile)
{
	velratio = infile.ReadItem(velratio, "VELOCITYRATIO");
	baseflow = infile.ReadItem(baseflow, "BASEFLOW");
	velcoef  = infile.ReadItem(velcoef,  "VELOCITYCOEF");
	flowexp  = infile.ReadItem(flowexp,  "FLOWEXP");
	timespan = timer->RemainingTime(0.0);
	dOtp     = timer->getOutputInterval();
	
	flowout   = 0.0;
	streamvel = 0.0;
	maxttimeInitial = 0.0; // added by Ara Ko 2017
	return;
}

/*****************************************************************************
**  
**  tFlowNet::MaxTravel()
**  
**  Calculates dimension of the output array
** 
*****************************************************************************/
int tFlowNet::MaxTravel() 
{
	flowboxes = static_cast<int>(ceil( (maxttime + timespan)/dOtp ));
	Cout<<"Hydrograph array size: \t\t"<<flowboxes<<endl; 
	if (maxttimeInitial == 0.0)	// added by Ara Ko 2017	
	   maxttimeInitial=maxttime; // added by Ara Ko 2017	
	return flowboxes; 
}

/*****************************************************************************
**  
**  tFlowNet::CalcSlopes()
**  
**  Calculates slope for each Edge  
**
*****************************************************************************/
void tFlowNet::CalcSlopes()
{
	tEdge *curedg;
	tMeshListIter<tEdge> i( gridPtr->getEdgeList() );
	double slp;
	
	// Loop through each pair of edges on the list
	for ( curedg = i.FirstP(); !( i.AtEnd() ); curedg = i.NextP() )
	{
		// Compute the slope and assign it to the current edge
		slp = ( curedg->getOrgZ() - curedg->getDestZ() )
		/ curedg->getLength();
		curedg->setSlope( slp );
		
		// Advance to the edge's complement, and assign it -slp
		curedg = i.NextP();
		curedg->setSlope( -slp );    
	}
	return;	
}

/*****************************************************************************
**  
**  tFlowNet::InitFlowDirs()
**  
**  Defines flow direction of flow for each Node  
**
*****************************************************************************/
void tFlowNet::InitFlowDirs() {
	
	tMeshListIter<tCNode> niter( gridPtr->getNodeList() );
	tCNode * curnode;
	tEdge  * flowedg;
	tListNode< tCNode > * nodeToMove;
	int ctr;
	int kMaxSpokes = 100; 
	
	// For every active (non-boundary) node, initialize it to flow to a
	// _non-boundary_ node (i.e., to a hillslope or stream node )
	
	curnode = niter.FirstP();
	while ( niter.IsActive() ) {
		// Start with the node's default edge
		flowedg = curnode->getEdg();
		
		// As long as the current edge is a no-flow edge, 
		// advance to the next one, counter-clockwise
		ctr = 0;
		while ( (!flowedg->FlowAllowed() || 
				 flowedg->getVEdgLen() < THRESH) && 
				niter.IsActive() ) {
			
			flowedg = flowedg->getCCWEdg();
			ctr++;
			if (ctr > kMaxSpokes) { //Make sure to prevent endless loops
				Cout << "\nError in InitFlowDirs(): Node " << curnode->getID()
				<< " surrounded by closed boundary nodes ..." << endl;
				Cout<<"\nNode boundary will be changed to kClosedBoundary..."<<endl;
				Cout<<curnode->getID()
					<<"\t"<<curnode->getX()<<"\t"<<curnode->getY()
					<<"\t"<<curnode->getZ()<<"\t"<<curnode->getBoundaryFlag()
					<<endl<<endl<<flush; 
				
				// Changing the boundary code and going to the next node
				curnode->setBoundaryFlag( kClosedBoundary );
				nodeToMove = niter.NodePtr();
				curnode = niter.NextP();
				flowedg = curnode->getEdg();
				ctr = 0;
				
				// Now, move that node to the end of the list
				gridPtr->getNodeList()->moveToBack( nodeToMove );
			}
		}
		curnode->setFlowEdg( flowedg ); //Flowedge for current node is set
		curnode = niter.NextP();
	}
	return;
}

/*****************************************************************************
**  
**  tFlowNet::FlowDirs()
**  
**  Defines flow direction of flow for each Node  
**
*****************************************************************************/
void tFlowNet::FlowDirs() 
{
	int ctr;
	double slp;        // steepest slope found so far
	int kLargeNegative  = -1000;
	int kMaxSpokes = 100;
	
	tMeshListIter<tCNode> i( gridPtr->getNodeList() );  
	tCNode *curnode;  // ptr to the current node
	tEdge * firstedg; // ptr to first edge
	tEdge * curedg;   // pointer to current edge
	tEdge * nbredg;   // steepest neighbouring edge so far
	
	Cout.setf( ios::fixed, ios::floatfield);
	
	// Find the connected edge with the steepest slope 
	curnode = i.FirstP();
	while ( i.IsActive() )  // DO for each non-boundary (active) node
	{
		firstedg = curnode->getFlowEdg();
		slp    = firstedg->getSlope();
		nbredg = firstedg;
		curedg = firstedg->getCCWEdg();
		ctr = 0;         
		
		// Check each of the various spokes, stopping when
		// we've gotten back to the beginning
		while ( curedg != firstedg )  {
			if (curedg->getSlope() > slp && 
				curedg->getDestinationPtrNC()->getBoundaryFlag() != kClosedBoundary &&
				curedg->getVEdgLen() > THRESH) {
				slp    = curedg->getSlope();
				nbredg = curedg;
			}
			curedg = curedg->getCCWEdg();
			ctr++;
			
			if (ctr > kMaxSpokes) { // Make sure to prevent endless loops
				Cout<<"\nError in FlowDirs(): Node "<<curnode->getID()
				<< " surrounded by closed boundary nodes" << endl;
			}
			/*
				if (curnode->getID() > 0 && curnode->getBoundaryFlag() == kStream) {
					Cout<<curedg->getDestinationPtrNC()->getID()
					<<"\t"<<curedg->getDestinationPtrNC()->getX()
					<<"\t"<<curedg->getDestinationPtrNC()->getY()
					<<"\t"<<curedg->getDestinationPtrNC()->getZ()
					<<"\t"<<curedg->getSlope()
					<<"\t"<<curedg->getDestinationPtrNC()->getBoundaryFlag()
					<<endl<<flush; }
			 */
		}
		curnode->setFlowEdg( nbredg ); // The steepest flowedge
		
		if ((slp <= 0.0) && (curnode->getBoundaryFlag() == kStream )) {
			//TellAboutNode(curnode);
		}
		
		// This is needed to check for pits
		if ( (slp > 0.0) && (curnode->getBoundaryFlag() != kClosedBoundary))
			curnode->setFloodStatus( kNotFlooded );
		else {
			curnode->setFloodStatus( kSink );
		}
		curnode = i.NextP();    
	}
	return;
}

//=========================================================================
//
//
//                  Section 3: tFlowNet Sorting Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  tFlowNet::SortNodesByNetOrder()
**  
**  This function sorts the list of nodes according to their order in the
**  network (upstream to downstream). The sorting algorithm is based on the 
**  "cascade" algorithm of Braun and Sambridge. 
**     The single-direction sorting algorithm works by initially assigning a 
**  tracer (like a packet of water) to each node. At each iteration, a tracer 
**  from each node issent downstream. Any nodes that have zero tracers left 
**  are moved to the bottom of the list (a FIFO stack), so that for example 
**  the very first node moved will be the first node on the list when the 
**  sorting is completed. The process continues until no unsorted nodes remain.
**    The multi-flow option was added to allow for multiple flow directions
**  and kinematic-wave routing. The algorithm is slightly different. At
**  each pass, the unsorted nodes are "de-flagged" by setting their
**  tracer variables to zero. Then each unsorted node flags _all_ of the
**  adjacent nodes that are downhill by setting their tracer to 1 (this
**  is accomplished through a call to ActivateSortTracer). Finally, any
**  unflagged nodes are moved to the back of the list, and the process
**  is repeated until all nodes have been sorted.
**
**  TODO: it is possible that the "flagging" method, as opposed to the
**  "tracer movement" method, is more efficient even for single-flow
**  routing. The added cost lies in unflagging the unsorted nodes at
**  each pass. However, the gain comes from the fact that you never have
**  multiple "tracers" at a node that need to be removed one by one.
**  The two methods should be tested and compared.
**
*****************************************************************************/
void tFlowNet::SortNodesByNetOrder()
{
	int nThisPass;                   // Number moved in current iteration
	int i;
	int done=0;
	
	tCNode * cn;
	tMeshList<tCNode> *nodeList = gridPtr->getNodeList();
	tMeshListIter<tCNode> listIter( nodeList );
	int nUnsortedNodes = nodeList->getActiveSize(); // Number not yet sorted
	
	// Assign initial tracers: use "qs" field which 
	// contains garbage at this stage
	for ( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() ) 
		cn->ActivateSortTracer();
	
	// Iterate: move tracers downstream and sort
	// until no nodes with tracers are left 
	do {
		// Send tracers downstream
		cn = listIter.FirstP();
		for ( i=1; i <= nUnsortedNodes; i++ ) {
			assert( cn != 0 );
			cn->MoveSortTracerDownstream();
			cn = listIter.NextP();
		}
		
		// Scan for nodes with no tracers, and move to the bottom of the list
		tListNode< tCNode > * nodeToMove;
		nThisPass = 0;
		done = TRUE;
		
		cn = listIter.FirstP();
		for ( i=1; i<=nUnsortedNodes; i++ ) {
			if ( cn->NoMoreTracers() ) {  // If no tracers, move to bottom of list
				nodeToMove = listIter.NodePtr();
				cn = listIter.NextP();
				nodeList->moveToActiveBack( nodeToMove );
				nThisPass++;
			}
			else {
				cn = listIter.NextP();
				done = FALSE;
			}
		}
		
		nUnsortedNodes -= nThisPass;
	} while ( !done );
	
	// Changed to make edge IDs consistent
	tEdge     * ce;
	tTriangle * ct;  
	tMeshListIter<tCNode>   niter( gridPtr->getNodeList() );
	tMeshListIter<tEdge>    eiter( gridPtr->getEdgeList() );
	tListIter<tTriangle>    titer( gridPtr->getTriList() );
	
	int nnodes = gridPtr->getNodeList()->getSize();
	int nedges = gridPtr->getEdgeList()->getSize();
	int ntri   = gridPtr->getTriList()->getSize();
	int id;
	
	Cout <<"\nRenumbering nodes, edges, triangles..."<<endl;
	for ( cn=niter.FirstP(), id=0; id<nnodes; cn=niter.NextP(), id++ )
        cn->setID( id );
	for ( ce=eiter.FirstP(), id=0; id<nedges; ce=eiter.NextP(), id++ )
        ce->setID( id );
	for ( ct=titer.FirstP(), id=0; id<ntri; ct=titer.NextP(), id++ )
        ct->setID( id );
	return;
}

//=========================================================================
//
//
//                  Section 4: tFlowNet Velocity Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  tFlowNet::setTravelVelocity()
**  
**  Sets travel velocities depending on the discharge at the oulet
**
*****************************************************************************/
void tFlowNet::setTravelVelocity(double curr_discharge) 
{
	// Get the node velocities first 
	if ( !curr_discharge )  { 
		// If it's zero, a baseflow value from the input
		// file is used to define the stream velocity     
		if ( !baseflow && flowexp ) {
			Cout<<"\nWarning: Baseflow is zero and thus the lower "
			<<"limit of stream velocity is undefined --> Set to 0.001"<<endl;
			baseflow = 0.001; 
		}
		flowout = baseflow; //velocity will correspond to baseflow
	}
	else 
		flowout = curr_discharge;
	
	streamvel = velcoef;
	hillvel   = streamvel/velratio;
	
	return;
}

/*****************************************************************************
** 
**  tFlowNet::initializeTravelTime()
**  
**  Needed to initialize some of the class members  
**
*****************************************************************************/
void tFlowNet::initializeTravelTime() 
{
	tCNode *cn;
	tCNode *ctimer;
	tEdge  *ce;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );
	tMeshListIter<tEdge>  edgIter( gridPtr->getEdgeList() );
	double hill, stream, tt;
	
	// Loop through the nodes and set velocity
	// Define maximum travel time        
	BasArea = 0.0;
	maxttime = 0.0;     		// SECONDS
	dist_hill_max = 0.0;   	// METERS
	dist_stream_max = 0.0; 	// METERS
	
	for ( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() ) {
		
		BasArea += cn->getVArea();

		//Initialize with distance from centroid to node
		cn->getCentroidX();
		cn->getCentroidY();
		hill = sqrt(pow(cn->getX()-cn->getCentroidX(),2.0)+pow(cn->getY()-cn->getCentroidY(),2.0));
		
		stream = tt = 0.0;
		
		ctimer = cn;
		assert(ctimer != 0);
		
		// If it is an hillslope, it flows at a hillslope velocity
		if (ctimer->getBoundaryFlag() == kNonBoundary) {
			while (ctimer->getBoundaryFlag() == kNonBoundary) {
				ce = ctimer->getFlowEdg(); //Get the steepest edge
				hill += ce->getLength();
				ctimer = ctimer->getDownstrmNbr();
			}
			// Set a stream node to which it contributes flow
			assert(ctimer != 0);
			cn->setStreamNode(ctimer);
		}
		
		// If it is a streamnode, it flows at a stream velocity
		// Go downstream to the outlet
		while ( ctimer->getDownstrmNbr() ) {
			ce = ctimer->getFlowEdg();  // Get the steepest edge
			if (ctimer->getBoundaryFlag() == kStream || 
				ctimer->getBoundaryFlag() == kOpenBoundary)
				stream += ce->getLength();            // METERS
			ctimer = ctimer->getDownstrmNbr();
		}
		
		cn->setHillPath(hill);     		// Set the HILLSLOPE path for node
		cn->setStreamPath(stream); 		// Set the STREAM path for node
		
		tt = hill/hillvel + stream/streamvel; // SECONDS 
		tt /= 3600.;                          // HOURS  
		
		// Set travel time: travel time is only suitable for HYDROLOGIC routing 
		// The KINEMATIC scheme uses independent calculation of travel time
		
		cn->setTTime( tt );        		// Set Travel Time in HOURS
		
		if (tt > maxttime) {  //MAX values of paths with maxttime
			maxttime = tt;
			dist_hill_max = hill;       	// METERS
			dist_stream_max = stream;   	// METERS
		}
	}
	Cout.setf( ios::fixed, ios::floatfield);
	Cout<<endl;
	Cout<<"\nFlow Characteristics: "<<endl;
	Cout<<"\nThe  Total  Basin  Area is   \t"<<BasArea<<" M^2"<<endl;
	Cout<<"Maximum travel time: \t\t"<<maxttime<<" hours"<<endl;
	Cout<<"Maximum hillslope path: \t"<<dist_hill_max<<" meters"<<endl;
	Cout<<"Maximum stream path: \t\t"<<dist_stream_max<<" meters"<<endl;
	Cout<<endl;
	
	return;
}

/*****************************************************************************
**  
**  initializeTravelTimeOnly()
**  
**  Difference with the preceeding function is that this one does  
**  not compute distances and deals only with the travel times
** 
*****************************************************************************/
void tFlowNet::initializeTravelTimeOnly() 
{
	tCNode *cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );
	double hill, stream, tt;
	
	// Loop through the nodes and set velocity
	// Define maximum travel time 
	maxttime = 0.0;        // HOURS
	dist_hill_max = 0.0;   // METERS
	dist_stream_max = 0.0; // METERS
	
	for ( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
	{
		hill   = cn->getHillPath();     // Get the HILLSLOPE path for node
		stream = cn->getStreamPath();   // Get the STREAM path for node
		tt = hill/hillvel + stream/streamvel;   //SECONDS
		tt /= 3600.0;                           //HOURS
		cn->setTTime( tt );             // Set Travel Time in HOURS
		
		if (tt > maxttime) {            // MAX values of paths with maxttime
			maxttime = tt;
			dist_hill_max = hill;       // METERS
			dist_stream_max = stream;   // METERS
		}
	}
	
	Cout<<"\nMaximum Travel Time: \t\t"<<maxttime<<" hours"<<endl;
	Cout<<"Maximum Hillslope path: \t"<<dist_hill_max<<" meters"<<endl;
	Cout<<"Maximum Stream path: \t\t"<<dist_stream_max<<" meters"<<endl;
	
	return;
}

/*****************************************************************************
**  
**  setMaxTravelTime()
**  
**  Needed to set Maximum Travel Time
**  Assumed that the travel velocities have been set by now!  
**
*****************************************************************************/
void tFlowNet::setMaxTravelTime() 
{
	// Update the limit time of simulation added Ara Ko in 2017
	//res->limit = flowboxes = static_cast<int>(ceil( (timer->RemainingTime(0.0) + maxttimeInitial)/dOtp)); 	//WR Debug restart, cause error in restart.
	
	maxttime = dist_hill_max/hillvel + dist_stream_max/streamvel; //SECONDS
	res->iimax = timer->getResStep(maxttime/(3600.0));
	
	// Set MAX index in hydrograph
	if (res->iimax > res->limit || res->iimax < 0) {
		Cout<<"\nError: setTravelTime(): iimax > limit,  iimax = "
        <<res->iimax<<"; limit = "<<res->limit
        <<" Exiting program..."<<endl<<flush;
		exit(1);
	}
	return;
}

//=========================================================================
//
//
//                  Section 5: tFlowNet Discharge and Flow Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  getCurrDischarge()
**  
**  Needed to set Travel Velocities
**
*****************************************************************************/
double tFlowNet::getCurrDischarge(int ihour)
{
	double curr_discharge;
	
	if (ihour > 0)
		curr_discharge = res->get_discharge(ihour);
	else
		curr_discharge = res->get_discharge(0);
	
	return curr_discharge;
}

/*****************************************************************************
** 
**  SurfaceFlow()
**  
**  First, defines current discharge at the basin outlet needed 
**         in stream velocitiy calculation
**  Second, define routine velocities at this moment
**  Third, setMaxTravelTime
**  Fourth, loop through the list of active nodes,
**         - if the model is non-linear, re-set travel time using
**           defined velocities and lengths of hillslope and stream path
**         - get runoff volume and store in appropriate array box  
**         - do the same with runoff types
**
*****************************************************************************/
void tFlowNet::SurfaceFlow()
{
	tCNode *cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );
	double ttime;           // Travel time for a current node, SECONDS
	double vRunoff = 0.0;   // Runoff volume, m^3
	double Area = 0.0;      // Voronoi cell area
	double AreaF;

	// SKY2008Snow from AJR2007
	double val;

	if (simCtrl->Verbose_label == 'Y') {
		Cout<<"\t->Surface flow simulation...\n"<<endl<<flush;
	}
	
	int ihour = (int)timer->getResStep(0.0);
	
	setTravelVelocity( getCurrDischarge(ihour-1) );
	
	setMaxTravelTime(); 
	
	// Runoff Storage in tFlowResults
	for ( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() ) {
		
		Area = cn->getVArea();                    // M^2
		AreaF = Area/BasArea;
		
		if (cn->getSrf() > 0.0) {
			
			if (flowexp > 0.0) { 
				cn->setTTime((cn->getHillPath()/hillvel + 
							  cn->getStreamPath()/streamvel)/3600.0);
			}
			
			ttime = cn->getTTime();                        // Travel time in HOURS
			vRunoff = cn->getSrf()*Area/1000.0;            // Srf, MM to M^3
			res->store_volume( ttime, vRunoff );           // PSEUDO-Routine function
			
			if ( cn->getHsrf() ) {
				vRunoff = cn->getHsrf()*Area/1000.0;         // MM to M^3
				res->store_volume_Type( ttime, vRunoff, 1 ); // PSEUDO-Routing function
			}
			if ( cn->getSbsrf() ) {
				vRunoff = cn->getSbsrf()*Area/1000.0;   
				res->store_volume_Type( ttime, vRunoff, 2 ); 
			}
			if ( cn->getPsrf() ) {
				vRunoff = cn->getPsrf()*Area/1000.0;   
				res->store_volume_Type( ttime, vRunoff, 3 ); 
			}
			if ( cn->getSatsrf() ) {
				vRunoff = cn->getSatsrf()*Area/1000.0;  
				res->store_volume_Type( ttime, vRunoff, 4 ); 
			}
		}
		 // Fix for converting MDGW sloped depth to vertical depth CJC2025
		// Get the slope correction factor for the current node.
		tEdge *flowEdge = cn->getFlowEdg();
		double cos_slope = cos(atan(flowEdge->getSlope()));
		if (cos_slope < 1E-9) cos_slope = 1.E-9;
		// Convert Nwt to a vertical depth.
		double nwt_vertical = cn->getNwtNew() / cos_slope;

		// Rainfall and Saturation Storage in tFlowResults
		res->store_rain(0.0, AreaF*cn->getRain());
		
		// Min/Max Rainfall and Fractional Rainfall 
		res->store_maxminrain(0.0, cn->getRain(), 0);
		if (cn->getRain() > 0.0)
			res->store_maxminrain(0.0, AreaF, 1);
		
		// Mean Soil Moisture (top 100 mm, Root and unsaturated zone) 
		// and Saturated Area
		res->store_saturation(0.0, AreaF*cn->getSoilMoistureSC(), 0);
		res->store_saturation(0.0, AreaF*cn->getRootMoistureSC(), 1);
		res->store_saturation(0.0, AreaF*cn->getSoilMoistureUNSC(), 2);
		if (cn->getSoilMoistureSC() >= 1)
			res->store_saturation(0.0, AreaF, 3);
		// Mean groundwater level
		res->store_saturation(0.0, AreaF*nwt_vertical, 4);
		//Mean Evapotranspiration
		res->store_saturation(0.0, AreaF*cn->getEvapoTrans(), 5);

		// SKY2008Snow from AJR2007
		//Mean SWE
		res->store_saturation(0.0, AreaF*(cn->getIceWE() + cn->getLiqWE()),6);//added by AJR 2007 @ NMT
		//Mean melt
		res->store_saturation(0.0, AreaF*cn->getLiqRouted(),7);//added by AJR 2007 @ NMT
		//Mean ST
		res->store_saturation(0.0, AreaF*cn->getSnTempC(),8);//added by AJR 2007 @ NMT
		//Mean DU
		res->store_saturation(0.0, AreaF*cn->getDU(),9);//added by AJR 2007 @ NMT
		//Mean sLHF
		res->store_saturation(0.0, AreaF*cn->getSnLHF(), 10);//added by AJR 2007 @ NMT
     		//Mean sSHF
		res->store_saturation(0.0, AreaF*cn->getSnSHF(), 11);//added by AJR 2007 @ NMT
     		//Mean sGHF
		res->store_saturation(0.0, AreaF*cn->getSnGHF(), 12);//added by AJR 2007 @ NMT
     		//Mean sPHF
		res->store_saturation(0.0, AreaF*cn->getSnPHF(), 13);//added by AJR 2007 @ NMT
     		//Mean sRLi
		res->store_saturation(0.0, AreaF*cn->getSnRLin(), 14);//added by AJR 2007 @ NMT
     		//Mean sRLo
		res->store_saturation(0.0, AreaF*cn->getSnRLout(), 15);//added by AJR 2007 @ NMT
     		//Mean sRSi
		res->store_saturation(0.0, AreaF*cn->getSnRSin(), 16);//added by AJR 2007 @ NMT
     		//Mean intSWE
		res->store_saturation(0.0, AreaF*cn->getIntSWE(),17);//added by AJR 2007 @ NMT
		//Mean intSub
		res->store_saturation(0.0, AreaF*cn->getIntSub(),18);//added by AJR 2007 @ NMT
		//Mean intUnl
		res->store_saturation(0.0, AreaF*cn->getIntSnUnload(),19);//added by AJR 2007 @ NMT
     		//SCA
		if ((cn->getIceWE() + cn->getLiqWE()) > 0.0)
			val = 1.0;
		else
			val = 0.0;
		res->store_saturation(0.0, AreaF*val,20);//added by AJR 2007 @ NMT -- SCA
			//Mean SnSub
		res->store_saturation(0.0, AreaF*cn->getSnSub(),21);// Calculated mean snowpack sublimation CJC2020 
			//Mean SnSub
		res->store_saturation(0.0, AreaF*cn->getSnEvap(),22);// Calculated mean snowpack sublimation CJC2020
		//Mean Qunsat
		res->store_saturation(0.0, AreaF*(cn->getQpout() - cn->getQpin()) * 1.E-6 / cn->getVArea(),24); // CJC 2025

        //ASM Percolation option
        if (percolationOption != 0)
            res->store_saturation(0.0, cn->getChannelPerc(),23);

	}
	return;
}

//=========================================================================
//
//
//                  Section 5: tFlowNet FillLakes Function
//
//
//=========================================================================

/*****************************************************************************
**
**  tFlowNet::FillLakes 
**
**  Finds drainage for closed depressions. The algorithm assumes
**  that sinks (nodes that are lower than any of their neighbors)
**  have already been identified during the flow directions
**  procedure. For each sink, the algorithm creates a list of
**  nodes in the current lake, which initially is just the sink
**  itself. It then iteratively looks for the lowest node on the
**  perimeter of the current lake. That node is checked to see 
**  whether it can be an outlet, meaning that one of its
**  neighbors is both lower than itself and is not already
**  flooded (or is an open boundary). If the low node is not an 
**  outlet, it is added to the current lake and the process 
**  is repeated. If it is an outlet, then all of the nodes on the 
**  current-lake list are identified draining it. The list is then
**  cleared, and the next sink is processed. If during the search
**  for the lowest node on the perimeter a flooded node is found
**  that isn't already part of the current lake (i.e., it was
**  flagged as a lake node when a previous sink was processed),
**  then it is simply added to the current-lake list --- in other
**  words, the "new" lake absorbs any "old" ones that are encountered.
**
**  Once an outlet has been found, flow directions for nodes in the
**  lake are resolved in order to create a contiguous path through
**  the lake.
**
**    Calls: FindLakeNodeOutlet
**    Called by: MakeFlow
**    Modifies:  flow direction and flood status flag of affected nodes
**
*****************************************************************************/
void tFlowNet::FillLakes()
{
	tCNode *cn, *cn2,        // Node on list: if a sink, then process
	*thenode,           // Node on lake perimeter
	*lowestNode,        // Lowest node on perimeter found so far
	*cln,               // current lake node
	*node;              // placeholder
	tMeshListIter< tCNode > nodIter( gridPtr->getNodeList() ); // node iterator
	tPtrList< tCNode > lakeList;                 // List of flooded nodes
	tPtrListIter< tCNode > lakeIter( lakeList ); // Iterator for lake list
	tEdge *ce;              
	double lowestElev;      // Lowest elevation found so far on lake perimeter
	int done;               // Flag indicating whether outlet has been found
	
	// Check each active node to see whether it is a sink
	for ( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() ) {
		
		if ( cn->getFloodStatus() == kSink )
		{
			// Create a new lake-list, initially containing just the sink node.
			lakeList.insertAtBack( cn );
			cn->setFloodStatus( kCurrentLake );
			
			// Iteratively search for an outlet along the perimeter of the lake 
			done = FALSE;
			do
			{
				lowestNode = lakeIter.FirstP();
				lowestElev = kVeryHigh; // Initialize lowest elev to very high val.
				
				// === Check the neighbors of every node on the lake-list ===
				for ( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
					  cln = lakeIter.NextP() )
				{
					ce = cln->getEdg();
					do
					{
						thenode = (tCNode *) ce->getDestinationPtrNC();
						// Is it a potential outlet (i.e. not flooded 
						// and not a basin boundary)?
						if (thenode->getBoundaryFlag() == kOpenBoundary ||
							(thenode->getFloodStatus()  == kNotFlooded &&
							 thenode->getBoundaryFlag() != kClosedBoundary))
						{ 
							// If it's the basin outlet OR
							// Is it lower than the lowest found so far?
							if (thenode->getBoundaryFlag() == kOpenBoundary || 
								thenode->getZ() < lowestElev) {
								lowestNode = thenode;
								lowestElev = thenode->getZ();
							}
						}
						// If it's a previous lake node or a sink, add it to the list
						// Added check for boundary flag
						else if ((thenode->getFloodStatus()  == kFlooded ||
								  thenode->getFloodStatus()  == kSink) && 
								 (thenode->getBoundaryFlag() == kStream ||
								  thenode->getBoundaryFlag() == kNonBoundary))  {
							lakeList.insertAtBack( thenode );
							thenode->setFloodStatus( kCurrentLake );
						}
					}  while ( ( ce=ce->getCCWEdg() ) != cln->getEdg() ); // END spokes
					
				} // End lakeList
				
				// Now we've found the lowest point on the perimeter. Now test
				// to see whether it's an outlet
				
				// 1.) If it's an open boundary, it's an outlet
				//     => Get out of the loop
				if ( lowestNode->getBoundaryFlag() == kOpenBoundary ) 
					done = TRUE;
				
				else  { 
					// 2.) it's also an outlet if it can drain to a "dry" location.
					//     Can lowestNode drain to a non-flooded location?
					if ( FindLakeNodeOutlet( lowestNode ) )
						done = TRUE;
					// no, it can't, so add it to the list and continue:
					else
					{
						if (lowestNode->getBoundaryFlag() == kStream ||
							lowestNode->getBoundaryFlag() == kNonBoundary) {                      
							lakeList.insertAtBack( lowestNode );
							lowestNode->setFloodStatus( kCurrentLake );
						}
					}
				}
				
				if ( lakeList.getSize() > gridPtr->getNodeList()->getActiveSize() )
				{
					Cout<<"\nLAKE LIST SIZE: "<<lakeList.getSize()<<endl;
					Cout<<"\nACTIVE SIZE: "<<gridPtr->getNodeList()->getActiveSize()
						<<endl<<endl;
					cln  = lakeIter.FirstP();
					node = lakeIter.NextP();
					Cout<<"\n----> THE FIRST ONE IN THE LAKE LIST: "<<endl;
					TellAboutNode(cln);
					Cout<<"\n----> THE SECOND ONE IN THE LAKE LIST: "<<endl;
					TellAboutNode(node);
				}
				assert(lakeList.getSize() <= gridPtr->getNodeList()->getActiveSize());
				
         } while ( !done );
			
			// Now we've found an outlet for the current lake.
			// This next bit of code assigns a flowsTo for each node so there's
			// a complete flow path toward the lake's outlet. This isn't strictly
			// necessary --- the nodes could all point directly to the outlet,
			// skipping anything in between --- but it prevents potential problems
			// in ordering the list by network order. This also works by pointing
			// each node toward the first neighboring node they happen to find
			// that has been flagged as having its flow direction resolved. 
			// Initially, the low node is thus flagged, and the algorithm repeats
			// until all the lake nodes are flagged as having a flow direction.
			// The algorithm isn't unique---there are many paths that could be
			// taken; this simply finds the most convenient one.
			
			lowestNode->setFloodStatus( kOutletFlag );
			
			// Test for error in mesh: if the lowestNode is a closed boundary, it
			// means no outlet can be found
			do
			{
				done = TRUE;  // assume done until proven otherwise
				for ( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
					  cln = lakeIter.NextP() )
				{
					if ( cln->getFloodStatus() != kOutletFlag )
					{
						done = FALSE;
						
						// Check each neighbor
						ce = cln->getEdg();
						do
						{
							node = (tCNode *) ce->getDestinationPtrNC();
							if ( node->getFloodStatus() == kOutletFlag )
							{  // found one!  
								cln->setFloodStatus( kOutletPreFlag );
								cln->setFlowEdg( ce );
								
							}
						} while ( cln->getFloodStatus() != kOutletFlag
								  && ( ce=ce->getCCWEdg() ) != cln->getEdg() );
					} // END if node not flagged as outlet
				} // END for each lake node
				
				// Now flag all the "preflagged" lake nodes as outlets
				for ( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
					  cln = lakeIter.NextP() )
					if ( cln->getFloodStatus() == kOutletPreFlag )
						cln->setFloodStatus( kOutletFlag );
				
			} while ( !done );
			lowestNode->setFloodStatus( kNotFlooded );
			
			// Finally, flag all of the nodes in it as "kFlooded" 
			// and clear the list so we can move on to the next sink
			for ( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
				  cln = lakeIter.NextP() )
				cln->setFloodStatus( kFlooded );
			lakeList.Flush();
      } // END if Sink
   } // END Active Nodes
	return;  
}

/*****************************************************************************
**
**  FindLakeNodeOutlet
**
**  This function is part of the lake-filling algorithm. It checks to see 
**  whether there is a valid outlet for the current node, and if so it 
**  assigns that outlet. An "outlet" essentially means a downhill neighbor 
**  that isn't already flooded to the level of the current node. The function 
**  performs basically the same operation as FlowDirs, but with stricter 
**  criteria. The criteria for a valid outlet are:
**
**  (1) It must be lower than the current node (slope > 0)
**  (2) It must not be part of the current lake (a lake can't outlet to itself)
**  (3) It must not be a closed boundary (_flowAllowed_ must be TRUE)
**  (4) If the outlet is itself part of a different lake, the water surface
**      elevation of that lake must be lower than the current node.
**
**  Returns: TRUE if a valid outlet is found, FALSE otherwise
**  Calls: (none)
**  Called by: FillLakes
**  Created: 6/97 GT
**  Updated: 12/19/97 SL; 1/15/98 gt bug fix (open boundary condition)
**
*****************************************************************************/
int tFlowNet::FindLakeNodeOutlet( tCNode *node )
{
	double maxslp = 0; // Maximum slope found so far
	tEdge * ce;        // Current edge
	tCNode *dn, *an;   // 'dn' - Potential outlet
					   // 'an' - Node ptr used to find outlet of a previously
					   // identified lake
	
	// Check all node's neighbors
	ce = node->getEdg();
	do
	{
		// If it passes this test, it's a valid outlet
		dn = (tCNode *) ce->getDestinationPtrNC();
		assert( dn != nullptr );
		
		// If the node checked IS Outlet - the outlet has been found
		if (dn->getBoundaryFlag() == kOpenBoundary) {
			maxslp = 99999.0;
			node->setFlowEdg( ce );
		}
		else if ( ce->getSlope() > maxslp &&
				  dn->getBoundaryFlag() != kClosedBoundary && 
				  dn->getFloodStatus()  != kCurrentLake )
		{
			// Handle a very special and rare case: if the "target" node dn is
			// part of a previous lake, it's still a valid exit as long as its
			// water surface elevation is lower than the current lake (whose
			// wse, assuming an outlet is found, would be equal to _node_'s
			// elevation). It can sometimes happen that the target lake's wse is
			// exactly equal in elevation to _node_, in which case
			// the point is not considered an outlet---if it were, infinite loops
			// could result
			if ( dn->getFloodStatus() == kFlooded )
			{
				// Iterate "downstream" through the "old" lake until reaching the
				// outlet, then test its elevation. If the elevation is exactly
				// equal to _node_, skip the rest and go on to the next iteration.
				an = dn;
				while ( an->getFloodStatus()  != kNotFlooded &&
						an->getBoundaryFlag() != kOpenBoundary )
                    an = an->getDownstrmNbr();
				if ( an->getZ() == node->getZ() && 
					 an->getBoundaryFlag() != kOpenBoundary) continue;
			}
			// Assign the new max slope and set the flow edge accordingly
			maxslp = ce->getSlope();
			node->setFlowEdg( ce );
		}
	}  while ( ( ce=ce->getCCWEdg() ) != node->getEdg() );
	
	return( maxslp > 0 );
}

/*****************************************************************************
**
**  SetBasinOutlet
**
**  Reset drainage area for all active nodes to zero
**  and find basin outlet node
**
*****************************************************************************/
void tFlowNet::SetBasinOutlet()
{
	tCNode *cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );
	
	// Finding basin outlet node and setting zero contributing area 
	
	for (cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP()) {
		if (cn->getBoundaryFlag() != kClosedBoundary)
			cn->setContrArea( 0.0 );
		
		if (cn->getBoundaryFlag() == kOpenBoundary) {
			OutletNode = cn; //Outlet does not have Voronoi area
			Cout<<"\nBasin outlet node determined..."<<endl<<flush;
		}
	}
	return;
}

/*****************************************************************************
**
**  DrainAreaVoronoi
**
**  Computes drainage area for each node by summing the Voronoi areas of all
**  nodes that drain to it, using the following algorithm:
**
**    FOR each active node
**      Cascade downstream, adding starting node's Voronoi area to each
**      ownstream node's drainage area, until an outlet or sink is reached
**
**    Note that each node's drainage area includes its own Voronoi area.
**
**    Modifies:  node contrbuting area
**
*****************************************************************************/
void tFlowNet::DrainAreaVoronoi()
{
	tCNode * cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );
	
	// Send voronoi area for each node to the node at 
	// the other end of the flowedge and downstream  
	Cout.setf( ios::fixed, ios::floatfield);
	
	for ( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
		RouteFlowArea( cn, cn->getVArea() );
	
	if (OutletNode->getContrArea() < THRESH)
		Cout<<"\nOutlet contributing area:\t"
			<<(int)(OutletNode->getContrArea())<<" m^2"<<endl;
	else 
		Cout<<"\nOutlet contributing area:\t"
			<<(OutletNode->getContrArea())*THRESH<<" km^2"<<endl;
	return;
}

/*****************************************************************************
**
**  RouteFlowArea
**
**  Starting from the current node 'cn', this routine increments 
**  the drainage area of the node and each node downstream by _addedArea_
**
*****************************************************************************/
void tFlowNet::RouteFlowArea( tCNode *cn, double addedArea )
{
	int niterations=0;  // Safety feature: prevents endless loops
	
	// As long as the current node is neither a boundary nor a sink, add
	// _addedArea_ to its total drain. area and advance to the next downstream
	
	do {
		cn->addContrArea( addedArea );
		cn = cn->getDownstrmNbr();
		niterations++;
		assert( niterations <= gridPtr->getNodeList()->getActiveSize() );
	} while ( cn != OutletNode);
	
	cn->addContrArea( addedArea ); //Add it to the Outlet node as well
	
	return;
}

/****************************************************************************
**
**  DeriveCurvature
**
**  Derives curvature for each element in the basin
**
*****************************************************************************/
void tFlowNet::DeriveCurvature() 
{
	double curv;
	tCNode *cnn, *cn;
	tEdge  *firstedg; 
	tEdge  *curedg;
	tMeshListIter< tCNode > niter( gridPtr->getNodeList() );
	
	for ( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() ) {
		curv = 0.0;
		firstedg = cn->getFlowEdg();
		curedg = firstedg->getCCWEdg();
		while (curedg != firstedg) {
			cnn = (tCNode*)curedg->getDestinationPtrNC();
			if (cnn->getBoundaryFlag() != kClosedBoundary &&
				cnn->getBoundaryFlag() != kOpenBoundary) {
				
				if (cnn->getFlowEdg()->getSlope() < 0) {
					if (simCtrl->Verbose_label == 'Y') {
						Cout<<"\nWarning: Negative slope: tFlowNet::DeriveCurvature()\n"
						<<"ID = "<<cnn->getID()<<endl;
					}
				}
				
				if (cnn->getFlowEdg()->getDestinationPtrNC() == (tNode*)cn )
					curv += cnn->getFlowEdg()->getSlope();
			}
			curedg = curedg->getCCWEdg();
		}
		curv -= firstedg->getSlope();
		cn->setCurvature(curv);
	}
	return;
}

/*****************************************************************************
**  
**  DeriveStreamReaches
**
**  Define stream reaches in the proper order according to their position
**  in the drainage network
**
*****************************************************************************/
void tFlowNet::DeriveStreamReaches(tInputFile &infile)
{
	int cnt, flag, ll; // niterations;
	tCNode *cn;
	tCNode *cmove, *cprev;
	
	tMeshListIter< tCNode > niter( gridPtr->getNodeList() );
	tPtrListIter< tCNode >  NodesIterO( NodesLstO );
	tPtrListIter< tCNode >  NodesIterH( NodesLstH );
	tListIter< int >        NNodesIter( NNodes );
	
	cerr.setf( ios::fixed, ios::floatfield);
	Cout.setf( ios::fixed, ios::floatfield);
	
	// Activating tracer and composing a list of river heads 
	
	for ( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
		// ...Tracer value...
		cn->ActivateSortTracer(); // tracer is assigned to '1' for each
	
	// The stream nodes must be ORDERED - that's why  
	// routine 'SortStreamNodes()' is called before.  
	// If the routine work correctly, then we should be
	// able to extract stream reaches quite correctly  
	// in terms of their hierarchical position in the drainage network.  
	
	ll = 0;
	for (cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP()) {
		
		// If the node is a stream node and it has not been
		// included in any stream reach yet (it also can be
		// an outlet for some reach), take it as a new origin 
		
		if (cn->getBoundaryFlag() == kStream && !(cn->NoMoreTracers())) {
			NodesLstH.insertAtBack( cn ); // Channel Head Node
			cprev = cmove = cn;
			flag = 0;
			cnt = 1;
			
			// Go downstream until you reach a confluence/outlet 
			// node, which we define in the code
			while ( !flag ) {
				cprev = cmove; // Always keeps track of the previous node
				cmove = cmove->getDownstrmNbr(); // Downstream node...
				cmove->DeactivateTracer();
				flag = IsConfluence(cmove, cprev);
				if (flag)
					cmove->ActivateSortTracer(); // Re-set tracer for confluence
				cnt++;
			}
			
			//Cout<<"SID = "<<ll<<"\t# of NODES = "<<cnt<<endl;
			NodesLstO.insertAtBack( cmove ); // Channel Outlet (basin O included)
			NNodes.insertAtBack( cnt );      //# of nodes in a stream reach
			ll++;
		}
	}
	// Now we have stream reaches between the points of confluence 
	// and origin. The first list (NodesLstH) contains the origin node. 
	// It is either starts at the river head OR at the node of confluence
	// Therefore, this node is always unique. The second list (NodesLstO)
	// contains outlet points, which are either the points of confluence 
	// for different stream reaches or the Outlet. There are various 
	// reaches having the same outlet node. This feature will be used to 
	// extract the hierarchical order of river network.
	
#ifdef PARALLEL_TRIBS
   // If running in parallel, only the Master will write this file
   if (tParallel::isMaster())
#endif
	PrintArcInfoLinks( infile ); 
	
	Cout<<"\nDerived stream reaches..."<<endl;
	return;
}

/*****************************************************************************
**  
**  PrintArcInfoLinks
**
**  Prints to a file coordinates of stream links derived in the 
**  'DeriveStreamReaches' routine in ArcInfo input format
**  Use 'generate' command to obtain the coverage
**  
*****************************************************************************/
void tFlowNet::PrintArcInfoLinks(tInputFile &infile) 
{
	int cnt; // flag;
	char fullName[kMaxNameSize+20];
	ofstream ControlOut;
	
	tCNode *cn;
	tCNode *cmove, *cprev;
	tPtrListIter< tCNode >  NodesIterO( NodesLstO );
	tPtrListIter< tCNode >  NodesIterH( NodesLstH );
	tListIter< int >        NNodesIter( NNodes );
	
	infile.ReadItem(fullName, "OUTFILENAME" ); // basename 
	strcat( fullName, "_reach");
	
	ControlOut.open(fullName);
	if ( !ControlOut.good() ) {
		cerr<<fullName<<" did not successfully open"<<endl<<flush;
		exit(2);
	}
	ControlOut.setf( ios::fixed, ios::floatfield);
	
	NodesIterO.First(); 
	for (cn=NodesIterH.FirstP(), NodesIterO.First(), NNodesIter.First(), cnt=1; 
		 !(NodesIterH.AtEnd()); 
		 cn=NodesIterH.NextP(),  NodesIterO.Next(),  NNodesIter.Next(),  cnt++) {
		
		cprev = NodesIterO.DatPtr(); // Point to the current Outlet
		cmove = cn;                  // Assign to the current one
		
		ControlOut<<cnt<<endl;  
		ControlOut<<cmove->getX()<<","<<cmove->getY()<<endl;
		
		while (cmove != cprev ) {
			cmove = cmove->getDownstrmNbr();
			ControlOut<<cmove->getX()<<","<<cmove->getY()<<endl;
		}
		ControlOut<<"END"<<endl;
	}
	ControlOut<<"END"<<endl;
	ControlOut.close();
	
	return;
}

/*****************************************************************************
**  
**  WeightedShortestPath()
**
*****************************************************************************/

#define STEDGWEIGHT 2.5
#define SETTLED    3
#define INSTACK    2
#define RECURS     25

void tFlowNet::WeightedShortestPath()
{
	int flag = 1;
	int niterations = 0;
	tCNode *cn;
	
	tMeshListIter<tCNode>  niter( gridPtr->getNodeList() );
	tPtrListIter< tCNode > HeadsIter( HeadsLst );
	tPtrList< tCNode >     NodesLst; 
	tPtrListIter< tCNode > NodesIter( NodesLst );
	tPtrList< tEdge >      EdgeLst; 
	tPtrListIter< tEdge >  EdgeIter( EdgeLst );
	
	Cout<<"\nConnect stream nodes by weighted short path..."<<endl<<flush;
	
	// Set tracer to '1' and use the contributing area variable 
	// to find the shortest path (summing weights for edge length) 
	
	for ( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() ) {
		cn->DeactivateTracer();
		cn->setContrArea(1.0E+9); // Inf at the beginning
	}
	
	// "Settle" the Outlet node
	OutletNode->setTracer(SETTLED);
	OutletNode->setContrArea(0.0);
	
	// Start adding stream nodes to the least
	AddUnsettledNeighbors(OutletNode, NodesLst, EdgeLst);
	
	// The first 'while' loop is intended not to miss  
	// any possible disjoints in stream network of the 
	// basin. There might be many and so the routine  
	// searches for them. The second 'while' is for  
	// correction of flow directions between disjoints
	
	while ( flag ) { 
		// Loop through the nodes of the list 
		// Update the list as we move up in the basin
		NodesIter.FirstP();
		while ( !NodesLst.isEmpty() ) {
			
			// Take the current node
			cn = NodesIter.DatPtr();
			
			AddUnsettledNeighbors(cn, NodesLst, EdgeLst);
			
			// Now, remove the current node from the list  
			NodesLst.removeFromFront();
			EdgeLst.removeFromFront();
			
			// If list is not empty, take the next node
			if (!NodesLst.isEmpty())
				NodesIter.FirstP();
			
			assert( niterations <= gridPtr->getNodeList()->getActiveSize());
			niterations++;
		}
		NodesLst.Flush();
		EdgeLst.Flush();
		
		// Now, define the basin stream heads among "settled" nodes
		
		if (simCtrl->Verbose_label == 'Y') {
			Cout<<"\n\nInitial Stream Head List\nID\tX\tY\tZ\tLength\tSlope\tBND"<<endl;
		}
		for ( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() ) {
			if (cn->getBoundaryFlag() == kStream && cn->getTracer() == SETTLED) {
				if ( IsStreamHead(cn) ) {
					HeadsLst.insertAtBack( cn );
					if (simCtrl->Verbose_label == 'Y')
						TellAboutNode(cn);
				}
			}
		}
		
		// Once we have some stream heads, let's see if we have missed 
		// some stream nodes due to the disjoints in the TIN 
		flag = 0;
		for (cn=HeadsIter.FirstP(); !(HeadsIter.AtEnd()); cn=HeadsIter.NextP()) {
			//Cout<<"\n\n # Checking stream heads:"<<endl;
			//TellAboutNode(cn);
			flag += FindStreamDisjoints(cn, 1, NodesLst);
		}
		
		// If disjoint nodes have been found Flush the stream head list, it would
		// be in wrong order otherwise 
		if (flag)
			HeadsLst.Flush();
		
		// If disjoint nodes have NOT been found by the above procedure 
		// (i.e. searching at a defined stream head), they still may exist, 
		// so use another procedure to figure out if they indeed exist 
		else {
			flag = FindConfluenceDisjoints(NodesLst);
		}
	}

    // setup conditional to check number of stream nodes
    int num_stream_nodes = 0;
    for ( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() ) {
        if (cn->getBoundaryFlag() == kStream) {
            num_stream_nodes++;
        }
    }

	// Find STREAM nodes that have not been passed by during the  
	// preceding operations and make them 'HILLSLOPE' nodes 
    int count = 0;
	for ( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() ) {
		if (cn->getBoundaryFlag() == kStream && cn->getTracer() < INSTACK ) {
			count++;
            cn->setBoundaryFlag( kNonBoundary );
		}
	}

    // Print warning if no stream nodes exist after above for loop.
    if(count == num_stream_nodes){
        cerr<<"Warning!!!\nAll stream nodes were converted to interior nodes.\n "
              "If this is unintentional, increase density of stream nodes"<<endl;
    }
	
	// Due to some peculiarities, we need to re-check stream heads again
	for (cn=HeadsIter.FirstP(); !(HeadsIter.AtEnd()); cn=HeadsIter.NextP()) {
		AddUnsettledNeighbors(cn, NodesLst, EdgeLst);
	}
	HeadsLst.Flush();
	
	// Re-define the basin stream heads one more time 
	for ( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() ) {
		if (cn->getBoundaryFlag() == kStream && cn->getTracer() == SETTLED) {
			if ( IsStreamHead(cn) )
				HeadsLst.insertAtBack( cn );
		}
	}
	
	// Eliminate tributaries that have only two nodes: a stream head
	// and an outlet. Make stream head as a hillslope node
	tPtrListNode < tCNode > * nodeToMove;
	niterations = HeadsLst.getSize();
	cn = HeadsIter.FirstP();
	for ( flag=1; flag <= niterations; flag++ ) {
		if ( IsToEliminate( cn ) ) {
			cn->setBoundaryFlag( kNonBoundary );
			nodeToMove = HeadsIter.NodePtr();
			cn = HeadsIter.NextP();
			HeadsLst.moveToFront(nodeToMove);
			HeadsLst.removeFromFront();
		}
		else 
			cn = HeadsIter.NextP();
	}
	
	// Set contributing area variable to back to '0.0' 
	for ( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() ) {
		cn->setContrArea(0.0);
	}
	
	// Displaying info about stream heads... 
	// Deactivate tracers for stream heads 
	if (simCtrl->Verbose_label == 'Y') {
		Cout<<"\nNumber of stream heads :\t"<<HeadsLst.getSize()<<endl;
		Cout<<"\nFinal Stream Head List\nID\tX\tY\tZ\tLength\tSlope\tBND"<<endl;
		for (cn=HeadsIter.FirstP(); !(HeadsIter.AtEnd()); cn=HeadsIter.NextP())
			TellAboutNode(cn);
	}

	return;
}

/*****************************************************************************
**  
**  AddUnsettledNeighbors()
**  
**  The function considers neighboring to 'cn' nodes and sets the shortest 
**  path. It first take an edge to a node that has been "settled" already 
**  and searches for the other "settled" nodes to set the flowedge. All
**  "unsettled" nodes are added to the general list of nodes that are to 
**  be checked by the calling function
**
*****************************************************************************/
void tFlowNet::AddUnsettledNeighbors(tCNode *cn,
					 tPtrList<tCNode> &NodesLst,
					 tPtrList<tEdge> &EdgeLst)
{
	int cnt = 0;
	double ttt = 0;
	tCNode *cnn;
	tEdge  *firstedg;   // pointer to first edge
	tEdge  *curedg;     // pointer to current edge
	
	ttt = cn->getContrArea();
	
	// 'FlowEdge' should point to a downstream node
	//  in some optimized manner --> after FillLake()
	if (cn->getBoundaryFlag() == kNonBoundary ||
		cn->getBoundaryFlag() == kStream)
		firstedg = cn->getFlowEdg();
	else
		firstedg = cn->getEdg();
	
	// First, find any downslope "settled" node
	if (cn != OutletNode) {
		cnn = (tCNode*)firstedg->getDestinationPtrNC();
		while (cnn->getTracer() != SETTLED && cnt < 100) {
			firstedg = firstedg->getCCWEdg();
			cnn = (tCNode*)firstedg->getDestinationPtrNC();
			cnt++;
		}
		if (cnt >= 100) { 
			Cout<<">>> Error in tFlowNet:: NO Settled neighbors found! "
			<<"Exiting..."<<endl;
			exit(2);
		}
	}
	// Start checking the nodes... 
	cnt = CheckNeighbor( cn, firstedg, NodesLst, EdgeLst );
	
	curedg = firstedg->getCCWEdg();
	while (curedg != firstedg) {
		cnt += CheckNeighbor( cn, curedg, NodesLst, EdgeLst );
		curedg = curedg->getCCWEdg();
	}

	// Make all the necessary checks
	if (cn->getBoundaryFlag() == kNonBoundary ||
		cn->getBoundaryFlag() == kStream) {
		if (cn->getContrArea() >= ttt && cn->getTracer() != SETTLED ) {
			Cout<<">>> tFlowNet: AddUnsettledNeighbors: "
			<<"The path variable has not changed! Exiting..."<<endl;
			exit(2);
		}
		// Once the flow direction is chosen, mark this node as "settled"
		else
			cn->setTracer(SETTLED);
	}
	return;
}

/*****************************************************************************
**  
**  CheckNeighbor()
**  
**  The function takes as arguments a ptr to a current node 'cn', an edge
**  'curedge' that originates at 'cn' and a node list NodesLst. The routine
**  analyzes the node located on another end of 'curedg'. Depending to what 
**  its tracer equals to, it:
**      - Does nothing with it (tracer == INSTACK);
**      - Checks if 'cnn' is suitable for flowing into (tracer == SETTLED);  
**      - Adds the node to the list of nodes to be analyzed ((tracer < INSTACK);
**
*****************************************************************************/
int tFlowNet::CheckNeighbor(tCNode *cn, tEdge *curedg,
				tPtrList<tCNode> &NodesLst,
				tPtrList<tEdge> &EdgeLst)
{
	int cnt = 0;
	double tempo;
	tCNode *cnn;
	
	cnn = (tCNode*)curedg->getDestinationPtrNC();
	
	// Check STREAM nodes
	if (cnn->getBoundaryFlag() == kStream || 
		cnn->getBoundaryFlag() == kOpenBoundary) {
		
		// If this is a "Settled" node --> compare its weight
		// plus the weight of the connecting edge to the 
		// current weight of the node
		if (cnn->getTracer() == SETTLED) { 
			tempo = cnn->getContrArea();
			tempo += ComputeEdgeWeight( curedg, STEDGWEIGHT );
			
			// If we change the flow direction for 'cn'
			// we need to check and update all its neighbors
			// if 'cn' is a better choice for them to flow to 
			// --> call 'UpdatePathVariable'
			if ( tempo < cn->getContrArea() ) {
				cn->setContrArea( tempo );
				cn->setFlowEdg( curedg );
				UpdatePathVariable( cn );
			}
			
			// If current weight assigned to 'cn' is 
			// less than the one assigned to 'cnn'
			// ---> change the flow direction for 'cnn'
			else {
				tempo = cn->getContrArea();
				tempo += ComputeEdgeWeight( curedg->FindComplement(), STEDGWEIGHT );
				
				if ( tempo < cnn->getContrArea() ) {
					cnn->setContrArea( tempo );
					cnn->setFlowEdg( curedg->FindComplement() );
					UpdatePathVariable( cnn );
				}
			}
		}
		// If the stream node is "Unsettled" --> put it in the stack 
		else if (cnn->getTracer() < INSTACK) {
			NodesLst.insertAtBack( cnn );
			EdgeLst.insertAtBack( curedg );
			cnn->setTracer(INSTACK);
			cnt++;
		}
	}
	// For the break points we also need to check HILLSLOPE nodes 
	else if (cnn->getBoundaryFlag() == kNonBoundary ) {
		;
	}
	return cnt;
}

/*****************************************************************************
**  
**  ComputeEdgeWeight()
**  
**  Computes the flow edge weight for the algorithm of shortest path
**
*****************************************************************************/
double tFlowNet::ComputeEdgeWeight( tEdge *curedg, double weight ) 
{
	double tempo;
	tempo = pow((curedg->getLength()/100.),weight);
	return tempo;
}

/*****************************************************************************
**  
**  IsStreamHead()
**  
**  Checks if node 'cn' is a stream head. For that it checks all neighbors of
**  'cn', if there is any that has a flow edge pointing to 'cn' --> 'cn' is 
**  NOT a stream head (return '0'). Otherwise, it is a stream head (return '1')
**  
*****************************************************************************/
int tFlowNet::IsStreamHead(tCNode *cn)
{
	int cnt = 0;
	tCNode *cnn;
	tEdge  *firstedg;   // pointer to first edge
	tEdge  *curedg;     // pointer to current edge
	
	if (cn->getBoundaryFlag() != kStream) 
		return cnt;
	
	// The stream node which is the destination node of current edge
	// may have FLOWEDGE pointed to the 'cn' node (the node being tested)
	
	firstedg = cn->getFlowEdg();
	cnn = (tCNode*)firstedg->getDestinationPtrNC();
	if (cnn->getBoundaryFlag() == kStream ) {
		if (cnn->getFlowEdg()->getDestinationPtrNC() == (tNode*)cn )
			cnt++;
	}
	
	curedg = firstedg->getCCWEdg();
	while (curedg != firstedg) {
		cnn = (tCNode*)curedg->getDestinationPtrNC();
		if (cnn->getBoundaryFlag() == kStream ) {
			if (cnn->getFlowEdg()->getDestinationPtrNC() == (tNode*)cn )
				cnt++;
		}
		curedg = curedg->getCCWEdg();
	}
	if (cnt > 0)
		cnt = 0;
	else if (cnt == 0)
		cnt = 1;
	
	return cnt;
}

/*****************************************************************************
**  
**  IsStreamHead()
**  
**  The function checks if the tributary starting at the stream head node 
**  'chead' has to be eliminated (the length of the tributary is 2 nodes:
**  the head 'chead' and the outlet which is the confluence node).  
**  
*****************************************************************************/
int tFlowNet::IsToEliminate(tCNode *chead)
{
	int cnt = 0;
	tCNode *cn, *cnn;
	tEdge  *firstedg;   // pointer to first edge
	tEdge  *curedg;     // pointer to current edge
	
	if (chead->getBoundaryFlag() != kStream) 
		return 0;
	
	// Get the downstream node of a stream head node
	cn = chead->getDownstrmNbr();
	
	// The stream node which is the destination node of an edge
	// may have FLOWEDGE pointed to the 'cn' node (the node being tested)
	
	if ( cn != OutletNode ) {
		firstedg = cn->getFlowEdg();
		cnn = (tCNode*)firstedg->getDestinationPtrNC();
		if (cnn->getBoundaryFlag() == kStream && cnn != chead) {
			if (cnn->getFlowEdg()->getDestinationPtrNC() == (tNode*)cn )
				cnt++;
		}
		
		curedg = firstedg->getCCWEdg();
		while (curedg != firstedg) {
			cnn = (tCNode*)curedg->getDestinationPtrNC();
			if (cnn->getBoundaryFlag() == kStream && cnn != chead) {
				if (cnn->getFlowEdg()->getDestinationPtrNC() == (tNode*)cn )
					cnt++;
			}
			curedg = curedg->getCCWEdg();
		}
	}
	return cnt;
}

/*****************************************************************************
**  
**  FindConfluenceDisjoints()
**  
**  Disjoints may exist at the confluence points so here we use another
**  procedure to figure out if they indeed exist. Return '0' if there are 
**  no disjoint tributaries at confluences, OR if they can not be found
**  by the existing algorithm (the like reason is the wrong finding of the 
**  node for search of the confluence)
**  
*****************************************************************************/
int tFlowNet::FindConfluenceDisjoints(tPtrList<tCNode> &NodesLst)
{
	int cnt = 1;
	int flag = 0;
	int niterations = 0;
	tCNode *cn, *cmove;
	
	tMeshListIter<tCNode> niter( gridPtr->getNodeList() );
	
	// Loop through the nodelist until all the possibilities are checked, i.e.
	// loop until you find stream head that may lead us to a confluence, if it
	// does - get out and work with this tributary, if not - continue searching 
	
	while (cnt > 0 && flag == 0) {
		
		cnt  = 0;
		flag = 0;
		
		cn = niter.FirstP();
		while ( niter.IsActive() && cnt == 0 ) {
			// Consider only nodes with (tracer == 0)
			if (cn->getBoundaryFlag() == kStream && cn->NoMoreTracers()) {
				if ( IsStreamHead(cn) ) {
					cnt++;
				}
			}
			cmove = cn;
			cn = niter.NextP();
		}
		
		// If an unsettled stream head has been found --> we need 
		// to find where this tributary joins the main network of 
		// settled nodes. For that, go down stream of 'cmove' until   
		// a settled node is reached. Figure out how to join them. 
		
		if ( cnt ) {
			niterations = 0;
			do { 
				cn = cmove;
				cmove->setTracer(-1); // <-- Assign tracer to '-1'
				cmove = cmove->getDownstrmNbr();
				
				// Get out of the loop if you have found a settled node 
				// OR have reached the basin outlet 
				if ((cmove->getTracer() == SETTLED && 
					 cmove->getBoundaryFlag() == kStream) || cmove == OutletNode)
					niterations = -9999;
				niterations++;
				assert( niterations < gridPtr->getNodeList()->getActiveSize() );
			} while ( niterations > 0 );
			
			// We now have found the outlet node for the tributary. In the 
			// vicinity of this node we should look for the node that is a 
			// true tributary outlet (the first one may be fictitious outlet) 
			// BUT: we will assume that 'cmove' IS the true outlet!
			// It is a strong assumption but in most cases would be true. 
			// If it is not - search for the true outlet among 'cmove' neighbors 
			// (most likely it would be upstream of 'cmove')
			
			// NOTE: Please note, that if no stream disjoints have been 
			// found in the vicinity of 'cmove', the algorithm will return 
			// (flag = 0) and therefore the tributary will not be ever found
			
			flag = FindStreamDisjoints(cmove, 1, NodesLst);
		}
	}
	return flag;
}

/*****************************************************************************
**  
**  FindStreamDisjoints()
**  
**  The function is intended to define if a current stream head node 'cn' is
**  actually a disjoint node in stream network.  For that, the function looks
**  at the hillslope neighbors of 'cn': if some of these nodes have connection
**  to a stream node that has "unsettled" value of tracer (< INSTACK) then the 
**  node 'cn' is a disjoint node. Correspondingly, an optimum path to the found
**  "unsettled" stream node is defined through a hillslope node (which becomes
**  'stream') and it is stored in the general stack of nodes 'NodesLst'.
**  
**  NOTE: - So far, it is assumed that there might be only ONE hillslope node 
**  that disjoints stream network. The function might be designed for recursive
**  calls to look for a stream node through several disjoining hillslope nodes.
**  Variable 'times' is intended for that.
** 
*****************************************************************************/
int tFlowNet::FindStreamDisjoints(tCNode *cn, int times, 
				  tPtrList<tCNode> &NodesLst)
{
	int cnt = 0;
	double tempo = 0;
	double min = 1.0E+9;
	tCNode *cnn;
	tEdge  *firstedg;   // pointer to first edge
	tEdge  *curedg;     // pointer to current edge
	
	tPtrList< tEdge >     EdgeLst1; 
	tPtrListIter< tEdge > EdgeIter1( EdgeLst1 );
	tPtrList< tEdge >     EdgeLst2; 
	tPtrListIter< tEdge > EdgeIter2( EdgeLst2 );
	
	if (times > 0) {
		times--;
		if ((cn->getBoundaryFlag() == kStream &&
			 cn->getTracer() == SETTLED) || cn == OutletNode) {
			
			if (cn == OutletNode)
				firstedg = cn->getEdg();
			else
				firstedg = cn->getFlowEdg();
			cnn = (tCNode*)firstedg->getDestinationPtrNC();
			
			// Check the node
			cnt += IsStreamDisjoint(firstedg, cnn, EdgeLst1, EdgeLst2);
			
			curedg = firstedg->getCCWEdg();
			while (curedg != firstedg) {
				cnn = (tCNode*)curedg->getDestinationPtrNC();
				cnt += IsStreamDisjoint(curedg, cnn, EdgeLst1, EdgeLst2);
				curedg = curedg->getCCWEdg();
			}
			
			// Now, if (cnt > 0), then some nodes have been found
			if (cnt > 0) {
				for (EdgeIter1.First(),EdgeIter2.First(); 
					 !(EdgeIter1.AtEnd()); 
					 EdgeIter1.Next(), EdgeIter2.Next() ) {
					
					// Compute path weight
					tempo  = ComputeEdgeWeight(EdgeIter1.DatPtr(), 1);
					tempo += ComputeEdgeWeight(EdgeIter2.DatPtr(), 1);
					
					// Choose, if it provides min path
					// Store the value and the edges
					if (tempo < min) {
						min = tempo;
						firstedg = EdgeIter1.DatPtr();
						curedg   = EdgeIter2.DatPtr();
					}
				}
				// By now, the min path has been chosen
				// Do appropriate assignments:
				
				cnn = (tCNode*)firstedg->getDestinationPtrNC();
				// 1.) Change its boundary flag
				cnn->setBoundaryFlag( kStream );
				
				// 2.) Assign the flow edge - complimentary to 'firstedg'
				cnn->setFlowEdg( firstedg->FindComplement() );
				
				// 3.) Set tracer to "settled"
				cnn->setTracer(SETTLED);
				
				// 4.) Assign the weight
				tempo = cn->getContrArea();
				tempo += ComputeEdgeWeight(firstedg, STEDGWEIGHT);
				cnn->setContrArea( tempo );
				
				// 5.) Put the found stream node in stack
				cnn = (tCNode*)curedg->getDestinationPtrNC();
				NodesLst.insertAtBack( cnn );
				
				// 6.) Clean up temporary stacks
				EdgeLst1.Flush();
				EdgeLst2.Flush();
      }
			// else --> for recursive calls (& times)
    }
  }
	return cnt;
}

/*****************************************************************************
**  
**  IsStreamDisjoint()
**  
**  From a current stream head (in the calling function), analyze a hillslope
**  node 'cn', a destination node of edge 'edg'. If 'cn' has connection 
**  to "unsettled" stream nodes, put two edges in the stacks, return cnt > 0
**  
*****************************************************************************/
int tFlowNet::IsStreamDisjoint(tEdge *edg, tCNode *cn, 
				   tPtrList<tEdge> &EdgeLst1,
				   tPtrList<tEdge> &EdgeLst2)
{
	int cnt = 0;
	int ntimes;
	double tempo;
	tCNode *cnn;
	tEdge  *firstedg;   // pointer to first edge
	tEdge  *curedg;     // pointer to current edge
	
	// If 'cn' is not a hillslope node, we do not consider it at all
	if (cn->getBoundaryFlag() == kNonBoundary) {
		
		// Now, consider edge by edge if there is 
		// any stream node that has not been assigned 
		// the flow direction in the weighted shortest path algorithm 
		firstedg = cn->getEdg();
		cnn = (tCNode*)firstedg->getDestinationPtrNC();
		if (cnn->getTracer() < INSTACK && cnn->getBoundaryFlag() == kStream) {
			
			// Now, check if this node has positive elevation tendency,
			// i.e. minimum elevation of stream nodes increases
			ntimes = 1;
			tempo = CompElevationTendency(cn, cnn, &ntimes, RECURS);
			tempo /= ntimes;
			
			// IF the tendency is higher than the node elevation - fine
			// ELSE - we've found a node on another tributary 
			if (tempo > cnn->getZ()) {
				EdgeLst1.insertAtBack( edg );        // to Hillslope node
				EdgeLst2.insertAtBack( firstedg );   // to Stream node 
				cnt++;
			}
		}
		
		curedg = firstedg->getCCWEdg();
		while (curedg != firstedg) {
			cnn = (tCNode*)curedg->getDestinationPtrNC();
			if (cnn->getTracer() < INSTACK && cnn->getBoundaryFlag() == kStream) {
				ntimes = 1;
				tempo = CompElevationTendency(cn, cnn, &ntimes, RECURS);
				tempo /= ntimes;
				
				if (tempo > cnn->getZ()) {
					EdgeLst1.insertAtBack( edg );
					EdgeLst2.insertAtBack( curedg );
					cnt++;
				}
			}
			curedg = curedg->getCCWEdg();
		}
	}
	return cnt;
}

/*****************************************************************************
**  
**  IsConfluence
**
**  Checks if stream node 'cn' is a confluence node. Node 'cup' defines the 
**  upstream node for 'cn' of a stream link for which junction node is sought 
**
**  - Returns '1' if 'cn' is a confluence or outlet
**  - Returns '0' otherwise
**
*****************************************************************************/
int tFlowNet::IsConfluence(tCNode *cn, tCNode *cup) 
{
	int flag = 0;
	int niterations = 0;
	tCNode *cnn, *citer;
	tEdge  *firstedg; 
	tEdge  *curedg;
	
	// This flowedge has been checked before
	// It points to a downstream STREAM node
	if ( cn->getBoundaryFlag() == kOpenBoundary ) //Outlet!
		flag = 1;  // Outlet Node -- Exit!
	else {
		firstedg = cn->getFlowEdg();
		curedg = firstedg->getCCWEdg();
		while ((curedg != firstedg) && (!flag)) {
			
			// CHECK ## 1 : First, we look only at stream nodes, != upper node
			cnn = (tCNode*)curedg->getDestinationPtrNC();
			if (cnn != cup &&
				cnn->getBoundaryFlag() == kStream) {
				
				// CHECK ## 2 : Slope condition may not work (ideally, 
				// upstream node from the tributary should be upslope); 
				// Check the contributing area condition
				if (cn->getContrArea() > cnn->getContrArea()) {
					
					// CHECK ## 3 : The stream node which is the destination node 
					// of 'curedg' must have FLOWEDGE pointed to 'cn' node (test node)
					if (cnn->getFlowEdg()->getDestinationPtrNC() == (tNode*)cn ) {
						
						// CHECK ## 4 : Last check. We have to make sure that 'cnn' 
						// is not one of the downstream nodes which streamflow from 
						// 'cn' would pass by on its way to Outlet (this would 
						// actually duplicate Check ## 2)
						flag = 1;
						citer = cn;
						do { 
							citer = citer->getDownstrmNbr();
							if (citer == cnn)
								flag = 0;
							niterations++;
							assert( niterations < gridPtr->getNodeList()->getActiveSize() );
						} while ( niterations < 30 && citer != OutletNode && flag == 1 );
					}
				}
			}
			curedg = curedg->getCCWEdg();
		}
	}
	return flag;
}

/*****************************************************************************
**  
**  CompElevationTendency()
**  
**  Computes average elevation of neighboring to 'cn_check' nodes which
**  are not shared with the 'cn' node ('cn' is usually upstream of 'cn_check'
**  and we therefore try to figure out if 'cn_check' leads us downslope)
**
*****************************************************************************/
double tFlowNet::CompElevationTendency(tCNode *cn, tCNode *cn_check, 
					   int *cnt, int flag)
{
	double Elev = 1.0E+6;
	double tempo = 0;
	tCNode *cnn, *cmm;
	tEdge  *firstedg;   // pointer to first edge
	tEdge  *curedg;     // pointer to current edge
	tPtrList< tCNode >     DnLst; 
	tPtrListIter< tCNode > DnIter( DnLst );
	
	tList< int >       Tracer;   // Tracer value in a stream node
	tListIter< int >   TracerIter( Tracer );
	
	if (flag > 0) { 
		flag--;
		
		firstedg = cn_check->getEdg();
		cnn = (tCNode*)firstedg->getDestinationPtrNC();
		
		// If 'cnn' is not shared with 'cn' 
		// use its elevation for computation
		if (cnn->getBoundaryFlag() == kOpenBoundary) {
			Elev = flag-10;
			flag = 0;
		}
		
		// 1.) If it is not an upstream node
		// 2.) If it is not connected to the upstream node
		// 3.) If it is a _Stream_ node
		// 4.) We have not previously stored it in the stack
		if (cnn != cn && !IsConnected(cn, cnn) && 
			cnn->getBoundaryFlag() == kStream && 
			cnn->getTracer() < 1) {
			Elev = cnn->getZ();
			DnLst.insertAtBack( cnn );
			Tracer.insertAtBack( cnn->getTracer() ); 
			cnn->setTracer(1);  // Set tracer to a "stored" value
		}
		curedg = firstedg->getCCWEdg();
		while (curedg != firstedg) {
			cnn = (tCNode*)curedg->getDestinationPtrNC();
			// If it's an Outlet - make it the lowest
			if (cnn->getBoundaryFlag() == kOpenBoundary) {
				Elev = flag-10;
				flag = 0;
			}
			// If 'cnn' is not connected to the upper node 'cn'
			// --> check its elevation and add to the list
			// of nodes that are to be checked
			if (cnn != cn && !IsConnected(cn, cnn) && 
				cnn->getBoundaryFlag() == kStream &&
				cnn->getTracer() < 1) {
				if (cnn->getZ() < Elev)
					Elev = cnn->getZ();
				DnLst.insertAtBack( cnn );
				Tracer.insertAtBack( cnn->getTracer() ); 
				
				// Set tracer to a "stored" value
				cnn->setTracer(1);
			}
			curedg = curedg->getCCWEdg();
		}
		
		// Now, for each of the nodes and compute its tendency
		// Estimate the average minimum tendency and add it to Elev
		
		double minElv = 0;
		
		if (DnLst.getSize() > 0) {
			for (cmm=DnIter.FirstP(); !(DnIter.AtEnd()); cmm=DnIter.NextP()) {
				tempo = CompElevationTendency(cn_check, cmm, cnt, flag);
				if (tempo != 1.0E+6) {
					minElv += tempo;  // if (tempo < minElv)
					*cnt = *cnt + 1;
				}
			}
			// Use "average minimum" value 
			if ( minElv > 0 && minElv < 1.0E+6 ) {
				// 1.) minElv /= cnt;  // Elev = (Elev + minElv)/2.;
				// 2.) Elev = (Elev + minElv)/(cnt+1);
				Elev += minElv;
			}
			// Re-set tracer values to what they were...
			for (cmm = DnIter.FirstP(), TracerIter.First(); 
				 !(DnIter.AtEnd());
				 cmm = DnIter.NextP(), TracerIter.Next())
				cmm->setTracer( TracerIter.DatRef() );
		}
		else
			flag = 0;
		
		// Release memory
		DnLst.Flush();
		Tracer.Flush();
	}
	return Elev;
}

/*****************************************************************************
**  
**  IsConnected()
**  
**  The routine checks if the node 'cn' has a common edge with 'cn_check'
**  Returns '1' if so, '0' otherwise
**
*****************************************************************************/
int tFlowNet::IsConnected(tCNode *cn, tCNode *cn_check) 
{
	int    flaggss = 0;
	tCNode *cnn;
	tEdge  *firstedg;   // pointer to first edge
	tEdge  *curedg;     // pointer to current edge
	
	firstedg = cn->getFlowEdg();
	cnn = (tCNode*)firstedg->getDestinationPtrNC();
	if (cnn == cn_check)
		return 1;
	else {
		curedg = firstedg->getCCWEdg();
		while (curedg != firstedg && flaggss == 0) {
			cnn = (tCNode*)curedg->getDestinationPtrNC();
			if (cnn == cn_check)
				flaggss = 1;    
			curedg = curedg->getCCWEdg();
		}
	}
	return flaggss;
}

/*****************************************************************************
**  
**  FindAngle()
**  
**  The routine finds angle for the nodes 'cn', 'cn1', and 'cn2'
**  It is assumed that they form vectors:
**           cn ---------- cn1
**             \ ) <--  alpha
**              \
**               \
**                \
**                cn2
**  Returned value is in degrees
**
*****************************************************************************/
double tFlowNet::FindAngle(tCNode *cn, tCNode *cn1, tCNode *cn2)
{
	double d1, d2, x1, y1, x2, y2, alpha;
	tArray<double> xy(2), xy1(2), xy2(2);
	
	xy2 = cn2->get2DCoords();
	xy1 = cn1->get2DCoords();
	xy  = cn->get2DCoords();
	
	x1  = xy1[0]-xy[0]; //<= THESE ARE COORDINATES
	y1  = xy1[1]-xy[1]; //<= OF THE VECTORS
	x2  = xy2[0]-xy[0]; //<= REQUIRED FOR FURTHER
	y2  = xy2[1]-xy[1]; //<= CALCULATIONS
	
	d1 = FindDistance(xy[0], xy[1], xy1[0], xy1[1]);
	d2 = FindDistance(xy[0], xy[1], xy2[0], xy2[1]);
	
	// Compute the angle between vectors cn-cn1 and cn-cn2
	// Angle is in degrees
	alpha = acos((x1*x2 + y1*y2)/(d1*d2))*180/(4*atan(1.0));  //Added decimal
	
	return alpha;
}

/***************************************************************************
**
**  FindDistance( )
**
**  Finds distance between two nodes located at (x1,y1) and (x2,y2)
**
***************************************************************************/
double tFlowNet::FindDistance(double x1, double y1, double x2, double y2) 
{
	return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
}

/***************************************************************************
**
**  TellAboutNode()
**
**  Prints node info
**
***************************************************************************/
void tFlowNet::TellAboutNode(tCNode *cn)
{
	if (simCtrl->Verbose_label == 'Y') {
		if (cn->getBoundaryFlag() == kNonBoundary ||
			cn->getBoundaryFlag() == kStream)
			Cout<<cn->getID()
				<<"\t"<<cn->getX()
				<<"\t"<<cn->getY()
				<<"\t"<<cn->getZ()
				<<"\t"<<cn->getFlowEdg()->getLength()
				<<"\t"<<cn->getFlowEdg()->getSlope()
				<<"\t"<<cn->getBoundaryFlag()<<endl<<flush;
		else 
			Cout<<cn->getID()
				<<"\t"<<cn->getX()
				<<"\t"<<cn->getY()
				<<"\t"<<cn->getZ()
				<<"\t-999\t-999"
				<<"\t"<<cn->getBoundaryFlag()<<endl<<flush; 
	}
	return;
}

/*****************************************************************************
**  
**  SortStreamNodes
**
**  Similar to SortNodesByStreamOrder but deals only with stream nodes.
**  Reason:  for some cases we have to re-adjust flow directions for stream
**           nodes and therefore it is necessary to ensure that the order 
**           of computations is correct
**  Algorithm: - Send all nodes to the back of the list such all hillslope 
**               nodes are preceding
**             - Sort stream nodes according to the relationship "drains to"
**             - Re-enumerate indices of the nodes
** 
*****************************************************************************/
void tFlowNet::SortStreamNodes() 
{
	int i;
	int nThisPass, nPassed; // Number moved in current iteration & in total
	int nStreams = 0;
	
	tCNode * cn;
	tMeshList<tCNode> *nodeList = gridPtr->getNodeList();
	tMeshListIter<tCNode> listIter( nodeList );
	tListNode< tCNode > * nodeToMove;
	
	int nUnsortedNodes = nodeList->getActiveSize(); // Number not yet sorted
	
	// First, move all the downstream nodes to the end of list,
	// following FIFO principle. Assign initial tracers: use "qs" 
	// field, which contains garbage at this stage.
	cn = listIter.FirstP();
	for ( i=1; i <= nUnsortedNodes; i++ ) {
		cn->ActivateSortTracer();  // <-- assign it to one
		
		// If stream, move to bottom of the list
		if (cn->getBoundaryFlag() == kStream) {
			nStreams++;
			nodeToMove = listIter.NodePtr();
			cn = listIter.NextP();
			nodeList->moveToActiveBack( nodeToMove );
		}
		else
			cn = listIter.NextP();
	}
	
	// Iterate: move tracers downstream and sort 
	// until no nodes with tracers are left
	nPassed = 0;
	do {
		// --- Send tracers downstream ---
		cn = listIter.FirstP();
		for ( i=1; i <= nUnsortedNodes; i++ ) {
			assert( cn!=0 );
			if (cn->getBoundaryFlag() == kStream)
				cn->MoveSortTracerDownstream();
			cn = listIter.NextP();
		}
		
		// Scan for any nodes that have no tracers, and move them 
		// to the bottom of the list. This simply means that there 
		// were no any nodes for which they would lie downstream! 
		nThisPass = 0;
		cn = listIter.FirstP();
		for ( i=1; i<=nUnsortedNodes; i++ ) {
			if (cn->getBoundaryFlag() == kStream && cn->NoMoreTracers() ) {
				// --- If no tracers, move to bottom of list ---
				nodeToMove = listIter.NodePtr();
				cn = listIter.NextP();
				nodeList->moveToActiveBack( nodeToMove );
				nThisPass++;
			}
			else 
				cn = listIter.NextP();
		}
		nPassed += nThisPass;
		nUnsortedNodes -= nThisPass;
	} while ( nPassed < nStreams );
	
	// Changed to make all IDs consistent
	tEdge     * ce;
	tTriangle * ct;  
	tMeshListIter<tCNode> niter( gridPtr->getNodeList() );
	tMeshListIter<tEdge>  eiter( gridPtr->getEdgeList() );
	tListIter<tTriangle>  titer( gridPtr->getTriList() );
	
	int nnodes = gridPtr->getNodeList()->getSize();
	int nedges = gridPtr->getEdgeList()->getSize();
	int ntri   = gridPtr->getTriList()->getSize();
	int id;
	
	for ( cn=niter.FirstP(), id=0; id<nnodes; cn=niter.NextP(), id++ )
		cn->setID( id );
	for ( ce=eiter.FirstP(), id=0; id<nedges; ce=eiter.NextP(), id++ )
		ce->setID( id );
	for ( ct=titer.FirstP(), id=0; id<ntri; ct=titer.NextP(), id++ )
		ct->setID( id );
	Cout<<"\nSorting stream nodes completed..."<<endl<<flush;
	return;
}

/*****************************************************************************
**  
**  ComputeDistanceToStream()
**  
**  Computes 'hillpath' for hillslope nodes, i.e. computes distance and 
**  initializes pointer to the nearest stream node 
**
*****************************************************************************/
void tFlowNet::ComputeDistanceToStream() 
{
	Cout<<"\nComputing distances to streams..."<<endl;
	
	tCNode *cn;
	tCNode *ctimer;
	tEdge  *ce;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );
	tMeshListIter<tEdge>  edgIter( gridPtr->getEdgeList() );
	double hill, tt;
	
	BasArea = 0.0;
	maxttime = 0.0;        // note: the travel times are given in SECONDS
	dist_hill_max = 0.0;   // METERS
	dist_stream_max = 0.0; // METERS
	
	for (cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP()) {
		if (cn->getBoundaryFlag() != kClosedBoundary)
			cn->setContrArea( 0.0 );
		
		if (cn->getBoundaryFlag() == kOpenBoundary) {
			OutletNode = cn; // Outlet does not have Voronoi area
			Cout<<"->> The outlet node has been set up!"<<endl<<flush;
		}
	}
	
	// Loop through the nodes and set velocity 
	// Assign pointer to the nearest stream node 
	// Define maximum hillslope travel time also 
	for ( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
	{
		BasArea += cn->getVArea();
	
		//Initialize with distance from centroid to node
		cn->getCentroidX();
		cn->getCentroidY();
		hill = sqrt(pow(cn->getX()-cn->getCentroidX(),2.0)+pow(cn->getY()-cn->getCentroidY(),2.0));
		
		tt = 0.0;
		ctimer = cn;
		
		// Go to all downstream nodes to the outlet
		while ( ctimer->getBoundaryFlag() == kNonBoundary) {
			ce = ctimer->getFlowEdg(); // Get the steepest flowedge
			hill += ce->getLength();
			ctimer = ctimer->getDownstrmNbr();
		}
		
		cn->setHillPath(hill);     // Set the HILLSLOPE path for node
		cn->setStreamNode(ctimer); // Set a stream node to which it contributes
		
		tt = hill/hillvel;  // SECONDS
		tt /= 3600.;        // HOURS 
		cn->setTTime( tt ); // Set Travel Time, HOURS <- To the STREAM Node
		
		if (tt > maxttime) {  // MAX values of path should correspond to maxttime
			maxttime = tt;
			dist_hill_max = hill; // METERS
		}
	}
	
	Cout.setf( ios::fixed, ios::floatfield);
	Cout<<"->> The Total Basin Area is "<<BasArea<<" M^2"<<endl;
	Cout<<"->> The Maximum Travel Time is "<<maxttime<<" HOURS"<<endl;
	Cout<<"->> The Maximum Hillslope path is\t"<<dist_hill_max<<" METERS"<<endl;
	Cout<<"tFlowNet: ComputeDistanceToStream() finished..."<<endl<<flush;
	return;
}

/*****************************************************************************
**  
**  UpdatePathVariable()
**  
**  The function updates the value of path variable for any neighbors of
**  'cn' that are "settled" and has flowedge to 'cn'
**  
*****************************************************************************/
void tFlowNet::UpdatePathVariable(tCNode *cn)
{
	double tempo;
	tCNode *cnn;
	tEdge  *firstedg;   // pointer to first edge
	tEdge  *curedg;     // pointer to current edge
	
	if (cn->getBoundaryFlag() != kStream) 
		return;
	
	// If for any of the _settled_ neighbors 'cn'
	// is a better choice to flow to ->
	// then we need to update its weight 
	firstedg = cn->getFlowEdg();
	curedg = firstedg->getCCWEdg();
	while (curedg != firstedg) {
		cnn = (tCNode*)curedg->getDestinationPtrNC();
		
		if (cnn->getBoundaryFlag() == kStream && cnn->getTracer() == SETTLED) {
			tempo = cn->getContrArea();
			tempo += ComputeEdgeWeight( curedg->FindComplement(), STEDGWEIGHT );
			if ( tempo < cnn->getContrArea() ) {
				cnn->setContrArea( tempo );
				cnn->setFlowEdg( curedg->FindComplement() );
				UpdatePathVariable( cnn );
			}
		}
		curedg = curedg->getCCWEdg();
	}
	return;
}

#undef STEDGWEIGHT
#undef SETTLED
#undef INSTACK
#undef RECURS

/*************************************************************************
**
**  CheckVDrainageWidths
**
**  The function loops through the active node list in search of the 
**  nodes that have zero (or close to zero) voronoi edge length in the
**  steepest drainage direction. It calls a function that computes a 
**  pseudo width instead defined as the area of voronoi sector
**  located between the CW & CCW neighbors of the flow edge divided by 
**  the flow edge length
**
*************************************************************************/
void tFlowNet::CheckVDrainageWidths() 
{
	tCNode *cn;
	tMeshListIter<tCNode> nodIter( gridPtr->getNodeList() );
	
	Cout<<"tFlowNet: Checking Voronoi drainage widths..."<<endl;
	
	for ( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
		FixVoronoiEdgeWidth( cn );
	
	return;
}

/*************************************************************************
**
**  FixVoronoiEdgeWidth
**
**  The function computes a pseudo width for a node 'cn'. The pseudo width
**  is defined as the area of voronoi sector located between the CW & CCW 
**  neighbors of the flow edge divided by the flow edge length
**
*************************************************************************/
void tFlowNet::FixVoronoiEdgeWidth(tCNode *cn) 
{
	int i, j, flag, flagg, NV;
	int changeWidth;
	double sArea, wWidth;
	double **polyg;
	
	tEdge  *flowedg, *ccwedg, *cwedg;
	tEdge  *firstedg, *curedg;

	tArray<double> centroid(2);		// Centroid of voronoi polygon
	tArray<double> vv_flow(2);		// Voronoi vertex for flow edge
	tArray<double> vv_ccw1(2);		// First CCW voronoi vertex
	tArray<double> vv_ccw3(2);		// Second CCW voronoi vertex
	tArray<double> vv_cw2(2);			// First CW voronoi vertex

	tArray<double> orig_ccw(2), dest_ccw(2);		// CCW edge orig and dest points
	tArray<double> orig_cw(2), dest_cw(2);			// CW edge orig and dest points
	tArray<double> orig_flow(2), dest_flow(2);	// Flow edge orig and dest pts

	tArray<double> intersect_ccw(2), intersect_cw(2);
	tArray<double> vv_cur(2), vv_tmp(2);
	
	cout<<setprecision(6);
	
	flag = flagg = 0;
	wWidth = sArea = 0;
	changeWidth = 0;
	
	// At most NV vertices are considered when fixing an edge width
	NV = 5;
	polyg = new double* [2];
	for (i=0; i < 2; i++) {
		polyg[i] = new double[NV];
		assert(polyg[i] != 0);
	}

	// Centroid information
	centroid = cn->get2DCoords();

	// Flow edge information
	flowedg = cn->getFlowEdg();
	vv_flow = flowedg->getRVtx();
	orig_flow = flowedg->getOriginPtrNC()->get2DCoords();
	dest_flow = flowedg->getDestinationPtrNC()->get2DCoords();

	// CCW edge information
	ccwedg = flowedg->getCCWEdg();
	orig_ccw = ccwedg->getOriginPtrNC()->get2DCoords();
	dest_ccw = ccwedg->getDestinationPtrNC()->get2DCoords();

	// Find the CW neighbor of flow edge by iterating around
	// Save voronoi vertices around flow edge for later
	curedg = flowedg->getCCWEdg();
	while ( curedg != flowedg )  {
		
		vv_cur = curedg->getRVtx();
		
		if (vv_cur[0] != vv_flow[0] || vv_cur[1] != vv_flow[1]) {
			// Eventually this will be the first CW voronoi vertex
			vv_cw2 = vv_cur; 

			// First time through this loop
			if (flag == 0) {
				// Save first distinct CCW voronoi vertex for later
				vv_ccw1 = vv_cur; 
				flag++;
			}

			// Second time through this loop
			else if (flagg == 0) {
				// Save the second distinct CCW voronoi vertex
				if (vv_cur[0] != vv_ccw1[0] || vv_cur[1] != vv_ccw1[1]) {
					vv_ccw3 = vv_cur; 
					flagg++;
				}
			}
		}
		// When loop exits this will be the CW edge
		cwedg = curedg;
		curedg = curedg->getCCWEdg();
	}
	
	// CW edge information
	orig_cw = cwedg->getOriginPtrNC()->get2DCoords();
	dest_cw = cwedg->getDestinationPtrNC()->get2DCoords();
	
	flag = flagg = 0;
	
	// ==========================================================
	//
	// Fix flow width of approximately zero
	//
	// ==========================================================
	if (flowedg->getVEdgLen() < THRESH) {
		
		// Check to see if the vertex vv_flow falls in the 
		// sector bounded by 'ccw_edg' & 'cw_edg'. If it 
		// does not -> the Voronoi geometry is complicated
		// and we just ignore correction of the Voronoi width
		
		if (IsInTriangle(centroid, dest_cw,   dest_flow, vv_flow[0], vv_flow[1]) ||
				IsInTriangle(centroid, dest_flow, dest_ccw,  vv_flow[0], vv_flow[1]) ||
				IsInTriangle(centroid, dest_cw,   dest_flow, vv_ccw1[0], vv_ccw1[1]) ||
				IsInTriangle(centroid, dest_flow, dest_ccw,  vv_ccw1[0], vv_ccw1[1])) {
			
			// Collect the 3,4 or 5 points which will define the sector area
			// Always the centroid and either 1 or 2 points on CW and CCW sides
			j = 0;

			// First vertex in the sector area (sArea) is the flow voronoi vertex
			for (i=0; i < 2; i++)
				polyg[i][j] = vv_flow[i];
			j++;
			
			// Now we have two lines that are defined by vertices 
			// vv_flow & vv_ccw1, vv_flow & vv_cw2, and two edges ccwedg & cwedg 
			// - neighbors of flowedg. We need to define where the
			// lines cross the edges and then find an area of a 
			// corresponding polygon. 
			
			// Find intersection of the first and second voronois with the CCW edge
			// If the intersection is in the triangle add one point to section area
			// otherwise add two points to section area
			intersect_ccw = FindIntersectionCoords(vv_flow, vv_ccw1, orig_ccw, dest_ccw);
			
			if (IsInTriangle(centroid, dest_flow, dest_ccw, 
					 intersect_ccw[0], intersect_ccw[1]) ) {

				// If that intersection is inside the triangle of centroid to dest node
				// of the flow edge to the dest node of the CCW edge
				if (intersect_ccw != vv_flow) {

					// And the intersection is not the flow voronoi vertex
					// then the intersection is the second point in the sector area
					for (i=0; i < 2; i++)
						polyg[i][j] = intersect_ccw[i];
					j++;
				}
			}
			else {
				// Otherwise the second voronoi vertex is the second point in the
				// sector area and we look for a third point
				for (i=0; i < 2; i++)
					polyg[i][j] = vv_ccw1[i];
				j++;
				
				// Find intersection of second and third voronoi vertices with CCW edge
				intersect_ccw = FindIntersectionCoords(vv_ccw1, vv_ccw3, orig_ccw, dest_ccw);

				if (IsInTriangle(centroid, dest_flow, dest_ccw,
						intersect_ccw[0], intersect_ccw[1])) {

					// If intersection falls in the triangle it is the third sector point
					for (i=0; i < 2; i++)
						polyg[i][j] = intersect_ccw[i];
					j++;
				}

				else if (IsInTriangle(centroid, dest_flow, dest_ccw, vv_ccw3[0], vv_ccw3[1])) {

					// Otherwise if the second voronoi vertex falls in triange it is
					// the third sector point
					for (i=0; i < 2; i++)
						polyg[i][j] = vv_ccw3[i];
					j++;
				}
			}
			// Now add the centroid as the fourth sector point
			for (i=0; i < 2; i++)
				polyg[i][j] = centroid[i];
			j++;
			
			// Wrap around and deal with the far side of the CW triangle
			// Find intersection of the first voronoi and the first CW voronoi
			// with the CW edge
			intersect_cw = FindIntersectionCoords(vv_flow, vv_cw2, orig_cw, dest_cw); 

			if (IsInTriangle(centroid, dest_cw, dest_flow,
					intersect_cw[0], intersect_cw[1]) ) {

				// If that intersection is inside the CW triangle
				if (intersect_cw != vv_flow) {

					// And the intersection is not the flow voronoi vertex
					// then the intersection is the fifth sector point
					for (i=0; i < 2; i++)
						polyg[i][j] = intersect_cw[i];
					j++;
				}
			}
			else {
				// If the intersection is outside the CW triangle
				if (vv_cw2 != vv_flow) {

					// And the last voronoi vertex is not the first voronoi vertex
					// then the last voronoi vertex is the fifth sector point
					for (i=0; i < 2; i++)
						polyg[i][j] = vv_cw2[i];
					j++;
				}
			}
			sArea = polygonArea(polyg, j);
			changeWidth++; // Flag is needed to skip the error portion of this method
		}
		else {
			if (simCtrl->Verbose_label == 'Y') {
				Cout<<"\ttFlowNet: Can not re-define zero flow width! "<<endl;
				Cout<<"\txy0 node is not within the nedg1 - nedg2 sector"<<endl;
				TellAboutNode(cn);
				Cout<<endl;
			}
		}
	}

	// ==========================================================
	//
	// Fix flow widths of size greater than 0 (all others)
	//
	// ==========================================================
	else {
		
		// Check to see if the voronoi vertex associated with the flow edge or the
		// voronoi vertex ccw from that fall in the sector bounded by ccwedg-cwedg
		// If they do not the Voronoi geometry is complicated and we skip the fix
		if (IsInTriangle(centroid, dest_cw,   dest_flow, vv_flow[0], vv_flow[1]) ||
			  IsInTriangle(centroid, dest_flow, dest_ccw,  vv_flow[0], vv_flow[1]) ||
			  IsInTriangle(centroid, dest_cw,   dest_flow, vv_ccw1[0], vv_ccw1[1]) ||
			  IsInTriangle(centroid, dest_flow, dest_ccw,  vv_ccw1[0], vv_ccw1[1])) {
			
			// Collect the 3,4 or 5 points which will define the sector area
			// Always the centroid and either 1 or 2 points on CW and CCW sides
			j = 0;
			
			// CW edge intersects line from flow voronoi vertex to next ccw vertex
			if (!AreSegmentsParallel(vv_flow, vv_ccw1, orig_cw, dest_cw)) {
				
				// Find intersection of vv_flow to vv_ccw1 with CW edge
				vv_tmp = FindIntersectionCoords(vv_flow, vv_ccw1, orig_cw, dest_cw);

				if ( IsBetweenEndPnts(vv_flow, vv_ccw1, orig_cw, dest_cw,
							vv_tmp[0], vv_tmp[1]) ) {
					// Intersection is within minimum bounding box of 4 endpoints
					flag = 1;
				} else {
					vv_tmp = vv_flow;
				}
			} else {
				vv_tmp = vv_flow;
			}

			for (i=0; i < 2; i++)
				polyg[i][j] = vv_tmp[i];
			j++;
			
			// CCW edge intersects line from flow voronoi vertex to next ccw vertex
			if (!AreSegmentsParallel(vv_flow, vv_ccw1, orig_ccw, dest_ccw)) {
				
				// Find intersection of vv_flow to vv_ccw1 with CCW edge
				vv_tmp  = FindIntersectionCoords(vv_flow, vv_ccw1, orig_ccw, dest_ccw);

				if ( IsBetweenEndPnts(vv_flow, vv_ccw1, orig_ccw, dest_ccw,
							vv_tmp[0], vv_tmp[1]) ) {
					// Intersection is within minimum bounding box of 4 endpoints
					flagg = 1;
				} else {
					vv_tmp = vv_ccw1;
				}
			} else {
				vv_tmp = vv_ccw1;
			}

			for (i=0; i < 2; i++)
				polyg[i][j] = vv_tmp[i];
			j++;
			
			// We didn't find a good intersection point on the CCW side
			// We'll add two points instead of one to the sector area polygon
			if (flagg == 0) {

				// Find intersection of CCW edge segment with next pair voronoi verts
				intersect_ccw = FindIntersectionCoords(vv_ccw1, vv_ccw3, orig_ccw, dest_ccw);

				// Add the intersection point to the sector area polygon
				for (i=0; i < 2; i++)
					polyg[i][j] = intersect_ccw[i];
				j++;
			}

			// Add the centroid to the sector area polygon
			for (i=0; i < 2; i++)
				polyg[i][j] = centroid[i];
			j++;

			// We've wrapped around to the CW side again so pick up the extra point
			// if the initial intersection point was not usable
			if (flag == 0) {

				// Find intersection of CW edge segment with flow vv and cw vv
				intersect_cw = FindIntersectionCoords(vv_flow, vv_cw2, orig_cw, dest_cw);
				// Add intersection to sector area polygon
				for (i=0; i < 2; i++)
					polyg[i][j] = intersect_cw[i];
				j++;
			}

			// Sector polygon was defined with 3 points if both intersections were
			// good or with 5 points if both were bad
			sArea = polygonArea(polyg, j);
			changeWidth++; // Flag is needed to skip the error portion of this method
		}
		// If neither vv_flow or vv_ccw1 is within the sector of CCW to CW then we
		// aren't changing the width (changeWidth flag is still 0) and we won't
		// execute the next block of code
	}
	
	// ==========================================================
	//
	// If we DID check the flow width for a current Voronoi cell,
	// then proceed with the following:
	//
	// ==========================================================
	
	// Just check of a situation when Voronoi cell is somewhat strange...
	// If we enter this code we will not change the voronoi width
	if (cn->getID() == -999 || ((changeWidth > 0) && (sArea > 1.0E+6)) ) {
		
		firstedg = cn->getFlowEdg();
		curedg = firstedg->getCCWEdg();
		
		Cout<<"\nID = \n"<<cn->getID();
		Cout <<"\twidth 1 = "<<firstedg->getVEdgLen()*1000.0;
		Cout <<"\twidth 2 = "<<curedg->getVEdgLen()*1000.0;

		vv_cur = firstedg->getRVtx();
		Cout<<"Voronoi Vertex = [ "<<vv_cur[0]<<" " << vv_cur[1] << "]" << endl;

		while ( curedg != firstedg )  {
			vv_cur = curedg->getRVtx();
			Cout<<"Voronoi Vertex = [ " <<vv_cur[0]<<" " << vv_cur[1] << "]" << endl;
			curedg = curedg->getCCWEdg();
		}
		
		curedg = firstedg;
		Cout<<"Slope: " << curedg->getSlope()<<"\t";
		TellAboutNode((tCNode *)curedg->getDestinationPtrNC());
		Cout<<"       -9999999\t";
		TellAboutNode(cn);
		curedg = firstedg->getCCWEdg();

		while ( curedg != firstedg )  {
			Cout<<"Slope: " << curedg->getSlope()<<"\t";
			TellAboutNode((tCNode *)curedg->getDestinationPtrNC());
			Cout<<"       -9999999\t";
			TellAboutNode(cn);
			curedg = curedg->getCCWEdg();
		}
		
		Cout<<"\nvv_flow = "<<vv_flow[0]<<", "<<vv_flow[1]<<endl;
		Cout<<"vv_ccw1 = "<<vv_ccw1[0]<<", "<<vv_ccw1[1]<<endl;
		Cout<<"vv_cw2 = "<<vv_cw2[0]<<", "<<vv_cw2[1]<<endl;
		Cout<<"vv_ccw3 = "<<vv_ccw3[0]<<", "<<vv_ccw3[1]<<endl;
		
		Cout<<"intersect_ccw = "<<intersect_ccw[0]<<", "<<intersect_ccw[1]<<endl;
		Cout<<"centroid = "<<centroid[0]<<", "<<centroid[1]<<endl;
		Cout<<"intersect_cw = "<<intersect_cw[0]<<", "<<intersect_cw[1]<<endl;
		
		if (changeWidth > 0) {
			NV = j;
			Cout<<"changeWidth = "<<changeWidth<<"; flagg = "<<flagg<<"; j = "<<j<<endl;
			for (i=0; i < 2; i++) {
				for (j=0; j < NV; j++)
					Cout<<polyg[i][j]<<" ";
				Cout<<endl;
			}
		}
		changeWidth = 0;
	}
	
	if ( changeWidth > 0 ) {
		if (flowedg->getLength() > 0.0)
			wWidth = sArea/flowedg->getLength();
		else 
			wWidth = -99999; // To ensure we dont get value
		
		if (flowedg->getVEdgLen() < wWidth) {
			// Output comparative information:
			//Cout<<"Node = " << cn->getID() << "\tsArea = "<<sArea
			//    <<"\tWIDTH BEFORE =  "<<flowedg->getVEdgLen()
			//    <<"\tWIDTH AFTER = "<<wWidth<<endl;
			flowedg->setVEdgLen(wWidth);
		}
	}

	for (i=0; i < 2; i++) 
		delete [] polyg[i];
	delete [] polyg;

	return;
}

/***************************************************************************
**
** Function:  polygonArea
** Arguments: *poly[2] - coordinates of NOT closed polygon
**            - poly[0][] - X coordinates
**            - poly[1][] - Y coordinates
** Objective: calculates polygon area: polygon is 
**            considered to be NON-closed
** Return value: - area of a polygon
** Algorithm: [O'Rourke], page 21 
**
***************************************************************************/
double tFlowNet::polygonArea(double **poly, int L) 
{
	double sum = 0.0;
	for (int i=0; i < L; i++) {
		if (i == (L-1))
			sum += (poly[0][i]*poly[1][0] - poly[1][i]*poly[0][0]);
		else 
			sum += (poly[0][i]*poly[1][i+1] - poly[1][i]*poly[0][i+1]);
	}
	return(fabs(sum/2.));
}

/*************************************************************************
**
**  IsBetweenEndPnts
**
**  The function checks if the point (x,y) falls between the endpoints 
**  of two line segments defined by xy1(x1,y1) - xy2(x2,y2) and 
**  xy3(xx1,yy1) - xy4(xx2,yy2).  Returns '1' if yes, '0' otherwise
**
*************************************************************************/
int tFlowNet::IsBetweenEndPnts(tArray< double > &xy1, tArray< double > &xy2,
							   tArray< double > &xy3, tArray< double > &xy4,
							   double x, double y)
{
	int result;
	double minx1, minx2, maxx1, maxx2;
	double miny1, miny2, maxy1, maxy2;
	double x1,x2,xx1,xx2,y1,y2,yy1,yy2;
	
	x1  = xy1[0];
	x2  = xy2[0];
	xx1 = xy3[0];
	xx2 = xy4[0];
	y1  = xy1[1];
	y2  = xy2[1];
	yy1 = xy3[1];
	yy2 = xy4[1];
	
	if (x1 < x2) {
		minx1 = x1;
		maxx1 = x2;
	}
	else { 
		minx1 = x2;  // <-- Can be a point if x1 = x2;
		maxx1 = x1;
	}
	if (xx1 < xx2) {
		minx2 = xx1;
		maxx2 = xx2;
	}
	else { 
		minx2 = xx2;  // <-- Can be a point if xx1 = xx2;
		maxx2 = xx1;
	}
	if (y1 < y2) {
		miny1 = y1;
		maxy1 = y2;
	}
	else { 
		miny1 = y2;  // <-- Can be a point if y1 = y2;
		maxy1 = y1;
	}
	if (yy1 < yy2) {
		miny2 = yy1;
		maxy2 = yy2;
	}
	else { 
		miny2 = yy2;  // <-- Can be a point if yy1 = yy2;
		maxy2 = yy1;
	}
	// In order to deal with the numerical issues...
	if (fabs(minx1-maxx1) < THRESH) {
		minx1 -= THRESH;
		maxx1 += THRESH;
	}
	if (fabs(minx2 - maxx2) < THRESH) {
		minx2 -= THRESH;
		maxx2 += THRESH;
	}
	if (fabs(miny1 - maxy1) < THRESH) {
		miny1 -= THRESH;
		maxy1 += THRESH;
	}
	if (fabs(miny2 - maxy2) < THRESH) {
		miny2 -= THRESH;
		maxy2 += THRESH;
	}
	
	if ((x >= minx1 && x <= maxx1) &&
		(x >= minx2 && x <= maxx2) &&
		(y >= miny1 && y <= maxy1) &&
		(y >= miny2 && y <= maxy2))
		result = 1;
	else
		result = 0;
	return result;
}

/*************************************************************************
**
**  AreSegmentsParallel
**
**  The function checks if the two segments defined by xy1(x1,y1)-xy2(x2,y2) 
**  and xy3(xx1,yy1)-xy4(xx2,yy2) are parallel. Returns '1' if yes, '0' 
**  otherwise
**
*************************************************************************/
int tFlowNet::AreSegmentsParallel(tArray< double > &xy1, tArray< double > &xy2,
								  tArray< double > &xy3, tArray< double > &xy4)
{
	int result;
	double dxa, dxb, dya, dyb, det;
	double x1,x2,xx1,xx2,y1,y2,yy1,yy2;
	x1  = xy1[0];
	x2  = xy2[0];
	xx1 = xy3[0];
	xx2 = xy4[0];
	y1  = xy1[1];
	y2  = xy2[1];
	yy1 = xy3[1];
	yy2 = xy4[1];
	
	// Check if  segments exist
	if ( ((x2 == x1) && (y1 == y2)) || 
		 ( (xx1 == xx2) && (yy1 == yy2)) ) {
		Cout<<"Segment doesn't exist:"<<endl;
		Cout<<"X1 = "<<x1<<"; Y1 = "<<y1<<"; X2 = "<<x2<<"; Y2 = "<<y2
			<<"XX1 = "<<xx1<<"; YY1 = "<<yy1
			<<"; XX2 = "<<xx2<<"; YY2 = "<<yy2<<endl;
		result = 0;
	}
	else {
		dxa = x2 - x1;
		dxb = xx2 - xx1;
		dya = y2 - y1;
		dyb = yy2 - yy1;
		// Check if the lines are not parallel
		det = dya*dxb-dyb*dxa;
		if (det !=0 )
			result = 0;
		else
			result = 1;
	}
	return result;
}

/***************************************************************************
**
**  IsInTriangle()
**
**  Defines if a point (x,y) falls in the triangle. The algorithm exploits
**  the fact that the 3 triangle points are always in counter-clockwise
**  order, so that the point is contained within a given triangle (p0,p1,p2)
**  if and only if the point lies to the left of vectors p0->p1, p1->p2,
**  and p2->p0. Here's how it works:
** 
***************************************************************************/
int tFlowNet::IsInTriangle(tArray< double > &xyp1,
                           tArray< double > &xyp2,
                           tArray< double > &xyp3,
						   double x, double y) 
{
	int k;
	double a, b, c;
	
	k = 1;
	for (int i=0; (i<3)&&(k>0) ; i++) {
		
		if (i == 0) {
			a = (xyp1[1] - y) * (xyp2[0] - x);
			b = (xyp1[0] - x) * (xyp2[1] - y);
			c = a - b;
		}
		else if (i == 1) {
			a = (xyp2[1] - y) * (xyp3[0] - x);
			b = (xyp2[0] - x) * (xyp3[1] - y);
			c = a - b;
		}
		else if (i == 2) {
			a = (xyp3[1] - y) * (xyp1[0] - x);
			b = (xyp3[0] - x) * (xyp1[1] - y);
			c = a - b;
		}
		
		if ( c > 0.0 )     // <--- Not to the LEFT
			k = -1;
		else { 
			if ( c == 0.0 )  // <--- on the BND
				;
		}
	}
	if (k > 0)
		return 1;
	else 
		return 0;
}
	
/***************************************************************************
**  
** tFlowNet::SetReachInformation() Function
**  
** Set the reach ID for each node
**
***************************************************************************/
 
void tFlowNet::SetReachInformation()
{
	// Assign reaches for stream nodes
	tCNode *creach, *coutlet, *chead, *cn;
	tPtrListIter< tCNode > HeadIter( NodesLstH );
	tPtrListIter< tCNode > OutletIter( NodesLstO );
	tMeshListIter<tCNode> niter( gridPtr->getNodeList() );

	int reach;
	for (chead = HeadIter.FirstP(), coutlet = OutletIter.FirstP(), reach = 0;
			!(HeadIter.AtEnd());
			chead = HeadIter.NextP(), coutlet = OutletIter.NextP(), reach++) {

		creach = chead;
		while (creach != coutlet) {
			creach->setReach(reach);
			creach = creach->getDownstrmNbr();
		} 

		// Set final reach
		if (reach == NodesLstH.getSize() - 1)
			coutlet->setReach(reach);
	}

	// Set all the other node reaches using the stream nodes
	for ( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() ) {
		if (cn->getBoundaryFlag() != kStream) {
			reach = cn->getStreamNode()->getReach();
			cn->setReach(reach);
		}
	}
}

/***************************************************************************
**  
** tFlowNet::writeRestart() Function
**  
** Called from tSimulator during simulation loop
**
***************************************************************************/
 
void tFlowNet::writeRestart(fstream & rStr) const
{
  BinaryWrite(rStr, flowboxes);
  BinaryWrite(rStr, hillvel);
  BinaryWrite(rStr, streamvel);
  BinaryWrite(rStr, velratio);
  BinaryWrite(rStr, velcoef);
  BinaryWrite(rStr, flowexp);
  BinaryWrite(rStr, baseflow);
  BinaryWrite(rStr, dOtp); 
  BinaryWrite(rStr, timespan);
  BinaryWrite(rStr, flowout);
  BinaryWrite(rStr, maxttime);
  BinaryWrite(rStr, dist_hill_max);
  BinaryWrite(rStr, dist_stream_max);
  BinaryWrite(rStr, BasArea);
 
  res->writeRestart(rStr);
} 

/***************************************************************************
**
** tFlowNet::readRestart() Function
**  
***************************************************************************/
                                        
void tFlowNet::readRestart(fstream & rStr)
{   
  BinaryRead(rStr, flowboxes);
  BinaryRead(rStr, hillvel);
  BinaryRead(rStr, streamvel);
  BinaryRead(rStr, velratio);
  BinaryRead(rStr, velcoef);
  BinaryRead(rStr, flowexp);
  BinaryRead(rStr, baseflow);
  BinaryRead(rStr, dOtp); 
  BinaryRead(rStr, timespan);
  BinaryRead(rStr, flowout);
  BinaryRead(rStr, maxttime);
  BinaryRead(rStr, dist_hill_max); 
  BinaryRead(rStr, dist_stream_max);
  BinaryRead(rStr, BasArea);
                                         
  res->readRestart(rStr);
}   

//=========================================================================
//
//
//                           End of tFlowNet 
//
//
//=========================================================================
