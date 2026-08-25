/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/

/***************************************************************************
**
**  tMesh.cpp: Functions for class tMesh (see tMesh.h) based on CHILD
**             routines for TIN Mesh Generation
**
***************************************************************************/

#include <cstdint>
#include <unordered_map>
#include "src/tMesh/tMesh.h"
#include "src/Headers/globalIO.h"
#include "src/tCNode/tCNode.h"

#ifdef PARALLEL_TRIBS
#include "src/tParallel/tParallel.h"
#endif

//=========================================================================
//
//
//                  Section 0: tIdArray Class
//
//
//=========================================================================

// tIdArray: lookup table per Id for a tList

//Class Definition
template< class T >
class tIdArray{
	tArray< T* > e_;
public:
	tIdArray(tList< T >& List);
	T* operator[]( int subscript ) const {
		// return a value and not a reference, hence the "const".
		return e_[subscript];
	}
};

//Class Constructor
template< class T >
tIdArray< T >::tIdArray(tList< T >& List) :
e_(List.getSize())
{
	tListIter< T > Iter( List );
	T *c;
	for( c=Iter.FirstP(); !(Iter.AtEnd()); c=Iter.NextP() )
		e_[c->getID()] = c;
}


//=========================================================================
//
//
//                  Section 1: Templated Global Functions
//
//
//=========================================================================

/***************************************************************************
**
**  Next3Delaunay
**
**  global function; determines whether nbr node currently pointed to
**  by iterator and the next two in the nbr list form a Delaunay triangle.
**
**  Inputs:  nbrList -- list of pointers to nodes
**           nbrIter -- iterator for this list
**  Returns: 1 if the next 3 nodes on the list are Delaunay, 0 otherwise
**  Called by: tMesh::RepairMesh
**
***************************************************************************/

template< class tSubNode >
int Next3Delaunay( tPtrList< tSubNode > &nbrList,
                   tPtrListIter< tSubNode > &nbrIter )
{
	static int ncalls = 0;
	ncalls++;
	tSubNode *cn, *nbrnd;
	
	//assert( (&nbrList != 0) && (&nbrIter != 0) ); //WR--09192023:reference cannot be bound to dereferenced null pointer in well-defined C++ code; comparison may be assumed to always evaluate to true
	
	nbrnd = nbrIter.DatPtr();
	tPtrListIter< tSubNode > nbrIterCopy( nbrList );
	int i = nbrIter.Where();
	nbrIterCopy.Get(i);
	
	tArray< double > p0( nbrIterCopy.DatPtr()->get2DCoords() );
	tArray< double > p1( nbrIterCopy.NextP()->get2DCoords() );
	tArray< double > p2( nbrIterCopy.NextP()->get2DCoords() );
	
	// If points aren't counter-clockwise, we know it's not Delaunay
	if( !PointsCCW( p0, p1, p2 ) ) return 0;
	
	tArray< double > ptest;
	cn = nbrIterCopy.NextP();  // Move to next point in the ring
	while( cn != nbrnd )       // Keep testing 'til we're back to p0
	{
		ptest = cn->get2DCoords();
		if( !TriPasses( ptest, p0, p1, p2 ) ){
			return 0;
		}
		cn = nbrIterCopy.NextP();  // Next point in ring
	}
	return 1;
}

/***************************************************************************
**
**  PointAndNext2Delaunay
**
**  Global function that determines whether nbr node currently pointed to
**  by iterator and the next two in the nbr list form a Delaunay triangle.
**  Similar to Next3Delaunay but p2 is an arbitrary node (testNode) rather
**  than one of the neighbor list nodes.
**
**  Inputs:  testNode -- a node to check (this is "p2")
**           nbrList -- list of pointers to nodes
**           nbrIter -- iterator for this list
**  Returns: 1 if they are Delaunay, 0 otherwise
**
***************************************************************************/

template< class tSubNode >
int PointAndNext2Delaunay( tSubNode &testNode, tPtrList< tSubNode > &nbrList,
                           tPtrListIter< tSubNode > &nbrIter )
{
	assert( (&nbrList != 0) && (&nbrIter != 0) && (&testNode != 0) );
	
	tPtrListIter< tSubNode > nbrIterCopy( nbrList );
	int i = nbrIter.Where();
	nbrIterCopy.Get( i );
	assert( nbrIterCopy.DatPtr() == nbrIter.DatPtr() );
	
	tArray< double > p0( nbrIterCopy.DatPtr()->get2DCoords() );
	assert( nbrIterCopy.Next() );
	tArray< double > p1( nbrIterCopy.DatPtr()->get2DCoords() );
	tArray< double > p2( testNode.get2DCoords() );
	
	// If the points aren't CCW then we know it's not Delaunay
	if( !PointsCCW( p0, p1, p2 ) ) return 0;
	
	// Otherwise, call TriPasses to compare
	tArray< double > ptest;
	assert( nbrIterCopy.Next() );
	while( nbrIterCopy.DatPtr() != nbrIter.DatPtr() ){
		ptest = nbrIterCopy.DatPtr()->get2DCoords();
		if( !TriPasses( ptest, p0, p1, p2 ) ){
			return 0;
		}
		assert( nbrIterCopy.Next() );
	}
	return 1;
}

//=========================================================================
//
//
//             Section 2: tMesh Default Constructors and Destructor
//
//
//=========================================================================

// DEFAULT CONSTRUCTOR 
template< class tSubNode >   
tMesh< tSubNode >:: tMesh() 
{
	nnodes = nedges = ntri = seed = 0;
	mSearchOriginTriPtr=0;
	miNextNodeID = miNextEdgID = miNextTriID = 0;
	// layerflag = FALSE; (Layering off in tRIBS)
}

template< class tSubNode >    
tMesh< tSubNode >:: tMesh(SimulationControl *simCtrPtr) :
	nnodes(0),
	nedges(0), 
	ntri(0), 
	nodeList(),
	seed(0),
	miNextNodeID(0),
	miNextEdgID(0),
	miNextTriID(0),   
	mSearchOriginTriPtr(0)
{
	simCtrl = simCtrPtr;
}

// COPY CONSTRUCTOR
template< class tSubNode >
tMesh<tSubNode>::tMesh( tMesh *originalMesh )
{
	nnodes = originalMesh->nnodes;
	nedges = originalMesh->nedges;
	ntri = originalMesh->ntri;
	nodeList = originalMesh->nodeList;
	edgeList = originalMesh->edgeList;
	triList = originalMesh->triList;
	seed = originalMesh->seed;
	// layerflag = originalMesh->layerflag; (Layering off in tRIBS)
	miNextNodeID = originalMesh->miNextNodeID;
	miNextEdgID = originalMesh->miNextEdgID;
	miNextTriID = originalMesh->miNextTriID;   
	mSearchOriginTriPtr=0;
}

// DESTRUCTOR 
template< class tSubNode >
tMesh< tSubNode >:: ~tMesh() {
	Cout << "tMesh Object has been destroyed..." << endl;
}                    

//=========================================================================
//
//
//                Section 3: tMesh( infile ) Constructor
//
//
//=========================================================================

/**************************************************************************
**
**   tMesh( infile ): Reads from infile whether it is to reconstruct a mesh
**                    from input, construct a mesh from a list of points,
**                    or other options including from Arc/Info files.
**
**   Calls: tInputFile::ReadItem, MakeMeshFromInputData( infile ),
**      MakeMeshFromPoints( infile ), MakeRandomPointsFromArcGrid( infile ), 
**      MakeHexMeshFromArcGrid( infile ), MakePointFromFileArcInfo( infile ), 
**      MakePointFromFileArcInfoGen( infile ), MakeMeshFromScratch( infile ),
**	MakeLayersFromInputData( infile )
**               
**   Options in MESHINPUT: 
**	1	Create mesh by reading data files (*.edges, *.nodes, *.tri)
**      2       Create new mesh from list of points (*.points)
**      3       Create random mesh from regular arc ascii grid
**      4       Create hex mesh from regular arc ascii grid
**      5       Create points file from Arc/Info Ungeneratetin (*.net)
**      6       Create points file from Arc/Info Ungeneratetin (*.lin & *.pnt)
**      7 	Create mesh from scratch using parameters
**
**************************************************************************/

template< class tSubNode >
tMesh< tSubNode >::
tMesh( SimulationControl *simCtrPtr, tInputFile &infile ) :
nnodes(0),
nedges(0), 
ntri(0), 
nodeList(),
seed(0),
miNextNodeID(0),
miNextEdgID(0),
miNextTriID(0),   
mSearchOriginTriPtr(0)
{
	simCtrl = simCtrPtr;
	
	// Read in MESHINPUT Option
	int read = infile.ReadItem( read, "OPTMESHINPUT" );  
	
	if( read < 1 || read > 2 ){
		cerr << "\nInvalid mesh input option requested.";
		cerr << "\nValid options for reading mesh input are:\n"
			<< "  1 -- read mesh from input data files\n"
			<< "  2 -- create mesh from points file using Triangulator\n\n";
		exit(1);
	}
	
	if( read == 1 ) {
		Cout<<"\n\nPart 2: Creating Mesh from Existing Mesh (Option 1)"<<endl;
		Cout<<"---------------------------------------------------"<<endl;      
		MakeMeshFromInputData( infile );
		
		// Layering Off in tRIBS 
		// ---------------------
		// int lay = infile.ReadItem( lay, "OPTREADLAYER" );
		// if( lay == 1 ) MakeLayersFromInputData( infile );
		Cout<<"\nMakeMeshFromInputData Successful Using Option 1"<<endl<<flush;
	}
	
	else if( read == 2){
		Cout<<"\n\nPart 2: Creating Mesh from Points File (Option 2)"<<endl;
		Cout<<"---------------------------------------------------"<<endl;
		MakeMeshFromTriangulator( infile ); 
		Cout<<"\nMakeMeshFromTriangulator Successful Using Option 2"<<endl<<flush;
	}
	
	// Layering Flag off in tRIBS
	// ----------------------------
	// int lflag = infile.ReadItem( lflag, "OPTINTERPLAYER" );
	// if(lflag > 0) layerflag = TRUE;
	// else layerflag = FALSE;
	
	// Find Geometric Center of domain:
	
	double cx = 0.0;
	double cy = 0.0;
	double sumarea = 0.0;
	double carea;
	tMeshListIter< tSubNode > nI( getNodeList() );
	tNode* cn;
	
	for( cn = nI.FirstP(); !nI.AtEnd(); cn = nI.NextP() ){
		carea = cn->getVArea();
		cx += cn->getX() * carea;
		cy += cn->getY() * carea; 
		sumarea += carea;
	}
	
	assert( sumarea>0.0 );
	cx /= sumarea;
	cy /= sumarea;
	
	// Find triangle in which these coordinates lie 
	// Designate it the search origin:
	
	mSearchOriginTriPtr = LocateTriangle( cx, cy );
	
}

//=========================================================================
//
//
//                  Section 4: tMesh:: MakeMeshFromInputData( )
//
//
//=========================================================================


/**************************************************************************
**
**   tMesh::MakeMeshFromInputData
**
**   Constructs tListInputData object and makes mesh from data in that object.
**                    
**   Calls: tListInputData( infile ), UpdateMesh(), CheckMeshConsistency()
**   Inputs: infile -- main input file from which various items are read
**                     
**************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromInputData( tInputFile &infile )
{
#ifdef PARALLEL_TRIBS
	tParallel::barrier();
#endif
	int i;
	tListInputData< tSubNode > input( infile );
	seed = 0;
	nnodes = input.x.getSize();
	nedges = input.orgid.getSize();
	ntri = input.p0.getSize();
	
	assert( nnodes > 0 );
	assert( nedges > 0 );
	assert( ntri > 0 );
	
	// Create the node list by creating a temporary node and iteratively
	// (1) assigning it values from the input data and (2) inserting it
	// onto the back of the node list.
	
	Cout << "\nCreating node list..." << endl;
	
	tSubNode tempnode( infile );
	int bound;
	for( i = 0; i< nnodes; i++ )
	{
		tempnode.set3DCoords( input.x[i], input.y[i], input.z[i] );
		tempnode.setID( i );
		bound = input.boundflag[i];
		assert( bound >= 0 && bound <= 3 );
		
		tempnode.setBoundaryFlag( bound );
		if( (bound == 0) || (bound==3) )
			nodeList.insertAtActiveBack( tempnode );
		else if( bound == kOpenBoundary )
			nodeList.insertAtBoundFront( tempnode );
		else
			nodeList.insertAtBack( tempnode );       //kClosedBoundary
	}
	
	//Initialize tIdArray Object
	
	const tIdArray< tSubNode > NodeTable(nodeList); // for fast lookup per ID
	
	// Create and initialize the edge list by creating two temporary edges
	// and then iteratively assigning values to the pair and inserting them 
	// onto the back of the edgeList
	
	Cout << "\nCreating edge list..." << endl;
	{
		//tMeshListIter< tSubNode > nodIter( nodeList );
		tEdge tempedge1, tempedge2;
		int obnd, dbnd;
		for( miNextEdgID = 0; miNextEdgID < nedges-1; miNextEdgID+=2 ){
			// Assign values: ID, origin and destination pointers
			tempedge1.setID( miNextEdgID );
			tempedge2.setID( miNextEdgID + 1 );
			{
				tSubNode *nodPtr1 = NodeTable[ input.orgid[miNextEdgID] ];
				tempedge1.setOriginPtr( nodPtr1 );
				tempedge2.setDestinationPtr( nodPtr1 );
				obnd = (*nodPtr1).getBoundaryFlag();
			}
			{
				tSubNode *nodPtr2 = NodeTable[ input.destid[miNextEdgID] ];
				tempedge1.setDestinationPtr( nodPtr2 );
				tempedge2.setOriginPtr( nodPtr2 );
				dbnd = (*nodPtr2).getBoundaryFlag();
			}
			
			// set the "flowallowed" status (FALSE if either endpoint is a
			// closed boundary, or both are open boundaries) 
			// and insert edge pair onto the list --- active
			// part of list if flow is allowed, inactive if not
			
			if( obnd == kClosedBoundary || dbnd == kClosedBoundary
				|| (obnd==kOpenBoundary && dbnd==kOpenBoundary) )
			{ 
				tempedge1.setFlowAllowed( 0 );
				tempedge2.setFlowAllowed( 0 );
				edgeList.insertAtBack( tempedge1 );
				edgeList.insertAtBack( tempedge2 );
			}
			else
			{
				tempedge1.setFlowAllowed( 1 );
				tempedge2.setFlowAllowed( 1 );
				edgeList.insertAtActiveBack( tempedge1 );
				edgeList.insertAtActiveBack( tempedge2 );
			}
		}
	}
	
	// Set up the lists of edges (spokes) connected to each node
	Cout << "\nSetting up spoke lists..." << endl;
	const tIdArray< tEdge > EdgeTable(edgeList); // for fast lookup per ID
	
	// set up the lists of edges (spokes) connected to each node
	// (GT added code to also assign the 1st edge to "edg" as an alternative
	// to spokelist implementation)
	
	{
		tMeshListIter< tSubNode > nodIter( nodeList );
		assert( nodIter.First() );
		do
		{
			tSubNode * curnode;
			curnode = nodIter.DatPtr();
			const int e1 = input.edgid[curnode->getID()];  
			
			tEdge *edgPtr = EdgeTable[e1];
			curnode->insertBackSpokeList( edgPtr );
			curnode->setEdg( edgPtr );
			
			int ne;
			for( ne = input.nextid[e1]; ne != e1; ne = input.nextid[ne] )
			{
				if( ne>=nedges )
				{
					cerr << "Warning: edge " << e1 
					<< " has non-existant ccw edge "
					<< ne << endl;
					cerr << "This is likely to be a problem in the edge input file"
						<< endl;
				}
				tEdge *edgPtr = EdgeTable[ne];
				curnode->insertBackSpokeList( edgPtr );
			}
		}
		while( nodIter.Next() );
	}
	
	// Assign ccwedg connectivity that tells each edge about its neighbor
	// immediately counterclockwise
	
	Cout << "\nSetting up CCW edges..." << endl;
	
	{
		tMeshListIter< tEdge > edgIter( edgeList );
		tMeshListIter< tSubNode > nodIter( nodeList );
		tEdge * curedg, * ccwedg;
		int ccwedgid;
		tMeshListIter<tEdge> ccwIter( edgeList ); // 2nd iter for performance
		for( i=0, curedg=edgIter.FirstP(); i<nedges; i++, curedg=edgIter.NextP() )
		{
			ccwedgid = input.nextid[i];
			ccwedg = EdgeTable[ccwedgid]; //test
			curedg->setCCWEdg( ccwedg );
		}
	}
	
	Cout << "\nSetting up triangle connectivity..." << endl;
	
	{
		tMeshListIter< tEdge > edgIter( edgeList );
		tMeshListIter< tSubNode > nodIter( nodeList );
		for ( i=0; i<ntri; i++ )
		{
			tTriangle newtri;
			newtri.setID( i );
			{
				newtri.setPPtr( 0, NodeTable[ input.p0[i] ] );
				newtri.setPPtr( 1, NodeTable[ input.p1[i] ] );
				newtri.setPPtr( 2, NodeTable[ input.p2[i] ] );
			}
			{
				newtri.setEPtr( 0, EdgeTable[ input.e0[i] ] );
				newtri.setEPtr( 1, EdgeTable[ input.e1[i] ] );
				newtri.setEPtr( 2, EdgeTable[ input.e2[i] ] );
			}
			triList.insertAtBack( newtri );
		}
		const tIdArray< tTriangle > TriTable(triList); // for fast lookup per ID
		
		tListIter< tTriangle > triIter( triList );
		tTriangle * ct, * nbrtri;
		for( i=0, ct=triIter.FirstP(); i<ntri; ct=triIter.NextP(), i++ )
		{
			nbrtri = ( input.t0[i]>=0 ) ? TriTable[ input.t0[i] ] : 0;
			ct->setTPtr( 0, nbrtri );
			nbrtri = ( input.t1[i]>=0 ) ? TriTable[ input.t1[i] ] : 0;
			ct->setTPtr( 1, nbrtri );
			nbrtri = ( input.t2[i]>=0 ) ? TriTable[ input.t2[i] ] : 0;
			ct->setTPtr( 2, nbrtri );
		}
	}
	
	Cout<<"\nTesting Mesh..."<<endl;
	UpdateMesh();
	CheckMeshConsistency();
	
}

//=========================================================================
//
//
//                  Section 5: tMesh:: MakeMeshFromPoints( )
//
//
//=========================================================================

/**************************************************************************
**
**   tMesh::MakeMeshFromTriangulator( infile )
**
**   Similar to tMesh::MakeMeshFromPoints but uses Tipper's triangulation
**   algorithm.
**
**   Created: 07/2002, Arnaud Desitter, Greg Tucker, Oxford
**   Modified: 08/2002, MIT
**
**************************************************************************/

// edge numbering translation
static inline int e_t2c(int ei, bool o){ // Tipper to child
	return o? 2*ei : 2*ei+1;
}
static inline int e_t2c(const oriented_edge &oe){
	return e_t2c(oe.e(), oe.o());
}

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromTriangulator( tInputFile &infile ){
#ifdef PARALLEL_TRIBS
	tParallel::barrier();
#endif
	int i, numpts;                      // no. of points in mesh
	tArray<double> x, y, z;          // arrays of x, y, and z coordinates
	tArray<int> bnd;                 // array of boundary codes 
	char pointFilenm[kMaxNameLength+kMaxExt];            // name of file containing (x,y,z,b) data
	ifstream pointfile;              // the file (stream) itself
	
	tMeshListIter< tSubNode > nodIter( nodeList );   //Node List
	tSubNode tempnode( infile );  // temporary node used to create node list
	
	//Read Points
	infile.ReadItem( pointFilenm, "POINTFILENAME" );
	pointfile.open( pointFilenm );
	if( !pointfile.good() ){
		cout << "\nPoint file name: '" << pointFilenm << "' not found\n";
		exit(1);    
	}
	
	Cout<<"\nReading in '"<<pointFilenm<<"' points file..."<<endl;
	pointfile >> numpts;
	x.setSize( numpts );
	y.setSize( numpts );
	z.setSize( numpts );
	bnd.setSize( numpts );
	
	//Read point file, make Nodelist 
	for( i=0; i<numpts; i++ ){
		if( pointfile.eof() )
			cout << "\nReached end-of-file while reading points.\n" ;
		pointfile >> x[i] >> y[i] >> z[i] >> bnd[i];
		tempnode.set3DCoords( x[i], y[i], z[i]);
		tempnode.setBoundaryFlag( bnd[i] );
		if( bnd[i]<0 || bnd[i]>3 ){
			cout << "\nInvalid boundary code.\n"<<endl;
			cout << "\n\nExiting Program..."<<endl;
			exit(2);
		}
		tempnode.setID( i );
		
		if(bnd[i]==kNonBoundary || bnd[i]==kStream)
			nodeList.insertAtActiveBack( tempnode );
		else if( bnd[i]== kOpenBoundary )
			nodeList.insertAtBoundFront( tempnode );
		else
			nodeList.insertAtBack( tempnode );
		unsortList.insertAtBack( tempnode );
	}
	
	pointfile.close();
	
	nnodes = nodeList.getSize();
	point *p = new point[nnodes];   // for Tipper triangulator 
	{
		tMeshListIter< tSubNode > nodIter(nodeList);
		tSubNode* cn;
		int inode = 0;
		for( cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP()){
			p[inode].x = cn->getX();
			p[inode].y = cn->getY();
			p[inode].id = cn->getID();
			++inode;
		}
	}
	
	const tIdArray< tSubNode > NodeTable(nodeList); // for fast lookup per ID
	
	// call mesh generator based on Tipper's method
	Cout << "\nComputing triangulation..." << flush;
	int nedgesl;
	int nelem;
	edge* edges(0);
	elem* elems(0);
	tt_sort_triangulate(nnodes, p, &nedgesl, &edges, &nelem, &elems);
	
	// set sizes
	nedges = 2*nedgesl;
	ntri = nelem;
	
	// Create and initialize the edge list by creating two temporary edges
	// (which are complementary, ie share the same endpoints) and then
	// iteratively assigning values to the pair and inserting them onto the
	// back of the edgeList
	
	Cout << "\nCreating edge list..." << endl<<flush;
	{
		for( int iedge = 0; iedge < nedgesl; ++iedge ) {
			tEdge tempedge1, tempedge2;
			int obnd, dbnd;
			
			// Assign values: ID, origin and destination pointers
			tempedge1.setID( e_t2c(iedge,true) );
			tempedge2.setID( e_t2c(iedge,false) );
			{
				tSubNode *nodPtr1 = NodeTable[p[edges[iedge].from].id];
				tempedge1.setOriginPtr( nodPtr1 );
				tempedge2.setDestinationPtr( nodPtr1 );
				obnd = (*nodPtr1).getBoundaryFlag();
			}
			{
				tSubNode *nodPtr2 = NodeTable[p[edges[iedge].to].id];
				tempedge1.setDestinationPtr( nodPtr2 );
				tempedge2.setOriginPtr( nodPtr2 );
				dbnd = (*nodPtr2).getBoundaryFlag();
			}
			
			// set the "flowallowed" status (FALSE if either endpoint is a
			// closed boundary, or both are open boundaries) 
			// and insert edge pair onto the list --- active
			// part of list if flow is allowed, inactive if not
			if( obnd == kClosedBoundary || dbnd == kClosedBoundary
				|| (obnd==kOpenBoundary && dbnd==kOpenBoundary) )
			{
				tempedge1.setFlowAllowed( 0 );
				tempedge2.setFlowAllowed( 0 );
				edgeList.insertAtBack( tempedge1 );
				edgeList.insertAtBack( tempedge2 );
			}
			else
			{
				tempedge1.setFlowAllowed( 1 );
				tempedge2.setFlowAllowed( 1 );
				edgeList.insertAtActiveBack( tempedge1 );
				edgeList.insertAtActiveBack( tempedge2 );
			}
		}
	}
	const tIdArray< tEdge > EdgeTable(edgeList); // for fast lookup per ID
	
	// set up the lists of edges (spokes) connected to each node
	Cout << "\nSetting up spoke lists..." << endl<<flush;
	{
		// connectivity point - sorted point
		tArray< int > p2sp(nnodes);
		for(int inodes=0;inodes!=nnodes;++inodes){
			p2sp[p[inodes].id] = inodes;
		}
		
		tMeshListIter< tSubNode > nodIter(nodeList);
		oriented_edge *oedge;
		tt_build_spoke(nnodes, nedgesl, edges, &oedge);
		
		tSubNode * curnode;
		assert( nodIter.First() );
		do
		{
			// first spoke
			curnode = nodIter.DatPtr();
			{
				const int e1 = e_t2c(oedge[p2sp[curnode->getID()]]);
				tEdge *edgPtr = EdgeTable[e1];
				curnode->insertBackSpokeList( edgPtr );
				curnode->setEdg( edgPtr );
			}
			// build rest of spoke list
			const oriented_edge& oe_ref = oedge[p2sp[curnode->getID()]];
			oriented_edge ccw_from = oe_ref.ccw_edge_around_from(edges);
			while( ccw_from.e() != oe_ref.e()) {
				assert(ccw_from.e() < nedgesl);
				const int ne = e_t2c(ccw_from);
				tEdge *edgPtr = EdgeTable[ne];
				curnode->insertBackSpokeList( edgPtr );
				ccw_from = ccw_from.ccw_edge_around_from(edges);
			}
		}
		while( nodIter.Next() );
		delete [] oedge;
	}
	
	// Assign ccwedg connectivity (that is, tell each edge about its neighbor
	// immediately counterclockwise)
	Cout << "\nSetting up CCW edges..." << endl<<flush;
	{
		int iedge;
		tEdge * curedg, * ccwedg;
		tMeshListIter< tEdge > edgIter( edgeList );
		for( iedge=0, curedg=edgIter.FirstP(); iedge<nedgesl; ++iedge)
		{
		{
			const oriented_edge e1(iedge,true);
			const oriented_edge ccw_from = e1.ccw_edge_around_from(edges);
			const int ccwedgid = e_t2c(ccw_from);
			ccwedg = EdgeTable[ccwedgid];
			curedg->setCCWEdg( ccwedg );
		}
			curedg = edgIter.NextP();
			{
				const oriented_edge e2(iedge,false);
				const oriented_edge ccw_to = e2.ccw_edge_around_from(edges);
				const int ccwedgid = e_t2c(ccw_to);
				ccwedg = EdgeTable[ccwedgid];
				curedg->setCCWEdg( ccwedg );
			}
			curedg = edgIter.NextP(); 
		}
	}
	
	Cout << "\nSetting up triangle connectivity..." <<endl<< flush;
	{
		int ielem;
		for ( ielem=0; ielem<nelem; ++ielem ) {
			
			tTriangle newtri;
			newtri.setID( ielem );
			
			{
				newtri.setPPtr( 0, NodeTable[p[elems[ielem].p1].id] );
				newtri.setPPtr( 1, NodeTable[p[elems[ielem].p2].id] );
				newtri.setPPtr( 2, NodeTable[p[elems[ielem].p3].id] );
			}
			{
				newtri.setEPtr( 0, EdgeTable[e_t2c(elems[ielem].e1, 
												   elems[ielem].eo1)] );
				newtri.setEPtr( 1, EdgeTable[e_t2c(elems[ielem].e2,
												   elems[ielem].eo2)] );
				newtri.setEPtr( 2, EdgeTable[e_t2c(elems[ielem].e3,
												   elems[ielem].eo3)] );
			}
			triList.insertAtBack( newtri );
		}
		const tIdArray< tTriangle > TriTable(triList); // for fast lookup per ID
		
		tTriangle * ct, * nbrtri;
		tListIter< tTriangle > triIter( triList );
		for( ielem=0, ct=triIter.FirstP(); ielem<nelem; ct=triIter.NextP(), 
			 ++ielem ) {
			nbrtri = ( elems[ielem].t1>=0 ) ? TriTable[ elems[ielem].t1 ] : 0;
			ct->setTPtr( 0, nbrtri );
			nbrtri = ( elems[ielem].t2>=0 ) ? TriTable[ elems[ielem].t2 ] : 0;
			ct->setTPtr( 1, nbrtri );
			nbrtri = ( elems[ielem].t3>=0 ) ? TriTable[ elems[ielem].t3 ] : 0;
			ct->setTPtr( 2, nbrtri );
		}
	}   
	
	// deallocation of Tipper triangulator data structures
	delete [] edges;
	delete [] elems;
	delete [] p;
	
	// assertions
	assert( edgeList.getSize() == 2*nedgesl );
	assert( triList.getSize() == nelem );
	
	UpdateMesh(); //calls CheckMeshConsistency()  
	CheckMeshConsistency();                    
	
}

//=========================================================================
//
//
//                  Section 9: tMesh:: CheckMeshConsistency( )
// 					ChangePointOrder( )
//
//=========================================================================

/*****************************************************************************
**
**  tMesh::CheckMeshConsistency
**
**  Performs a series of tests to make sure the mesh connectivity is correct.
**  Should be called immediately after reading in a user-defined mesh 
**
**  The consistency checks include the following:
**
**  1) Each edge:
**     - Has valid origin and destination pointers
**     - Has a valid counter-clockwise edge, which shares the same origin but
**       not the same destination
**     - Is paired with its complement in the list
**
**  2) Each node:
**     - Points to a valid edge which has the node as its origin
**     - If the node is not a boundary, it has at least one neighbor that
**       is not a closed boundary (unless boundaryCheckFlag is FALSE).
**     - Has a consistent spoke list (ie, you can go around the spokes and
**       get back to where you started)
**
**  3) Each triangle:
**     - Has 3 valid points and edges
**     - Each edge Ei has Pi as its origin and P((i+2)%3) as its
**       destination
**     - If an opposite triangle Ti exists, points P((i+1)%3) and
**       P((i+2)%3) are the same as points PO((n+2)%3) and PO((n+1)%3) in
**       the opposite triangle, where PO denotes a point in the opposite
**       triangle and n is the vertex ID (0, 1, or 2) of the point in the
**       opposite triangle that is opposite from the shared face.
**     - If an opposite triange Ti does not exist, points P((i+1)%3) and
**       and P((i+2)%3) should both be boundary points.
**
**      Parameters:  boundaryCheckFlag -- defaults to TRUE; if FALSE,
**                                        node connection to open node or
**                                        open boundary isn't tested
**
*****************************************************************************/

template<class tSubNode>
void tMesh< tSubNode >::
CheckMeshConsistency( int boundaryCheckFlag )
{   
	tMeshListIter<tSubNode> nodIter( nodeList );
	tMeshListIter<tEdge> edgIter( edgeList );
	tListIter<tTriangle> triIter( triList );
	tPtrListIter< tEdge > sIter;
	tNode * cn, * org, * dest;
	tEdge * ce, * cne, * ccwedg;
	tTriangle * ct, * optr;
	int boundary_check_ok, i, nvop;
	int kMaxSpokes = 100;
	
	// Edges: make sure complementary pairs are together in the list
	// (each pair Ei and Ei+1, for i=0,2,4,...nedges-1, should have the same
	// endpoints but the opposite orientation)
	
	for( ce=edgIter.FirstP(); !(edgIter.AtEnd()); ce=edgIter.NextP() ){
		cne = edgIter.NextP();
		if( ce->getOriginPtrNC() != cne->getDestinationPtrNC()
			|| ce->getDestinationPtrNC() != cne->getOriginPtrNC() ){
			cerr << "EDGE #" << ce->getID()
			<< " must be followed by its complement in the list\n";
			goto error;
		}
		
	}
	
	// Edges: check for valid origin, destination, and ccwedg
	for( ce=edgIter.FirstP(); !(edgIter.AtEnd()); ce=edgIter.NextP() ){
		if( !(org=ce->getOriginPtrNC() ) ){
			cerr << "EDGE #" << ce->getID()
			<< " does not have a valid origin point\n";
			goto error;
		}
		if( !(dest=ce->getDestinationPtrNC() ) ){
			cerr << "EDGE #" << ce->getID()
			<< " does not have a valid destination point\n";
			goto error;
		}
		if( !(ccwedg=ce->getCCWEdg() ) ){
			cerr << "EDGE #" << ce->getID()
			<< " does not point to a valid counter-clockwise edge\n";
			goto error;
		}
		if( ccwedg->getOriginPtrNC()!=org ){
			cerr << "EDGE #" << ce->getID()
			<< " points to a CCW edge with a different origin\n";
			goto error;
		}
		if( ccwedg->getDestinationPtrNC()==dest ){
			cerr << "EDGE #" << ce->getID()
			<< " points to a CCW edge with the same destination\n";
			goto error;
		}
		if( org==dest ){
			cerr << "EDGE #" << ce->getID()
			<< " has the same origin and destination nodes\n";
			goto error;
		}   
	}
	
	// Nodes: check for valid edg pointer, spoke connectivity, and connection
	// to at least one non-boundary or open boundary node
	
	for( cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP() ){
		// edg pointer
		if( !(ce = cn->getEdg()) ){
			cerr << "NODE #" << cn->getID()
			<< " does not point to a valid edge\n";
			goto error;
		}
		if( ce->getOriginPtrNC()!=cn ){
			cerr << "NODE #" << cn->getID()
			<< " points to an edge that has a different origin\n";
			goto error;
		}
		
		boundary_check_ok = ( cn->getBoundaryFlag()==kNonBoundary &&
							  boundaryCheckFlag ) ? 0 : 1;
		i = 0;
		// Loop around the spokes until we're back at the beginning
		do{
			
			if( ce->getDestinationPtrNC()->getBoundaryFlag()!=kClosedBoundary )
				boundary_check_ok = 1;  // OK, there's at least one open nbr
			i++;
			if( i>kMaxSpokes ){
				cerr << "NODE #" << cn->getID()
				<< ": infinite loop in spoke connectivity\n";
				goto error;
			}
			
			// Make sure node is the origin --- and not the destination
			if( ce->getOriginPtrNC()!=cn ){
				cerr << "EDGE #" << ce->getID()
				<< " is in the spoke chain of NODE " << cn->getID()
				<< " but does not have the node as an origin\n";
				goto error;
			}
			if( ce->getDestinationPtrNC()==cn ){
				cerr << "EDGE #" << ce->getID()
				<< " is in the spoke chain of NODE " << cn->getID()
				<< " but has the node as its destination\n";
				goto error;
			}   
			
		} while( (ce=ce->getCCWEdg())!=cn->getEdg() );
		
		if( !boundary_check_ok ){ 
			tArray<double> x;
			x = cn->get2DCoords();
			cerr << "NODE #" << cn->getID()
				<<" ( "<<x[0]<< " , "<<x[1]<<" )"
				<< " is surrounded by closed boundary nodes\n";
			goto error;
		}
		
		//make sure node coords are consistent with edge endpoint coords:
		sIter.Reset( cn->getSpokeListNC() );
		for( ce = sIter.FirstP(); !(sIter.AtEnd()); ce = sIter.NextP() ){
			if( ce->getOriginPtrNC()->getX() != cn->getX() ||
				ce->getOriginPtrNC()->getY() != cn->getY() ){
				cerr << "NODE #" << cn->getID()
				<< " coords don't match spoke origin coords\n";
				goto error;
			}
		}
		
	}
	
	// Triangles: check for valid points and connectivity
	
	for( ct=triIter.FirstP(); !(triIter.AtEnd()); ct=triIter.NextP() ){
		for( i=0; i<=2; i++ ){
			// Valid point i?
			if( !(cn=ct->pPtr(i)) ){
				cerr << "TRIANGLE #" << ct->getID()
				<< " has an invalid point " << i << endl;
				goto error;
			}
			// Valid edge i?
			if( !(ce=ct->ePtr(i)) ){
				cerr << "TRIANGLE #" << ct->getID()
				<< " has an invalid edge " << i << endl;
				goto error;
			}
			// Edge and point consistency
			if( ce->getOriginPtrNC()!=cn ){
				cerr << "TRIANGLE #" << ct->getID()
				<< ": edge " << i << " does not have point " << i
				<< " as origin\n";
				goto error;
			}
			// changed from (i+1) to (i+2) for "right-hand" format
			if( ce->getDestinationPtrNC()!=ct->pPtr((i+2)%3) ){
				cerr << "TRIANGLE #" << ct->getID()
				<< ": edge " << i << " does not have point " << (i+1)%3
				<< " as destination\n";
				goto error;
			}
			// Opposite triangle: if it exists, check common points
			if( (optr = ct->tPtr(i)) ){
				nvop = optr->nVOp(ct); // Num (0,1,2) of opposite vertex in optr
				if( nvop < 3 ){
					if( ct->pPtr((i+1)%3) != optr->pPtr((nvop+2)%3)
						|| ct->pPtr((i+2)%3) != optr->pPtr((nvop+1)%3) ){
						cerr << "TRIANGLE #" << ct->getID()
						<< ": opposite triangle " << i << " does not share nodes "
						<< (ct->pPtr((i+1)%3))->getID() << " and "
						<< (ct->pPtr((i+2)%3))->getID() << endl;
						goto error;
					}
				}
				else{
					cerr << "TRIANGLE #" << ct->getID()
                    << ": opposite triangle " << i << ", triangle #"
                    << optr->getID() << ",does not have current tri as neighbor\n";
					goto error;
				}
			}
			// If no opposite triangle, make sure it really is a boundary
			else{
				if( (ct->pPtr((i+1)%3))->getBoundaryFlag()==kNonBoundary
					|| (ct->pPtr((i+2)%3))->getBoundaryFlag()==kNonBoundary )
				{
					cerr << "TRIANGLE #" << ct->getID()
                    << ": there is no neighboring triangle opposite node "
                    << cn->getID() << " but one (or both) of the other nodes "
                    << "is a non-boundary point."
					<<"\nX = "<<cn->getX()<<"\tY = "<<cn->getY()
					<<"\tZ = "<<cn->getZ()<<endl;
					goto error;
				}
			}       
		}
	}
	return;
	
error:
		cerr<<"Error in mesh consistency.";
	
}

/***************************************************************************
**
**  tMesh::CheckMeshConsistency(tInputFile)
**          see also previous routine  CheckMeshconsistency()
**
**  Perform the same checking that the previous CheckMeshconsitency()
**  In the case of real data from dems, the constructed mesh may have
**  some points nodes surrounded by closed boundary nodes. That comes
**  from the succession of the points in the input file, those points
**  have to be at the end of the file read by MakeMeshFromPoints and 
**  may be at the second run be removed.
**
****************************************************************************/

template<class tSubNode>
void tMesh< tSubNode >::CheckMeshConsistency( tInputFile &infile, int boundaryCheckFlag ){
	tMeshListIter<tSubNode> nodIter( nodeList );
	tMeshListIter<tEdge> edgIter( edgeList );
	tListIter<tTriangle> triIter( triList );
	tPtrListIter< tEdge > sIter;
	tNode * cn, * org, * dest;
	tEdge * ce, * cne, * ccwedg;
	tTriangle * ct, * optr;
	int boundary_check_ok, i, nvop;
	int kMaxSpokes = 100;
	tList<double> xy;
	
	for( ce=edgIter.FirstP(); !(edgIter.AtEnd()); ce=edgIter.NextP() ){
		cne = edgIter.NextP();
		if( ce->getOriginPtrNC() != cne->getDestinationPtrNC()
			|| ce->getDestinationPtrNC() != cne->getOriginPtrNC() ){
			cerr << "EDGE #" << ce->getID()
			<< " must be followed by its complement in the list\n";
		}
	}
	
	// Edges: check for valid origin, destination, and ccwedg
	
	for( ce=edgIter.FirstP(); !(edgIter.AtEnd()); ce=edgIter.NextP() ){
		if( !(org=ce->getOriginPtrNC() ) ){
			cerr << "EDGE #" << ce->getID()
			<< " does not have a valid origin point\n";
		}
		if( !(dest=ce->getDestinationPtrNC() ) ){
			cerr << "EDGE #" << ce->getID()
			<< " does not have a valid destination point\n";
		}
		if( !(ccwedg=ce->getCCWEdg() ) ){
			cerr << "EDGE #" << ce->getID()
			<< " does not point to a valid counter-clockwise edge\n";
		}
		if( ccwedg->getOriginPtrNC()!=org ){
			cerr << "EDGE #" << ce->getID()
			<< " points to a CCW edge with a different origin\n";
		}
		if( ccwedg->getDestinationPtrNC()==dest ){
			cerr << "EDGE #" << ce->getID()
			<< " points to a CCW edge with the same destination\n";
		}
		if( org==dest ){
			cerr << "EDGE #" << ce->getID()
			<< " has the same origin and destination nodes\n";
		}   
	}
	
	// Triangles: check for valid points and connectivity
	
	for( ct=triIter.FirstP(); !(triIter.AtEnd()); ct=triIter.NextP() ){
		for( i=0; i<=2; i++ ){
			// Valid point i?
			if( !(cn=ct->pPtr(i)) ){
				cerr << "TRIANGLE #" << ct->getID()
				<< " has an invalid point " << i << endl;
			}
			if( !(ce=ct->ePtr(i)) ){
				cerr << "TRIANGLE #" << ct->getID()
				<< " has an invalid edge " << i << endl;
			}
			if( ce->getOriginPtrNC()!=cn ){
				cerr << "TRIANGLE #" << ct->getID()
				<< ": edge " << i << " does not have point " << i
				<< " as origin\n";
			}
			if( ce->getDestinationPtrNC()!=ct->pPtr((i+2)%3) ){
				cerr << "TRIANGLE #" << ct->getID()
				<< ": edge " << i << " does not have point " << (i+1)%3
				<< " as destination\n";
			}
			if( (optr = ct->tPtr(i)) ){
				nvop = optr->nVOp(ct); // Num (0,1,2) of opposite vertex in optr
				if( nvop < 3 ){
					if( ct->pPtr((i+1)%3) != optr->pPtr((nvop+2)%3)
						|| ct->pPtr((i+2)%3) != optr->pPtr((nvop+1)%3) ){
						cerr << "TRIANGLE #" << ct->getID()
						<< ": opposite triangle " << i << " does not share nodes "
						<< (ct->pPtr((i+1)%3))->getID() << " and "
						<< (ct->pPtr((i+2)%3))->getID() << endl;
					}
				}
				else{
					cerr << "TRIANGLE #" << ct->getID()
                    << ": opposite triangle " << i << ", triangle #"
                    << optr->getID() << ", does not have current tri as neighbor\n";
				}
			}
			else{
				if( (ct->pPtr((i+1)%3))->getBoundaryFlag()==kNonBoundary
					|| (ct->pPtr((i+2)%3))->getBoundaryFlag()==kNonBoundary ){
					cerr << "TRIANGLE #" << ct->getID()
                    << ": there is no neighboring triangle opposite node "
                    << cn->getID() << " but one (or both) of the other nodes "
                    << "is a non-boundary point\n";
				}
			}       
		}
	}
	
	// Nodes: check for edges, boundary
	
	for( cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP() ){
		if( !(ce = cn->getEdg()) ){
			cerr << "NODE #" << cn->getID()
			<< " does not point to a valid edge\n";
		}
		if( ce->getOriginPtrNC()!=cn ){
			cerr << "NODE #" << cn->getID()
			<< " points to an edge that has a different origin\n";
		}
		
		boundary_check_ok = ( cn->getBoundaryFlag()==kNonBoundary &&
							  boundaryCheckFlag ) ? 0 : 1;
		i = 0;
		
		// Loop around the spokes until we're back at the beginning
		do{ 
			if( ce->getDestinationPtrNC()->getBoundaryFlag()!=kClosedBoundary )
				boundary_check_ok = 1;  // OK, there's at least one open nbr
			i++;
			if( i>kMaxSpokes ){
				cerr << "NODE #" << cn->getID()
				<< ": infinite loop in spoke connectivity"<<endl;}
			
			// Make sure node is the origin --- and not the destination
			if( ce->getOriginPtrNC()!=cn ){
				cerr << "EDGE #" << ce->getID()
				<< " is in the spoke chain of NODE " << cn->getID()
				<< " but does not have the node as an origin\n";
			}
			if( ce->getDestinationPtrNC()==cn ){
				cerr << "EDGE #" << ce->getID()
				<< " is in the spoke chain of NODE " << cn->getID()
				<< " but has the node as its destination\n";
			}           
		}while( (ce=ce->getCCWEdg())!=cn->getEdg() );
		
		if( !boundary_check_ok ){ 
			tArray<double> x;
			x= cn->get2DCoords();
			xy.insertAtBack( x[0] );
			xy.insertAtBack( x[1] );
			cout << "NODE #" << cn->getID()
				<<" ( "<<x[0]<< " , "<<x[1]<<" )"
				<< " is surrounded by closed boundary nodes\n";
		}
		
		//make sure node coords are consistent with edge endpoint coords:
		sIter.Reset( cn->getSpokeListNC() );
		for( ce = sIter.FirstP(); !(sIter.AtEnd()); ce = sIter.NextP() ){
			if( ce->getOriginPtrNC()->getX() != cn->getX() ||
				ce->getOriginPtrNC()->getY() != cn->getY() ){
				cerr << "NODE #" << cn->getID()
				<< " coords don't match spoke origin coords\n";
			}
		}     
	}
	
	if(!xy.isEmpty()){ 
		//xy.makeCircular();
		int check = ChangePointOrder(infile,xy); 
		if(check!=1)
			cout<<"Error in returning ChangePointOrder..."<<endl;
	}
	return;
}

/*****************************************************************************
**
**  tMesh: ChangePointOrder( tInputFile &, tList<double>)
**
**	Function created to change the order in a points file for those
**      interior points exclusively connected to exterior nodes. This function
**      becomes unnecessary if an inner ring is included in the TIN Mesh.
**       
*****************************************************************************/

template<class tSubNode>
int tMesh< tSubNode >::ChangePointOrder( tInputFile &infile, tList<double> XY ){
	tMeshListIter<tSubNode> nodIter( nodeList );
	tArray<double> p1,p2,p3,nod;
	tArray<int> p4;
	int i,b,k,n,nt;
	tNode *cn;
	char pointFilenm[80];
	
	Cout<<"\nChanging Point Order..."<<endl;
	
	//opening the file
	infile.ReadItem( pointFilenm, "POINTFILENAME" );
	ofstream outfile(strcat(pointFilenm,".corr"));
	
	//write in the file the number of points
	nt = nodeList.getSize();
	nod.setSize(3);
	
	//nb of points to be moved at the end of the file
	n = XY.getSize(); 
	p1.setSize( n );
	p2.setSize( n );
	p3.setSize( n );
	p4.setSize( n );
	
	outfile<<nt-n/2<<endl;  
	int flag;
	k = 0;
	for( cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP() ){  
		flag=0;
		nod = cn->get3DCoords();
		b = cn->getBoundaryFlag();
		if (XY.getSize()!=0){
			double * x=XY.FirstP();
			for (i=0; i<XY.getSize();i=i+2){
				if ( (nod[0]!=*x) || (nod[1]!=*(XY.NextP())) )  
					x=XY.NextP();
				else{
					flag=1;
					p1[k]=nod[0]; p2[k]=nod[1]; p3[k]=nod[2]; p4[k]=b;
					k++;
					double v;
					if(XY.getSize()>2){     
						x=XY.NextP();
						XY.removePrev(v,XY.getCurrentItem());
						XY.removePrev(v,XY.getCurrentItem());
						i=XY.getSize();
					}
					else{
						int a =XY.removeFromFront(v);
						a=XY.removeFromFront(v);
					}
				}
			}
		}
		if(flag!=1) 
			outfile<<setprecision(14)<<nod[0]<<" "<<nod[1]<<" "<<nod[2]<<" "<<b<<endl;
    } 
	
	outfile.close();
	return 1;
}

//=========================================================================
//
//
//             Section 11: tMesh Functions using MeshElements
//
//
//=========================================================================

template< class tSubNode >
void tMesh< tSubNode >::Print(){
	triList.print();
	nodeList.print();
	edgeList.print();
}

template< class tSubNode >
void tMesh< tSubNode >::MakeCCWEdges(){
	tMeshListIter< tSubNode > nodIter( nodeList );
	tSubNode *cn;
	for( cn = nodIter.FirstP(); !( nodIter.AtEnd() ); cn = nodIter.NextP() ){
		cn->makeCCWEdges();
	}
}

/*****************************************************************************
**
**  tMesh::setVoronoiVertices
**
**  Each Delaunay triangle is associated with an intersection between
**  three Voronoi cells, called a Voronoi vertex. These Voronoi vertices
**  are used in computing the area of each Voronoi cell. The Voronoi
**  vertex associated with each triangle is the circumcenter of the
**  triangle. This routine finds the Voronoi vertex associated with
**  each triangle by finding the triangle's circumcenter. 
**
**    Assumes: correct triangulation with valid edge pointers in each tri.
**    Data members modified: none
**    Other objects modified: Voronoi vertices set for each tEdge
**    Modifications:
**     - reverted to earlier triangle-based computation, from an edge-based
**       computation that takes 3x as long because NE = 3NT. In so doing,
**       the definition of the Voronoi vertex stored in a tEdge is changed
**       to "left-hand", meaning the V. vertex associated with the edge's
**       lefthand triangle (the vertex itself may or may not lie to the left
							**       of the edge). 1/98 GT
**     - also moved circumcenter computation into a tTriangle mbr fn.
**     - copied function to tMesh member from tStreamNet member, gt 3/98.
**       Other fns now use "right-hand" definition; this fn may have to
**       be changed.
**
*****************************************************************************/

template <class tSubNode>
void tMesh<tSubNode>::setVoronoiVertices(){
	tArray< double > xy;
	tListIter< tTriangle > triIter( triList );
	tTriangle * ct;
	
	// Find the Voronoi vertex associated with each Delaunay triangle
	
	for( ct = triIter.FirstP(); !(triIter.AtEnd()); ct = triIter.NextP() ){
		xy = ct->FindCircumcenter();    
		
		// Assign the Voronoi point as the left-hand point of the three edges 
		// associated with the current triangle
		
		ct->ePtr(0)->setRVtx( xy );
		ct->ePtr(1)->setRVtx( xy );
		ct->ePtr(2)->setRVtx( xy );
	}
}

/**************************************************************************
**
**  tMesh::CalcVoronoiEdgeLengths
**
**  Updates the length of the Voronoi cell edge associated with each
**  triangle edge. Because complementary edges are stored pairwise on
**  the edge list, we can save computation time by only computing the
**  vedglen once for the first of the pair, then assigning it to the
**  second. For boundary triangle edges, the corresponding Voronoi edge
**  is infinitely long, so the calculation is only done for interior
**  (active) edges.
**
**************************************************************************/

template <class tSubNode>
void tMesh<tSubNode>::CalcVoronoiEdgeLengths()
{
    tEdge *ce;
    double vedglen;
    tMeshListIter<tEdge> edgIter( edgeList );
	
    for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() ){
		vedglen = ce->CalcVEdgLen();     // Compute Voronoi edge length
		ce = edgIter.NextP();            // Advance to complement edge and
		ce->setVEdgLen( vedglen );       // assign the same edge length.
    }
}

/**************************************************************************
**
**  tMesh::CalcVAreas
**
**  Computes Voronoi area for each active (non-boundary) node in the
**  mesh (Voronoi area is only defined for interior nodes). Accomplishes
**  this by calling ComputeVoronoiArea for each node. 
**
**************************************************************************/

template <class tSubNode>
void tMesh<tSubNode>::CalcVAreas()
{
	tSubNode* curnode;
	tMeshListIter< tSubNode > nodIter( nodeList );
	
	for( curnode = nodIter.FirstP(); nodIter.IsActive();
		 curnode = nodIter.NextP() ){
		
		curnode->ComputeVoronoiArea(); 
	} 
}

//=========================================================================
//
//
//         Section 12: tMesh Functions related to Adding/Deleting Nodes
//
//
//=========================================================================

/**************************************************************************
**
**  tMesh::DeleteNode( tListNode<tSubNode> *, int =1 )
**    (see DeleteNode( tSubNode *, int =1 ) below)
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
DeleteNode( tListNode< tSubNode > *nodPtr, int repairFlag )
{
	if( !DeleteNode( nodPtr->getDataPtrNC(), repairFlag ) ) return 0;
	return 1;   
}

/**************************************************************************
**
**  tMesh::DeleteNode( tSubNode *, int =1 )
**
**  Deletes a node from the mesh. This is done by first calling
**  ExtricateNode to detach the node by removing its edges and their
**  associated triangles, then removing the node from the nodeList.
**  Normally, RepairMesh is then called to retriangulate the "hole" left
**  behind in the mesh. (However, if the node was on the hull of the
**  mesh there's no "hole" to fix --- the caller is assumed to be smart
**  enough to recognize this case and let us know about it by setting
**  repairFlag to kNoRepair. This is the case, for example, when deleting
**  the nodes that form a "supertriangle" as in MakeMeshFromPoints).
**
**  Once the mesh is repaired, the nodes are renumbered and as a safety
**  measure for debugging/testing purposes, UpdateMesh is called.
**
**  Data mbrs modified:  nnodes, nedges, and ntri are updated;
**                       the node is deleted from nodeList; other edges &
**                       triangles are removed and/or modified by
**                       ExtricateNode and RepairMesh (qv)
**  Calls:  tMesh::ExtricateNode, tMesh::RepairMesh, plus utility member
**               functions of tNode, tMeshList, etc. 
**  Returns:  error code: 0 if ExtricateNode or RepairMesh fails,
**            		  1 otherwise.
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
DeleteNode( tSubNode *node, int repairFlag ){
	tPtrList< tSubNode > nbrList;
	tListNode< tSubNode > *nodPtr;
	tMeshListIter< tSubNode > nodIter( nodeList );
	nodIter.Get( node->getID() );
	tSubNode nodeVal;
	
	nodPtr = nodIter.NodePtr();
	if( !( ExtricateNode( node, nbrList ) ) ) return 0;
	
	//tRIBS compatability for Stream Cells
	
	if((node->getBoundaryFlag()==2) || (node->getBoundaryFlag()==1)){
		nodeList.moveToBack( nodPtr );
		nodeList.removeFromBack( nodeVal );
	}
	else{
		nodeList.moveToFront( nodPtr );
		nodeList.removeFromFront( nodeVal );
	}
	
	nnodes = nodeList.getSize();
	nedges = edgeList.getSize();
	ntri = triList.getSize();
	
	if( repairFlag ){
		nbrList.makeCircular();
		if( !RepairMesh( nbrList ) ) return 0;
	}
	
	//reset node id's
	assert( nodIter.First() );
	miNextNodeID = 0;
	do{
		nodIter.DatRef().setID( miNextNodeID );
		miNextNodeID++;
	}
	while( nodIter.Next() );
	
	if( repairFlag ) UpdateMesh();
	return 1;
}


/**************************************************************************
**
**  tMesh::ExtricateNode
**
**  Detaches a node from the mesh by deleting all of its edges (which in
																**  turn removes the affected triangles). Returns a list of the node's
**  former neighbors by modifying the nbrList input parameter. Also
**  returns a code that indicates failure if the node still has a non-empty
**  spoke list after edge deletion.
**
**  Data mbrs modified:  nnodes; edges and triangles are removed from
**                       edgeList and triList
**  Calls:  tMesh::DeleteEdge and utility member functions of tNode,
**               tPtrList, tPtrListIter
**  Output:  list of node's (former) neighbors, in nbrList
**  Returns:  1 if all edges successfully deleted, 0 if not
**  Calls: DeleteEdge
**  Assumes:  
**  Notes:
**  Created: SL fall, '97
**  Modifications: if node is a closed boundary, any of its neighbors that
**             are non-boundaries are switched to closed boundaries, so
**             that nodes along the edge of the domain (including nodes of
														**             a "supertriangle" used in MakeMeshFromPoints) may be removed
**             without causing errors, GT 4/98
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
ExtricateNode( tSubNode *node, tPtrList< tSubNode > &nbrList )
{
	tPtrListIter< tEdge > spokIter( node->getSpokeListNC() );
	tEdge edgeVal1, edgeVal2, *ce;
	tSubNode *nbrPtr;
	
	for( ce = spokIter.FirstP(); !(spokIter.AtEnd()); ce = spokIter.FirstP() ){
		nbrPtr = ( tSubNode * ) ce->getDestinationPtrNC();
		nbrList.insertAtBack( nbrPtr );
		if( node->getBoundaryFlag()                      // If node is a bdy make
			&& nbrPtr->getBoundaryFlag()==kNonBoundary ) // sure nbrs are also
		{                                                // boundaries.
			nbrPtr->ConvertToClosedBoundary();
			nodeList.moveToBack( nbrPtr );
		}
		if( !DeleteEdge( ce ) ) return 0;
	}  
	
	if( node->getSpokeList().isEmpty() ) return 1;
	return 0;
}

/**************************************************************************
**
**  tMesh::DeleteEdge
**
**  Deletes an edge from the mesh, returning 1 if deletion succeeds and
**  0 if not. Starts by calling ExtricateEdge to detach the edge from
**  the other mesh elements. This function actually deletes two directed
**  edges: edgePtr and its complement.
**
**  Inputs:  edgePtr -- ptr to the edge to be deleted
**  Returns:  1 if successful, 0 if not
**  Calls: ExtricateEdge 
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
DeleteEdge( tEdge * edgePtr )
{
	tEdge edgeVal1, edgeVal2;
	
	// Detach the edge from other mesh elements
	if( !ExtricateEdge( edgePtr ) ) return 0;
	
	if( edgePtr->getBoundaryFlag() ){
		if( !( edgeList.removeFromBack( edgeVal1 ) ) ) return 0;
		if( !( edgeList.removeFromBack( edgeVal2 ) ) ) return 0;
	}
	else{
		if( !( edgeList.removeFromFront( edgeVal1 ) ) ) return 0;
		if( !( edgeList.removeFromFront( edgeVal2 ) ) ) return 0;
	}
	
	//if( &edgeVal1 == 0 || &edgeVal2 == 0 ) return 0;//WR--09192023:  comparison of address of 'edgeVal*' equal to a null pointer is always false
	return 1;
}


/**************************************************************************
**
**  tMesh::ExtricateEdge
**
**  Here we detach an edge and its complement from the surrounding mesh 
**  elements prior to deletion. Adjacent triangle(s) are also deleted
**  via a call to DeleteTriangle. Calls the virtual node function
**  WarnSpokeLeaving to signal the affected nodes to take appropriate
**  action. (Appropriate action might depend on the application; that's
**  why it is a virtual function that can be handled by any descendents
**  of tNode). The two complementary edges are then placed at the back
**  of the edge list, where DeleteEdge can conveniently find them.
**
**  Inputs:  edgePtr -- ptr to the edge to be deleted
**  Returns: 1 if successful, 0 otherwise
**  Calls: DeleteTriangle, <tSubNode>::WarnSpokeLeaving
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
ExtricateEdge( tEdge * edgePtr )
{
	assert( edgePtr != 0 );
	tEdge *tempedgePtr=0, *ce, *cce, *spk;
	tEdge *ceccw, *cceccw;
	tMeshListIter< tEdge > edgIter( edgeList );
	tPtrListIter< tEdge > spokIter;
	tPtrList< tEdge > *spkLPtr;
	tListNode< tEdge > *listnodePtr;
	tTriangle triVal1, triVal2;
	tArray< tTriangle * > triPtrArr(2);
	
	
	ce = edgIter.GetP( edgePtr->getID() );  
	spkLPtr = &( ce->getOriginPtrNC()->getSpokeListNC() );
	spokIter.Reset( *spkLPtr );
	for( spk = spokIter.FirstP(); spk != ce && !( spokIter.AtEnd() ); spk = spokIter.NextP() );
	if( spk == ce ){
		spk = spokIter.NextP();
		spkLPtr->removePrev( tempedgePtr, spokIter.NodePtr() );
	}
	
	// Find the triangle that points to the edge
	triPtrArr[0] = TriWithEdgePtr( edgePtr ); 
	
	// Find the edge's complement
	listnodePtr = edgIter.NodePtr();
	assert( listnodePtr != 0 );
	
	if( edgePtr->getID()%2 == 0 ) cce = edgIter.NextP();
	else if( edgePtr->getID()%2 == 1 ) cce = edgIter.PrevP();
	else return 0; //NB: why whould this ever occur??
	
	// Find the triangle that points to the edges complement
	triPtrArr[1] = TriWithEdgePtr( cce );
	
	if( triPtrArr[0] != 0 )
		if( !DeleteTriangle( triPtrArr[0] ) ) return 0;
	if( triPtrArr[1] != 0 )
		if( !DeleteTriangle( triPtrArr[1] ) ) return 0;
	
	// Update complement's origin's spokelist
	spkLPtr = &(cce->getOriginPtrNC()->getSpokeListNC());
	spokIter.Reset( *spkLPtr );
	for( spk = spokIter.FirstP(); spk != cce && !( spokIter.AtEnd() );
		 spk = spokIter.NextP() );
	if( spk == cce ){
		spk = spokIter.NextP();
		spkLPtr->removePrev( tempedgePtr, spokIter.NodePtr() );
	}
	
	tSubNode * nodece = (tSubNode *) ce->getOriginPtrNC();
	nodece->WarnSpokeLeaving( ce ); 
	tSubNode * nodecce = (tSubNode *) cce->getOriginPtrNC();
	nodecce->WarnSpokeLeaving( cce );
	
	//Take care of the edges who had as thier ccwedge ce or cce
	
	ceccw=ce->getCCWEdg();
	tempedgePtr=ceccw;
	do{
		tempedgePtr=tempedgePtr->getCCWEdg();
	}while(tempedgePtr->getCCWEdg() != ce);
	
	//Set tempedgeptrs ccwedge to ceccw
	tempedgePtr->setCCWEdg( ceccw);
	
	cceccw=cce->getCCWEdg();
	tempedgePtr=cceccw;
	do{
		tempedgePtr=tempedgePtr->getCCWEdg();
	}while(tempedgePtr->getCCWEdg() != cce);
	
	//Set tempedgeptrs ccwedge to cceccw
	tempedgePtr->setCCWEdg(cceccw);
	
	if( ce->getBoundaryFlag() )
	{
		//move edges to back of list
		edgeList.moveToBack( listnodePtr );
		edgeList.moveToBack( edgIter.NodePtr() );
	}
	else
	{
		//move edges to front of list
		edgeList.moveToFront( edgIter.NodePtr() );
		edgeList.moveToFront( listnodePtr );
	}
	nedges-=2;
	return 1;
}


/***************************************************************************
**
**  tMesh::LocateTriangle
**
**  Locates the triangle in which point (x,y) falls. The algorithm exploits
**  the fact that the 3 triangle points are always in counter-clockwise
**  order, so that the point is contained within a given triangle (p0,p1,p2)
**  if and only if the point lies to the left of vectors p0->p1, p1->p2,
**  and p2->p0. Here's how it works:
**   1 - start with a given triangle (currently first on the list, but a
**       smarter initial guess could be used -- TODO)
**   2 - lv is the number of successful left-hand checks found so far:
**       initialize it to zero
**   3 - check whether (x,y) lies to the left of p(lv)->p((lv+1)%3)
**   4 - if so, increment lv by one (ie, move on to the next vector)
**   5 - if not, (x,y) is to the right of the current face, so move to
**       the triangle that lies opposite that face and reset lv to zero
**   6 - continue steps 3-5 until lv==3, which means that we've found
**       our triangle.
**   7 - so far, a point "on the line", i.e., colinear w/ two of the
**       three points, still passes; that's OK unless that line is on
**       the boundary, so we need to check
**
**  Input: x, y -- coordinates of the point
**  Modifies: (nothing)
**  Returns: a pointer to the triangle that contains (x,y)
**  Assumes: the point is contained within one of the current triangles
**
***************************************************************************/

template< class tSubNode >
tTriangle * tMesh< tSubNode >::
LocateTriangle( double x, double y )
{
	int n, lv=0;
	tListIter< tTriangle > triIter( triList );  
	tTriangle *lt = ( mSearchOriginTriPtr != nullptr ) ? mSearchOriginTriPtr
												: triIter.FirstP(); //Updated to new c++ standards
	double a, b, c;
	int online = -1;
	tArray< double > xy1, xy2;
	
	for (n=0 ;(lv!=3)&&(lt); n++){
		xy1 = lt->pPtr(lv)->get2DCoords();
		xy2 = lt->pPtr( (lv+1)%3 )->get2DCoords();
		a = (xy1[1] - y) * (xy2[0] - x);
		b = (xy1[0] - x) * (xy2[1] - y);
		c = a - b;
		
		if ( c > 0.0 ){
			lt=lt->tPtr( (lv+2)%3 );
			lv=0;
			online = -1;
		}
		else{
			if( c == 0.0 ) online = lv;
			lv++;
		}
		
		assert( n < 3*ntri );
	}
	
	if( online != -1 )
		if( lt->pPtr(online)->getBoundaryFlag() != kNonBoundary &&
			lt->pPtr( (online+1)%3 )->getBoundaryFlag() != kNonBoundary )
			return 0;
	
	return(lt);
}


/**************************************************************************
**
**  tMesh::LocateNewTriangle
**
**  Called by: AddNodeAt
**
**************************************************************************/

template< class tSubNode >
tTriangle * tMesh< tSubNode >::
LocateNewTriangle( double x, double y )
{
	int n, lv=0;
	tListIter< tTriangle > triIter( triList ); 
	tTriangle *lt = triIter.FirstP();
	tSubNode *p1, *p2;
	
	tArray< double > xy1, xy2;
	
	for (n=0 ;(lv!=3)&&(lt); n++){
		p1 = (tSubNode *) lt->pPtr(lv);
		
		xy1 = p1->get2DCoords();
		p2 = (tSubNode *) lt->pPtr( (lv+1)%3 );
		
		xy2 = p2->get2DCoords();
		if ( ( (xy1[1] - y) * (xy2[0] - x) ) > ( (xy1[0] - x) * (xy2[1] - y)) ){
			lt=lt->tPtr( (lv+2)%3 );
			lv=0;
		}
		else {lv++;}  
	}
	return(lt);
}


/**************************************************************************
**
**  tMesh::TriWithEdgePtr
**
**  Finds and returns the triangle that points to edgPtr as one of its
**  clockwise-oriented edges.
**
**************************************************************************/

template< class tSubNode >
tTriangle *tMesh< tSubNode >::
TriWithEdgePtr( tEdge *edgPtr ){
	assert( edgPtr != 0 );
	tTriangle *ct;
	tListIter< tTriangle > triIter( triList );
	
	for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
		if( ct != 0 ) //TODO: is this test nec? why wd it be zero?
			if( ct->ePtr(0) == edgPtr ||
				ct->ePtr(1) == edgPtr ||
				ct->ePtr(2) == edgPtr ) return ct;
	return 0;
}

/**************************************************************************
**
**  tMesh::DeleteTriangle
**
**  Deletes a triangle from the mesh. Starts off with a call to 
**  ExtricateTriangle to detach the triangle from other mesh elements,
**  after which the triangle is at the front of the triangle list,
**  from whence it is then deleted.
**
**  Inputs:  triPtr -- ptr to the triangle to be deleted
**  Returns:  1 if successful, 0 if not
**  Calls: ExtricateTriangle
**  Called by: DeleteEdge, AddNode, AddNodeAt
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
DeleteTriangle( tTriangle * triPtr ){ 
	tTriangle triVal;
	
	if( !ExtricateTriangle( triPtr ) ) return 0;
	
	if( !(triList.removeFromFront(triVal) ) ){
		cerr << "DeleteTriangle(): triList.removeFromFront( triPtr ) failed\n";
		return 0;
	}
	
	//if( &triVal == 0 ) // Produces warning on ALPHA -VIVA
	//  return 0;
	
	return 1;
}

/**************************************************************************
**
**  tMesh::ExtricateTriangle
**
**  Detaches a triangle from surrounding mesh elements and places it at
**  the head of the triangle list, where it can be easily deleted by
**  DeleteTriangle.
**
**  Inputs: triPtr -- ptr to the triangle to be extricated
**  Returns: 1 if successful, 0 if not
**  Called by: DeleteTriangle
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
ExtricateTriangle( tTriangle *triPtr )
{
	tListIter< tTriangle > triIter( triList );
	tTriangle *ct;
	
	// Find the triangle on the list
	for( ct = triIter.FirstP(); ct != triPtr && !( triIter.AtEnd() );
		 ct = triIter.NextP() );
	if( ( triIter.AtEnd() ) ) return 0;
	
	int i, j;
	for( i=0; i<3; i++ ) for( j=0; j<3; j++ )
		if( triPtr->tPtr(i) != 0 )
			if( triPtr->tPtr(i)->tPtr(j) == triPtr ) 
				triPtr->tPtr(i)->setTPtr( j, 0 );
	
	// Move the triangle to the head of the list where it can be deleted
	triList.moveToFront( triIter.NodePtr() );
	
	ntri--;
	
	if( triPtr == mSearchOriginTriPtr ){
		mSearchOriginTriPtr = 0;
		for( i=0; i<3; ++i )
			if( triPtr->tPtr(i) != 0 )
				mSearchOriginTriPtr = triPtr->tPtr(i);
	}
	
	return 1;
}


/**************************************************************************
**
**  tMesh::RepairMesh
**
**  This function repairs the "hole" in the mesh that is created when
**  a node is deleted. Essentially, this function stiches the hole back
**  together by adding edges and triangles as needed, preserving
**  Delaunay-ness. The nodes around the hole are stored in the input
**  parameter nbrList. As each new triangle is added, its "interior"
**  point is removed from the neighbor list. The process of stitching
**  proceeds iteratively until the hole itself is a Delaunay triangle.
**
**  For each set of 3 successive counter-clockwise points along the rim
**  of the whole, the function calls Next3Delaunay to compare the
**  potential triangle p0, p1, p2 with other potential triangles
**  p0, p1, ptest (where ptest is each of the other nodes along the rim).
**  When a Delaunay triangle is found, AddEdgeAndMakeTriangle is called
**  to create the necessary edge and triangle objects, and the interior
**  node is removed from the neighbor list. (or at least that's what
**  gt is able to deduce...)
**
**  Inputs: nbrList -- list of nodes surrounding the "hole"
**  Returns: 1 if successful, 0 if not
**  Calls: Next3Delaunay, AddEdgeAndMakeTriangle, MakeTriangle
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
RepairMesh( tPtrList< tSubNode > &nbrList )
{
	//assert( &nbrList != 0 );//WR--09192023: reference cannot be bound to dereferenced null pointer in well-defined C++ code; comparison may be assumed to always evaluate to true
	if( nbrList.getSize() < 3 ) return 0;
	tSubNode * meshnodePtr = 0;
	nbrList.makeCircular();
	tPtrListIter< tSubNode > nbrIter( nbrList );
	
	// Keep stitching until only 3 nodes are left
	while( nbrList.getSize() > 3 ){
		if( Next3Delaunay( nbrList, nbrIter ) ) //checks for ccw and Del.
		{
			AddEdgeAndMakeTriangle( nbrList, nbrIter );
			//remove "closed off" pt
			nbrList.removeNext( meshnodePtr, nbrIter.NodePtr() );
		}
		
		nbrIter.Next();                    //step forward once in nbrList
	}
	assert( nbrList.getSize() == 3 );
	assert( ntri == triList.getSize() );
	assert( nedges == edgeList.getSize() );
	assert( nnodes == nodeList.getSize() );       //make sure numbers are right
	MakeTriangle( nbrList, nbrIter );             //make final triangle
	
	return 1;
}

/**************************************************************************
**
**  tMesh::AddEdge
**
**  Function to add edge pair between two nodes. Resets edge IDs.
**
**  Inputs: three nodes; edge is added between first two, and third
**   should be CCW 3rd member of triangle.
**  Returns: 1 if successful, 0 if not
**
**  Created: SL fall, '97
**  Modified: SL 10/98--routine sometimes failed when node1 (or node2) had
**   had edges to neither node2 (or node1) nor node3; to fix, replaced
**   the "assert( !( spokIter.AtEnd() ) )"'s with new algorithm to find
**   where new spoke should be inserted: finds where the sequence of 3 spoke
**   unit vectors, including the new one in the middle, are CCW; calls new
**   global function, tArray< double > UnitVector( tEdge* ).
**   - GT 1/99 -- to avoid compiler warning, now stores output of 
**     UnitVector calls in arrays p1, p2, p3, which are then sent as
**     arguments to PointsCCW.
**   - GT 2/99 -- added calls to WelcomeCCWNeighbor and AttachNewSpoke
**     to update CCW edge connectivity
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
AddEdge( tSubNode *node1, tSubNode *node2, tSubNode *node3 ) 
{
	assert( node1 != 0 && node2 != 0 && node3 != 0 );
	
	int flowflag = 1;                    // Boundary code for new edges
	tEdge tempEdge1, tempEdge2;          // The new edges
	tEdge *ce, *le;
	tMeshListIter< tEdge > edgIter( edgeList );
	tMeshListIter< tSubNode > nodIter( nodeList );
	tPtrListIter< tEdge > spokIter;
	tArray<double> p1, p2, p3;           // Used to store output of UnitVector
	
	// Set origin and destination nodes and find boundary status
	
	tempEdge1.setOriginPtr( node1 );                 //set edge1 ORG
	tempEdge2.setDestinationPtr( node1 );            //set edge2 DEST
	if( node1->getBoundaryFlag() == kClosedBoundary ) 
		flowflag = 0;
	tempEdge2.setOriginPtr( node2 );                 //set edge2 ORG
	tempEdge1.setDestinationPtr( node2 );            //set edge1 DEST
	if( node2->getBoundaryFlag() == kClosedBoundary ) 
		flowflag = 0;
	if( node1->getBoundaryFlag()==kOpenBoundary       // Also no-flow if both
		&& node2->getBoundaryFlag()==kOpenBoundary )  //  nodes are open bnds
		flowflag = 0;
	
	// Set boundary status and ID
	
	tempEdge1.setID( miNextEdgID );               //set edge1 ID
	miNextEdgID++;
	tempEdge2.setID( miNextEdgID );               //set edge2 ID
	miNextEdgID++;
	tempEdge1.setFlowAllowed( flowflag );         //set edge1 FLOWALLOWED
	tempEdge2.setFlowAllowed( flowflag );         //set edge2 FLOWALLOWED
	
	// Place new edge pair on the list: 
	// active back if not a boundary edge, back otherwise
	
	if( flowflag == 1 ){
		edgeList.insertAtActiveBack( tempEdge1 );    //put edge1 active in list
		edgeList.insertAtActiveBack( tempEdge2 );    //put edge2 active in list
		le = edgIter.LastActiveP();                  //set edgIter to lastactive
	}
	else{
		edgeList.insertAtBack( tempEdge1 );          //put edge1 in list
		edgeList.insertAtBack( tempEdge2 );          //put edge2 in list
		le = edgIter.LastP();                        //set edgIter to last
	}
	
	// Add pointers to the new edges to nodes' spokeLists
	// Three possible cases: (1) there aren't any spokes currently attached,
	// so just put the new one at the front of the list and make it circ'r;
	// (2) there is only one spoke, so it doesn't matter where we attach
	// (3) there is already >1 spoke (the general case)
	
	spokIter.Reset( node2->getSpokeListNC() );
    
	if( node2->getSpokeListNC().isEmpty() ){
		node2->insertFrontSpokeList( le);
		node2->getSpokeListNC().makeCircular();
		node2->AttachFirstSpoke( le ); // gt added to update ccwedg 2/99
	}
	else if( spokIter.ReportNextP() == spokIter.DatPtr() ){
		node2->insertFrontSpokeList( le);
		ce = node2->getEdg();  // these 2 lines added by gt 2/99
		assert( ce!=0 );
		ce->WelcomeCCWNeighbor( le );
	}
	else // general case: figure out where to attach spoke
	{
		for( ce = spokIter.FirstP();
			 ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
			 ce = spokIter.NextP() );
		
		if( spokIter.AtEnd() )
		{
			for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
			{
				p1 = UnitVector( ce );
				p2 = UnitVector( le );
				p3 = UnitVector( spokIter.ReportNextP() );
				if( PointsCCW( p1, p2, p3 ) )
					break;
			}
		}
		node2->getSpokeListNC().insertAtNext( le,
											  spokIter.NodePtr() ); 
		assert( ce!=0 );
		ce->WelcomeCCWNeighbor( le );
	}
	spokIter.Reset( node1->getSpokeListNC() );
	le = edgIter.PrevP();                     //step backward once in edgeList
    
	if( node1->getSpokeListNC().isEmpty() ){
		node1->insertFrontSpokeList( le );
		node1->getSpokeListNC().makeCircular();
		node1->AttachFirstSpoke( le ); // Tell node it's getting a spoke
	}
	else if( spokIter.ReportNextP() == spokIter.DatPtr() )
	{
		node1->insertFrontSpokeList( le );
		ce = node1->getEdg();  // these 2 lines added by gt 2/99
		assert( ce!=0 );
		ce->WelcomeCCWNeighbor( le ); // Tell node it has a new neighbor
	}
	else
	{
		for( ce = spokIter.FirstP();
			 ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
			 ce = spokIter.NextP() );
		if( spokIter.AtEnd() )
		{
			for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
			{
				p1 = UnitVector( ce );
				p2 = UnitVector( le );
				p3 = UnitVector( spokIter.ReportNextP() );
				if( PointsCCW( p1, p2, p3 ) )
				{
					spokIter.Next();
					break;
				}
			}
		}
		node1->getSpokeListNC().insertAtPrev( le,
											  spokIter.NodePtr() );
		assert( ce!=0 );
		ce->WelcomeCCWNeighbor( le );  // Tell node it has a new neighbor!
	}
	
	nedges+=2;
	
	// Reset edge id's
	for( ce = edgIter.FirstP(), miNextEdgID = 0; !( edgIter.AtEnd() ); 
		 ce = edgIter.NextP(), miNextEdgID++ ){
		ce->setID( miNextEdgID );
	}
	
	return 1;
}


/**************************************************************************
**
**  tMesh::AddEdgeAndMakeTriangle 
**
**  Function to add the "missing" edge and make
**  the triangle. Formerly more complicated than AddEdge() and
**  MakeTriangle(); now simply calls these functions.
**
**  Inputs: a tPtrList<tSubNode> of nodes in triangle; a tPtrListIter
**   object, the iterator of the latter list; edge is added between node
**   currently pointed to by iterator and the node-after-next. List should
**   be circular.
**  Calls: AddEdge(), MakeTriangle()
**  Created: SL fall, '97
**  Modified: SL 10/98 to call AddEdge() and MakeTriangle()
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
AddEdgeAndMakeTriangle( tPtrList< tSubNode > &nbrList,
                        tPtrListIter< tSubNode > &nbrIter )
{
	//assert( (&nbrList != 0) && (&nbrIter != 0) ); //WR--09192023: reference cannot be bound to dereferenced null pointer in well-defined C++ code; comparison may be assumed to always evaluate to true
	
	tSubNode *cn, *cnn, *cnnn;
	tPtrList< tSubNode > tmpList;
	tPtrListIter< tSubNode > tI( tmpList );
	
	cn = nbrIter.DatPtr();
	cnn = nbrIter.NextP();
	cnnn = nbrIter.ReportNextP();
	nbrIter.Prev();
	if( !AddEdge( cnnn, cn, cnn ) ) return 0;
	tmpList.insertAtBack( cn );
	tmpList.insertAtBack( cnn );
	tmpList.insertAtBack( cnnn );
	tmpList.makeCircular();
	if( !MakeTriangle( tmpList, tI ) ) return 0;
	return 1;
}


/**************************************************************************
**  
**  tMesh::MakeTriangle
**
**   Function to make triangle and add it to mesh; called
**   when all necessary nodes and edges already exist, i.e., the triangle
**   exists geometrically but not as a "triangle" member of the data
**   structure. Checks to make sure points are CCW. Resets triangle IDs at
**   end (necessary?). This function is relatively messy and complicated
**   but is extensively commented below.
**
**  Inputs: a tPtrList<tSubNode> of nodes in triangle; a tPtrListIter
**   object, the iterator of the latter list; edge is added between node
**   currently pointed to by iterator and the node-after-next. List must
**   contain three, and only three, members and be circular.
**  Created: SL fall, '97
**
**  Modifications:
**   - mSearchOriginTriPtr is reset to point to the newly added triangle
**     in an attempt to speed up triangle searches, especially during
**     mesh creation. GT 1/2000
**
**************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
MakeTriangle( tPtrList< tSubNode > &nbrList,
              tPtrListIter< tSubNode > &nbrIter )
{
	//assert( (&nbrList != 0) && (&nbrIter != 0) ); //WR-09192023: warning: reference cannot be bound to dereferenced null pointer in well-defined C++ code; comparison may be assumed to always evaluate to true
	assert( nbrList.getSize() == 3 );
	int i, j;
	
	tTriangle *nbrtriPtr;
	tSubNode *cn, *cnn, *cnnn;
	tEdge *ce;
	tTriangle *ct;
	tListIter< tTriangle > triIter( triList );
	tMeshListIter< tEdge > edgIter( edgeList );
	tPtrListIter< tEdge > spokIter;
	assert( nbrList.getSize() == 3 );
	
	cn = nbrIter.FirstP();      // cn, cnn, and cnnn are the 3 nodes in the tri
	cnn = nbrIter.NextP();
	cnnn = nbrIter.NextP();
	nbrIter.Next();
	tArray< double > p0( cn->get2DCoords() ), p1( cnn->get2DCoords() ),
		p2( cnnn->get2DCoords() );
	
	// Create the new triangle and insert a pointer to it on the list.
	// Here, the triangle constructor takes care of setting pointers to
	// the 3 vertices and 3 edges. The neighboring triangle pointers are
	// initialized to zero.
	
	triList.insertAtBack( tTriangle( miNextTriID++, cn, cnn, cnnn ) );//put 
		
		ct = triIter.LastP();                  //set triIter to last
		assert( cn == ct->pPtr(0) );           //make sure we're where we think we are
		
		// To speed up future searches in calls to LocateTriangle, assign the
		// starting triangle, mSearchOriginTriPtr, to the our new triangle.
		// The idea here is that there's a good chance that the next point
		// to be added will be close to the current location. (added 1/2000)
		
		mSearchOriginTriPtr = ct;
		
		// Now we assign the neighbor triangle pointers. The loop successively
		// gets the spokelist for (p0,p1,p2) and sets cn to the next ccw point
		// (p1,p2,p0). It then finds the edge (spoke) that joins the two points
		// (p0->p1, p1->p2, p2->p0). These are the edges that are shared with
		// neighboring triangles (t2,t0,t1) and are pointed to by the neighboring
		// triangles. This means that in order to find neighboring triangle t2,
		// we need to find the triangle that points to edge (p0->p1), and so on.
		// In general, t((j+2)%3) is the triangle that points to edge
		// p(j)->p((j+1)%3).
		
		nbrtriPtr = 0;
		cn = nbrIter.FirstP();
		for( j=0; j<3; j++ )
		{
			// get spokelist for p(j) and advance to p(j+1)
			spokIter.Reset( cn->getSpokeListNC() );
			cn = nbrIter.NextP();               //step forward once in nbrList
			
			// Find edge ce that connects p(j)->p(j+1)
			for( ce = spokIter.FirstP();
				 ce->getDestinationPtrNC() != cn && !( spokIter.AtEnd() );
				 ce = spokIter.NextP() );
			assert( !( spokIter.AtEnd() ) );
			
			if( !( TriWithEdgePtr( ce ) != nbrtriPtr || nbrtriPtr == 0 ) )
			{
				p0 = cn->get2DCoords();
				p1 = cnn->get2DCoords();
				p2 = cnnn->get2DCoords();
				
				if( PointsCCW( p0, p1, p2 ) )
					cerr << "\nWarning: In MakeTriangle()";
				else cerr << "tri not CCW: " << nbrtriPtr->getID() << endl;
			}
			
			// Find the triangle, if any, that shares (points to) this edge
			// and assign it as the neighbor triangle t((j+2)%3).
			
			nbrtriPtr = TriWithEdgePtr( ce );
			
			ct->setTPtr( (j+2)%3, nbrtriPtr );      //set tri TRI ptr (j+2)%3
			
			if( nbrtriPtr != 0 ){
				for( i=0; i<3; i++ ){
					assert( nbrtriPtr->ePtr(i) != 0 );
					assert( ce != 0 );
					if( nbrtriPtr->ePtr(i) == ce ) break;
				}
				assert( i < 3 );
				nbrtriPtr->setTPtr( (i+1)%3, ct );  //set NBR TRI ptr to tri
			}
		}   
		ntri++;
		
		for( ct = triIter.FirstP(), miNextTriID=0; !( triIter.AtEnd() );
			 ct = triIter.NextP(), miNextTriID++ )
		{
			ct->setID( miNextTriID );
		}
		
		return 1;
}

/**************************************************************************
**
**   tMesh::AddNode ( tSubNode nodeRef& )
**
**   Adds a new node with the properties of nodRef to the mesh.
**
**   Calls: tMesh::LocateTriangle, tMesh::DeleteTriangle, tMesh::AddEdge,
**            tMesh::AddEdgeAndMakeTriangle, tMesh::MakeTriangle,
**            tMesh::CheckForFlip; various member functions of tNode,
**            tMeshList, tMeshListIter, tPtrList, etc. Also tLNode
**            functions (TODO: this needs to be removed somehow),
**            and temporarily, tMesh::UpdateMesh
**   Parameters: nodeRef -- reference to node to be added (really,
**                          duplicated)
**   Returns:  
**   Assumes:
**   Created: SL fall, '97
**   Modifications:
**        - 4/98: node is no longer assumed to be a non-boundary (GT)
**        - 7/98: changed return type from int (0 or 1) to ptr to
**                the new node (GT)
**        -10/98: if node is open boundary,
**                added with tMeshList::insertAtBoundFront() (SL)
**        -5/99: removed unreferenced vars tedg1, tedg3 (GT)
**
**************************************************************************/

#define kLargeNumber 1000000000
template< class tSubNode >
tSubNode * tMesh< tSubNode >::
AddNode( tSubNode &nodeRef, int updatemesh, double time )
{
	int i, ctr;
	tTriangle *tri;
	tSubNode *cn;
	tArray< double > xyz( nodeRef.get3DCoords() );
	tMeshListIter< tSubNode > nodIter( nodeList );
	//assert( &nodeRef != 0 ); //WR--09192023:  reference cannot be bound to dereferenced null pointer in well-defined C++ code; comparison may be assumed to always evaluate to true
	
	tri = LocateTriangle( xyz[0], xyz[1] );
	if( tri == 0 ){
		cout<<"\ntMesh::AddNode(): coords out of bounds: "<<xyz[0]<<" "<<xyz[1]<<endl;
		return 0;
	}
	// Layering off in tRibs 
	// ---------------------
	// if( layerflag && time > 0 ) 
	//    nodeRef.LayerInterpolation( tri, xyz[0], xyz[1], time );
	
	// The next two statements are both unneccessary, and awkward-SMR
	//nodeRef.setID( miNextNodeID );  
	
	miNextNodeID++;
	
	if((nodeRef.getBoundaryFlag()==kNonBoundary) || (nodeRef.getBoundaryFlag()== kStream)){
		nodeList.insertAtActiveBack( nodeRef );
	}
	else if( nodeRef.getBoundaryFlag() == kOpenBoundary ){
		nodeList.insertAtBoundFront( nodeRef );
	}
	else{
		nodeList.insertAtBack( nodeRef );
	}
	
	unsortList.insertAtBack( nodeRef );
	
	nnodes++;
	
	// Retrieve a pointer to the new node and flush its spoke list
	if( (nodeRef.getBoundaryFlag() == kNonBoundary) || (nodeRef.getBoundaryFlag()== kStream))
		cn = nodIter.LastActiveP();
	else if( nodeRef.getBoundaryFlag() == kOpenBoundary )
		cn = nodIter.FirstBoundaryP();
	else
		cn = nodIter.LastP();
	assert( cn!=0 );
	cn->getSpokeListNC().Flush();
	
	//Make ptr list of triangle's vertices
	tPtrList< tSubNode > bndyList;
	tSubNode *tmpPtr;
	for( i=0; i<3; i++ )
	{
		tmpPtr = (tSubNode *) tri->pPtr(i);
		bndyList.insertAtBack( tmpPtr );
	}
	bndyList.makeCircular();
	
	
	// Delete the triangle in which the node falls
	i = DeleteTriangle( tri );
	assert( i != 0 );  //if ( !DeleteTriangle( tri ) ) return 0;
	
	// Make 3 new triangles
	
	tPtrListIter< tSubNode > bndyIter( bndyList );
	tSubNode *node3 = bndyIter.FirstP();     // p0 in original triangle
	tSubNode *node2 = cn;                    // new node
	tSubNode *node1 = bndyIter.NextP();      // p1 in orig triangle
	tSubNode *node4 = bndyIter.NextP();      // p2 in orig triangle
	tArray< double > p1( node1->get2DCoords() ),
		p2( node2->get2DCoords() ), p3( node3->get2DCoords() ),
		p4( node4->get2DCoords() );
	
	//if( xyz.getSize() == 3){
	//if(!PointsCCW(p3,p1,p2) || !PointsCCW(p2,p1,p4) || !PointsCCW(p2,p4,p3))
	//cout << "new tri not CCW" << endl; }
	//else{
	//if( node1->Meanders() ) p1 = node1->getNew2DCoords();
	//if( node2->Meanders() ) p2 = node2->getNew2DCoords();
	//if( node3->Meanders() ) p3 = node3->getNew2DCoords();
	//if( node4->Meanders() ) p4 = node4->getNew2DCoords();  
	//if(!PointsCCW(p3,p1,p2) || !PointsCCW(p2,p1,p4) || !PointsCCW(p2,p4,p3))
	//cout << "new tri not CCW" << endl;
	//}
	
	// Here's how the following works. Let the old triangle vertices be A,B,C
	// and the new node N. The task is to create 3 new triangles ABN, NBC, and
	// NCA, and 3 new edge-pairs AN, BN, and CN.
	// First, edge pair BN is added. Then AEMT is called to create triangle
	// ABN and edge pair AN. AEMT is called again to create tri NBC and edge
	// pair CN. With all the edge pairs created, it remains only to call
	// MakeTriangle to create tri NCA.
	
	assert( node1 != 0 && node2 != 0 && node3 != 0 );
	AddEdge( node1, node2, node3 );  //add edge between node1 and node2
	tPtrList< tSubNode > tmpList;
	tmpList.insertAtBack( node3 );  // ABN
	tmpList.insertAtBack( node1 );
	tmpList.insertAtBack( node2 );
	tPtrListIter< tSubNode > tmpIter( tmpList );
	AddEdgeAndMakeTriangle( tmpList, tmpIter );
	tmpList.Flush();
	tmpList.insertAtBack( node2 );  // NBC
	tmpList.insertAtBack( node1 );
	tmpList.insertAtBack( node4 );
	tmpIter.First();
	AddEdgeAndMakeTriangle( tmpList, tmpIter );
	tmpList.Flush();
	tmpList.insertAtBack( node2 );  // NCA
	tmpList.insertAtBack( node4 );
	tmpList.insertAtBack( node3 );
	tmpList.makeCircular();
	tmpIter.First();
	MakeTriangle( tmpList, tmpIter );
	
	//hasn't changed yet, put 3 resulting triangles in ptr list
	//cout << "Putting tri's on list\n" << flush;
	
	if( xyz.getSize() == 3 )
	{
		tPtrList< tTriangle > triptrList;
		tListIter< tTriangle > triIter( triList );
		tPtrListIter< tTriangle > triptrIter( triptrList );
		tTriangle *ct;
		
		triptrList.insertAtBack( triIter.LastP() );
		triptrList.insertAtBack( triIter.PrevP() );
		triptrList.insertAtBack( triIter.PrevP() );
		
		//check list for flips; if flip, put new triangles at end of list
		int flip = 1;
		ctr = 0;
		while( !( triptrList.isEmpty() ) ){
			ctr++;
			if( ctr > kLargeNumber ) // Make sure to prevent endless loops
			{                       
				cerr << "Mesh error: adding node " << node2->getID()
				<< " flip checking forever" << endl;
				cerr << "Bailing out of AddNode()";
			}
			ct = triptrIter.FirstP();
			
			for( i=0; i<3; i++ )
			{
				if( ct->tPtr(i) != 0 )
				{
					if( CheckForFlip( ct, i, flip ) )
					{
						triptrList.insertAtBack( triIter.LastP() );
						triptrList.insertAtBack( triIter.PrevP() );
						break;
					}
				}
			}
			
			triptrList.removeFromFront( ct );
		}
	}
	
	node2->makeCCWEdges();
	node2->InitializeNode();
	
	if( updatemesh ) UpdateMesh();
	
	return node2;  // Return ptr to new node
}


/**************************************************************************
**
**  tMesh::AddNodeAt
**
**   add a node with referenced coordinates to mesh;
**   this fn duplicates functionality of AddNode
**
**  Created: SL fall, '97
**  Modified: NG summer, '98 to deal with layer interpolation
**
**************************************************************************/

template< class tSubNode >
tSubNode *tMesh< tSubNode >::
AddNodeAt( tArray< double > &xyz, double time )
{
	//assert( &xyz != 0 ); //WR--09192023: reference cannot be bound to dereferenced null pointer in well-defined C++ code; comparison may be assumed to always evaluate to true
	tTriangle *tri;
	if( xyz.getSize() == 3 ) tri = LocateTriangle( xyz[0], xyz[1] );
	else tri = LocateNewTriangle( xyz[0], xyz[1] );
	if( tri == 0 )      return 0;
	
	int i, ctr;
	tMeshListIter< tSubNode > nodIter( nodeList );
	tSubNode tempNode, *cn;
	tempNode.set3DCoords( xyz[0], xyz[1], xyz[2]  );
	
	// Layering off in tRIBS
	// ---------------------
	//if( layerflag && time > 0.0) tempNode.LayerInterpolation( tri, xyz[0], xyz[1], time );
	
	if( xyz.getSize() != 3 ) tempNode.setNew2DCoords( xyz[0], xyz[1] );
	tempNode.setBoundaryFlag( 0 );
	
	// Assign ID to the new node and insert it at the back of the active
	// portion of the node list (NOTE: node is assumed NOT to be a boundary)
	
	tempNode.setID( miNextNodeID );
	miNextNodeID++;
	Cout << miNextNodeID << endl;
	nodeList.insertAtActiveBack( tempNode );
	assert( nodeList.getSize() == nnodes + 1 );
	nnodes++;
	
	//Make ptr list of triangle's vertices:
	tPtrList< tSubNode > bndyList;
	tSubNode *tmpPtr;
	for( i=0; i<3; i++ ){
		tmpPtr = (tSubNode *) tri->pPtr(i);
		bndyList.insertAtBack( tmpPtr );
	}
	bndyList.makeCircular();
	
	if ( !DeleteTriangle( tri ) ) return 0;
	
	// Make 3 new triangles
	tPtrListIter< tSubNode > bndyIter( bndyList );
	tSubNode *node3 = bndyIter.FirstP();
	tSubNode *node2 = nodIter.LastActiveP();
	tSubNode *node1 = bndyIter.NextP();
	tSubNode *node4 = bndyIter.NextP();
	tArray< double > p1( node1->get2DCoords() ),
		p2( node2->get2DCoords() ), p3( node3->get2DCoords() ),
		p4( node4->get2DCoords() );
	if( xyz.getSize() == 3)
	{
		if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
			cout << "new tri not CCW" << endl;
	}
	else
	{
		// Meandering off in tRIBS
		// -----------------------
		//if( node1->Meanders() ) p1 = node1->getNew2DCoords();
		//if( node2->Meanders() ) p2 = node2->getNew2DCoords();
		//if( node3->Meanders() ) p3 = node3->getNew2DCoords();
		//if( node4->Meanders() ) p4 = node4->getNew2DCoords();  
		
		if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
			cout << "new tri not CCW" << endl;
	}
	
	assert( node1 != 0 && node2 != 0 && node3 != 0 );
	AddEdge( node1, node2, node3 );  //add edge between node1 and node2
	tPtrList< tSubNode > tmpList;
	tmpList.insertAtBack( node3 );
	tmpList.insertAtBack( node1 );
	tmpList.insertAtBack( node2 );
	tPtrListIter< tSubNode > tmpIter( tmpList );
	AddEdgeAndMakeTriangle( tmpList, tmpIter );
	tmpList.Flush();
	tmpList.insertAtBack( node2 );
	tmpList.insertAtBack( node1 );
	tmpList.insertAtBack( node4 );
	tmpIter.First();
	AddEdgeAndMakeTriangle( tmpList, tmpIter );
	tmpList.Flush();
	tmpList.insertAtBack( node2 );
	tmpList.insertAtBack( node4 );
	tmpList.insertAtBack( node3 );
	tmpList.makeCircular();
	tmpIter.First();
	MakeTriangle( tmpList, tmpIter );
	
	// Put 3 resulting triangles in ptr list
	if( xyz.getSize() == 3 )
	{
		tPtrList< tTriangle > triptrList;
		tListIter< tTriangle > triIter( triList );
		tPtrListIter< tTriangle > triptrIter( triptrList );
		tTriangle *ct;
		triptrList.insertAtBack( triIter.LastP() );
		triptrList.insertAtBack( triIter.PrevP() );
		triptrList.insertAtBack( triIter.PrevP() );
		
		//Check list for flips; if flip, put new triangles at end of list
		int flip = 1;
		ctr = 0;
		while( !( triptrList.isEmpty() ) )
		{
			ctr++;
			if( ctr > kLargeNumber ) // Make sure to prevent endless loops
			{
				cerr << "Mesh error: adding node " << node2->getID()
				<< " flip checking forever"
				<< endl;
				cerr << "Bailing out of AddNodeAt()";
			}
			ct = triptrIter.FirstP();
			for( i=0; i<3; i++ )
			{
				if( ct->tPtr(i) != 0 )
				{
					if( CheckForFlip( ct, i, flip ) )
					{
						triptrList.insertAtBack( triIter.LastP() );
						triptrList.insertAtBack( triIter.PrevP() );
						break;
					}
				}
			}
			triptrList.removeFromFront( ct );
		}
	}
	
	// Reset node id's
	Cout << "reset ids\n";
	for( cn = nodIter.FirstP(), miNextNodeID=0; !( nodIter.AtEnd() ); cn = nodIter.NextP(), miNextNodeID++ ){
		cn->setID( miNextNodeID );
	}
	
	node2->makeCCWEdges();
	node2->InitializeNode();
	
	UpdateMesh();
	
	tEdge *ce, *fe;
	fe = node2->getFlowEdg();
	ce = fe;
	
	int hlp=0;
	do{
		ce=ce->getCCWEdg();
		hlp++;
	}while(ce != fe );
	
	return node2;
}
#undef kLargeNumber


//=========================================================================
//
//
//                 Section 12: tMesh Get Functions
//
//
//=========================================================================

/**************************************************************************
**
**  tMesh "get" functions
**
**************************************************************************/

template <class tSubNode>
tMeshList<tSubNode> * tMesh<tSubNode>::
getUnsortList() {return &unsortList;}

template <class tSubNode>
tList< tTriangle > * tMesh<tSubNode>::   
getTriList() {return &triList;}


/**************************************************************************
**
**  tMesh::getEdgeComplement
**
**  Returns the complement of _edge_ (i.e., the edge that shares the same
**  endpoints but points in the opposite direction). To find the complement,
**  it exploits the fact that complementary pairs of edges are stored 
**  together on the edge list, with the first of each pair having an 
**  even-numbered ID and the second having an odd-numbered ID.
**
**  Modifications: gt replaced 2nd IF with ELSE to avoid compiler warning
**
**************************************************************************/

template< class tSubNode >
tEdge *tMesh< tSubNode >::
getEdgeComplement( tEdge *edge )
{
	tMeshListIter< tEdge > edgIter( edgeList );
	int edgid = edge->getID();
	
	assert( edgIter.Get( edgid ) );
	edgIter.Get( edgid ); 
	if( edgid%2 == 0 ) return edgIter.GetP( edgid + 1 );
	else return edgIter.GetP( edgid - 1 );
}


//=========================================================================
//
//
//                 Section 13: Updating and Flipping Mesh
//
//
//=========================================================================

/**************************************************************************
**
**  tMesh::UpdateMesh
**
**  Updates mesh geometry:
**   - computes edge lengths
**   - finds Voronoi vertices
**   - computes Voronoi edge lengths
**   - computes Voronoi areas for interior (active) nodes
**   - updates CCW-edge connectivity
**
**  Note that the call to CheckMeshConsistency is for debugging
**  purposes and should be removed prior to release.
**
**  Calls: MakeCCWEdges(), setVoronoiVertices(), CalcVoronoiEdgeLengths(),
**   CalcVAreas(), CheckMeshConsistency()
**  Assumes: nodes have been properly triangulated
**  Created: SL fall, '97
**
**************************************************************************/

template <class tSubNode>
void tMesh<tSubNode>::
UpdateMesh()
{ 
	tMeshListIter<tEdge> elist( edgeList );
	tEdge * curedg = 0;
	double len;
	
	// Edge lengths
	curedg = elist.FirstP();
	do{
		len = curedg->CalcLength();
		if(len <= 0.0){
			cout<<"Point Destin X = " <<curedg->getDestinationPtr()->getX();
			cout<<"\nPoint Destin Y = " <<curedg->getDestinationPtr()->getY();
			cout<<"\nPoint Destin Z = " <<curedg->getDestinationPtr()->getZ();
			cout<<"\nPoint Destin B = " <<curedg->getDestinationPtr()->getBoundaryFlag();
			cout<<"\n\nPoint Origin X = "<<curedg->getOriginPtr()->getX();
			cout<<"\nPoint Origin Y = "<<curedg->getOriginPtr()->getY();
			cout<<"\nPoint Origin Z = "<<curedg->getOriginPtr()->getZ();
			cout<<"\nPoint Origin B = "<<curedg->getOriginPtr()->getBoundaryFlag();
			cout<<"\n\n";
		}
		assert( len>0.0 );
		curedg = elist.NextP();
		assert( curedg != nullptr ); // failure = complementary edges not consecutive, compiler error indicates comparison between pointer and zero, so replaced with null pter -WR
		curedg->setLength( len );
	} while( (curedg=elist.NextP()) != nullptr);//TODO: is this correct or semantic error? WR warning: using the result of an assignment as a condition without parentheses [-Wparentheses]
	
	MakeCCWEdges();
	
	setVoronoiVertices();
	CalcVoronoiEdgeLengths();
	CalcVAreas();
}

/*****************************************************************************
**
**  tMesh::CheckForFlip
**
**  Checks whether edge between two triangles should be
**  flipped; may either check, flip, and report, or just check and report.
**  Checks whether the present angle or the possible angle
**  is greater. Greater angle wins. Also uses flip variable
**  to determine whether to use newx, newy, or x, y.
**
**      Inputs: tri -- ptr to the triangle to be tested
**              nv -- the number of the vertex opposite the edge that
**                    might be flipped (0, 1, or 2)
**              flip -- flag indicating whether we want to actually flip
**                      the edge if needed (TRUE) or simply test the flip
**                      condition for a point that is about to be moved to
**                      a new position (FALSE)
**      Returns: 1 if flip is needed, 0 otherwise
**      Modifies: edge may be flipped
**      Called by: AddNode, AddNodeAt, CheckLocallyDelaunay,
**                 tStreamMeander::CheckBrokenFlowedge
**      Calls: PointsCCW, FlipEdge, TriPasses                                             
**                                                  
*****************************************************************************/

template< class tSubNode >
int tMesh< tSubNode >::
CheckForFlip( tTriangle * tri, int nv, int flip )
{
	if( tri == 0 ){
		cout << "CheckForFlip: tri == 0" << endl;
		return 0;
	}
	assert( nv < 3 );
	
	tSubNode *node0, *node1, *node2, *node3;
	node0 = ( tSubNode * ) tri->pPtr(nv);
	node1 = ( tSubNode * ) tri->pPtr((nv+1)%3);
	node2 = ( tSubNode * ) tri->pPtr((nv+2)%3);
	
	tTriangle *triop = tri->tPtr(nv);
	int nvop = triop->nVOp( tri );
	node3 = ( tSubNode * ) triop->pPtr( nvop );
	tArray< double > ptest( node3->get2DCoords() ), p0( node0->get2DCoords() ),
		p1( node1->get2DCoords() ), p2( node2->get2DCoords() );
	
	// Meandering Off in tRIBS
	// -----------------------
	//if( !flip ){
	//   if( node0->Meanders() ) p0 = node0->getNew2DCoords();
	//   if( node1->Meanders() ) p1 = node1->getNew2DCoords();
	//   if( node2->Meanders() ) p2 = node2->getNew2DCoords();
	//   if( node3->Meanders() ) ptest = node3->getNew2DCoords();
	//}
	
	// If p0-p1-p2 passes the test, no flip is necessary
	if( TriPasses( ptest, p0, p1, p2 ) ) return 0;
	
	// Otherwise, a flip is needed, provided that the new triangles are
	// counter-clockwise and that the node isn't a moving node 
	
	if( flip ){
		if( !PointsCCW( p0, p1, ptest ) || !PointsCCW( p0, ptest, p2 ) ) 
			return 0;
		FlipEdge( tri, triop, nv, nvop );
		
	}
	return 1;
}


/******************************************************************
**
**  tMesh::FlipEdge
**
**  Flips the edge pair between two adjacent triangle to
**  re-establish Delaunay-ness.
**
**  Note on notation in flip edge:
**
**                d
**               /|\
**       tri->  / | \ <-triop
**             /  |  \
**            a   |   c
**             \  |  /
**              \ | /
**               \|/
**                b
**        Edge bd will be removed
**        and an edge ac will be made.
**        nbrList contains the points a, b, c, d
**
**    Inputs:  tri, triop -- the triangles sharing the edge to be
**                           flipped
**             nv -- the number of tri's vertex (0, 1 or 2) opposite
**                   the edge (ie, point a)
**             nvop -- the number of triop's vertex (0, 1 or 2)
**                     opposite the edge (ie, point c)
**    Calls: DeleteEdge, AddEdgeAndMakeTriangle, MakeTriangle
**    Called by: CheckForFlip, CheckTriEdgeIntersect
**
*******************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
FlipEdge( tTriangle * tri, tTriangle * triop ,int nv, int nvop )
{
	
	tSubNode *cn = 0;
	tPtrList< tSubNode > nbrList;
	
	// Place the four vertices of the two triangles on a list
	nbrList.insertAtBack( (tSubNode *) tri->pPtr(nv) );
	nbrList.insertAtBack( (tSubNode *) tri->pPtr((nv+1)%3) );
	nbrList.insertAtBack( (tSubNode *) triop->pPtr( nvop ) );
	nbrList.insertAtBack( (tSubNode *) tri->pPtr((nv+2)%3) );
	nbrList.makeCircular();
	
	// Delete the edge pair between the triangles, along with the tri's
	
	DeleteEdge( tri->ePtr( (nv+2)%3 ) );  // Changed for right-hand data struc
	
	// Recreate the triangles and the edges in their new orientation
	tPtrListIter< tSubNode > nbrIter( nbrList );
	AddEdgeAndMakeTriangle( nbrList, nbrIter );
	nbrIter.First();
	nbrList.removeNext( cn, nbrIter.NodePtr() );
	MakeTriangle( nbrList, nbrIter );
	
}


//=========================================================================
//
//
//                  Section 14: Meandering Functions on Mesh
//
//
//=========================================================================

/*****************************************************************************
**
**  tMesh::CheckLocallyDelaunay
**
**  Updates the triangulation after moving some points.
**  Only uses x and y values, which have already been updated in
**  MoveNodes (frmr PreApply).
**  MoveNodes SHOULD BE CALLED BEFORE THIS FUNCTION IS CALLED
**
**  The logic here is somewhat complicated. Here is GT's understanding
**  of it (Stephen, can you confirm?):
**
**  1. We create a list of triangles that have at least one vertex that has
**     moved (triPtrList) and which therefore might no longer be
**     Delaunay.
**  2. For each of these, we do a flip check across each face. Before
**     doing so, however, we find the triangle on the triPtrList, if any,
**     that comes just before this neighboring triangle. If the edge between
**     the triangles gets flipped, both the triangles will be deleted and
**     recreated on the master triangle list; thus, we will need to delete
**     both affected triangles from triPtrList and re-add the new ones.
**  3. If a flip occurs, remove the opposite triangle pointer from the
**     list if needed in order to prevent a dangling pointer. The two
**     affected triangles will have been replaced by new triangles which
**     are now at the back of the master triangle list; add these two to
**     the triPtrList to be rechecked, and break out of the vertex loop.
**  4. Remove the triangle in question from the head of the triPtrList
**     (regardless of whether it was flipped or not; if it was, its a 
		**     dangling pointer; if not, it is Delaunay and we no longer need
		**     worry about it)
**  5. Continue until there are no more triangles to be checked.
**
**        
*****************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
CheckLocallyDelaunay()
{
	tTriangle *at;
	tPtrList< tTriangle > triPtrList;
	tPtrListIter< tTriangle > triPtrIter( triPtrList );
	tListIter< tTriangle > triIter( triList );
	int i, change;
	tArray< int > npop(3);
	tSubNode *nodPtr;
	int flip = 1;
	
	// Search through tri list to find triangles with at least one
	// moving vertex, and put these on triPtrList, put each triangle into the stack
	
	for( at = triIter.FirstP(); !( triIter.AtEnd() ); at = triIter.NextP() ){
		change = FALSE;
		for( i = 0; i < 3; i++ )
		{
			nodPtr = ( tSubNode * ) at->pPtr(i);
			//if( nodPtr->Meanders() ) change = TRUE;
		}
		if( change ) triPtrList.insertAtBack( at );
	}
	
	// Check list for flips; if flip, put new triangles at end of list
	
	tPtrListIter< tTriangle > duptriPtrIter( triPtrList );
	tTriangle *tn, *tp;
	while( !( triPtrList.isEmpty() ) ) {
		at = triPtrIter.FirstP();
		for( i=0; i<3; i++ )
		{
			// If a neighboring triangle exists across this face, check for flip
			if( at->tPtr(i) != 0 )
			{
				tp = at->tPtr(i);
				for( tn = duptriPtrIter.FirstP();
					 duptriPtrIter.ReportNextP() != tp &&
                     !( duptriPtrIter.AtEnd() );
					 tn = duptriPtrIter.NextP() );
				tn = 0;
				if( !( duptriPtrIter.AtEnd() ) ){
					tn = duptriPtrIter.ReportNextP();
				}
				
				// Check triangle _at_ for a flip across face opposite vertex i,
				// and do the flip if needed
				if( CheckForFlip( at, i, flip ) )
				{
					
					if( tn != 0 )
						triPtrList.removeNext( tn, duptriPtrIter.NodePtr() );
					
					tn = triIter.LastP(); 
					
					triPtrList.insertAtBack( tn );
					tn = triIter.PrevP();
					
					triPtrList.insertAtBack( tn );
					break;
				}
			}
		}
		
		triPtrList.removeFromFront( at );
	}
}

/*****************************************************************************
**
**  tMesh::CheckTriEdgeIntersect
**
**        This function implements node movement.
**        We want to know if the moving point has passed beyond the polygon
**        defined by its spoke edges; if it has, then we will have edges
**        intersecting one another. In the case where the point has simply
**        passed into one of the 'opposite' triangles, then we can just do a
**        flip operation. In the other case, the remedial action is much more
**        complicated, so we just delete the point and add it again.
**
**
*****************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
CheckTriEdgeIntersect()
{
	int i, j, nv, nvopp;
	int flipped = TRUE;
	int crossed;
	tSubNode *subnodePtr, tempNode, newNode;  
	tEdge * cedg, *ce;
	tTriangle * ct, * ctop, *rmtri;
	tListIter< tTriangle > triIter( triList );
	
	tMeshListIter< tEdge > edgIter( edgeList );
	tMeshListIter< tSubNode > nodIter( nodeList );
	tMeshListIter< tEdge > xedgIter( edgeList );
	tPtrListIter< tEdge > spokIter;
	tMeshList< tSubNode > tmpNodeList;
	tMeshListIter< tSubNode > tmpIter( tmpNodeList );
	tArray< double > p0, p1, p2, xy, xyz, xy1, xy2;
	tSubNode *cn;
	tPtrList< tTriangle > triptrList;
	tPtrListNode< tTriangle > *tpListNode;
	tPtrListIter< tTriangle > tpIter( triptrList );
	
	//check for triangles with edges which intersect (an)other edge(s)
	
	while( flipped )
	{
		flipped = FALSE;
		
		// Make a list of triangles containing at least one moving vertex
		for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() ){
			for( i=0; i<3; i++ ){
				cn = (tSubNode *) ct->pPtr(i);
				//if( cn->Meanders() ) break;
			}
			if( i!=3 ) triptrList.insertAtBack( ct );
		}
		
		for( ct = tpIter.FirstP(); !(triptrList.isEmpty());
			 triptrList.removeFromFront( ct ), ct = tpIter.FirstP() ){ 
			if( !NewTriCCW( ct ) ){
				flipped = TRUE;
				for( i=0, j=0; i<3; i++ ){
					if( ct->pPtr(i)->getBoundaryFlag() != kNonBoundary ) j++;
				}
				if( j > 1 ){
					for( i=0, j=0; i<3; i++ ){
						subnodePtr = (tSubNode *) ct->pPtr(i);
						subnodePtr->RevertToOldCoords();
					}
				}
				else{   
					crossed = FALSE;
					for( i=0; i<3; i++ ){
						cn = (tSubNode *) ct->pPtr(i);
						if( cn->Meanders() ){
							cedg = ct->ePtr( (i+2)%3 );
							spokIter.Reset( cn->getSpokeListNC() );
							for( ce = spokIter.FirstP(); !( spokIter.AtEnd() );
								 ce = spokIter.NextP() ){
								if( Intersect( ce, cedg ) ){
									if( ct->tPtr(i) == 0 ) {
										subnodePtr = (tSubNode *) ct->pPtr(i);
										subnodePtr->RevertToOldCoords();
									}
									else{
										crossed = TRUE;
										ctop = ct->tPtr(i);
										xy = cn->getNew2DCoords();
										if( NewTriCCW( ctop ) && InNewTri( xy, ctop ) ){
											for( rmtri = tpIter.FirstP();
												 tpIter.ReportNextP() != ctop && !(tpIter.AtEnd());
												 rmtri = tpIter.NextP() );
											if( !(tpIter.AtEnd()) ) {
												tpListNode = tpIter.NodePtr();
												triptrList.removeNext( rmtri, tpListNode );
											}                           
											nv = ct->nVOp( ctop );
											nvopp = ctop->nVOp( ct );
											FlipEdge( ct, ctop, nv, nvopp );
											rmtri = triIter.LastP();
											triptrList.insertAtBack( rmtri );
											rmtri = triIter.PrevP();
											triptrList.insertAtBack( rmtri );
										}
										else{
											if( LocateTriangle( xy[0], xy[1] ) != 0 ){
												for( ce = spokIter.FirstP(); !(spokIter.AtEnd());
													 ce = spokIter.NextP() )
												{
													rmtri = TriWithEdgePtr( ce );
													for(tpIter.FirstP();
														tpIter.ReportNextP() != rmtri &&
														!(tpIter.AtEnd());
														tpIter.NextP() );
													if( !(tpIter.AtEnd()) ) {
														tpListNode = tpIter.NodePtr();
														triptrList.removeNext( rmtri, tpListNode );
													}
												}
												//delete the node;
												xyz = cn->getNew3DCoords();
												tmpNodeList.insertAtBack( *cn );
												DeleteNode( cn, kRepairMesh );
											}
											else{
												subnodePtr = (tSubNode *) ct->pPtr(i);
												subnodePtr->RevertToOldCoords();
											}
										}
									}
									break;
								}
							}
						}
						if( crossed ) break;
					}
				}      
			}
		}
	}
	
	// Update coordinates of moving nodes.
	
	// Meandering off
	// --------------
	// for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() )
	// if ( cn->Meanders() ) cn->UpdateCoords();
	
	for( cn = tmpIter.FirstP(); !(tmpIter.AtEnd()); cn = tmpIter.NextP() )
	{
		//if ( cn->Meanders() ) cn->UpdateCoords();
		cn->getSpokeListNC().Flush();
		cn = AddNode( *cn );
		assert( cn!=0 );
	}
}


/*****************************************************************************
**
**  tMesh::MoveNodes 
**
**  Once the new coordinates for moving nodes have been established, this
**  function is called to update the node coordinates, modify the 
**  triangulation as needed, and update the mesh geometry (Voronoi areas,
**  edge lengths, etc) through a series of calls to helper functions.
**  
**  Interpolation is performed on nodes with layering (3D vertical 
**  component) here. TODO: make interpolation general, perhaps by
**  defining a virtual tNode function called "AlertNodeMoving" 
**
**      Inputs: time -- simulation time (for layer updating)
**      Data members updated: Mesh elements & their geometry
**      Called by:  called outside of tMesh by routines that compute
**                  node movement (e.g., stream meandering, as implemented
**                  by tStreamMeander)
**      Calls: CheckTriEdgeIntersect, CheckLocallyDelaunay, UpdateMesh,
**             LocateTriangle, tLNode::LayerInterpolation
**      Created: SL       
**
*****************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MoveNodes( double time )
{
	//tSubNode * cn;  
	//tMeshListIter< tSubNode > nodIter( nodeList );
	
	// Layering off in tRIBS
	// ----------------------
	// if( layerflag && time > 0.0 ) {   
	//   tTriangle *tri;
	//   tArray<double> newxy(2);
	//   for(cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP()){
	//      newxy=cn->getNew2DCoords();
	//      if( (cn->getX()!=newxy[0]) || (cn->getY()!=newxy[1]) ){
	//         tri = LocateTriangle( newxy[0], newxy[1] );
	//         cn->LayerInterpolation( tri, newxy[0], newxy[1], time );   
	//      } } }
	
	//check for triangles with edges which intersect (an)other edge(s)
	//calls tCNode::UpdateCoords() for each node
	
	CheckTriEdgeIntersect();
	CheckLocallyDelaunay();
	UpdateMesh();
	CheckMeshConsistency(); 
}

/*****************************************************************************
**
**  tMesh::AddNodesAround
**
**  Densifies the mesh in the vicinity of a given node (centerNode) by
**  adding new nodes at the coordinates of the centerNode's Voronoi
**  vertices.
**
**  Properties of each node are initially those of the centerNode, except
**  z which is computed using interpolation by getVoronoiVertexXYZList.
**
**      Inputs: centerNode -- the node around which to add new nodes
**              time -- simulation time (for layer updating)
**      Data members updated: Mesh elements & their geometry
**      Called by:  called outside of tMesh by routines that handle
**                  adaptive meshing
**      Calls: AddNode, UpdateMesh, tNode::getVoronoiVertexXYZList
**      Created: GT, for dynamic mesh updating, Feb 2000  
**
*****************************************************************************/

template<class tSubNode>
void tMesh< tSubNode >::
AddNodesAround( tSubNode * centerNode, double time )
{
	tList< Point3D > vvtxlist;  
	tListIter< Point3D > vtxiter( vvtxlist );
	
	assert( centerNode!=0 );
	
	centerNode->getVoronoiVertexXYZList( &vvtxlist );
	tCNode tmpnode = *centerNode;  // New node to be added -- passed to AddNode                                 
	Point3D *xyz;  		  // Coordinates of current vertex
	
	for( xyz=vtxiter.FirstP(); !(vtxiter.AtEnd()); xyz=vtxiter.NextP() ){
		tmpnode.set3DCoords( xyz->x, xyz->y, xyz->z );  // Assign to tmpnode
		AddNode( tmpnode, FALSE, time );  	      // Add the node
	}
	UpdateMesh();
	
}

//=========================================================================
//
//
//                  Section 15: Mesh Data Debugging Functions
//
//
//=========================================================================

#ifndef NDEBUG

/*****************************************************************************
**
**      DumpEdges(), DumpSpokes(), DumpTriangles(), DumpNodes(): debugging
**         routines which simply write out information pertaining to the mesh;
**      DumpNodes() calls DumpSpokes for each node;
**      DumpSpokes() takes a pointer to a node as an argument.
**
*****************************************************************************/

template<class tSubNode>
void tMesh<tSubNode>::
DumpEdges()
{
	tMeshListIter< tEdge > edgIter( edgeList );
	tEdge *ce;
	tTriangle *ct;
	int tid;
	for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() ){
		ct = TriWithEdgePtr( ce );
		tid = ( ct != 0 ) ? ct->getID() : -1;
		cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
			<< " to " << ce->getDestinationPtrNC()->getID() << "; in tri "
			<< tid << " (flw " << ce->getBoundaryFlag() << ")" << endl;
	}
}

template<class tSubNode>
void tMesh<tSubNode>::
DumpSpokes( tSubNode *cn )
{
	tEdge *ce;
	tPtrListIter< tEdge > spokIter( cn->getSpokeListNC() );
	Cout << "node " << cn->getID() << " with spoke edges " << endl;
	for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() ){
		Cout << "   " << ce->getID()
		<< " from node " << ce->getOriginPtrNC()->getID()
		<< " to " << ce->getDestinationPtrNC()->getID() << endl;
	}
}

template<class tSubNode>
void tMesh<tSubNode>::
DumpTriangles()
{
	tListIter< tTriangle > triIter( triList );
	tTriangle *ct, *nt;
	int tid0, tid1, tid2;
	Cout << "triangles:" << endl;
	for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
	{
		nt = ct->tPtr(0);
		tid0 = ( nt != 0 ) ? nt->getID() : -1;
		nt = ct->tPtr(1);
		tid1 = ( nt != 0 ) ? nt->getID() : -1;
		nt = ct->tPtr(2);
		tid2 = ( nt != 0 ) ? nt->getID() : -1;
		Cout << ct->getID() << " with vertex nodes "
			<< ct->pPtr(0)->getID() << ", "
			<< ct->pPtr(1)->getID() << ", and "
			<< ct->pPtr(2)->getID() << "; edges "
			<< ct->ePtr(0)->getID() << ", "
			<< ct->ePtr(1)->getID() << ", and "
			<< ct->ePtr(2)->getID() << "; nbr triangles "
			<< tid0 << ", "
			<< tid1 << ", and "
			<< tid2 << endl;
	}
}

template<class tSubNode>
void tMesh<tSubNode>::
DumpNodes(){
	tMeshListIter< tSubNode > nodIter( nodeList );
	tSubNode *cn;
	Cout << "nodes: " << endl;
	for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() ){
		Cout << " at " << cn->getX() << ", " << cn->getY() << ", " << cn->getZ()
		<< "; bndy: " << cn->getBoundaryFlag() << "; ";
		DumpSpokes( cn );
	}
}

template<class tSubNode>
void tMesh<tSubNode>::TellAboutNode(tSubNode *cn){
    cout<<cn->getID()
	<<"\t"<<cn->getX()
	<<"\t"<<cn->getY()
	<<"\t"<<cn->getZ()
	<<"\t"<<cn->getBoundaryFlag()<<endl<<flush;
	return;
}

#endif
/***************************************************************************
**
** tMesh::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/


template<class tSubNode>
void tMesh<tSubNode>::writeRestart(ostream & rStr)
{
  tMeshListIter< tSubNode > nodIter( nodeList );
  tSubNode *cn;
  for (cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP())
    cn->writeRestart(rStr);
}


/***************************************************************************
**
** tMesh::readRestart() Function
** For each node restart state is read in
**
***************************************************************************/


template<class tSubNode>
void tMesh<tSubNode>::readRestart(istream & rStr)
{
  tMeshListIter< tSubNode > nodIter( nodeList );
  tSubNode *cn;
  for (cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP())
    cn->readRestart(rStr);
}

template<class tSubNode>
void tMesh<tSubNode>::readRestartGlobal(istream & rStr, int totalNodes)
{
  // Build a local nodeID → node pointer map so we can safely look up
  // nodes by ID without risking out-of-bounds access on NodeTable.
  unordered_map<int, tSubNode*> nodeMap;
  tMeshListIter<tSubNode> nI(nodeList);
  for (tSubNode* cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP())
    nodeMap[cn->getID()] = cn;

  for (int i = 0; i < totalNodes; i++) {
    int64_t nodeID;
    BinaryRead(rStr, nodeID);
    auto it = nodeMap.find(static_cast<int>(nodeID));
    if (it != nodeMap.end())
      it->second->readRestartBody(rStr);
    else
      tCNode::skipRestartBody(rStr);
  }
}


//=========================================================================
//
//
//                           End of tMesh.cpp
//
//
//=========================================================================
