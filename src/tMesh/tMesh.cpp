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
**  tMesh.cpp: Functions for class tMesh (see tMesh.h) based on CHILD
**             routines for TIN Mesh Generation
**
***************************************************************************/

#include "src/tMesh/tMesh.h"
#include "src/Headers/globalIO.h"
#include <iostream>

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
	
	if( read < 1 || read > 9 ){
		cerr << "\nInvalid mesh input option requested.";
		cerr << "\nValid options for reading mesh input are:\n"
			<< "  1 -- read mesh from input data files\n"
			<< "  2 -- create mesh from a list of (x,y,z,b) points\n"
			<< "  3 -- create random mesh from ArcGrid ascii output\n"
			<< "  4 -- create hexagonal mesh from ArcGrid ascii output\n"
			<< "  5 -- create mesh from ArcInfo file.net via points\n"
			<< "  6 -- create mesh from ArcInfo file.pnt and .lin via points\n"
			<< "  7 -- create mesh from scratch\n\n"
			<< "  8 -- create mesh from points file using Triangulator\n\n"
			<< "  9 -- create mesh from MeshBuilder files\n\n";
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
	
	else if( read == 2 ){
		Cout<<"\n\nPart 2: Creating Mesh from Points File (Option 2)"<<endl;
		Cout<<"---------------------------------------------------"<<endl;
		MakeMeshFromPoints( infile ); 
		Cout<<"\nMakeMeshFromPoints Successful Using Option 2"<<endl<<flush;
	}
	
	else if( read == 3 ){
		Cout<<"\n\nPart 2: Creating Mesh from ArcGrid - Random (Option 3)"<<endl;
		Cout<<"---------------------------------------------------"<<endl;
		MakeRandomPointsFromArcGrid( infile ); 
		Cout<<"\nMakeRandomPointsFromArcGrid Successful Using Option 3"<<endl<<flush;
	}
	
	else if( read == 4 ){
		Cout<<"\n\nPart 2: Creating Mesh from ArGrid - Hexagonal (Option 4)"<<endl;
		Cout<<"---------------------------------------------------"<<endl;
		MakeHexMeshFromArcGrid( infile );
		Cout<<"\nMakeHexMeshFromArcGrid Successful Using Option 4"<<endl<<flush;
	}
	
	else if( read == 5 ){
		Cout<<"\n\nPart 2: Creating Points from ArcInfo *net (Option 5)"<<endl; 
		Cout<<"---------------------------------------------------"<<endl;
		MakePointFromFileArcInfo( infile ); 
		
		Cout<<"\n----------------------------------------------"<<endl<<flush;
		Cout<<"\nMakePointFromFileArcInfo Successful Using Option 5...."<<endl<<flush;
		Cout<<"\t1. Check the points file for consistent boundary codes"<<endl<<flush;
		Cout<<"\t2. Rerun tRIBS with Option 2 using the points file"<<endl<<flush;
		Cout<<"\n---------------------------------------------"<<endl<<flush;
		exit(1);
	} 
	
	else if( read == 6 ){
		Cout<<"\n\nPart 2: Creating Points from ArcInfo *pnt (Option 6)"<<endl; 
		Cout<<"---------------------------------------------------"<<endl;
		MakePointFromFileArcInfoGen( infile ); 
		
		Cout<<"\n----------------------------------------------"<<endl<<flush;
		Cout<<"\nMakePointFromFileArcInfo Successfull Using Option 6...."<<endl<<flush;
		Cout<<"\t1. Check the points file for consistent boundary codes"<<endl<<flush;
		Cout<<"\t2. Rerun tRIBS with Option 2 using the points file"<<endl<<flush;
		Cout<<"\n----------------------------------------------"<<endl<<flush;
		exit(2);
	} 
	
	else if( read == 7){
		Cout<<"\n\nPart 2: Creating Points from Scratch (Option 7)"<<endl; 
		Cout<<"---------------------------------------------------"<<endl;
		MakeMeshFromScratch( infile ); 
		Cout<<"\nMakeMeshFromScratch Successful Using Option 7"<<endl<<flush;
	}
	
	else if( read == 8){
		Cout<<"\n\nPart 2: Creating Mesh from Points File (Option 8)"<<endl;
		Cout<<"---------------------------------------------------"<<endl;
		MakeMeshFromTriangulator( infile ); 
		Cout<<"\nMakeMeshFromTriangulator Successful Using Option 8"<<endl<<flush;
	}

	else if( read == 9){
		Cout<<"\n\nPart 2: Creating Mesh from MeshBuilder Files (Option 9)"<<endl;
		Cout<<"---------------------------------------------------"<<endl;
		MakeMeshFromMeshBuilder( infile ); 
		Cout<<"\nMakeMeshFromMeshBuilder Successful Using Option 9"<<endl<<flush;
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
**   tMesh::MakeMeshFromPoints
**
**   Constructs a mesh from a given set of (x,y,z,b) points.
**
**   The format of the file containing points is:
**        NP
**        X0 Y0 Z0 B0
**        X1 Y1 Z1 B1 ...etc.
**   where NP is the number of points in the file, X, Y, and Z are
**   x, y, and z coords, and B is the boundary code.
**
**   Calls: tInputFile::ReadItem, MakeCCWEdges(),
**          UpdateMesh(), CheckMeshConsistency()
**   Parameters: infile -- main parameter input file
**
**************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromPoints( tInputFile &infile )
{
	int i;                           // loop counter
	int numpts;                      // no. of points in mesh
	tArray<double> x, y, z;          // arrays of x, y, and z coordinates
	tArray<int> bnd;                 // array of boundary codes 
	char pointFilenm[80];            // name of file containing (x,y,z,b) data
	ifstream pointfile;              // the file (stream) itself
	double minx = 1e12, miny = 1e12, // minimum x and y coords
		maxx = 0, maxy=0,            // maximum x and y coords 
		dx, dy;                      // max width and height of region
	tSubNode tempnode( infile ),     // temporary node used in creating new pts
		*stp1, *stp2, *stp3;         // supertriangle vertices 
	
	
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
	
	for( i=0; i<numpts; i++ ){
		if( pointfile.eof() )
			cout << "\nReached end-of-file while reading points.\n" ;
		pointfile >> x[i] >> y[i] >> z[i] >> bnd[i];
		if( bnd[i]<0 || bnd[i]>3 ){
			cout << "\nInvalid boundary code.\n"<<endl;
			cout << "\n\nExiting Program..."<<endl;
			exit(2);
		}
		if( x[i]<minx ) minx = x[i];
		if( x[i]>maxx ) maxx = x[i];
		if( y[i]<miny ) miny = y[i];
		if( y[i]>maxy ) maxy = y[i];    
	}
	
	pointfile.close();
	dx = maxx - minx;
	dy = maxy - miny;
	
	//for(int j=0;j<numpts;j++){
	//  cout<<"\nx = "<<x[j]<<" y = "<<y[j]<<" z = "<<z[j]<<" b = "<<bnd[j];
	//}
	
	// Create the 3 nodes that form the supertriangle and place them on the
	// node list in counter-clockwise order. (Note that the base and height
	// of the supertriangle are 7 times the width and height, respectively,
	// of the rectangle that encloses the points.) Assigning the IDs allows 
	// us to retrieve and delete these nodes when we're done creating the mesh.
	
	Cout << "\nCreating supertriangle..."<< endl;
	
	tempnode.set3DCoords( minx-3*dx, miny-3*dy, 0.0 );
	tempnode.setBoundaryFlag( kClosedBoundary );
	tempnode.setID( -1 );
	nodeList.insertAtBack( tempnode );
	tempnode.set3DCoords( maxx+3*dx, miny-3*dy, 0.0 );
	tempnode.setID( -2 );
	nodeList.insertAtBack( tempnode );
	tempnode.set3DCoords( minx+0.5*dx, maxy+3*dy, 0.0 );
	tempnode.setID( -3 );
	nodeList.insertAtBack( tempnode );
	
	// Set # of nodes, edges, and triangles
	nnodes = 3;
	nedges = ntri = 0;
	
	// Create the edges that connect the supertriangle vertices and place
	// them on the edge list.
	
	tMeshListIter<tSubNode> nodIter( nodeList );
	stp1 = nodIter.FirstP();
	stp2 = nodIter.NextP();
	stp3 = nodIter.NextP();
	AddEdge( stp1, stp2, stp3 );  // edges 1->2 and 2->1
	AddEdge( stp2, stp3, stp1 );  // edges 2->3 and 3->2
	AddEdge( stp3, stp1, stp2 );  // edges 3->1 and 1->3
	
	// Set up the triangle itself and place it on the list. To do this, we
	// just set up a list of pointers to the three nodes in the super tri
	// and pass the list (along with an iterator) to MakeTriangle.
	
	tPtrList<tSubNode> supertriptlist;
	supertriptlist.insertAtBack( stp1 );
	supertriptlist.insertAtBack( stp2 );
	supertriptlist.insertAtBack( stp3 );
	supertriptlist.makeCircular();
	tPtrListIter<tSubNode> stpIter( supertriptlist );
	MakeTriangle( supertriptlist, stpIter );
	
	// Now add the points one by one to construct the mesh.
	// Added kNoUpdate based on Lancaster 
	
	Cout<<"\nAdding new nodes from points file..."<<endl;
	for( i=0; i<numpts; i++ ){
		tempnode.setID( i );
		tempnode.set3DCoords( x[i], y[i], z[i] );
		tempnode.setBoundaryFlag( bnd[i] );
		if(simCtrl->Verbose_label == 'Y'){
			//cout<<"Adding Node: "<<i+1<<" ..."<<endl<<flush;
		}
		//AddNode( tempnode, kNoUpdateMesh );  //Possibly switch based on Lancaster
		AddNode( tempnode );
	}
	
	// Remove supertriangle by deleting its vertices.
	
	Cout<<"\nDeleting supertriangle..."<<endl;
	DeleteNode( stp1, kNoRepair);  
	DeleteNode( stp2, kNoRepair);
	DeleteNode( stp3, kNoRepair);
	
	// Update Voronoi areas, edge lengths and test the consistency of mesh  
	Cout<<"\nTesting Mesh..."<<endl;
	UpdateMesh();
	CheckMeshConsistency( infile );  //Possibly Removed based on Lancaster
	supertriptlist.Flush();
}


//=========================================================================
//
//
//                  Section 5a: tMesh:: MakeMeshFromTriangulator()
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
//                  Section 5b: tMesh:: MakeMeshFromMeshBuilder()
//
//
//=========================================================================
                                                                                
                                                                                
/**************************************************************************
**
**   tMesh::MakeMeshFromMeshBuilder( infile )
**
**   Similar to tMesh::MakeMeshFromTriangulator using Tipper's triangulation
**   algorithm within the separate MeshBuilder code.  This mesh generator
**   is for large data sets.
**
**   Reads mesh which was already checked and enhanced with voronoi
**   information from data files.  Reads in basic tNode information
**   but creates a tCNode.
**
**************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromMeshBuilder( tInputFile &infile )
{
	int id, boundary, firstEdgeID, flood, tracer, flowEdgeID, streamNodeID;
	int origID, destID, ccwID, origBoundary, destBoundary, reach;
	double hillpath, traveltime, srf, hsrf, psrf, satsrf, sbsrf;
	double SoilMoistureSC, SoilMoistureUNSC, RootMoistureSC;
	double EvapoTranspiration, ContrArea, NwtNew, Rain, Curvature, streampath;
	double x, y, z, varea, rvtx[2], length, slope, vedglen;
	fstream nodeStr, edgeStr;

	// Read the node file setting all variables except the pointer to first edge
	// Since that can't be found until edges are read
	//
	nodeStr.open("nodes.meshb", ios::in | ios::binary);
	BinaryRead(nodeStr, nnodes);
	assert( nnodes > 0 );
	tSubNode curnode;

	for (int i = 0; i < nnodes; i++) {
		BinaryRead(nodeStr, id);
		BinaryRead(nodeStr, boundary);
		BinaryRead(nodeStr, x);
		BinaryRead(nodeStr, y);
		BinaryRead(nodeStr, z);
		BinaryRead(nodeStr, varea);
		BinaryRead(nodeStr, firstEdgeID);
		BinaryRead(nodeStr, flood);
		BinaryRead(nodeStr, tracer);
		BinaryRead(nodeStr, hillpath);
		BinaryRead(nodeStr, traveltime);
		BinaryRead(nodeStr, srf);
		BinaryRead(nodeStr, hsrf);
		BinaryRead(nodeStr, psrf);
		BinaryRead(nodeStr, satsrf);
		BinaryRead(nodeStr, sbsrf);
		BinaryRead(nodeStr, EvapoTranspiration);
		BinaryRead(nodeStr, SoilMoistureSC);
		BinaryRead(nodeStr, SoilMoistureUNSC);
		BinaryRead(nodeStr, RootMoistureSC);
		BinaryRead(nodeStr, ContrArea);
		BinaryRead(nodeStr, NwtNew);
		BinaryRead(nodeStr, Rain);
		BinaryRead(nodeStr, Curvature);
		BinaryRead(nodeStr, streampath);
		BinaryRead(nodeStr, reach);
		BinaryRead(nodeStr, flowEdgeID);
		BinaryRead(nodeStr, streamNodeID);

		curnode.setID(id);
		curnode.setBoundaryFlag(boundary);
		curnode.setX(x);
		curnode.setY(y);
		curnode.setZ(z);
		curnode.setVArea(varea);
		curnode.setVArea_Rcp(1.0 / varea);
		curnode.setFloodStatus(flood);
		curnode.setTracer(tracer);
		curnode.setHillPath(hillpath);
		curnode.setTTime(traveltime);
		curnode.setsrf(srf);
		curnode.sethsrf(hsrf);
		curnode.setpsrf(psrf);
		curnode.setsatsrf(satsrf);
		curnode.setsbsrf(sbsrf);
		curnode.setEvapoTrans(EvapoTranspiration);
		curnode.setSoilMoistureSC(SoilMoistureSC);
		curnode.setSoilMoistureUNSC(SoilMoistureUNSC);
		curnode.setRootMoistureSC(RootMoistureSC);
		curnode.setContrArea(ContrArea);
		curnode.setNwtNew(NwtNew);
		curnode.setRain(Rain);
		curnode.setCurvature(Curvature);
		curnode.setStreamPath(streampath);
		curnode.setReach(reach);

		if ((boundary == 0) || (boundary==3)) {
			nodeList.insertAtActiveBack(curnode);
		} else if (boundary == kOpenBoundary) {
			nodeList.insertAtBoundFront(curnode);
		} else {
			nodeList.insertAtBack(curnode);       //kClosedBoundary
		}
	}
	nodeStr.close();

	// Create the lookup table of nodes indexed by node id
	// Used when assigning node pointers to edge origin and destination
	//
	NodeTable = new tSubNode*[nnodes];
	tSubNode* cn;
	tMeshListIter< tSubNode > nodIter( nodeList );
	for (cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP())
		NodeTable[cn->getID()] = cn;

	// Read the edge file setting all variables except ccw edge which can't
	// be found until all edges are read once
	//
	edgeStr.open("edges.meshb", ios::in | ios::binary);
	edgeStr.read((char*) &nedges, sizeof(int));
	assert( nedges > 0 );
	tEdge curedge;
	tArray<double> RVtx(2);

	for (int i = 0; i < nedges; i++) {
		BinaryRead(edgeStr, id);
		BinaryRead(edgeStr, rvtx[0]);
		BinaryRead(edgeStr, rvtx[1]);
		BinaryRead(edgeStr, length);
		BinaryRead(edgeStr, slope);
		BinaryRead(edgeStr, vedglen);
		BinaryRead(edgeStr, origID);
		BinaryRead(edgeStr, destID);
		BinaryRead(edgeStr, ccwID);

		curedge.setID(id);
		RVtx[0] = rvtx[0];
		RVtx[1] = rvtx[1];
		curedge.setRVtx(RVtx);
		curedge.setLength(length);
		curedge.setSlope(slope);
		curedge.setVEdgLen(vedglen);

		// Look up the node origin and destination for the given id
		tSubNode *origNode = NodeTable[origID];
		curedge.setOriginPtr(origNode);
		origBoundary = origNode->getBoundaryFlag();
       
		tSubNode *destNode = NodeTable[destID];
		curedge.setDestinationPtr(destNode);
		destBoundary = destNode->getBoundaryFlag();

		// set the "flowallowed" status (FALSE if either endpoint is a
		// closed boundary, or both are open boundaries)
		// and insert edge pair onto the list --- active
		// part of list if flow is allowed, inactive if not

		if (origBoundary == kClosedBoundary || 
		    destBoundary == kClosedBoundary || 
		    (origBoundary == kOpenBoundary && destBoundary == kOpenBoundary) ) {
			curedge.setFlowAllowed(0);
			edgeList.insertAtBack(curedge);
		}
		else {
			curedge.setFlowAllowed(1);
			edgeList.insertAtActiveBack(curedge);
		}
	}
	edgeStr.close();

	// Create the lookup table of edges indexed by edge id
	tEdge** EdgeTable = new tEdge*[nedges];
	tEdge* ce;
	tMeshListIter< tEdge > edgIter( edgeList );
	for (ce = edgIter.FirstP(); !(edgIter.AtEnd()); ce = edgIter.NextP())
		EdgeTable[ce->getID()] = ce;

	// Reread nodes file to get the id of the first edge which we now
	// can get a pointer to because of EdgeTable
	//
	nodeStr.open("nodes.meshb", ios::in | ios::binary);
	BinaryRead(nodeStr, nnodes);

	for (int i = 0; i < nnodes; i++) {
		BinaryRead(nodeStr, id);
		BinaryRead(nodeStr, boundary);
		BinaryRead(nodeStr, x);
		BinaryRead(nodeStr, y);
		BinaryRead(nodeStr, z);
		BinaryRead(nodeStr, varea);
		BinaryRead(nodeStr, firstEdgeID);

		BinaryRead(nodeStr, flood);
		BinaryRead(nodeStr, tracer);
		BinaryRead(nodeStr, hillpath);
		BinaryRead(nodeStr, traveltime);
		BinaryRead(nodeStr, srf);
		BinaryRead(nodeStr, hsrf);
		BinaryRead(nodeStr, psrf);
		BinaryRead(nodeStr, satsrf);
		BinaryRead(nodeStr, sbsrf);
		BinaryRead(nodeStr, EvapoTranspiration);
		BinaryRead(nodeStr, SoilMoistureSC);
		BinaryRead(nodeStr, SoilMoistureUNSC);
		BinaryRead(nodeStr, RootMoistureSC);
		BinaryRead(nodeStr, ContrArea);
		BinaryRead(nodeStr, NwtNew);
		BinaryRead(nodeStr, Rain);
		BinaryRead(nodeStr, Curvature);
		BinaryRead(nodeStr, streampath);
		BinaryRead(nodeStr, reach);
		BinaryRead(nodeStr, flowEdgeID);
		BinaryRead(nodeStr, streamNodeID);
                                                                                
		// Index of the first spoke coming off this node
		NodeTable[id]->setEdg(EdgeTable[firstEdgeID]);

		// Index of flow edge
		if (flowEdgeID >= 0)
			NodeTable[id]->setFlowEdg(EdgeTable[flowEdgeID]);

		// Index of stream node
		if (streamNodeID >= 0)
			NodeTable[id]->setStreamNode(NodeTable[streamNodeID]);
	}
	nodeStr.close();

	// Reread edge file to get the pointer to the counter clockwise edge
	// which we can not get a pointer to because of EdgeTable
	edgeStr.open("edges.meshb", ios::in | ios::binary);
	BinaryRead(edgeStr, nedges);
	assert( nedges > 0 );

	tEdge *curedg, *ccwedg;

	for (int i = 0; i < nedges; i++) {
		BinaryRead(edgeStr, id);
		BinaryRead(edgeStr, rvtx[0]);
		BinaryRead(edgeStr, rvtx[1]);
		BinaryRead(edgeStr, length);
		BinaryRead(edgeStr, slope);
		BinaryRead(edgeStr, vedglen);
		BinaryRead(edgeStr, origID);
		BinaryRead(edgeStr, destID);
		BinaryRead(edgeStr, ccwID);

		if (ccwID >= 0 && ccwID < nedges)
			EdgeTable[id]->setCCWEdg(EdgeTable[ccwID]);
		else
			cout << "Edge " << id << " has bad ccwedge "
			     << ccwID << " should be 0 to " << nedges << endl;
	}
	edgeStr.close();

	// Can't delete NodeTable because tFlowNet and tGraph must look up pointers
	//delete [] NodeTable;
	delete [] EdgeTable;
}

//=========================================================================
//
//
//                  Section 6: tMesh:: MakePointFromFileArcInfo( )
//                                     MakePointFromFileArcInfoGen( )
//
//=========================================================================

/**************************************************************************** 
**   
**    tMesh::MakePointFromFileArcInfo
**
**    Routine to convert ArcInfo files .net (NODES, EDGES, TRIANGLES)
**    describing a TIN  (generated by ungeneratetin function) to a 
**    file of points ( x, y, z, boundary) readable by the function 
**    MakeMeshFromPoints
**   
**    calls: tInputFile::ReadItem
**    Created: 10/2000 VT
**
**
******************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakePointFromFileArcInfo( tInputFile &infile ){
	tArray<double> n, x, y, z;        // arrays of x, y, and z  coordinates and n
	double minx,maxx,miny,maxy,minz;
	minx=1e12; miny=1e12;minz=1e12;
	maxx=0; maxy=0;
	tArray<int> nbbound;		    // array of number of boundary points
	char arcFilenm[kName];               // name of file containing data from ArcInfo
	char pointFilenm[kName];             // name of output file containing points
	char pointText[kName];
	ifstream filenet;
	
	Cout<<"\nReading Arc/Info *.net base name..."<<endl;    
	infile.ReadItem( arcFilenm, "ARCINFOFILENAME" );
	
	Cout<<"\nReading Point file base name..."<<endl;
	infile.ReadItem( pointFilenm, "POINTFILENAME" );
	ofstream outfile (pointFilenm);
	
	strcpy(pointText, pointFilenm);
	strcat(pointText, ".txt");
	ofstream outtextfile (pointText);
	
	// Open the file .net from Arcinfo
	
	filenet.open(strcat(arcFilenm,".net"));
	if( !filenet.good()){
		cerr << "\nError: File name: '" << arcFilenm << "' not found.";
		cerr << "\nExiting program...";
		exit(1);
	}
	
	char oneline[90];
	
	//Read the first word and check it is "NODES"
	filenet >> oneline;  
	if (strcasecmp(oneline, "NODES"))
		cerr<<"\nFirst line should begin with the word NODES \n";
	
	//Find the numbers of nodes 
	filenet >> oneline;
	int nbpnts=0; 
	while (strcasecmp(oneline, "EDGES")){
		char  u[80],v[80],w[80];
		nbpnts++;
		filenet>>u>>v>>w;
		filenet>>oneline;   
	}
	filenet.close();
	
	//Open again the file to extract de coordinates
	filenet.open( arcFilenm );
	filenet >> oneline; 
	n.setSize( nbpnts );
	x.setSize( nbpnts );
	y.setSize( nbpnts ); 
	z.setSize( nbpnts );
	nbbound.setSize( nbpnts );
	
	//Fill the arrays (x,y,z) with the coordinates of the nodes  
	for(int  i=0; i<nbpnts; i++ ){
		if(filenet.eof())
			cerr << "Reached end-of-file while reading points.";
		filenet >> n[i] >> x[i] >> y[i] >> z[i];    
		if( x[i]<minx ) minx = x[i];
		if( x[i]>maxx ) maxx = x[i];
		if( y[i]<miny ) miny = y[i];
		if( y[i]>maxy ) maxy = y[i];
		if( z[i]<minz ) minz = z[i];
	}
	
	filenet >> oneline;
	if(strcasecmp(oneline, "EDGES"))
		cerr<<"cannot find the EDGES \n";
	
	// Store in the array nbbound the # of the boundary nodes 
	// in file.net from ARC, the edges have boundary codes:
	//  boundary= 0, interior=1
	
	int  p1, p2, bd, ct, ict;
	filenet >> oneline;
	ict = 0; ct=0;
	int *parray = new int[nbpnts]; 
	while (strcasecmp(oneline, "TRIANGLES")){
		filenet >> p1>> p2>> bd;
		int cflag, cflag2;
		cflag = 0;cflag2 = 0;
		for(int j=0;j<ct;j++){           //parray stores nodes with valid edges
			if(parray[j]== p1) cflag=1;
			if(parray[j]== p2) cflag2=1;
		}
		if(cflag == 0) {parray[ct] = p1; ct++;}
		if(cflag2 == 0) {parray[ct] = p2; ct++;}
		if (bd ==  0){ 
			int flag, flag2;
			flag = 0;flag2 = 0;
			for (int k=0; k<ict; k++){ 
				if (nbbound[k]==p1) flag =1;
				if (nbbound[k]==p2) flag2 =1;
			}
			if (flag == 0) {nbbound[ict] = p1; ict++;}
			if (flag2 == 0) {nbbound[ict] = p2; ict++;}
		}
		filenet >> oneline;  	 
	}
	
	int nb= ict;
	int cc = ct;
	filenet.close();
	
	// Write the new file with number of nodes and for each node, 
	// the coordinates x,y,z and the boundary code => 0:interior, 1 boundary
	
	outfile<<cc<<endl;
	outtextfile<<"x\t"<<"y\t"<<"z\t"<<"b\n";
	
	for(int ikt=0; ikt<nbpnts; ikt++ )
	{ 
		for(int k=0;k<cc;k++){  
			if(n[ikt]==parray[k]){
				outfile<< setprecision(14)<< x[ikt]<< " ";
				outfile<< setprecision(14)<< y[ikt]<< " ";
				outfile<< z[ikt]<< " ";
				outtextfile<< setprecision(14)<< x[ikt]<< "\t";
				outtextfile<< setprecision(14)<< y[ikt]<< "\t";
				outtextfile<< z[ikt]<< "\t";
				int  flag {};
				for(int j=0; j<nb; j++ ){
					flag =0;
					if (n[ikt]==nbbound[j]){       
						if (z[ikt]!= minz) flag =1;
						else flag =2;   //outlet
						j=nb;
					}
				}
				outfile <<flag << endl;  
				outtextfile<<flag<<endl;
				k=cc;
			}
		}
	}
	outfile.close();
	outtextfile.close();
	delete [] parray; //TODO delete[]? WR 'delete' applied to a pointer that was allocated with 'new[]'; did you mean 'delete[]'? [-Wmismatched-new-delete]
}

/***************************************************************************** 
**    
**   tMesh::MakePointFromFileArcInfoGen
**
**    Routine to convert ArcInfo files .lin and .pnt (LINES and POINTS)
**    describing a TIN  (generated by ungeneratetin function using generate)
**    to a file of points ( x, y, z, boundary) readable by the function 
**    MakeMeshFromPoints. Used in cases where the stream network must
**    be inputed as a separate boundary code.
**
**
**    calls: tInputFile::ReadItem
**
******************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakePointFromFileArcInfoGen( tInputFile &infile )
{
	long double *xstream, *ystream, *zstream;  
	long double *xbound, *ybound, *zbound;
	long double *xinterior, *yinterior, *zinterior;
	int *bstream, *bbound, *binterior;
	char arcFilenm[kName];                    // base name of Arc/Info 
	char pointFilenm[kName];                  // base name of output point file 
	char lineFile1[kName], lineFile2[kName];     // names of Arc/Info *.lin
	char pntFile1[kName], pntFile2[kName];       // names of Arc/Info *.pnt
	char oneline[25];
	char pointText[kName];
	ifstream filelin, filelin2, filelin3, filelin4, filepnt, filepnt2; 
	int i, j, k;
	
	Cout<<"\nReading Arc/Info *.pnt and *.lin base name..."<<endl;    
	infile.ReadItem( arcFilenm, "ARCINFOFILENAME" ); 
	
	Cout<<"\nReading Point file base name..."<<endl;
	infile.ReadItem( pointFilenm, "POINTFILENAME" );
	ofstream outfile (pointFilenm);
	
	strcpy(pointText, pointFilenm);
	strcat(pointText, ".txt");
	ofstream outtextfile (pointText);
	
	// Open the files .lin and .pnt and concatenate extensions 
	
	strcpy(lineFile1, arcFilenm);
	strcpy(pntFile1, arcFilenm);
	strcpy(lineFile2, arcFilenm);
	strcpy(pntFile2, arcFilenm); 
	
	//Extract number of interior points: Code 1 (Mass Points)
	
	filepnt.open(strcat(pntFile2,".pnt2"));
	if( !filepnt.good()){
		cerr << "File name: " << arcFilenm << ".pnt2 not found!" << "\n";  
		cerr << "Corresponds to Double Ring TIN with Stream Network"<<endl<<flush;
		exit(2);
	}
	
	int nbpnts=0; 
	filepnt >> oneline; 
	while(strcmp(oneline, "END")!=0){
		char *x = new char[25]; 
		char *y = new char[25];
		char *z = new char[25];
		nbpnts++;
		filepnt>>x>>y>>z;
		filepnt>>oneline; 
		delete [] x;
        delete [] y;
        delete [] z;
	}
	filepnt.close();
	
	//Extract number of boundary points: Code 7 (HardClip)
	
	filelin.open(strcat(lineFile2,".lin2"));
	if( !filelin.good()){
		cerr << "File name: " << arcFilenm << ".lin2 not found" << "\n";
		cerr << "Corresponds to Double Ring TIN with Stream Network 1"<<endl<<flush;
		exit(1);
	}
	
	filelin >> oneline;  
	while(strcmp(oneline,"7")!=0){
		filelin>>oneline; }
	
	int numbound = 0;
	while(strcmp(oneline,"END")!=0){
		char *x = new char[25];
		char *y = new char[25];
		char *z = new char[25];
		filelin>>x>>y>>z;
		numbound++; 
		for(i=0;i<25;i++)
        {
            oneline[i] = x[i];
        }
        delete [] x;
        delete [] y;
        delete [] z;
	}
	numbound--;      
	filelin.close();
	
	//Extract Boundary Point Coordinates: Code 7 (HardClip)
	
	filelin2.open(lineFile2);
	if( !filelin2.good()){
		cerr << "File name: " << arcFilenm << ".lin2 not found!" << "\n";
		cerr << "Corresponds to Double Ring TIN with Stream Network 2"<<endl<<flush;
		exit(1);
	}
	
	filelin2>>oneline;
	while(strcmp(oneline,"7")!=0){
		filelin2>>oneline; }
	
	char **pbx, **pby, **pbz;
	pbx = new char*[nbpnts];   //Hold boundary node x,y,z as char
	pby = new char*[nbpnts];
	pbz = new char*[nbpnts];
	
	//Take care of repeting bound point (first and last)
	numbound--;
	for(i=0;i<numbound;i++){  
		pbx[i] = new char[25];
		pby[i] = new char[25]; 
		pbz[i] = new char[25];
		filelin2>>pbx[i]>>pby[i]>>pbz[i]; 
	}
	filelin2.close();
	
	//Extract Stream Point Coordinates: Code 3 (Hardline)
	
	//Read and save lines that include inner ring only
	
	filelin3.open(strcat(lineFile1,".lin1"));
	if( !filelin3.good() ){
		cout << "File name: " << arcFilenm << ".lin1" << " not found!\n";
		cout << "Corresponds to Double Ring TIN without Stream Network 1"<<endl<<flush;
		cout <<"\nExiting Program..."<<endl;
		exit(1);
	}
	
	char **pdx, **pdy, **pdz;
	char **rdx, **rdy, **rdz;
	pdx = new char*[nbpnts];   //Hold inner ring node x,y,z as char
	pdy = new char*[nbpnts];
	pdz = new char*[nbpnts];
	rdx = new char*[nbpnts];   //Read inner ring node x,y,z as char
	rdy = new char*[nbpnts];
	rdz = new char*[nbpnts];
	int numring = 0;
	int cdt = 0; int cdta = 0;
	int flagd = 0;
	
	filelin3 >> oneline;  
	while(strcmp(oneline,"7")!=0){ 
		for(j=0;j<2;j++){  
			flagd = 0;    
			rdx[cdta] = new char[25]; rdy[cdta] = new char[25]; rdz[cdta] = new char[25]; 
			filelin3>>rdx[cdta]>>rdy[cdta]>>rdz[cdta];
			for(i=0;i<cdt;i++){
				if(strcmp(rdx[cdta],pdx[i])==0&&strcmp(rdy[cdta],pdy[i])==0&&strcmp(rdz[cdta],pdz[i])==0){ 
					flagd = 1;
				}
			} 
			if(flagd==0){
				pdx[cdt] = new char[25]; pdy[cdt] = new char[25]; pdz[cdt] = new char[25]; 
				pdx[cdt] = rdx[cdta]; pdy[cdt] = rdy[cdta]; pdz[cdt] = rdz[cdta];
				cdt++;} 
			cdta++;
		}
		filelin3 >> oneline;
		filelin3 >> oneline;
	}
	numring = cdt;   //number of double ring nodes
	filelin3.close();
	
	//Save only stream lines by eliminating inner ring lines, boundary nodes
	
	filelin4.open(lineFile2);
	if( !filelin4.good()){
		cout << "File name: " << arcFilenm << ".lin2" << "not found!\n";
		cout << "Corresponds to Double Ring TIN with Stream Network"<<endl<<flush;
		cout <<"\nExiting Program..."<<endl;
		exit(1);
	}
	
	filelin4 >> oneline;  
	int numstream = 0; 
	int flag, flag2, flagb; 
	int ct = 0; int cta = 0;
	char **psx, **psy, **psz;
	char **rx, **ry, **rz;
	psx = new char*[nbpnts];   //Hold stream node x,y,z as char
	psy = new char*[nbpnts];
	psz = new char*[nbpnts];
	rx = new char*[nbpnts];    //Read stream node x,y,z as char
	ry = new char*[nbpnts];
	rz = new char*[nbpnts];
	
	while(strcmp(oneline,"7")!=0){
		for(j=0;j<2;j++){   
			flag=0; flagb=0; flagd=0;   
			rx[cta] = new char[25]; ry[cta] = new char[25]; rz[cta] = new char[25]; 
			filelin4>>rx[cta]>>ry[cta]>>rz[cta];
			for(i=0;i<ct;i++){
				if(strcmp(rx[cta],psx[i])==0&&strcmp(ry[cta],psy[i])==0&&strcmp(rz[cta],psz[i])==0){ 
					flag = 1;
				}
			}
			for(i=0;i<numbound;i++){
				if(strcmp(rx[cta],pbx[i])==0&&strcmp(ry[cta],pby[i])==0&&strcmp(rz[cta],pbz[i])==0){
					flagb = 1;
				}
			}
			for(i=0;i<numring;i++){
				if(strcmp(rx[cta],pdx[i])==0&&strcmp(ry[cta],pdy[i])==0&&strcmp(rz[cta],pdz[i])==0){
					flagd = 1;
				}
			}

			if(flag==0 && flagb==0 && flagd==0){
				psx[ct] = new char[25]; psy[ct] = new char[25]; psz[ct] = new char[25]; 
				psx[ct]=rx[cta]; psy[ct]=ry[cta]; psz[ct]=rz[cta];
				ct++;
            }
            cta++;
		}   
		filelin4 >> oneline;
		filelin4 >> oneline;
	}
	numstream = ct;   //number of stream nodes
	filelin4.close();
	
	// Extract Interior Point Coordinates from *.pnt 
	// Exclude Stream and Boundary Points 
	
	filepnt2.open(pntFile2); 
	if( !filepnt2.good()){
		cerr << "File name: " << arcFilenm << ".pnt" << "\n";  
		exit(3);
	}
	
	char **pix, **piy, **piz;
	char **rix, **riy, **riz;
	pix = new char*[nbpnts];   //Hold interior node x,y,z as char
	piy = new char*[nbpnts];
	piz = new char*[nbpnts];
	rix = new char*[nbpnts];   //Read interior node x,y,z as char
	riy = new char*[nbpnts];
	riz = new char*[nbpnts];
	int onept; 
	int cnt = 0;
	
	for(i=0;i<nbpnts;i++){
		rix[i] = new char[25]; riy[i] = new char[25]; riz[i] = new char[25];
		flag = 0; flag2 = 0;
		filepnt2 >> onept >> rix[i] >> riy[i] >> riz[i];
		for(j=0;j<numstream;j++){
			if(strcmp(rix[i],psx[j])==0&&strcmp(riy[i],psy[j])==0&&strcmp(riz[i],psz[j])==0){ 
				flag = 1;
			}
		}
		for(k=0;k<numbound;k++){
			if(strcmp(rix[i],pbx[k])==0&&strcmp(riy[i],pby[k])==0&&strcmp(riz[i],pbz[k])==0){ 
				flag2 = 1;
			}
		}
		if(flag==0 && flag2==0){
			pix[cnt] = new char[25]; piy[cnt] = new char[25]; piz[cnt] = new char[25];
			pix[cnt]=rix[i]; piy[cnt]=riy[i]; piz[cnt]=riz[i];
			cnt++;
		}
	}
	int numinterior = cnt;
	filepnt2.close();
	
	// Convert chars to doubles
	
	xstream = new long double[numstream];
	ystream = new long double[numstream];
	zstream = new long double[numstream];
	xbound = new long double[numbound];
	ybound = new long double[numbound];
	zbound = new long double[numbound];
	xinterior = new long double[numinterior];
	yinterior = new long double[numinterior];
	zinterior = new long double[numinterior];
	bstream = new int[numstream];
	bbound = new int[numbound];
	binterior = new int[numinterior];
	
	char tKD[] = "D"; char tKE[] = "E";
	char *cnumx = new char[20]; char *cnumy = new char[20]; char *cnumz = new char[20];
	char *enumx = new char[5]; char *enumy = new char[5]; char *enumz = new char[5];
	
	for(i=0;i<numstream;i++){
		cnumx = strtok(psx[i],tKD);
		enumx = strtok(NULL,tKD);
		xstream[i] = atof(cnumx)*pow((double)10.,atof(enumx));
		cnumy = strtok(psy[i],tKD);
		enumy = strtok(NULL,tKD);
		ystream[i] = atof(cnumy)*pow((double)10.,atof(enumy));
		cnumz = strtok(psz[i],tKE);
		enumz = strtok(NULL,tKE);
		zstream[i] = atof(cnumz)*pow((double)10.,atof(enumz));
		bstream[i] = 3;
	}
	
	double minz = 1e12;
	int minbound = 0;
	for(i=0;i<numbound;i++){
		cnumx = strtok(pbx[i],tKD);
		enumx = strtok(NULL,tKD);
		xbound[i] = atof(cnumx)*pow((double)10.,atof(enumx));
		cnumy = strtok(pby[i],tKD);
		enumy = strtok(NULL,tKD);
		ybound[i] = atof(cnumy)*pow((double)10.,atof(enumy));
		cnumz = strtok(pbz[i],tKE);
		enumz = strtok(NULL,tKE);
		zbound[i] = atof(cnumz)*pow((double)10.,atof(enumz));
		if( zbound[i]<minz ){
			minz = zbound[i];
			minbound = i;
		}
		bbound[i] = 1;
	}
	bbound[minbound] = 2;
	
	for(i=0;i<numinterior;i++){
		cnumx = strtok(pix[i],tKD);
		enumx = strtok(NULL,tKD);
		xinterior[i] = atof(cnumx)*pow((double)10.,atof(enumx));
		cnumy = strtok(piy[i],tKD);
		enumy = strtok(NULL,tKD);
		yinterior[i] = atof(cnumy)*pow((double)10.,atof(enumy));
		cnumz = strtok(piz[i],tKE);
		enumz = strtok(NULL,tKE);
		zinterior[i] = atof(cnumz)*pow((double)10.,atof(enumz));
		binterior[i] = 0;
	}
	
	// Write to *.points File and *.points.txt files
	
	int totalnum = numstream+numbound+numinterior;  
	outfile << totalnum <<endl;
	outtextfile<<"x\t"<<"y\t"<<"z\t"<<"b\n";
	
	for(i=0; i<numinterior;i++){
		outfile << setprecision(14) << xinterior[i] << " ";
		outfile << setprecision(14) << yinterior[i] << " ";
		outfile << setprecision(14) << zinterior[i] << " ";
		outfile << binterior[i] << endl;
		
		outtextfile << setprecision(14) << xinterior[i] << "\t";
		outtextfile << setprecision(14) << yinterior[i] << "\t";
		outtextfile << setprecision(14) << zinterior[i] << "\t";
		outtextfile << binterior[i] << endl;
	}
	
	
	for(i=0;i<numstream;i++){
		outfile << setprecision(14) << xstream[i] << " ";
		outfile << setprecision(14) << ystream[i] << " ";
		outfile << setprecision(14) << zstream[i] << " ";
		outfile << bstream[i] << endl;
		
		outtextfile << setprecision(14) << xstream[i] << "\t";
		outtextfile << setprecision(14) << ystream[i] << "\t";
		outtextfile << setprecision(14) << zstream[i] << "\t";
		outtextfile << bstream[i] << endl;
	}
	
	int skpone2 = 0;
	for(i=0;i<numbound;i++){
		if(zbound[i] == zbound[minbound]){
			skpone2 = 1;
		}
		if(skpone2 == 0){
			outfile << setprecision(14) << xbound[i] << " ";
			outfile << setprecision(14) << ybound[i] << " ";
			outfile << setprecision(14) << zbound[i] << " ";
			outfile << bbound[i] << endl;
			
			outtextfile << setprecision(14) << xbound[i] << "\t";
			outtextfile << setprecision(14) << ybound[i] << "\t";
			outtextfile << setprecision(14) << zbound[i] << "\t";
			outtextfile << bbound[i] << endl;
			
		}
		else if(skpone2 == 1&&i!=numbound-1){
			outfile << setprecision(14) << xbound[i+1] << " ";
			outfile << setprecision(14) << ybound[i+1] << " ";
			outfile << setprecision(14) << zbound[i+1] << " ";
			outfile << bbound[i+1] << endl;
			
			outtextfile << setprecision(14) << xbound[i+1] << "\t";
			outtextfile << setprecision(14) << ybound[i+1] << "\t";
			outtextfile << setprecision(14) << zbound[i+1] << "\t";
			outtextfile << bbound[i+1] << endl;
		}
	}
	
	outfile << setprecision(14) << xbound[minbound] << " ";
	outfile << setprecision(14) << ybound[minbound] << " ";
	outfile << setprecision(14) << zbound[minbound] << " ";
	outfile << setprecision(14) << bbound[minbound] << " ";
	
	outtextfile << setprecision(14) << xbound[minbound] << "\t";
	outtextfile << setprecision(14) << ybound[minbound] << "\t";
	outtextfile << setprecision(14) << zbound[minbound] << "\t";
	outtextfile << setprecision(14) << bbound[minbound] << "\n";
	
	outfile.close();
	outtextfile.close();

    // TODO - WR there is some odd behavior deleting variables here, mostly because of the pointers pointing to pointers, could update to make more clear and make sure is up to date with current practices
    //delete rx, ry, rz first since psx... etc point to them first
    for(i=0; i < cta; i++){
        delete [] rx[i];
        delete [] ry[i];
        delete [] rz[i];
    }

	for(i=0;i<numbound;i++){
		delete pbx[i];
        delete pby[i];
        delete pbz[i];
	}
	for(i=0;i<nbpnts;i++){
		delete rix[i];
        delete riy[i];
        delete riz[i];
	}
	for(i=0;i<cdta;i++){
		delete rdx[i];
        delete rdy[i];
        delete rdz[i];
	}

	delete [] psx; delete[]  psy; delete[] psz;
	delete [] rx; delete[] ry; delete[] rz;
	delete [] pbx; delete[] pby; delete[] pbz;
	delete [] pix; delete[] piy; delete[] piz; 
	delete [] rix; delete[] riy; delete[] riz; 
	delete [] pdx; delete[] pdy; delete[] pdz;
	delete [] rdx; delete[] rdy; delete[] rdz;
    delete xstream;
    delete ystream;
    delete zstream;
    delete xbound;
    delete ybound;
    delete zbound;
    delete xinterior;
    delete yinterior;
    delete zinterior;
    delete bstream;
    delete bbound;
    delete binterior;
// WR--09192023: replaced with above as left operand of comma operator has no effect
//    delete xstream, ystream, zstream, xbound, ybound, zbound; //
//	delete xinterior, yinterior, zinterior;
//	delete bstream, bbound, binterior;
//	delete[] cnumx;
//  delete[] cnumy;
//  delete[] cnumz;
//  delete[] enumx;
//  delete[] enumy;
//  delete[] enumz;
	
	return;
}

//=========================================================================
//
//
//                  Section 7: tMesh:: MakeMeshFromScratch( )
//
//
//=========================================================================

/**************************************************************************
**
**   tMesh::MakeMeshFromScratch( infile )
**
**   Formerly tMesh( infile ). Makes mesh from scratch; reads parameters
**   from input file to get mesh size, spacing, method of point placement.
**
**      Could probably be done more gracefully, but here's how it does it:
**        1) makes boundary nodes and edges between them;
**        2) makes triangulation with only boundary nodes;
**        3) adds the rest of the nodes one at a time, i.e., triangulation
**           redone for each node added.
**
**   Created: 2/11/98, SL
**   Calls: MakeCCWEdges(), UpdateMesh(), CheckMeshConsistency()
**   Parameters: xGrid, yGrid, boundType, mElev, ptPlace, delGrid, numPts,
**               upperZ, xout, yout
**
**   TYP_BOUND keyword used to read type of open boundary:
**     0 = corner outlet (lower left)
**     1 = one open side (lower)
**     2 = two opposite sides (upper and lower)
**     3 = all sides
**     4 = specify outlet coordinates
**
**   OPT_PT_PLACE keyword used to read method of point placement:
**     0 = uniform grid;
**     1 = perturbed grid;
**     2 = random scatter;
**
**************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromScratch( tInputFile &infile )
{
	int i, j,                     // counters 
	n, nx, ny;                // no. of nodes along a side
	int numPts;                   // total no. of interior pts (if random)
	double dist;                  // current distance along boundary
	double delGrid,               // average spacing between nodes
		slope;
	double upperZ;
	double xout=0, yout=0;        // coordinates of user-specified outlet
	tArray< double > xyz(3);
	tSubNode tempnode( infile ),  // temporary node used to create node list
		*node0, *node1, *node2;
	tMeshListIter< tEdge > edgIter( edgeList );
	tMeshListIter< tSubNode > nodIter( nodeList );
	tPtrList< tSubNode > bndList;
	
	
	Cout<<"\nReading in Mesh Parameters..."<<endl;
	
	seed = infile.ReadItem( seed, "SEED" );
	double xGrid = infile.ReadItem( xGrid, "X_GRID_SIZE" );
	double yGrid = infile.ReadItem( yGrid, "Y_GRID_SIZE" );
	int boundType = infile.ReadItem( boundType, "TYP_BOUND" );
	
	int kSloped = 0;
	if(boundType == kOpenSide){
		kSloped = infile.ReadItem( kSloped, "SLOPED_SURF" );
		if(kSloped)
			upperZ = infile.ReadItem( upperZ, "UPPER_BOUND_Z" );
	}
	
	double mElev = infile.ReadItem( mElev, "MEAN_ELEV" );
	int ptPlace = infile.ReadItem( ptPlace, "OPT_PT_PLACE" );
	
	if( ptPlace == kUniformMesh || ptPlace == kPerturbedMesh ){
		delGrid = infile.ReadItem( delGrid, "GRID_SPACING" );
		if( delGrid >= xGrid || delGrid >= yGrid )
			cerr<< "Mesh point spacing must be smaller than total mesh width.";
	}
	else{
		numPts = infile.ReadItem( numPts, "NUM_PTS" );
		delGrid = sqrt( xGrid * yGrid / numPts );
	}
	
	// MAKE BOUNDARY
	
	Cout<<"\nMaking Boundary..."<<endl;
	if( boundType == kCornerOutlet ){
		miNextNodeID = 0;
		tempnode.setBoundaryFlag( kOpenBoundary );
		tempnode.set3DCoords( 0, 0, 0 );
		tempnode.setID( miNextNodeID );
		n = (int)(xGrid / delGrid); 
		tempnode.setBoundaryFlag( kOpenBoundary );
		nodeList.insertAtBack( tempnode );
		bndList.insertAtBack( nodIter.LastP() );
		tempnode.setBoundaryFlag( kClosedBoundary );
		for( i=1, miNextNodeID++; i<n; i++, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, 0, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		n = (int)(yGrid / delGrid);
		for( i=0; i<n; i++, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( xGrid, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		n = (int)(xGrid / delGrid);
		for( i=n; i>0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, yGrid, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		n = (int)(yGrid / delGrid); 
		for( i=n; i>0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( 0, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
	}
	else if( boundType == kOpenSide ){
		n = (int)(xGrid / delGrid);
		tempnode.setBoundaryFlag( kOpenBoundary );
		for( i=1, miNextNodeID=0; i<n; i++, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, 0, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		tempnode.setBoundaryFlag( kClosedBoundary );
		n = (int)(yGrid / delGrid);
		for( i=0; i<n; i++, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( xGrid, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		n = (int)(xGrid / delGrid);
		for( i=n; i>0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, yGrid, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		n = (int)(yGrid / delGrid); 
		for( i=n; i>=0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( 0, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
	}
	if( boundType == kOppositeSidesOpen ){
		upperZ = infile.ReadItem( upperZ, "UPPER_BOUND_Z" );
		n = (int)(xGrid / delGrid); 
		tempnode.setBoundaryFlag( kOpenBoundary );
		for( i=1, miNextNodeID=0; i<n; i++, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, 0, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		tempnode.setBoundaryFlag( kClosedBoundary );
		n = (int)(yGrid / delGrid);
		for( i=0; i<=n; i++, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( xGrid, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		tempnode.setBoundaryFlag( kOpenBoundary );
		n = (int)(xGrid / delGrid);
		for( i=n-1; i>0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, yGrid, upperZ );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBoundFront( tempnode );
			bndList.insertAtBack( nodIter.FirstBoundaryP() );
		}
		tempnode.setBoundaryFlag( kClosedBoundary );
		n = (int)(yGrid / delGrid);
		for( i=n; i>=0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( 0, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
	}
	else if( boundType == kAllSidesOpen ){
		miNextNodeID = 0;
		n = (int)(xGrid / delGrid);
		tempnode.setBoundaryFlag( kOpenBoundary );
		for( i=0; i<n; i++, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, 0, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		n = (int)(yGrid / delGrid);
		for( i=0; i<n; i++, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( xGrid, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		n = (int)(xGrid / delGrid);
		for( i=n; i>0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, yGrid, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
		n = (int)(yGrid / delGrid);
		for( i=n; i>0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( 0, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
		}
	}
	else if( boundType == kSpecifyOutlet ){
		// Read outlet coordinates from input file
		xout = infile.ReadItem( xout, "OUTLET_X_COORD" );
		yout = infile.ReadItem( yout, "OUTLET_Y_COORD" );
		
		// Create nodes for bottom (Y=0) boundary and place them on list
		n = (int)(xGrid / delGrid);
		tempnode.setBoundaryFlag( kClosedBoundary );
		for( i=0, miNextNodeID=0; i<n; i++, miNextNodeID++ ){
			// Assign node coords to tempnode and add tempnode to list
			dist = i * delGrid + 0.01 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, 0, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
			
			// If user wants outlet on this side between this and the next pt,
			// create the outlet now
			if( yout == 0 && xout > dist && xout < dist + delGrid ){
				tempnode.set3DCoords( xout, yout, 0 );
				tempnode.setBoundaryFlag( kOpenBoundary );
				miNextNodeID++;
				tempnode.setID( miNextNodeID );
				nodeList.insertAtBoundFront( tempnode );
				bndList.insertAtBack( nodIter.FirstBoundaryP() );
				tempnode.setBoundaryFlag( kClosedBoundary );
			}
		}
		
		// Create nodes for right (X=xGrid) boundary and place them on list
		n = (int)(yGrid / delGrid);
		for( i=0; i<n; i++, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( xGrid, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
			if( xout == xGrid && yout > dist && yout < dist + delGrid ){
				tempnode.set3DCoords( xout, yout, 0 );
				tempnode.setBoundaryFlag( kOpenBoundary );
				miNextNodeID++;
				tempnode.setID( miNextNodeID );
				nodeList.insertAtBoundFront( tempnode );
				bndList.insertAtBack( nodIter.FirstBoundaryP() );
				tempnode.setBoundaryFlag( kClosedBoundary );
			}
		}
		
		// Create nodes for top (Y=yGrid) boundary and place them on list
		n = (int)(xGrid / delGrid);
		for( i=n; i>0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( dist, yGrid, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
			if( yout == yGrid && xout < dist && xout > dist - delGrid ){
				tempnode.set3DCoords( xout, yout, 0 );
				tempnode.setBoundaryFlag( kOpenBoundary );
				miNextNodeID++;
				tempnode.setID( miNextNodeID );
				nodeList.insertAtBoundFront( tempnode );
				bndList.insertAtBack( nodIter.FirstBoundaryP() );
				tempnode.setBoundaryFlag( kClosedBoundary );
			}
		}
		
		// Create nodes for left (X=0) boundary and place them on list
		n = (int)(yGrid / delGrid);
		for( i=n; i>0; i--, miNextNodeID++ ){
			dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
			tempnode.set3DCoords( 0, dist, 0 );
			tempnode.setID( miNextNodeID );
			nodeList.insertAtBack( tempnode );
			bndList.insertAtBack( nodIter.LastP() );
			if( xout == 0 && yout < dist && yout > dist - delGrid ){
				tempnode.set3DCoords( xout, yout, 0 );
				tempnode.setBoundaryFlag( kOpenBoundary );
				miNextNodeID++;
				tempnode.setID( miNextNodeID );
				nodeList.insertAtBoundFront( tempnode );
				bndList.insertAtBack( nodIter.FirstBoundaryP() );
				tempnode.setBoundaryFlag( kClosedBoundary );
			}
		}
	}
	bndList.makeCircular();
	
	// Add edges
	Cout<<"\nAdding edges..."<<endl;
	
	tPtrListIter< tSubNode > bndIter( bndList );
	for( node0 = bndIter.FirstP(); !( bndIter.AtEnd() ); node0 = bndIter.NextP() ){
		node1 = bndIter.ReportNextP();
		node2 = bndIter.ReportPrevP();
		AddEdge( node0, node1, node2 );
	}
	nnodes = nodeList.getSize();
	nedges = edgeList.getSize();
	ntri = 0;
	
	int meshok = RepairMesh( bndList );
	assert( meshok );
	
	// Add nodes
	Cout<<"\nAdding nodes..."<<endl;
	
	tempnode.setBoundaryFlag( kNonBoundary );
	if( ptPlace == kUniformMesh || ptPlace == kPerturbedMesh ){
		nx = (int)(xGrid / delGrid);  // no. points in x direction
		ny = (int)(yGrid / delGrid);  // no. points in y direction
		for( i=1; i<nx; i++ ){
			for( j=1; j<ny; j++, miNextNodeID++ ){
				xyz[0] = i * delGrid - 0.25 * delGrid * (j%2)
                + 0.25 * delGrid * ((j+1)%2);
				xyz[1] = j * delGrid;
				if( ptPlace == kPerturbedMesh ){
					xyz[0] += 0.5 * delGrid * ( ran3( &seed ) - 0.5 );
					xyz[1] += 0.5 * delGrid * ( ran3( &seed ) - 0.5 );
				}
				xyz[2] = mElev + mElev * ( ran3( &seed ) - 0.5 );
				if( boundType == kOppositeSidesOpen ){
					slope = upperZ / yGrid;
					xyz[2] += slope * xyz[1] - mElev;
				}
				tempnode.set3DCoords( xyz[0], xyz[1], xyz[2] );
				tempnode.setID( miNextNodeID );
				AddNode( tempnode );
			}
		}
	}
	else if( ptPlace == kRandomMesh ){
		for( i=0; i<numPts; i++ ){
			// Randomize x,y, and z coordinates
			xyz[0] = ran3(&seed) * xGrid;
			xyz[1] = ran3(&seed) * yGrid;
			xyz[2] = mElev + mElev * ( ran3( &seed ) - 0.5 );
			if( xyz[0] != 0 && xyz[0] != xGrid && xyz[1] != 0 && xyz[1] != yGrid )
			{
				tempnode.set3DCoords( xyz[0], xyz[1], xyz[2] );
				tempnode.setID( miNextNodeID );
				AddNode( tempnode );
				miNextNodeID++;
			}
		}
	}
	
	if( boundType==kSpecifyOutlet && xout!=0 && yout!=0 ){
		tempnode.setBoundaryFlag( kOpenBoundary );
		tempnode.set3DCoords( xout, yout, 0 );
		tempnode.setID( miNextNodeID );
		miNextNodeID++;
		AddNode( tempnode );
	}
	
	MakeCCWEdges();
	
	Cout<<"\nTesting Mesh..."<<endl;
	UpdateMesh(); 
	CheckMeshConsistency();  
	
}

//=========================================================================
//
//
//                  Section 8: tMesh:: Make_PointsFromArcGrid( )
//
//
//=========================================================================

/**************************************************************************
**
**   tMesh::MakeRandomPointsFromArcGrid
**
**   Routine to make irregular mesh from regular Arc grid by randomly 
**   sampling the grid such that the average resolution of the irregular
**   mesh is equal to that of the Arc grid.
**
**   Calls: tInputFile::ReadItem, MakeCCWEdges(),
**          UpdateMesh(), CheckMeshConsistency()
**   Parameters: infile -- main parameter input file
**
**************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakeRandomPointsFromArcGrid( tInputFile &infile )
{
	int i, j;                        // loop counter
	int numrows, numcols, numpts;    // no. of rows, cols, points in mesh
	int delgrid,                     // arc grid cell size (side, meters) 
		nodata;                      // value signifying and no. w/ "no data"
	int numActNbrs;
	double x, y;                     // x, y coordinates
	tArray<int> bnd;                 // array of boundary codes 
	char arcgridFilenm[80];          // name of file containing arc grid data
	ifstream gridfile;               // the file (stream) itself
	double minx = 1e12, miny = 1e12, // minimum x and y coords
		maxx = 0, maxy=0,            // maximum x and y coords 
		minz= 1e12,	            // min. and max. elevs.
		dx, dy,                      // max width and height of region (meters)
		di, dj,                      // width and height of region in pixels
		xgen, ygen,                  // randomly generated x and y val's
		zinterp;                     // interp'd elev.
	double mindist;
	double delx, dely, dist;
	void srand48();
	void drand48();
	
	tSubNode tempnode( infile ),     // temporary node used in creating new pts
		*stp1, *stp2, *stp3;         // supertriangle vertices
	
	char dumhead[3];
	//tSubNode *cn, *minzPtr;
    auto *cn = new tSubNode();  // Initialize with a dynamically allocated object
    auto *minzPtr = new tSubNode();
    tEdge* ce;
	tMeshListIter< tSubNode > nI( nodeList );
	tPtrList<tSubNode> supertriptlist, deletelist;
	tPtrListIter<tSubNode> stpIter( supertriptlist ), dI( deletelist );
	tPtrListIter< tEdge > sI;
	
	infile.ReadItem( arcgridFilenm, "ARCGRIDFILENAME" );
	gridfile.open( arcgridFilenm );
	if( !gridfile.good() ){
		cerr << "Arc grid file name: '" << arcgridFilenm << "'\n";
		cerr << "I can't find a file by this name.";
	}
	
	// Read Header 
	gridfile >> dumhead >> numcols >> dumhead >> numrows >> dumhead 
		>> minx >> dumhead >> miny >> dumhead >> delgrid 
		>> dumhead >> nodata;
	
	maxx = minx + ( numcols - 1 ) * delgrid;
	maxy = miny + ( numrows - 1 ) * delgrid;
	dx = maxx - minx;
	dy = maxy - miny;
	di = numcols;
	dj = numrows;
	
	// Create matrix from input file using tMatrix
	
	tMatrix< double > elev( numrows, numcols );
	for( j=0; j<numrows; j++ ){
		for( i=0; i<numcols; i++ ){
			if( gridfile.eof() )
				cerr<<"Reached end-of-file while reading points.";
			gridfile >> elev( j, i );
		}
	}
	gridfile.close();
	
	// Create the 3 nodes that form the supertriangle and place them on the
	// node list in counter-clockwise order. (Note that the base and height
	// of the supertriangle are 5 times the width and height, respectively,
	// of the rectangle that encloses the points.) 
	
	tempnode.set3DCoords( minx-2*dx, miny-2*dy, 0.0 );
	tempnode.setBoundaryFlag( kClosedBoundary );
	tempnode.setID( -3 );
	nodeList.insertAtBack( tempnode );
	tempnode.set3DCoords( maxx+2*dx, miny-2*dy, 0.0 );
	tempnode.setID( -2 );
	nodeList.insertAtBack( tempnode );
	tempnode.set3DCoords( minx+0.5*dx, maxy+2*dy, 0.0 );
	tempnode.setID( -1 );
	nodeList.insertAtBack( tempnode );
	
	// set # of nodes, edges, and triangles
	nnodes = 3;
	nedges = ntri = 0;
	
	// Create the edges that connect the supertriangle vertices and place
	// them on the edge list.
	
	stp1 = nI.FirstP();
	stp2 = nI.NextP();
	stp3 = nI.NextP();
	AddEdge( stp1, stp2, stp3 );  // edges 1->2 and 2->1
	AddEdge( stp2, stp3, stp1 );  // edges 2->3 and 3->2
	AddEdge( stp3, stp1, stp2 );  // edges 3->1 and 1->3
	
	// Set up the triangle itself and place it on the list. To do this, we
	// just set up a list of pointers to the three nodes in the super tri
	// and pass the list (along with an iterator) to MakeTriangle.
	
	supertriptlist.insertAtBack( stp1 );
	supertriptlist.insertAtBack( stp2 );
	supertriptlist.insertAtBack( stp3 );
	supertriptlist.makeCircular();
	MakeTriangle( supertriptlist, stpIter );
	
	Cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: ";
	Cout << nedges << " NT: " << ntri << endl;
	Cout << "Begin Grid Interpolation\n";
	
	// Read and initialize seed for random number generation
	
	seed = infile.ReadItem( seed, "SEED" );
	srand( seed );
	numpts = numcols * numrows;
	tempnode.setBoundaryFlag( kNonBoundary );
	
	mindist = delgrid / 10.0;
	for( i=0; i<numpts; ++i ){
		xgen = std::rand() * (di - 1.0);
		ygen = std::rand() * (dj - 1.0);
		
		zinterp = InterpSquareGrid( xgen, ygen, elev, nodata );
		
		// Reset to "real" coords and add node
		
		if( zinterp != nodata )
			tempnode.setBoundaryFlag( kNonBoundary );
		else tempnode.setBoundaryFlag( kClosedBoundary );
		x = xgen * delgrid + minx;
		y = ygen * delgrid + miny;
		tempnode.set3DCoords( x, y, zinterp );
		
		// Check whether node is closer than delgrid/10 to another node
		
		cn = nI.FirstP();
		dist = delgrid;
		while( dist > mindist && !( nI.AtEnd() ) ){
			delx = cn->getX() - tempnode.getX();
			dely = cn->getY() - tempnode.getY();
			dist = sqrt( delx * delx + dely * dely );
			cn = nI.NextP();
		}
		if( dist > mindist ){
			cn = AddNode( tempnode, /*updatemesh =*/ 0 );
			if( zinterp != nodata && zinterp < minz ){
				minz = zinterp;
				minzPtr = cn;
			}
		}
		else --i; 
	}
	
	Cout << "Finished interpolation";
	Cout << "\n2 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: ";
	Cout << nedges << " NT: " << ntri << endl;
	
	// Make lowest node outlet (open boundary)
	
	minzPtr->setBoundaryFlag( kOpenBoundary );
	nI.Get( minzPtr->getID() );
	nodeList.moveToBoundFront( nI.NodePtr() );  
	
	Cout << "created open boundary outlet: " << nI.FirstBoundaryP()->getID() << "\n";
	Cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: ";
	Cout << nedges << " NT: " << ntri << endl;
	
	// Remove closed boundary nodes that do not have non-boundary nbrs
	
	for( cn = nI.FirstBoundaryP(); !( nI.AtEnd() ); cn = nI.NextP() ){
		numActNbrs = 0;
		sI.Reset( cn->getSpokeListNC() );
		for( ce = sI.FirstP(); !( sI.AtEnd() ); ce = sI.NextP() )
			if( ce->getDestinationPtr()->getBoundaryFlag() == kNonBoundary )
				++numActNbrs;
		if( numActNbrs ) cn->setZ( minz );
		else deletelist.insertAtBack( cn ); 
	}
	for( cn = dI.FirstP(); !( dI.AtEnd() ); cn = dI.FirstP() ){
		DeleteNode( cn, kNoRepair );
		deletelist.removeFromFront( cn );
	}
	
	Cout << "deleted superfluous boundary nodes\n";
	Cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: ";
	Cout << nedges << " NT: " << ntri << endl;
	
	UpdateMesh();
}

/**************************************************************************
**
**   tMesh::MakeHexMeshFromArcGrid
**
**   Routine to make irregular mesh from regular Arc grid by randomly 
**   sampling the grid such that the average resolution of the irregular
**   mesh is equal to that of the Arc grid. Designed to read from a grid 
**   containing points either within an already isolated basin or
**   containing "no data".
**
**   Calls: tInputFile::ReadItem, MakeCCWEdges(),
**          UpdateMesh(), CheckMeshConsistency()
**   Parameters: infile -- main parameter input file
**
**************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakeHexMeshFromArcGrid( tInputFile &infile )
{
	int keepgoing;
	int i, j;                        // loop counter
	int numrows, numcols;            // no. of rows, cols, points in mesh
	int delgrid,                     // arc grid cell size (side, meters) 
		nodata;                      // value signifying and no. w/ "no data"
	int numActNbrs;
	double x, y;                     // x, y coordinates
	tArray<int> bnd;                 // array of boundary codes 
	char arcgridFilenm[80];          // name of file containing arc grid data
	ifstream gridfile;               // the file (stream) itself
	double minx = 1e12, miny = 1e12, // minimum x and y coords
		maxx = 0, maxy=0,            // maximum x and y coords 
		minz= 1e12,                  // min. and max. elevs
		dx, dy,                      // max width and height of region (meters)
		di, dj,                      // width and height of region in pixels
		xgen, ygen,                  // randomly generated x and y val's
		zinterp;                     // interp'd elev.
	tSubNode tempnode( infile ),     // temporary node used in creating new pts
		*stp1, *stp2, *stp3;         // supertriangle vertices
	
	char dumhead[3];
	tSubNode *cn, *minzPtr;
	tEdge* ce;
	tMeshListIter< tSubNode > nI( nodeList );
	tPtrList<tSubNode> supertriptlist, deletelist;
	tPtrListIter<tSubNode> stpIter( supertriptlist ), dI( deletelist );
	tPtrListIter< tEdge > sI;
	
	infile.ReadItem( arcgridFilenm, "ARCGRIDFILENAME" );
	gridfile.open( arcgridFilenm );
	if( !gridfile.good() ){
		cerr << "Arc grid file name: '" << arcgridFilenm << "'\n";
		cerr << "I can't find a file by this name.";
	}
	
	// Read Header 
	gridfile >> dumhead >> numcols >> dumhead >> numrows >> dumhead 
		>> minx >> dumhead >> miny >> dumhead >> delgrid 
		>> dumhead >> nodata;
	
	maxx = minx + ( numcols - 1 ) * delgrid;
	maxy = miny + ( numrows - 1 ) * delgrid;
	dx = maxx - minx;
	dy = maxy - miny;
	di = numcols;
	dj = numrows;
	
	// Create matrix from input file
	
	tMatrix< double > elev( numrows, numcols );
	for( j=0; j<numrows; j++ ){
		for( i=0; i<numcols; i++ ){
			if( gridfile.eof() )
				cerr<<"Reached end-of-file while reading points.";
			gridfile >> elev( j, i );
		}
	}
	gridfile.close();
	
	// Create the 3 nodes that form the supertriangle and place them on the
	// node list in counter-clockwise order.
	
	tempnode.set3DCoords( minx-2*dx, miny-2*dy, 0.0 );
	tempnode.setBoundaryFlag( kClosedBoundary );
	tempnode.setID( -3 );
	nodeList.insertAtBack( tempnode );
	tempnode.set3DCoords( maxx+2*dx, miny-2*dy, 0.0 );
	tempnode.setID( -2 );
	nodeList.insertAtBack( tempnode );
	tempnode.set3DCoords( minx+0.5*dx, maxy+2*dy, 0.0 );
	tempnode.setID( -1 );
	nodeList.insertAtBack( tempnode );
	
	// Set # of nodes, edges, and triangles
	nnodes = 3;
	nedges = ntri = 0;
	
	// Create the edges 
	
	stp1 = nI.FirstP();
	stp2 = nI.NextP();
	stp3 = nI.NextP();
	AddEdge( stp1, stp2, stp3 );  // edges 1->2 and 2->1
	AddEdge( stp2, stp3, stp1 );  // edges 2->3 and 3->2
	AddEdge( stp3, stp1, stp2 );  // edges 3->1 and 1->3
	
	// Set up the triangle itself and place it on the list. 
	
	supertriptlist.insertAtBack( stp1 );
	supertriptlist.insertAtBack( stp2 );
	supertriptlist.insertAtBack( stp3 );
	supertriptlist.makeCircular();
	MakeTriangle( supertriptlist, stpIter );
	Cout << "formed supertriangle\n";
	Cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: ";  
	Cout << nedges << " NT: " << ntri << endl;
	
	// Place points in hexagonal mesh, i.e., equilateral triangles
	
	xgen = ygen = 0.0;
	j = 0;
	miNextNodeID = 0;
	keepgoing = 1;
	while( keepgoing ){
		zinterp = InterpSquareGrid( xgen, ygen, elev, nodata );
		
		// Reset to "real" coords and add node:
		if( zinterp != nodata )
			tempnode.setBoundaryFlag( kNonBoundary );
		else tempnode.setBoundaryFlag( kClosedBoundary );
		x = xgen * delgrid + minx;
		y = ygen * delgrid + miny;
		tempnode.setID( miNextNodeID );
		miNextNodeID++;
		tempnode.set3DCoords( x, y, zinterp );
		cn = AddNode( tempnode, /*updatemesh =*/ 0 );
		if( zinterp != nodata && zinterp < minz ){
			minz = zinterp;
			minzPtr = cn;
		}
		if( ygen <= dj - 1.732 ){
			// along x-dir:
			if( ( j%2 == 0 && xgen < di - 1.0 ) || ( j%2 == 1 && xgen < di - 2.0 ) )
				++xgen;
			// at end of row, start odd row:
			else {
				++j;
				xgen = 0.5 * (j%2);
				ygen += 0.866; // sqrt(3)/2
			}
		}
		else keepgoing = 0;
	}
	
	Cout << "finished interpolation:";
	Cout << "\n2 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: ";
	Cout << nedges << " NT: " << ntri << endl;
	
	// Make lowest node outlet (open boundary)
	
	minzPtr->setBoundaryFlag( kOpenBoundary );
	nI.Get( minzPtr->getID() );
	nodeList.moveToBoundFront( nI.NodePtr() ); 
	
	Cout << "created open boundary outlet: " << nI.FirstBoundaryP()->getID() << "\n";
	Cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: ";
	Cout  << nedges << " NT: " << ntri << endl;
	
	// Remove closed boundary nodes that do not have non-boundary nbrs
	
	for( cn = nI.FirstBoundaryP(); !( nI.AtEnd() ); cn = nI.NextP() ){
		numActNbrs = 0;
		sI.Reset( cn->getSpokeListNC() );
		for( ce = sI.FirstP(); !( sI.AtEnd() ); ce = sI.NextP() )
			if( ce->getDestinationPtr()->getBoundaryFlag() == kNonBoundary )
				++numActNbrs;
		if( numActNbrs ) cn->setZ( minz );
		else deletelist.insertAtBack( cn ); 
	}
	
	for( cn = dI.FirstP(); !( dI.AtEnd() ); cn = dI.FirstP() ){
		DeleteNode( cn, kNoRepair );
		deletelist.removeFromFront( cn );
	}
	
	Cout << "deleted superfluous boundary nodes\n";
	Cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: ";
	Cout << nedges << " NT: " << ntri << endl;
	
	UpdateMesh();
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
void tMesh<tSubNode>::writeRestart(fstream & rStr)
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
void tMesh<tSubNode>::readRestart(fstream & rStr)
{
  tMeshListIter< tSubNode > nodIter( nodeList );
  tSubNode *cn;
  for (cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP())
    cn->readRestart(rStr);
}


//=========================================================================
//
//
//                           End of tMesh.cpp
//
//
//=========================================================================
