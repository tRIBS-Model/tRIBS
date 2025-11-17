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
**  tResample.cpp:   Functions for tResample class (see tResample.h)
**
***************************************************************************/

#include "src/tRasTin/tResample.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//                  Section 1: tResample Constructors and Destructors
//
//
//=========================================================================

tResample::tResample(SimulationControl *simCtrPtr, tMesh<tCNode> * gridPtr)
{ 
	assert(gridPtr != nullptr);
	mew = gridPtr;
	simCtrl = simCtrPtr;
	
	// Total number of active Voronoi cells
	NVor = mew->getNodeList()->getActiveSize();
	
	// Contains the # of vertices for each Voronoi cell
	nPoints = new int [NVor];
	assert(nPoints != 0);
	for (int i=0; i < NVor; i++)
		nPoints[i] = -999;
	
	vXs = new double* [NVor];
	assert(vXs != 0);
	vYs = new double* [NVor];
	assert(vYs != 0);
	
	varFromGrid = new double [NVor];
	assert(varFromGrid != 0);
	
	varFromPoint = new int [NVor];
	assert(varFromPoint != 0);
	
	// Creates Dummy Cell
	eta = new vCell();
	
	Cout<<"\ntResample: Making boundary polygons..."<<endl<<endl;
	MakeBoundaryPolygons();
	
	Cout<<"\ntResample: Making interior polygons..."<<endl;
	MakeInteriorPolygons();

	
	NR=MR=0;
	gridIn = NULL;
	coorXG = NULL;
	coorYG = NULL;
	
	dR = dummy = 0.0;
	xllcR=yllcR=xulcR=yulcR=0.0;

}

tResample::~tResample()
{
	for (int i=0; i<NVor; i++) {
		delete [] vXs[i];       
		delete [] vYs[i];
	}
	
	delete [] vXs;
	delete [] vYs;
	
	delete [] nPoints;
	delete [] varFromGrid;
	delete [] varFromPoint;
	mew = NULL;
	if ( eta )
		delete eta;
	
	Cout<<"tResample Object has been destroyed..."<<endl<<flush;
	return;
}

void tResample::DestrtResample()
{
	for (int i=0; i<NR; i++)
		delete [] gridIn[i];
	delete [] gridIn;
	delete [] coorXG;
	delete [] coorYG;
	return;
}

//=========================================================================
//
//
//                  Section 2: tResample Functions
//
//
//=========================================================================

void tResample::allocMemory(double **xy, int n, int m)
{
	assert(n > 0 && m > 0);
	xy = new double* [n];
	assert(xy != 0);
	for (int i=0; i < n; i++) {
		xy[i] = new double[m];
		assert(xy[i] != 0);
	}
	return;
}

void tResample::printToFile(ofstream &Otp0, double **qq, int n, int m)
{
	for (int i=0; i < n; i++) {
		for (int j=0; j < m; j++)
			Otp0 <<qq[i][j]<<" ";
		Otp0 << endl;
	}
	return;
}

/***************************************************************************
**
**  tResample::MakeBoundaryPolygons()
**
**  Part of the calculated Voronoi area may lie outside the mesh domain. 
**  Sometimes the area can even blow up.
** 
**  This function gets rid of the outlying voronoi vertices by considering 
**  boundary edges and the problematic interior nodes associated with them.
**  
**  The algorithm is the following:
**      - take a boundary node, find two consecutive boundary edges for it; 
**      - make a list of non-boundary nodes which are connected to the  
**        current boundary node;
**      - test each of the nodes as a potentially problematic associating
**        it with one of the consecutive boundary edges (associating is 
**        done by constructing an equilateral triangle which has a 
**        boundary edge as its base side with the opposite node located
**        _inside_ of the domain; if a test interior node falls inside of 
**        this triangle, it is associated with the boundary edge); 
**      - if there is an outlying vertex (vertices), call the function  
**        'FixVoronoiPolygon()' which resolves the problem by cutting 
**        the outlying part and saving the info 'tCNode' VertsX & VertsY 
**        arrays
**
***************************************************************************/
void tResample::MakeBoundaryPolygons()
{ 
	int cnt=0;
	int i;
	int flag;
	double d, d1, d2, x1, y1, x2, y2, m1, m2, m11, m22;
	
	tCNode *cn;
	tCNode  *tan;
	tEdge  *curedg;   //pointer to current edge
	tEdge  *ce;
	tEdge  *be1, *be2;
	tMeshListIter<tCNode> niter ( mew->getNodeList() );
	
	tArray<double> xy(2), xy1(2), xy2(2);
	tArray<double> xyn1(2), xyn2(2), xynn(2), xynnn(2);
	
	tPtrList< tCNode >     NodesLst;
	tPtrListIter< tCNode > NodesIter( NodesLst );
	
	tPtrList< tEdge >     vedgList;        
	tPtrListIter< tEdge > vtxIter( vedgList ); 
	
	Cout.setf( ios::fixed, ios::floatfield);
	
	i=0;
	for (cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP()) {
		
		cnt = 0;
		
		// Check boundary status: search only for BND
		if ( cn->getBoundaryFlag() == kClosedBoundary || 
			 cn->getBoundaryFlag() == kOpenBoundary ) {
			
			NodesLst.Flush();
			
			ce = cn->getEdg();
			
			// Loop till you get an edge pointing to an interior node
			cnt = 0;
			while ((ce->getDestinationPtrNC()->getBoundaryFlag() == kClosedBoundary ||
					ce->getDestinationPtrNC()->getBoundaryFlag() == kOpenBoundary ) &&
				   (cnt < 50) ) {
				ce = ce->getCCWEdg();
				cnt++;
			}
			
			if (cnt >= 50 ) {
				if (simCtrl->Verbose_label == 'Y') {
					cout<<"\nIn tResample: Boundary Node ID = "<<cn->getID()
					<<" has only boundary spokes... going to the next one"<<endl;
				}
			}
			else {
				curedg = ce;
				
				// Now 'ce' is NOT a boundary edge, find a sequence of two CCW boundary edges
				
				flag = 1;
				do {
					assert( ce != nullptr);
					ce = ce->getCCWEdg();
					if ( ce->getDestinationPtrNC()->getBoundaryFlag() == kClosedBoundary ||
						 ce->getDestinationPtrNC()->getBoundaryFlag() == kOpenBoundary ) {
						if (flag) {
							be1 = ce;
							flag = 0;
						}
						else {
							be2 = ce; //be2 is CCW Boundary of be1 which is boundary 
									  // whichever first edge is found, it is the right one
						}
					}
					else { // Makes a stack of nodes connected to the current BND node
						assert((ce->getDestinationPtrNC()) != nullptr);
						NodesLst.insertAtBack( (tCNode*)(ce->getDestinationPtrNC()) );
					}
				} while ( ce != curedg );
				NodesLst.makeCircular(); // Now the LAST one points to the FIRST
				
				xy2 = be2->getDestinationPtrNC()->get2DCoords();
				xy1 = be1->getDestinationPtrNC()->get2DCoords();
				xy  = be1->getOriginPtrNC()->get2DCoords();
				x1  = xy1[0]-xy[0]; //These are coordinates of vectors be1 & be2
				y1  = xy1[1]-xy[1]; //required for any further calculation
				x2  = xy2[0]-xy[0]; //(x1,y1) & (x2,y2)
				y2  = xy2[1]-xy[1]; 
				d1 = be1->getLength();
				d2 = be2->getLength();
				
				// Compute the angle between be1 & be2
				// alpha = acos((x1*x2 + y1*y2)/(d1*d2))*180/(4*atan(1)); 
				// Angle in degrees
				// cout<<"ANGLE = "<<alpha<<endl;   
				
				// Now let's build a normal to the boundary 'be1' either in the
				// current point OR in the point dividing be1 in two halfs
				m1  = (xy1[0]+xy[0])/2.0;  //X Coordinate of the center of be1 
				m2  = (xy1[1]+xy[1])/2.0;  //Y Coordinate of the center of be1 
				
				m11  = (xy2[0]+xy[0])/2.0;  //X Coordinate of the center of be2
				m22  = (xy2[1]+xy[1])/2.0;  //Y Coordinate of the center of be2 
				
				x1  = xy1[0]-m1;   
				y1  = xy1[1]-m2;    
				
				x2  = xy[0]-m11;    
				y2  = xy[1]-m22;   
				
				d = 10.0*d1; // Arbitrary value: the length of normal
				xyn1 = FindNormal(d, x1, y1, m1, m2);
				
				d = 10.0*d2; // Arbitrary value: the length of normal
				xyn2 = FindNormal(d, x2, y2, m11, m22);
				
				// Let's make a triangle with Origin & Destination points of 'be1' 
				// at its base and the node of built perpendicular at the top.
				// The CCW sequence of points is: xy, xyn1, xy1. 
				// Loop node list to define if any of the nodes falls in triangle
				
				for ( tan=NodesIter.FirstP(); !(NodesIter.AtEnd()); tan=NodesIter.NextP() ) {
					xynn = tan->get2DCoords();
					
					// This node is associated with 'be1'
					if ( IsInTriangle(xy, xyn1, xy1, xynn[0], xynn[1]) ) { 
						
						// We need to check the edges FIRST, it's been shown that can form  
						// loops of Voronoi, which creates problems for running the model    
						
						if ((tan->bndEdge1 != be1) && (tan->bndEdge2 != be1)
							&& (tan->bndEdge1 != be1->FindComplement())
							&& (tan->bndEdge2 != be1->FindComplement()) ) { 
							
							FixVoronoiPolygon(cn, tan, xy, xy1);
							if (tan->bndEdge1 != nullptr) {
								if (!(tan->bndEdge2 != nullptr))  //Not assigned yet
									tan->bndEdge2 = be1;
							}  
							else
								tan->bndEdge1 = be1;
						}
					}
					
					// Let's check if the node falls in the second triangle formed for 'be2'. 
					// Reason is that for edges 'be1' and 'be2' having angle << 180o we may 
					// have Voronoi vertices outside for both edges for a single interior node
					
					// This node is associated with 'be2'
					if ( IsInTriangle(xy2, xyn2, xy, xynn[0], xynn[1]) ) { 
						
						if ((tan->bndEdge1 != be2) && (tan->bndEdge2 != be2)
							&& (tan->bndEdge1 != be2->FindComplement())
							&& (tan->bndEdge2 != be2->FindComplement()) ) { 
							
							FixVoronoiPolygon(cn, tan, xy2, xy);
							if (tan->bndEdge1 != nullptr) {
								if (!(tan->bndEdge2 != nullptr))
									tan->bndEdge2 = be2;
							}
							else
								tan->bndEdge1 = be2;
						}
					}
					vedgList.Flush();
				}
			} // ELSE statement: BND node has NO spokes to NON-BND nodes
		}
		i++;
	}
	Cout<<"\nTotal Number of Nodes: \t\t"<<mew->getNodeList()->getSize()<<endl<<flush;
	return;
}

/***************************************************************************
**
**  tResample::FindNormal()
**
**  The function returns coordinates 'xyn' of the normal length 'd'
**  with respect to vector having coordinates 'x1' & 'y1'. 'm1' and
**  'm2' are the coordinates of the origin point of the perpendicular
** 
***************************************************************************/
tArray< double > tResample::FindNormal(double d, double x1, double y1, 
									   double m1, double m2) 
{
	double xn, yn, tmp;
	tArray<double> xyn(2);
	
	tmp = x1*x1+y1*y1;
	if ( fabs(tmp) <= 1.0E-9) {
		cout<<"\ntResample::FindNormal found Zero Vector!..."<<endl;
		cout<<"Change the TIN Mesh"<<endl;
		cout<<"Exiting Program..."<<endl;
		exit(2);
	}
	
	yn = fabs(x1*d/sqrt( tmp ));
	if ( fabs(y1)<= 1.0E-9 )
		xn = 0;
	else 
		xn = fabs(sqrt( d*d - yn*yn ));
	
	if (x1 >= 0 && y1 >= 0) {     //Perpendicular is in IVth quarter
		if (yn)
			yn = -yn;
	}
	else if (x1 > 0 && y1 < 0) {  //Perpendicular is in IIIrd quarter
		xn = -xn;
		yn = -yn;
	}
	else if (x1 <= 0 && y1 <= 0) { //Perpendicular is in IInd quarter
		if (xn)
			xn = -xn;
	}
	
	//else if (x1 < 0 && y1 > 0)    // <=== Perpendicular is in Ist quarter
	
	xn += m1;
	yn += m2;
	xyn[0] = xn;
	xyn[1] = yn;  
	
	return xyn;
}

/***************************************************************************
**
**  tResample::IsInTriangle()
**
**  Defines if a point (x,y) falls in the triangle. The algorithm exploits
**  the fact that the 3 triangle points are always in counter-clockwise
**  order, so that the point is contained within a given triangle (p0,p1,p2)
**  if and only if the point lies to the left of vectors p0->p1, p1->p2,
**  and p2->p0. Here's how it works:
** 
***************************************************************************/
int tResample::IsInTriangle(tArray< double > &xyp1, tArray< double > &xyp2,
                            tArray< double > &xyp3, double x, double y)
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
		
		if ( c > 0.0 )  //Not to the LEFT
			k = -1;
		else { 
			if ( c == 0.0 )  //on the BND
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
**  tResample::FixVoronoiPolygon()
**
** The function takes as arguments pointer to a BND node 'cn' and 
** ptr to an interior node 'tan' which is connected to 'cn' through 
** an edge.
** NOTE: 'tan' may NOT necessarily be connected to 'cn' but you 
** will need to slightly modify the code:
**    1) Make sure that the BND edge starting at 'cn' and going to
**       some other bnd node is the APPROPRIATE edge for analysis. 
**       It is done through checking ALL nodes CONNECTED to 'cn' 
**       by constructing a triangle (see code before the function 
**       call).  If 'cn' and 'tan' are NOT connected through an edge,
**       the trick is how to find 'tan' that causes the problem.
**    2) Modify the line 'curedg=cn->EdgToNod(tan)->FindComplement()'
**       to something like 'curedg = tan->getFlowEdg()'
**       Everything else should still work fine.
** 
** I did not make the code more modular to account for these cases 
** because there is no any _fast_ way of finding node 'tan' when 
** it is not connected to 'cn'. Time is a factor.
** 
** o The node 'tan' is to be checked if some of its Voronoi vertices 
** are outside of the boundary edge defined by Origin point 'xy' 
** (first BND node) and Destination point 'xy1' (the other BND node). 
** 
** o The node 'tan' might already have had the correction procedure, 
** so if it did, the updated arrays of vertices are used which are 
** extracted from the 'tCNode' object data members.
** 
** NOTE: the code is also not modular in this sense but it can be  
** easily modified with the inclusion of templates. 
** 
** o Once the node 'tan' is checked and (possibly) updated, write 
** the arrays back to 'tCNode'. These data will be used later in the 
** tResample code
** 
** o Need to recheck the 'VoronoiArea' data-member in 'tNode'
**
** The code heavily relies on the global function 'PointsCCW()', so 
** if it iss not robust I do not guarantee that the algorithm 
** will work...
**
***************************************************************************/
void tResample::FixVoronoiPolygon(tCNode *cn, tCNode *tan, 
                                  tArray< double > &xy, 
				  tArray< double > &xy1) {

	int cnt=0;
	int ll;
	int flag;
	tEdge  *curedg;
	tEdge  *ce, *nne;
	
	tArray<double> tvtx(2), tt(2), tt1(2), tt2(2), dumm(2);
	
	tList< tArray< double > > vcL;             // list of vertex coordinates
	tListIter< tArray< double > > vcI( vcL );  // iterator for coord list
	
	tList< tArray< double > > NewL;            // list of vertex coordinates
	tListIter< tArray< double > > NewI( NewL );	// iterator for coord list
	
	// tan->allocVertArrays( NewL.getSize() ); // SKY2008Snow, AJR2008 -- Trial, REVISIT

	// Check to see if some of the vertices have been already
	// modified, compose a list of vertices  depending on this
	
	if (tan->nVerts <= 0) { //Vertices have not been modified
		
		// We need to start from the edge that has the right Voronoi    
		// vertex INSIDE of the domain.  We start from edge connecting nodes 
		// 'cn' and 'tan', then search for the appropriate edge. The edge
		// direction should be 'tan->xxx'. The reason for doing this is that  
		// two/or more Voronoi v-s can be outside relative the current edge
		// NOTE: Could we use nPoints[i] instead ? 
		//       The problem is how to relate i & tan ?
		curedg = cn->EdgToNod( tan )->FindComplement();
		assert(curedg != nullptr);
		ce = curedg;
		tvtx = curedg->getRVtx();
		cnt = 0;
		while ( (!PointsCCW( tvtx, xy1, xy )) && (cnt<50) ) {
			ce   = ce->getCCWEdg();
			tvtx = ce->getRVtx();
			cnt++;
		}
		assert(ce != nullptr);
		/*  Scheme of nodes/edges:         'tan'
			o
			'tvtx' O        /             
			/ 'ce' 
			L 
			xy1 o<-------------- o xy
			'be1'
			*/
		
		if (cnt >= 50) { 
			if (simCtrl->Verbose_label == 'Y') {
				cout<<"FixVoronoiPolygon: ERROR!!! The node does not have any Voronoi "
				<<"vertices inside of domain! ID = "<<tan->getID()<<endl<<flush;
				cout.setf( ios::fixed, ios::floatfield);
				cout<<"Coordinates:  ("<<tan->getX()<<","<<tan->getY()<<")"<<endl<<flush;
				cout<<"\tZ = "<<tan->getZ()<<endl<<flush;
				cout<<"\tBND Code: "<<tan->getBoundaryFlag()<<endl<<flush;
				cout<<"LEFT UNCHANGED!"<<endl<<flush;
			}
			return;
		}
		
		// Now, start composing a list of vertices. We know that the
		// first vertex is always INSIDE and this will help us later
		vcL.insertAtBack( tvtx ); //this node is INSIDE
		
		nne = ce->getCCWEdg();
		do {
			tt = tvtx;
			tvtx = nne->getRVtx();
			if (fabs(tvtx[0] - tt[0]) <= 1.0E-3 && 
				fabs(tvtx[1]- tt[1])  <= 1.0E-3) { // Don't duplicate!
				nne = nne->getCCWEdg();
			}
			else {
				vcL.insertAtBack( tvtx );
				nne = nne->getCCWEdg();
			}
		} while ( nne != ce );
	}
	
	// Some of the vertices have been modified before
	else { 
		// We need to check if the 1st point is outside
		flag = 0; 
		for (ll=0; ll<tan->nVerts; ll++) {
			tvtx[0] = tan->VertsX[ll];
			tvtx[1] = tan->VertsY[ll];
			
			// If first point(s) happens to be 
			// outside, then we'll need to fix that.
			if (ll == 0) {
				if ( !PointsCCW( tvtx, xy1, xy ) )
					flag = 1; 
			}
			vcL.insertAtBack( tvtx );
		}
		
		if (flag) {
			tvtx = *(vcI.FirstP());
			while ( (!PointsCCW( tvtx, xy1, xy )) && (cnt<50) ) {
				vcL.insertAtBack( tvtx );
				vcL.removeFromFront(tvtx);
				tvtx = *(vcI.FirstP());
			}
		}
	}
	vcL.makeCircular(); //Last points to first
	
	// Now let's do the actual work: see if any of the vertices 
	// is outside of the boundary. If so - adjust introducing new vertices
	for ( tvtx=*(vcI.FirstP()); !(vcI.AtEnd()); tvtx=*(vcI.NextP()) ) {
		flag = 0; //Indicator of how many vertices outside: 1 or 2
		
		// Check if the vertex is OUTSIDE
		if ( !PointsCCW( tvtx, xy1, xy ) ) {
			
			// As long as we always start from a point inside
			// we ensure that we always have *Previous* element
			tt = *(vcI.PrevP()); 
			vcI.Next(); // Go back to the current element again
			flag = 1;
			
			if (!(vcI.AtEnd())) {
				tt2 = *(vcI.NextP());
				dumm = tt2;
				while ( !PointsCCW( dumm, xy1, xy ) && !(vcI.AtEnd()) ) {
					flag++;
					tt2 = dumm; // always the last one which is OUTSIDE/INSIDE
					dumm = *(vcI.NextP()); // at the end, the one INSIDE
				}
			}
			// Take 1st in the list - INSIDE
			else
				tt2 = vcL.getFirst()->getDataRef();
			
			// 1ST INTERSECTION 
			tt1 = FindIntersectionCoords(tt, tvtx, xy1, xy);
			NewL.insertAtBack( tt1 );
			
			// 2ND INTERSECTION 
			if (flag == 1) { // 1 point is OUTSIDE
				tt1 = FindIntersectionCoords(tvtx, tt2, xy1, xy);
				NewL.insertAtBack( tt1 );
			}
			else if  (flag > 1) { //( >= 2 ) points are OUTSIDE
				tt1 = FindIntersectionCoords(tt2, dumm, xy1, xy);
				NewL.insertAtBack( tt1 );
			}
			
			if ( !(vcI.AtEnd()) )
				NewL.insertAtBack( dumm );
			
			if ( vcI.AtEnd() ) // Need to step back to finish the loop
				vcI.PrevFull();
		}
		else 
			NewL.insertAtBack( tvtx );
	}
	NewL.makeCircular(); //Last points to first
	
	
	// Record the changes so that next time we know 
	// which polygon vertices have been modified 
	
	if (tan->nVerts <= 0)  // Vertices have NOT been modified 
		tan->allocVertArrays( NewL.getSize() );
	
	else if (tan->nVerts == NewL.getSize())
		; // the array is of the same size 
	
	else {
		tan->deleteVertArrays();
		tan->allocVertArrays( NewL.getSize() );
	}
	ll = 0;
	for ( tvtx=*(NewI.FirstP()); !(NewI.AtEnd()); tvtx=*(NewI.NextP()) ) {
		tan->VertsX[ll] = tvtx[0];
		tan->VertsY[ll] = tvtx[1];
		ll++;
	}
	
	// Free memory
	vcL.Flush();
	NewL.Flush();
	
	return;
}

/***************************************************************************
**
**  tResample::MakeInteriorPolygons()
**
**  This function is called after the 'MakeBoundaryPolygons()' function.
**  It fills the arrays nPoints, vXs, vYs by reading the coordinates
**  of the voronoi vertices according to the right-hand-rule location
**  with respect to a current edge of the node.
**  It accounts for the adjustments made by the f-n 'MakeBoundaryPolygons()'
**
***************************************************************************/
void tResample::MakeInteriorPolygons() 
{ 
	int cnt=0;
	int i, iv, ll;
	double xCentroid, yCentroid, areaT;
	tCNode *cn;
	tEdge  *firstedg; //ptr to first edge
	tEdge  *curedg;   //pointer to current edge
	tMeshListIter<tCNode> niter ( mew->getNodeList() );
	tArray<double> xy(2);
	
	// Memory allocation and filling array 'nPoints'
	i=0;
	for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
		
		if (cn->nVerts <= 0) {   //Vertices have been modified
			cnt = 0;
			firstedg = cn->getFlowEdg();
			if (!firstedg)
				firstedg = cn->getEdg();
			assert(firstedg != nullptr);
			curedg = firstedg->getCCWEdg();
			cnt++;
			while (curedg != firstedg) {
				curedg = curedg->getCCWEdg();
				cnt++;
			}
		}
		else if (cn->nVerts > 0) //Vertices have not been modified
			cnt = cn->nVerts;
		
		nPoints[i] = cnt;
		vXs[i] = new double [cnt];
		assert(vXs[i] != 0);
		vYs[i] = new double [cnt];
		assert(vYs[i] != 0);
		
		if (cnt < 3 || cnt > 50) {    
			cout<<"\nError: tResample-> # of vertices < 3 OR > 50!!!"<<endl;
			cout<<"Node ID = "<<cn->getID()<<endl<<flush;
			cout<<"Coordinates:  ("<<cn->getX()<<","<<cn->getY()<<")"<<endl<<flush;
			cout<<"\tZ = "<<cn->getZ()<<endl<<flush;
			cout<<"\tBND Code: "<<cn->getBoundaryFlag()<<endl<<flush;
			cout<<"\tnVerts: "<<cn->nVerts<<endl<<flush;
			cout<<"Exiting Program..."<<endl;
			exit(2);
		}
		i++;
	}
	
	if (simCtrl->Verbose_label == 'Y') {
		cout<<"\tConstructor: NVor = "<<NVor<<"; I = "<<i
		<<"; <- Must be the same!"<<endl<<flush;
	}
	
	// Filling the ragged arrays 'vXs' and 'vYs'
	i = 0;
	for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
		
		if (cn->nVerts <= 0) { 
			iv = 0;
			firstedg = cn->getFlowEdg(); 
			
			xy = firstedg->getRVtx();
			vXs[i][iv] = xy[0];
			vYs[i][iv] = xy[1];
			iv++;
			curedg = firstedg->getCCWEdg();
			while (curedg != firstedg) {
				xy = curedg->getRVtx();
				if (xy[0] == vXs[i][iv-1] && xy[1] == vYs[i][iv-1]) {
					iv--;         //If points coincide
					nPoints[i]--; //just skip it...
				} 
				else { 
					vXs[i][iv] = xy[0];
					vYs[i][iv] = xy[1];
				}
				iv++;
				curedg = curedg->getCCWEdg();
			}
			
			if (iv < nPoints[i] || iv > nPoints[i]) { 
				cout<<"\nError: Constructor iv != nPoints[i]: iv = "<<iv
				<<"; nPoints[i] = "<<nPoints[i]<<endl<<flush;
				cout<<"Exiting Program..."<<endl;
				exit(2);
			}
		}
		
		else if (cn->nVerts > 0) {
			for (ll=0; ll<cn->nVerts; ll++) {
				vXs[i][ll] = cn->VertsX[ll];
				vYs[i][ll] = cn->VertsY[ll];
			}
			
			// We need to update Voronoi Area because it was
			// computed incorrectly in 'meshElements.cpp'
			eta-> initializeVCell(vXs[i], vYs[i], nPoints[i]);
			
			cnt = eta->polyCentroid(eta->VoronX, eta->VoronY, nPoints[i]+1, 
									&xCentroid, &yCentroid, &areaT);
			
			if (cnt > 0) { 
				cout<<"\nError: MakeInteriorPolygons in function 'polyCentroid'..."<<endl;
				cout<<"Node ID = "<<cn->getID()<<endl<<flush;
				cout<<"NODES ARE: vx = [";
				for (ll=0; ll<cn->nVerts; ll++) 
					cout<<cn->VertsX[ll]<<" ";
				cout<<"];\nNODES ARE: vy = [";
				for (ll=0; ll<cn->nVerts; ll++)
					cout<<cn->VertsY[ll]<<" ";
				cout<<"];"<<endl<<flush;
				cout<<"Exiting Program..."<<endl;
				exit(2);
			}
			cn->setVArea(areaT);
			cn->setVArea_Rcp(1.0/areaT);
			
			// Destroy temporary arrays in 'vCell'
			eta->DestrtvCell();
		}
		i++;
	}
	return;
}

/***************************************************************************
**
**  tResample::VerticesNoAccBndEff()
**
**  The algorithm below does not account for the possible effects at the
**  boundaries of the domain of interest, as opposed to the function 
**  'MakeBoundaryPolygons()' which creates corrected voronoi vertices.
**  The function fills the arrays nPoints, vXs, vYs by reading the 
**  coordinates of the voronoi vertices according to the right-hand-rule
**  location with respect to a current edge of the node.
**
***************************************************************************/
void tResample::VerticesNoAccBndEff() 
{ 
	int cnt=0;
	int i, iv;
	
	tCNode *cn;
	tEdge  *firstedg;
	tEdge  *curedg; 
	tMeshListIter<tCNode> niter ( mew->getNodeList() );
	tArray<double> xy(2);
	
	// Memory allocation and filling array 'nPoints'
	i=0;
	for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
		cnt = 0;
		firstedg = cn->getFlowEdg(); 
		curedg = firstedg->getCCWEdg();
		cnt++;
		while (curedg != firstedg) {
			curedg = curedg->getCCWEdg();
			cnt++;
		}
		
		nPoints[i] = cnt;
		vXs[i] = new double [cnt];
		assert(vXs[i] != 0);
		vYs[i] = new double [cnt];
		assert(vYs[i] != 0);
		
		if (cnt < 3 || cnt > 50) {
			cout<<"\nError: Number of voronoi vertices < 3 OR > 50..."<<endl;
			cout<<"Exiting Program..."<<endl;
			exit(2);
		}
		i++;
	}
	
	if (simCtrl->Verbose_label == 'Y') {
		cout<<"\tCONSTRUCTOR: NVor = "<<NVor<<"; I = "<<i
		<<"; <-- Must be the SAME!"<<endl<<flush;
	}
	
	// Filling the ragged arrays 'vXs' and 'vYs'
	i = 0;
	for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
		iv = 0;
		firstedg = cn->getFlowEdg(); 
		xy = firstedg->getRVtx();
		vXs[i][iv] = xy[0];
		vYs[i][iv] = xy[1];
		iv++;
		curedg = firstedg->getCCWEdg();
		
		while (curedg != firstedg) {
			xy = curedg->getRVtx();
			if (xy[0] == vXs[i][iv-1] && xy[1] == vYs[i][iv-1]) {
				iv--;         //If points coincide, skip it
				nPoints[i]--; 
			} 
			else { 
				vXs[i][iv] = xy[0];
				vYs[i][iv] = xy[1];
			}
			iv++;
			curedg = curedg->getCCWEdg();
		}
		if (iv < nPoints[i] || iv > nPoints[i]) { 
			cout<<"\nError: Constructor iv != nPoints[i]: iv = "<<iv
			<<"; nPoints[i] = "<<nPoints[i]<<endl<<flush;
			cout<<"Exiting Program..."<<endl;
			exit(2);
		}
		i++;
	}
	return;
}

/***************************************************************************
**
**  tResample::doIT(char *GridIn, int flag)
**
**  The function resamples ASCII grid in ArcInfo/View format located
**  in a file specified by name *GridIn with two possible options (flag):
**   '1'- to use weighted average
**   '2'- to use a discrete grid value 
**  The algorithm below HEAVILY relies on the assumption that node IDs
**  correspond to the order of nodes in the list created by 
**  tMeshListIter< tCNode > niter ( gridPtr->getNodeList() );
**  i.e. the first node must have ID = 0, second must have ID = 1, etc.
**  The algorithm has been added the capability of correction procedure
**  if the node does not lie in the proper grid area: a local stack of 
**  nodes is composed for which later  'FindNeighbValue( tCNode * )'
**  routine is used to solve the problem.
**
***************************************************************************/
double* tResample::doIt(char *GridIn, int flag) 
{
	tCNode *cn;
	tMeshListIter< tCNode > niter ( mew->getNodeList() );
	tPtrList< tCNode >     NodesLst;
	tPtrListIter< tCNode > NodesIter( NodesLst );
	NodesLst.Flush();

    if (simCtrl->Verbose_label == 'Y')
	    Cout<<"\nResampling "<<GridIn<<endl<<flush;
	
	// Reads Input Grid
	readInputGrid(GridIn);
	
	int i=0;
	for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
		
		// Initialize current Cell
		eta->initializeVCell(simCtrl, this, vXs[i], vYs[i], nPoints[i]);

		// Assign YmaxInd, YminInd, etc. to tCNode to increase efficiency
		// NOTE: NOT done yet
		
		// Get an appropriate grid value
		varFromGrid[i] = dummy;
		varFromGrid[i] = eta->convertToVoronoiFormat(flag);
		
    // Set resample index for nodes
    cn->setResIndex(i);
		
		// Destroy temporary arrays in 'vCell'
		eta->DestrtvCell();
		
		// Check if there is a problem with the value
		if (varFromGrid[i] == dummy)
			NodesLst.insertAtBack( (tCNode*)cn );
		
		i++;
	}
	
	// Run the correction procedure if needed
	for ( cn=NodesIter.FirstP(); !(NodesIter.AtEnd()); cn=NodesIter.NextP() )
		FindNeighbValue( cn );
	
	// Memory Management
	DestrtResample();
	NodesLst.Flush();
	
	return varFromGrid;
}

/***************************************************************************
**
**  tResample::FindNeighbValue(tCNode *cn )
**
**  The function finds a value for the current node *cn going first to 
**  its downstream neighbor, if it does not have the right value too,
**  it checks all local adjacent nodes... if this doesn't work, it goes
**  and checks adjacent nodes for the downstream node ... so on until 
**  it finds appropriate value.
**
***************************************************************************/
void tResample::FindNeighbValue(tCNode *cn ) 
{
	int ID;
	double value;
	tCNode *nbn;
	tEdge  *firstedg;
	tEdge  *curedg;
	
	firstedg = cn->getFlowEdg(); 
	assert(firstedg != nullptr);
	curedg = firstedg->getCCWEdg();
	
   ID = firstedg->getDestinationPtrNC()->getResIndex();
	value = varFromGrid[ID];
	while (value == dummy) {
		// Move to downstream nodes
		nbn = (tCNode *)firstedg->getDestinationPtrNC();
		while ((value == dummy) && (curedg != firstedg)) {
			if ( curedg->getDestinationPtrNC()->getBoundaryFlag() != kClosedBoundary &&
				 curedg->getDestinationPtrNC()->getBoundaryFlag() != kOpenBoundary ) { 

            ID = curedg->getDestinationPtrNC()->getResIndex();
				value = varFromGrid[ID];
			}
			curedg = curedg->getCCWEdg();
		}
		// The termination above was not due to this equality 
		if (value == dummy) { 
			firstedg = nbn->getFlowEdg();  
			curedg = firstedg->getCCWEdg();

         ID = firstedg->getDestinationPtrNC()->getResIndex();
			value = varFromGrid[ID];
		}
	}

        varFromGrid[cn->getResIndex()] = value;
	return;
}

/***************************************************************************
**
**  tResample::doIT(int *statID, double *XX, double *YY, int NN)
**
**  The function resamples values of point data (say gauge data) 
**  into the nodes of the mesh.  Station ID is specified in 'statID' array,
**  coordinates of the stations are in the arrays 'XX' and 'YY', the 
**  number of stations is given by NN.
**  Option '1'- weighted average is used in the algorithm
**
***************************************************************************/
int* tResample::doIt(int *statID, double *XX, double *YY, int NN) 
{
	tMeshListIter<tCNode> niter ( mew->getNodeList() );
	
	int i=0;
	for (tCNode *cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
		
		// Initialize current Cell
		eta->initializeVCell(vXs[i], vYs[i], nPoints[i]);
		
		// Get an appropriate grid value
		varFromPoint[i] = eta->convertPointData(statID, XX, YY, NN);

    // Set resample index
    cn->setResIndex(i);
		
		// Destroy temporary arrays in 'vCell'
		eta->DestrtvCell();
		
		i++;
	}
	return varFromPoint;
}

/***************************************************************************
**
**  tResample::readInputGrid(char *GridIn)
**
**  The function reads a grid in the standard ArcInfo/ArcView ASCII format
**  'GridIn' is the full pathname for the input grid file which is to be
**  specified by the calling function
**
***************************************************************************/
void tResample::readInputGrid(char *GridIn) 
{
	int i,j;
	double tempo, mott;
	char lineIn[300];
	char tmp[20];
	
	ifstream Inp0(GridIn);
	if (!Inp0) {
		cout <<"File "<<GridIn<<" not found!!!"<<endl;
		cout<<", tResample.cpp Error Location"<<endl;
		exit(2);
	}
	
	Inp0.getline(lineIn, 256);
	sscanf(lineIn, "%10s %d", tmp, &MR);
	
	Inp0.getline(lineIn, 256);
	sscanf(lineIn, "%10s %d", tmp, &NR);
	
	Inp0 >> tmp >> xllcR;
	Inp0 >> tmp >> yllcR;
	Inp0 >> tmp >> dR;
	Inp0 >> tmp >> dummy;
	
	coorYG = new double [NR+1];  //NR+1
	assert(coorYG != 0);
	coorXG = new double [MR+1];  //MR+1
	assert(coorXG != 0);
	
	for (i=0; i < NR+1; i++)
		coorYG[i] = 0.;
	for (j=0; j < MR+1; j++)
		coorXG[j] = 0.;
	
	gridIn = new double* [NR]; 
	assert(gridIn != 0);
	
	for (i=0; i < NR; i++)  {
		gridIn[i] = new double[MR];
		assert(gridIn[i] != 0);
		for (j=0; j < MR; j++)    {
			Inp0 >> tempo;
			gridIn[i][j] = tempo;
		}
	}
	
	// Get upper left corner coordinates
	xulcR = xllcR;
	yulcR = yllcR + NR*dR;
	
	tempo = yulcR;
	for (i=0; i < NR+1; i++)  {
		coorYG[i] = tempo;
		tempo -= dR;
	}
	mott  = xulcR;
	for (j=0; j < MR+1; j++)  {
		coorXG[j] = mott;
		mott += dR;
	}
	Inp0.close();
	return;
}

/***************************************************************************
**
**  tResample::In_Mrain_Name( )
**
***************************************************************************/
void tResample::In_Mrain_Name(char *pref, char *ext, 
                              int hour, int day, int month, int year) 
{
	snprintf(fileIn, sizeof(fileIn), "%s%02d%02d%02d%04d.%s", pref, month, day, year, hour, ext);
	Cout<<"File IN  --> "<<fileIn<<endl;
	return;
}

/***************************************************************************
**
**  tResample::Out_Mrain_Name( )
**
***************************************************************************/
void tResample::Out_Mrain_Name(char *prefix, char *mrf, int hour) 
{   
	int min=0;
	snprintf(fileOut, sizeof(fileOut), "%s.%04d_%02d.%s", prefix, hour, min, mrf);
	Cout<<"File OUT --> "<<fileOut<<endl;
	return;
}

/***************************************************************************
**
**  PrintEdgeInfo
**
**  An auxiliary function that prints out tEdge data 
** 
***************************************************************************/
void tResample::PrintEdgeInfo(tEdge *be)
{ 
	cout<<"-----------------------------------------------------"<<endl;
	cout<<"ORIGIN Node ID = "<<be->getOriginPtrNC()->getID()<<"  (" 
		<<be->getOriginPtrNC()->getX()<<","
		<<be->getOriginPtrNC()->getY()<<")"<<endl<<flush;
	cout<<"DESTIN Node ID = "<<be->getDestinationPtrNC()->getID()<<"  (" 
		<<be->getDestinationPtrNC()->getX()<<","
		<<be->getDestinationPtrNC()->getY()<<")"<<endl<<flush;
	cout<<endl;
	cout<<"-----------------------------------------------------"<<endl;
	return;
}


//=========================================================================
//
//
//                  Section 3: vCell Constructor and Destructor
//
//
//=========================================================================

vCell::vCell() {
	nv = -999;  
	inX    = NULL;
	inY    = NULL;
	VoronX = NULL;
	VoronY = NULL;
	base   = NULL; 
	
	xMax = -999.0;
	yMax = -999.0;
	xMin = 1.0E+37;
	yMin = 1.0E+37;
}

vCell::vCell(SimulationControl *simCtrPtr, tResample *tres, double *x, 
			 double *y, int n) 
{  
	simCtrl = simCtrPtr;
	initializeVCell(simCtrl, tres, x, y, n);
	
}

vCell::~vCell() 
{  
	base = NULL;
	DestrtvCell();
	return;
}

void vCell::DestrtvCell() 
{
	delete [] inX;
	delete [] inY;
	delete [] VoronX;
	delete [] VoronY;
	inX    = NULL;
	inY    = NULL;
	VoronX = NULL;
	VoronY = NULL;
	return;
}

//=========================================================================
//
//
//                  Section 4: vCell Functions
//
//
//=========================================================================

/***************************************************************************
**
**  vCell::initializeVCell(double *x, double *y, int n)
**
** Arguments: - x and y coordinates of the current voronoi cell vertices
**            - number of vertices
** Objective: To initialize 'vCell' object for ''Point'' resampling 
** Return value: void
** Algorithm: - assigns # of vertices
**            - copies their coordinates
**            - gets closed poygon
**
***************************************************************************/
void vCell::initializeVCell(double *x, double *y, int n) 
{
	nv = n+1;
	
	if ( (!(inX==NULL)) || (!(inY==NULL)) ||
		 (!(VoronX==NULL)) || (!(VoronX==NULL)) )
		DestrtvCell();
	
	if (n > 50) {
		cout<<"\nError: The number of vertices is larger than allowed..."<<endl;
		cout<<"Exiting Program..."<<endl;
		exit(2);
	}
	
	// Arrays inX & inY contain X-Y coordinates of an input cell
	inX = new double [n];
	assert(inX != 0);
	inY = new double [n];
	assert(inY != 0);
	
	VoronX = new double [nv];
	assert(VoronX != 0);
	VoronY = new double [nv];
	assert(VoronY != 0);
	
	XminInd = -1;
	YminInd = -1;
	XmaxInd = -1;
	YmaxInd = -1;
	
	getClosedPolygon(x, y, nv); //So that x,y(n) = x,y(0)
	
	return;
}

/***************************************************************************
**
**  vCell::initializeVCell( )
**
**  Function:  initializeVCell(double *x, double *y, int n)
**  Arguments: - x and y coordinates of the current voronoi cell vertices
**             - number of vertices
**  Objective: To initialize 'vCell' object for ''Grid'' resampling 
**  Return value: void
**  Algorithm: - assigns # of vertices
**           - copies their coordinates
**           - gets closed poygon
**	     - gets Xmax-min & Ymax-min indices of the bounding grid box
**
***************************************************************************/
void vCell::initializeVCell(SimulationControl *simCtrPtr, tResample *tres, 
							double *x, double *y, int n) 
{
	simCtrl = simCtrPtr;
	nv = n+1;     //Number of vertices of the CLOSED polygon
	
	base = tres;  //Pointer to tResample
	
	if ( (!(inX==NULL)) || (!(inY==NULL)) || 
		 (!(VoronX==NULL)) || (!(VoronY==NULL)) )
		DestrtvCell();
	
	if (n > 50) {
		cout<<"\nError: Number of vertices is larger than allowed..."<<endl<<flush;
		cout<<"Exiting Program..."<<endl;
		exit(2);
	}
	
	// Arrays inX & inY contain X-Y coordinates of the input cell
	inX = new double [n];
	assert(inX != 0);
	inY = new double [n];
	assert(inY != 0);
	
	VoronX = new double [nv];
	assert(VoronX != 0);
	VoronY = new double [nv];
	assert(VoronY != 0);
	
	xMax = -999.0;
	yMax = -999.0;
	xMin = 1.0E+37;
	yMin = 1.0E+37;
	
	getClosedPolygon(x, y, nv);
	
	findMaxMin();
	
	// Defining bounding box which contains current Voronoi cell.
	// Before calling these functions we need to define arrays of 
	// coordinates coorXG[] & coorYG[] and their dimesions MR & NR
	
	YmaxInd = getBoxYs(yMin);     // min & max are swapped here
	if (yMin == base->coorYG[YmaxInd-1]) //the point is ON the prec. row
		YmaxInd--;
	
	YminInd = getBoxYs(yMax);
	if (yMax != base->coorYG[YminInd]) //the point is NOT ON the curr. row
		YminInd--;
	
	XminInd = getBoxXs(xMin);
	if (xMin != base->coorXG[XminInd]) //the point is NOT ON the curr. row
		XminInd--;
	
	XmaxInd = getBoxXs(xMax);
	if (xMax == base->coorXG[XmaxInd-1]) //the point is on the preced. row
		XmaxInd--;  
	
	return;
}

/***************************************************************************
**
**  vCell::getClosedPolygon( )
**
**  Copies values from arrays of X-s and Y-s of vertices for non-closed
**  polygon contour to arrays which will contain closed contour coordinates
**  i.e. the first node in the array is the last one too.
**
***************************************************************************/
void vCell::getClosedPolygon(double *X, double *Y, int inN) 
{
	for (int i=0; i < inN; i++) { //To get a CLOSED polygon
		if (i == inN-1) {
			VoronX[i] = X[0];
			VoronY[i] = Y[0];
		}
		else {
			inX[i] = X[i];
			inY[i] = Y[i];
			VoronX[i] = X[i];
			VoronY[i] = Y[i];
		}
	}
	return;
}

/***************************************************************************
**
**  vCell::findMaxMin( )
**
**  Finds Max and Min in an array of Voronoi polygon vertices
**
***************************************************************************/
void vCell::findMaxMin() 
{
	for (int j=0; j < nv; j++) { //Loop in a single Voronoi cell
		if (VoronX[j] >= xMax)
			xMax = VoronX[j];
		else if (VoronX[j] <= xMin)
			xMin = VoronX[j];
		if (VoronY[j] >= yMax)
			yMax = VoronY[j];
		else if (VoronY[j] <= yMin)
			yMin = VoronY[j];
	}
	return;
}

/***************************************************************************
**
**  vCell::getBoxYs()
**
**  Function:  getBoxYs(double y)
**  Arguments: value in the Y direction bounding voronoi cell (y_max or y_min)
**             - INPUT grid dimensions must be set before: NR, MR
**  	   - array of node coordinates *coorYG[NR+1]* must be defined 
**  Objective: to find index of *coorYG* array  corresponding to bounding ROW
**  Return value: (ki+1) <-- corresponds to LOWER bound row (lies below/w point y)
**  Algorithm: - split the original array in two subarrays
**             - look if the point lies in one of the subarrays
**             - take the appropriate one... continue splitting...
**             - 100 m is assumed to be allowable distance to be out of the box
**  
***************************************************************************/
int vCell::getBoxYs(double y) 
{
	int offs  = 0;        //offset
	int range = base->NR; //number of ROWS in the input grid & max index
	int kki = range-1;
	int ki  = range-1;
	
	if (y < (base->coorYG[range]) || y > (base->coorYG[0])) {
		cout<<"\nError: 'y' coordinate of one of the voronoi cell vertices";
		cout<<"\n\t is out of the domain of the input grid!!!"<<endl;
		cout<<"\t y = "<<y<<";  YMinGrid\t= "<<base->coorYG[range]<<endl;
		cout<<"\t y = "<<y<<";  YMaxGrid\t= "<<base->coorYG[0]<<endl;
		cout<<"\tThe vertex is now assumed to fall on the grid boundary"<<endl; 
		
		if (y > base->coorXG[0]) // Falls out on the NORTH side of the grid 
			ki = 0;
		else 
			ki = range-1;
	}
	else {
		while (y > base->coorYG[ki] && y > base->coorYG[ki+1]) {
			kki = ki;
			range = (int)(ceil((double)range/2.));
			ki = offs + range;
			if (ki >= kki)
				ki = kki - range;
			if (ki > base->NR) {
				cout<<"\n\tERROR! Function 'getBoxYs' Check the algorithm!\n"<<endl;
				cout<<"Exiting Program..."<<endl;
				exit(2);
			}
			
			if (y < base->coorYG[ki] && y < base->coorYG[ki+1]) {
				offs = ki;
				ki = kki;
			}
		}
	}
	return(ki+1); //point 'y' always lies either ON or AT the preceding row
				  // or between, i.e. it is guaranteed that (coorYG[ki+1] <= y)
}               

/***************************************************************************
**
**  vCell::getBoxXs( )
**
**  Function:  getBoxXs(double x)
**  Arguments: value in the X direction bounding voronoi cell (x_max or x_min)
**             - INPUT grid dimensions must be set before: NR, MR
**  	   - array of node coordinates *coorXG[MR+1]* must be defined 
**  Objective: to find index of *coorXG* array  corresponding to bounding COLUMN
**  Return value: (ki+1) <-- corresponds to upper bound column (lies below/w point x)
**  Algorithm: - split the original array in two subarrays
**             - look if the point lies in one of the subarrays
**             - take the appropriate one... continue splitting...
**             - 100 m is assumed to be allowable distance to be out of the box
**
***************************************************************************/
int vCell::getBoxXs(double x) 
{
	int offs  = 0;        //offset
	int range = base->MR; //number of ROWS in the input grid 
	int kkj = range-1;
	int kj  = range-1;
	
	if (x > (base->coorXG[range]) || (x < base->coorXG[0])) {
		cout<<"\nError: 'x' coordinate of one of the voronoi cell vertices";
		cout<<"\n\t is out of the domain of the input grid domain!"<<endl;
		cout<<"\t x = "<<x<<";  XMinGrid\t= "<<base->coorXG[0]<<endl;
		cout<<"\t x = "<<x<<";  XMaxGrid\t= "<<base->coorXG[range]<<endl;
		cout<<"\tThe vertex is now assumed to fall on the grid boundary"<<endl;
		
		if (x < base->coorXG[0]) // Falls out on the WEST side of the grid 
			kj = 0;
		else 
			kj = range-1;
	}
	
	else {
		while (x < base->coorXG[kj] && x < base->coorXG[kj+1]) {
			kkj = kj;
			range = (int)(ceil((double)range/2.));
			kj = offs + range;
			if (kj >= kkj)
				kj = kkj - range;
			if (kj > base->MR) {
				cout<<"\n\tERROR! Function 'getBoxXs' Check the algorithm!\n"<<endl;
				cout<<"Exiting Program..."<<endl;
				exit(2);
			}
			if (x > base->coorXG[kj] && x > base->coorXG[kj+1]) {
				offs = kj;
				kj = kkj;
			}
		}
	}
	return(kj+1); //point 'x' always lies either ON or AT the preceding 
				  //column, i.e. it is guaranteed that (coorXG[kj+1] >= x)
}              

/***************************************************************************
**
**  vCell::convertPointData()
**
**  Function:  convertPointData (int *statID, double *XX, double *YY, int NN)
**  Arguments: 
**     - # of stations (points)
**     - station ID array
**     - arrays of station coordinates XX and YY
**  Objective:    to define appropriate index 
**  Return value: value of variable extracted from the input POINTS
**  Algorithm: - loops through the points (stations)
**             - finds the closest station
**  	   - returns its index
**  
***************************************************************************/
int vCell::convertPointData(int *statID, double *XX, double *YY, int NN) 
{
	int    index;
	int    ki = 0;
	double dist, areaT;
	double mindist = 1.0E+9;
	double xCentroid, yCentroid;
	
	for (int i=0; i < NN; i++) {
		
		ki = polyCentroid(VoronX, VoronY, nv, &xCentroid, &yCentroid, &areaT);
		if (ki > 0) { 
			cout<<"\nERROR in function 'polyCentroid'!!!"<<endl;
			cout<<"Exiting Program..."<<endl;
			exit(2);
		}
		
		dist = findDistance(xCentroid, yCentroid, XX[i], YY[i]);
		
		if (dist < mindist) {
			mindist = dist;
			index = statID[i];
		}
	}
	return index;
}

/***************************************************************************
**
**  vCell::findDistance( )
**
**  Objective: Finds distance between the two points whose coordinates
**  are specified by (x1,y1) and (x2,y2)
**
***************************************************************************/
double vCell::findDistance(double x1, double y1, double x2, double y2) 
{
	return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
}

/***************************************************************************
**
**  vCell::convertToVoronoiFormat(int flag)
**
**  Function:  convertToVoronoiFormat(int flag)
**  Arguments: * flag of the option should be used: 
**       '1'-to use weighted average
**       '2'- to use a discrete grid value 
**     - input GRID dimensions must be set before: NR, MR
**     - GRID values must be read 
**     - arrays of node coordinates *coorXG[MR+1]* *coorYG[NR+1]* must be defined 
**     - bounding coordinates of a current voronoi cell: yMin, yMax, xMin, xMax 
**  Objective: to obtain appropriate value for Voronoi cell
**  Return value: value of variable extracted from the input GRID
**  
***************************************************************************/
double vCell::convertToVoronoiFormat(int flag)
{
	
	// GMnSKY2008MLE
	//int i, ii, L, M, l, m;
	//int k, j, kk, tempSum;
   	//int tempInd1, tempInd2;
   	//int tempInd3, tempInd4;
	int i_row, ii, L, M, l, m;
	int k, i_col, kk, tempSum;
   	int skip, ind_poly;

   	int goodness(0),cnt(0);
	int ki=0; 
	int cntr=0;
	int sumInOrOut(0);

	// GMnSKY2008MLE
	int Vert_sumInOrOut(0);

	int dirXint, dirYint;
	int InOrOut[4]; //1 is In, 0 is out; Prior to 3/23/2023 InOrOut was creating and address error, or should have as it was InOrOut[3] -WR
					//initialize InOrOut to all 0s
	for (int i = 0; i < 4; i++)
		InOrOut[i] = 0;   
	
	// GMnSKY2008MLE
	int i_all, i_deall; // Giuseppe
	int flag_int; // Giuseppe
	std::vector<int> Vert_InOrOut(nv - 1); //1 is In, 0 is out	Giuseppe
	
	double area1, area2, areaT, value = 0.;
	double xCentroid, yCentroid;
	double y1NormVert, x1NormVert, y2NormVert, x2NormVert;
	
	// GMnSKY2008MLE
	//double angleCorn, angleVert1, angleVert2;
	//double cornerX, cornerY, cornAngle;
	double cornerX, cornerY;	
	double angleVert1Vert2, angleCornerVert1, angleCornerVert2; // Giuseppe
	
	double cornerX1, cornerY1;
	double baseMag, sMag, heightMag;
	
	// GMnSKY2008MLE
	//double x2Unit, y2Unit;
	//double xProj, yProj;
	//double trapX1, trapX2, trapY1, trapY2;	//used for intercepts
	//double trapX3, trapX4, trapY3, trapY4;	//used for intercept
	double X_inters(-999.0), Y_inters(-999.0); // Giuseppe
	double cornerXUL, cornerYUL; // Giuseppe
	double coord_int; // Giuseppe

	double areaTriPoly(0.0), areaTriCorn(0.0);
	double fracOfGrid(0.0), weightedAveSum(0.0);
	double dirOfRot(0.0);

	// GMnSKY2008MLE
	double TOL(1e-5);
	
	xCentroid = yCentroid = 0.0;
	
	double **polyXY_1, **polyXY_2;     // these are NOT closed polygons
	double **pxy_1, **pxy_2;           // these are NOT closed polygons
	
	// GMnSKY2008MLE
	//polyXY_1 = new double* [2];
	//polyXY_2 = new double* [2];
	//for (i=0; i < 2; i++) {
	//	polyXY_1[i] = new double[nv+2];
	//	assert(polyXY_1[i] != 0);
	//	polyXY_2[i] = new double[nv+2];
	//	assert(polyXY_2[i] != 0);
	//}
	polyXY_1 = 0;
	polyXY_2 = 0;
	pxy_1 = 0;
	pxy_2 = 0;
	
	// So, now the limiting boundaries have been defined
	// we should look HOW the voronoi cell is bounded...
	
	// CENTER OF MASS SITUATION
	if (flag == 2) {
		
		// Finding centroid of the polygon: VoronX & VoronY are CLOSED 
		ki = polyCentroid(VoronX, VoronY, nv, &xCentroid, &yCentroid, &areaT);
		
		if (ki > 0) { 
			cout<<"\nError in function 'polyCentroid'!"<<endl;
			exit(2);
		}
		
		L = XmaxInd;
		M = YmaxInd;
		findCorrespInd(xCentroid, yCentroid, &L, &M);
				
		// If the center of mass is outside of the domain
		// use vertices to define if it crosses at all... 
		ii = -1;
		while (((base->gridIn[M][L]) == base->dummy) && (ii < (nv-1)) ) {
			L = XmaxInd;
			M = YmaxInd;
			ii++;   // loops through the vertices
			findCorrespInd(VoronX[ii], VoronY[ii], &L, &M);
		}
		
		
		if (base->gridIn[M][L] != base->dummy)
			value = base->gridIn[M][L];
		else {
			cout<<"\n\n\tWarning! tResample::convertToVoronoiFormat: CASE #00\n"
			<<"The current cell is not within the grid domain!"<<endl;
			printCellData();
			value = base->dummy;
			cout<<"...accepted dummy value "<<base->dummy<<endl;
		}
	}
	
	// WEIGHTED AVERAGE SITUATIONS 
	else if (flag == 1) {
		
		// The cell is completely within 1 grid cell
		if (abs(YmaxInd-YminInd)==1 && abs(XmaxInd-XminInd)==1) {
			if (base->gridIn[YminInd][XminInd] != base->dummy)
				value = base->gridIn[YminInd][XminInd];
			else {
				cout<<"\n\n\tWarning! tResample::convertToVoronoiFormat: CASE #1\n"
				<<"The current cell is NOT within the grid domain!"<<endl;
				printCellData();
				value = base->dummy; 
				cout<<"...accepted dummy value "<<base->dummy<<endl;
			}
		}
		
		// The cell is intersected OX, i.e. lies in 2 grid cells
		else if (abs(YmaxInd-YminInd)==2 && abs(XmaxInd-XminInd)==1) { 
			if ((base->gridIn[YminInd][XminInd]   != base->dummy) && 
				(base->gridIn[YminInd+1][XminInd] != base->dummy) ) {
			
				// GMnSKY2008MLE
				// allocate memory for polyXY_1 and polyXY_2
				polyXY_1 = new double* [2];
				assert(polyXY_1 != 0);
				polyXY_2 = new double* [2];
				assert(polyXY_2 != 0);
				for (i_all=0; i_all < 2; i_all++) {
					polyXY_1[i_all] = 0;
					polyXY_1[i_all] = new double[nv+2];
					assert(polyXY_1[i_all] != 0);
					polyXY_2[i_all] = 0;
					polyXY_2[i_all] = new double[nv+2];
					assert(polyXY_2[i_all] != 0);
				}
		
				// GMnSKY2008MLE
				// initialize polyXY_1 and polyXY_2
				for (i_all=0; i_all < 2; i_all++) {
					for (k=0; k<(nv+2) ; k++)
					{
						polyXY_1[i_all][k] = 0;
						polyXY_1[i_all][k] = 0;
						polyXY_2[i_all][k] = 0;
						polyXY_2[i_all][k] = 0;
					}
				}

				// Get polygons by intersection
				defineTwoSubArrays(polyXY_1, polyXY_2, &L, &M, VoronX, VoronY, 
								   nv, base->coorYG[YmaxInd-1], 1);
				area1 = polygonArea(polyXY_1, L);
				area2 = polygonArea(polyXY_2, M);
				areaT = area1 + area2;
				value = (base->gridIn[YminInd][XminInd]*area1 +
						 base->gridIn[YminInd+1][XminInd]*area2)/areaT;
			
				// GMnSKY2008MLE
				// deallocate memory
				for (i_deall=0; i_deall < 2; i_deall++) {
					delete [] polyXY_1[i_deall];
					polyXY_1[i_deall] = 0;
					delete [] polyXY_2[i_deall];
					polyXY_2[i_deall] = 0;
				}
				delete [] polyXY_1;
				polyXY_1 = 0;
				delete [] polyXY_2;
				polyXY_2 = 0;
			
			}
			// Current Voronoi cell is partly outside of the grid domain
			// A value of an adjacent non-void grid cell is assigned
			else if ((base->gridIn[YminInd][XminInd]   == base->dummy) &&
					 (base->gridIn[YminInd+1][XminInd] != base->dummy) ) {
				value = base->gridIn[YminInd+1][XminInd];
			}
			// Current Voronoi cell is partly outside of the grid domain
			// A value of an adjacent non-void grid cell is assigned
			else if ((base->gridIn[YminInd][XminInd]   != base->dummy) &&
					 (base->gridIn[YminInd+1][XminInd] == base->dummy) ) {
				value = base->gridIn[YminInd][XminInd];
			}
			// Current Voronoi cell is not within the grid domain
			// A dummy value is accepted (will be modified later)
			else {
				printCellData();
				value = base->dummy; 
			}
		}
		
		// The cell is intersected by OY, i.e. lies in 2 grid cells   
		else if (abs(YmaxInd-YminInd)==1 && abs(XmaxInd-XminInd)==2) {
			if ((base->gridIn[YminInd][XminInd]   != base->dummy) && 
				(base->gridIn[YminInd][XminInd+1] != base->dummy) ) {
			
				// GMnSKY2008MLE
				// allocate memory for polyXY_1 and polyXY_2
				polyXY_1 = new double* [2];
				assert(polyXY_1 != 0);
				polyXY_2 = new double* [2];
				assert(polyXY_2 != 0);
				for (i_all=0; i_all < 2; i_all++) {
					polyXY_1[i_all] = 0;
					polyXY_1[i_all] = new double[nv+2];
					assert(polyXY_1[i_all] != 0);
					polyXY_2[i_all] = 0;
					polyXY_2[i_all] = new double[nv+2];
					assert(polyXY_2[i_all] != 0);
				}

				// GMnSKY2008MLE
				// initialize polyXY_1 and polyXY_2
				for (i_all=0; i_all < 2; i_all++) {
					for (k=0; k<(nv+2) ; k++)
					{
						polyXY_1[i_all][k] = 0;
						polyXY_1[i_all][k] = 0;
						polyXY_2[i_all][k] = 0;
						polyXY_2[i_all][k] = 0;
					}
				}

				// Get polygons by intersection
				defineTwoSubArrays(polyXY_1, polyXY_2, &L, &M, VoronX, VoronY,
								   nv, base->coorXG[XmaxInd-1], 2);   
				area1 = polygonArea(polyXY_1, L);
				area2 = polygonArea(polyXY_2, M);
				areaT = area1 + area2;
				value = (base->gridIn[YminInd][XminInd]*area2 +
						 base->gridIn[YminInd][XminInd+1]*area1)/areaT;
			
				// GMnSKY2008MLE
				// deallocate memory
				for (i_deall=0; i_deall < 2; i_deall++) {
					delete [] polyXY_1[i_deall];
					polyXY_1[i_deall] = 0;
					delete [] polyXY_2[i_deall];
					polyXY_2[i_deall] = 0;
				}
				delete [] polyXY_1;
				polyXY_1 = 0;
				delete [] polyXY_2;
				polyXY_2 = 0;			
			
			}
			// Current Voronoi cell is partly outside of the grid domain
			// A value of an adjacent non-void grid cell is assigned
			else if ((base->gridIn[YminInd][XminInd]   == base->dummy) && 
					 (base->gridIn[YminInd][XminInd+1] != base->dummy) ) {
				value = base->gridIn[YminInd][XminInd+1];
			}
			// Current Voronoi cell is partly outside of the grid domain
			// A value of an adjacent non-void grid cell is assigned
			else if ((base->gridIn[YminInd][XminInd]   != base->dummy) && 
					 (base->gridIn[YminInd][XminInd+1] == base->dummy) ) {
				value = base->gridIn[YminInd][XminInd];
			}
			// Current Voronoi cell is not within the grid domain
			// A dummy value is accepted (will be modified later)
			else {
				printCellData();
				value = base->dummy; 
			}
		}
		
		// The node is intersected by OX and OY, i.e. lies in 4 grid cells
		else if (abs(YmaxInd-YminInd)==2 && abs(XmaxInd-XminInd)==2) {
		
			// GMnSKY2008MLE
			// allocate memory for polyXY_1 and polyXY_2
			polyXY_1 = new double* [2];
			assert(polyXY_1 != 0);
			polyXY_2 = new double* [2];
			assert(polyXY_2 != 0);
			for (i_all=0; i_all < 2; i_all++) {
				polyXY_1[i_all] = 0;
				polyXY_1[i_all] = new double[nv+2];
				assert(polyXY_1[i_all] != 0);
				polyXY_2[i_all] = 0;
				polyXY_2[i_all] = new double[nv+2];
				assert(polyXY_2[i_all] != 0);
			}
			
			// GMnSKY2008MLE
			// initialize polyXY_1 and polyXY_2
			for (i_all=0; i_all < 2; i_all++) {
				for (k=0; k<(nv+2) ; k++)
				{
					polyXY_1[i_all][k] = 0;
					polyXY_1[i_all][k] = 0;
					polyXY_2[i_all][k] = 0;
					polyXY_2[i_all][k] = 0;
				}
			}

			// Get polygons by intersection
			defineTwoSubArrays(polyXY_1, polyXY_2, &L, &M, VoronX, VoronY, 
							   nv, base->coorYG[YmaxInd-1], 1);
			
			// To get CLOSED polygon for further calculations
			
			// GMnSKY2008MLE
			//for (i=0; i < 2; i++)  {  
			//	polyXY_1[i][L] = polyXY_1[i][0]; 
			//	polyXY_2[i][M] = polyXY_2[i][0]; 
			//}
			for (i_all=0; i_all < 2; i_all++)  {  
				polyXY_1[i_all][L] = polyXY_1[i_all][0]; 
				polyXY_2[i_all][M] = polyXY_2[i_all][0]; 
			}
			
			L++;
			M++;
			// To get areaT 
			ki = polyCentroid(VoronX, VoronY, nv, &xCentroid, &yCentroid, &areaT); 
			if (ki > 0) { 
				cout<<"\nError in function 'polyCentroid'!!!"<<endl;
				exit(2);
			}
			if (areaT < 0.0)
				areaT = fabs(areaT);
			
			// Part 1.)
			pxy_1 = new double* [2];

			// GMnSKY2008MLE
			assert(pxy_1 != 0);

			pxy_2 = new double* [2];

			// GMnSKY2008MLE
			assert(pxy_2 != 0);
			for (i_all=0; i_all < 2; i_all++) {
				pxy_1[i_all] = 0;
				pxy_1[i_all] = new double[L+2];
				assert(pxy_1[i_all] != 0);
				pxy_2[i_all] = 0;
				pxy_2[i_all] = new double[L+2];
				assert(pxy_2[i_all] != 0);
			}
			// initialize pxy_1 and pxy_2
			for (i_all=0; i_all < 2; i_all++) {
				for (k=0; k<(L+2) ; k++)
				{
					pxy_1[i_all][k] = 0;
					pxy_1[i_all][k] = 0;
					pxy_2[i_all][k] = 0;
					pxy_2[i_all][k] = 0;
				}
			}
			//for (i=0; i < 2; i++) {
			//	pxy_1[i] = new double[L+2];
			//	assert(pxy_1[i] != 0);
			//	pxy_2[i] = new double[L+2];
			//	assert(pxy_2[i] != 0);
			//}
			
			l = m = 0;
			area1 = area2 = 0.0;
			defineTwoSubArrays(pxy_1, pxy_2, &l, &m, polyXY_1[0], polyXY_1[1],
							   L, base->coorXG[XmaxInd-1], 2);
			
			if (l > 0) { // If the polygon has an area at all
				area1 = polygonArea(pxy_1, l);
				// Current Voronoi cell is partly outside of the grid domain
				// Looking at other grid nodes sharing the cell
				if (base->gridIn[YminInd][XminInd+1] == base->dummy) { 
					areaT -= area1;
					cntr++;
				}
				else 
					value += base->gridIn[YminInd][XminInd+1]*area1;
			}
			if (m > 0) { //If polygon has area at all...
				area2 = polygonArea(pxy_2, m);
				// Current Voronoi cell is partly outside of the grid domain
				// Looking at other grid nodes sharing the cell
				if (base->gridIn[YminInd][XminInd] == base->dummy) {
					areaT -= area2;
					cntr++;
				}
				else 
					value += base->gridIn[YminInd][XminInd]*area2;
			}
			
			// GMnSKY2008MLE
			//for (i=0; i < 2; i++) {
			//	delete [] pxy_1[i];
			//	delete [] pxy_2[i];
			//}
			for (i_deall=0; i_deall < 2; i_deall++) {
				delete [] pxy_1[i_deall];
				pxy_1[i_deall] = 0;
				delete [] pxy_2[i_deall];
				pxy_2[i_deall] = 0;
			}

			delete [] pxy_1;

			// GMnSKY2008MLE
			pxy_1 = 0;

			delete [] pxy_2;

			// GMnSKY2008MLE
			pxy_2 = 0;
			
			pxy_1 = new double* [2];

			// GMnSKY2008MLE
			assert(pxy_1 != 0);

			pxy_2 = new double* [2];

			// GMnSKY2008MLE
			assert(pxy_2 != 0);
			for (i_all=0; i_all < 2; i_all++) {
				pxy_1[i_all] = 0;
				pxy_1[i_all] = new double[M+2];
				assert(pxy_1[i_all] != 0);
				pxy_2[i_all] = 0;
				pxy_2[i_all] = new double[M+2];
				assert(pxy_2[i_all] != 0);
			}
			// initialize pxy_1 and pxy_2
			for (i_all=0; i_all < 2; i_all++) {
				for (k=0; k<(M+2) ; k++)
				{
					pxy_1[i_all][k] = 0;
					pxy_1[i_all][k] = 0;
					pxy_2[i_all][k] = 0;
					pxy_2[i_all][k] = 0;
				}
			}
			//for (i=0; i < 2; i++)     {
			//	pxy_1[i] = new double[M+2];
			//	assert(pxy_1[i] != 0);
			//	pxy_2[i] = new double[M+2];
			//	assert(pxy_2[i] != 0); 
			//}
			
			l = m = 0;
			area1 = area2 = 0.0;
			defineTwoSubArrays(pxy_1, pxy_2, &l, &m, polyXY_2[0], polyXY_2[1],
							   M, base->coorXG[XmaxInd-1], 2);
			
			if (l > 0) {
				area1 = polygonArea(pxy_1, l);
				// Current Voronoi cell is partly outside of the grid domain
				// Looking at other grid nodes sharing the cell
				if (base->gridIn[YminInd+1][XminInd+1] == base->dummy) {
					areaT -= area1;
					cntr++;
				}
				else 
					value += base->gridIn[YminInd+1][XminInd+1]*area1;
			}
			
			if (m > 0) {
				area2 = polygonArea(pxy_2, m);
				// Current Voronoi cell is partly outside of the grid domain
				// Looking at other grid nodes sharing the cell
				if (base->gridIn[YminInd+1][XminInd] == base->dummy) {
					areaT -= area2;
					cntr++;
				}
				else 
					value += base->gridIn[YminInd+1][XminInd]*area2;
			}
			
			if (cntr != 4)
				value /= areaT;
			// Current Voronoi cell is not within the grid domain
			// A dummy value is accepted (will be modified later)
			else {
				printCellData();
				value = base->dummy; 
			}
			
			// GMnSKY2008MLE
			//for (i=0; i < 2; i++) {
			//	delete [] pxy_1[i];
			//	delete [] pxy_2[i];
			//}
			//delete [] pxy_1;
			//delete [] pxy_2;
			for (i_deall=0; i_deall < 2; i_deall++) {
				delete [] pxy_1[i_deall];
				pxy_1[i_deall] = 0;
				delete [] pxy_2[i_deall];
				pxy_2[i_deall] = 0;
			}
			delete [] pxy_1;
			pxy_1 = 0;
			delete [] pxy_2;
			pxy_2 = 0;
			// deallocate memory
			for (i_deall=0; i_deall < 2; i_deall++) {
				delete [] polyXY_1[i_deall];
				polyXY_1[i_deall] = 0;
				delete [] polyXY_2[i_deall];
				polyXY_2[i_deall] = 0;
			}
			delete [] polyXY_1;
			polyXY_1 = 0;
			delete [] polyXY_2;
			polyXY_2 = 0;

		}

		
		// The cell is intersected by numerous grid elements
		
		// NOTE: NO weighted average algorithm had ever been implemented.
		// It is cumbersome and has less importance because this situation  
		// may only happen in groundwater data, which are inherently erroneous
		
		//--------------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------------
		
		// GMnSKY2008MLE
		//BEGIN OF EDITS BY GIUSEPPE MASCARO FOLLOWING PART OF THE FOLLOWING STRUCTURE OF ALEX RINEHART MAY 2008 @ NMT

		//BEGIN OF EDITS BY ALEX RINEHART JULY 17TH 2007 @ NMT	
		/*	Created September 2006 -- AJR
		The numerous intersection of grid cells is required for the sheltering
		schemes. We have to implement it. The case of >4 grid cells intersected
		SEEMS to imply that the cell area is < the polygon area. This algorithm
		is based on this assumption.
		
		Basically, how this algorithm works is that it loops through the part of
		the active (i.e., intersecting) grid and finds the percentage of each
		element overlapping the polygon. These ratios are then summed (this will
		not be one as it is by cell, not polygon) and the inner product of active
		part of 'gridIn' and the active part of the ratio grid. This means that we
		have to construct a new grid in the h-file.*/
		
		else if ((int(fabs(double(YmaxInd-YminInd))) > 2) || (int(fabs(double(XmaxInd-XminInd))) > 2)) {
			
			ki = polyCentroid(VoronX, VoronY, nv, &xCentroid, &yCentroid, &areaT);
			if (ki > 0) {
				cout<<"\nERROR in function 'polyCentroid'!!!"<< endl;
				exit(2);
			}
			if (areaT < 0)
				areaT = fabs(areaT);
			
			// GMnSKY2008MLE
			//place holders of indices
			//L = XmaxInd;	//x-holder
			//M = YmaxInd;	//y-holder			
			//initialize fracOfGrid to all 0s
			fracOfGrid = 0.0;
			weightedAveSum = 0.0;
			
			// GMnSKY2008MLE
			//loop through active grid area
			//for ( j = XminInd; j < XmaxInd; j++ ) {
			for ( i_col = XminInd; i_col < XmaxInd; i_col++ ) {	
				//for (i = YminInd; i < YmaxInd; i++ ) {
				for (i_row = YminInd; i_row < YmaxInd; i_row++ ) {	
					
					//find coordinate of upper left corner
					//cornerX = base->xulcR + i*base->dR;
					
					//check sign
					//cornerY = base->yulcR - j*base->dR;
					
					//centralize corner vector
					//	This change of coordinates allows a simple check of inside/outside for
					//	each of the cases.
					//cornerX -= xCentroid;
					//cornerY -= yCentroid;

					// GMnSKY2008MLE
					//initialize InOrOut and Vert_InOrOut to all 0s
					for (k = 0; k < 4; k++)
						InOrOut[k] = 0;   
					for (cnt = 0; cnt < nv-1; cnt++) {
						Vert_InOrOut[cnt] = 0;
					}

					// GMnSKY2008MLE
					// If the value in the gridcell is = dummy, then do not consider the
					// contribution of that gridcell
					if (base->gridIn[i_row][i_col] == base->dummy) {
						fracOfGrid += 0.0;
						weightedAveSum += 0.0*base->gridIn[i_row][i_col];	
					}
					else
					{
						//Find coordinate of upper left corner of the grid[i_row][i_col]
						//translated with respect to the centroid.
						//This change of coordinates allows a simple check of inside/outside for
						//each of the cases.
						cornerXUL = base->xulcR + i_col*base->dR - xCentroid;
						cornerYUL = base->yulcR - i_row*base->dR - yCentroid;	
						// initialize cornerX and cornerY
						cornerX = cornerXUL;
						cornerY = cornerYUL;
					
						//loop through the corners and find which are inside/outside of the vCell object
						for (k = 0; k < 4; k++) {
						
							// GMnSKY2008MLE
							//areaTriCorn = 0.0;						
							//---------------------------------------------------------------------------
							//find out if the corner is in the polygon
						
							goodness = 0;
							cnt = 0;
						
							//move around the corners of the grid cell
							if (k == 1)
								cornerY -= base->dR;
							else if (k == 2)
								cornerX += base->dR;
							else if (k == 3)
								cornerY += base->dR;
						
							
							// GMnSKY2008MLE
							//find angle from center
							//  Different relative locations require different compuatations 
							//  of global (not -pi/2 to pi/2) angles.
							//
							//  This is still a right-handed coordinate system
							//angleCorn = fabs(atan(cornerY/cornerX));
							//if ((cornerX > 0.0) && (cornerY > 0.0)) {
							//	angleCorn = angleCorn;}
							//else if ( (cornerX < 0.0) && (cornerY > 0.0) ){
							//	angleCorn = 3.1416 - angleCorn;}
							//else if ( (cornerX < 0.0) && (cornerY < 0.0) ){
							//	angleCorn = 3.1416 + angleCorn;}
							//else if ( (cornerX > 0.0) && (cornerY < 0.0) ){
							//	angleCorn = 2*3.1416 - angleCorn;}		
						
							//initialize first vertex
							x1NormVert = VoronX[cnt] - xCentroid;
							y1NormVert = VoronY[cnt] - yCentroid;

							// GMnSKY2008MLE
							skip = 0;
						
							while(goodness == 0) {

								// GMnSKY2008MLE
								// check if cnt is greater than the number of nodes
								if (cnt > (nv-2)){
								//	cout << "cnt is "<< cnt+1 << " larger than the number of nodes nv = "<< nv << endl;
									skip = 1;
									break;
								}
								else
								{
							
									//get next vertex over
									x2NormVert = VoronX[cnt+1] - xCentroid;
									y2NormVert = VoronY[cnt+1] - yCentroid;
							
									// GMnSKY2008MLE
									//find angles
									//angleVert1 = fabs(atan(x1NormVert/y1NormVert));
									//angleVert2 = fabs(atan(x2NormVert/y2NormVert));
									//if ((x1NormVert > 0.0) && (y1NormVert > 0.0)){
									//	angleVert1 = angleVert1;}
									//else if ( (x1NormVert < 0.0) && (y1NormVert > 0.0) ){
									//	angleVert1 = 3.1416 - angleVert1;}
									//else if ( (x1NormVert < 0.0) && (y1NormVert < 0.0) ){
									//	angleVert1 = 3.1416 + angleVert1;}
									//else if ( (x1NormVert > 0.0) && (y1NormVert < 0.0) ){
									//	angleVert1 = 2*3.1416 - angleVert1;}
									//if ((x2NormVert > 0.0) && (y2NormVert > 0.0)){
									//	angleVert2 = angleVert2;}
									//else if ( (x2NormVert < 0.0) && (y2NormVert > 0.0) ){
									//	angleVert2 = 3.1416 - angleVert2;}
									//else if ( (x2NormVert < 0.0) && (y2NormVert < 0.0) ){
									//	angleVert2 = 3.1416 + angleVert2;}
									//else if ( (x2NormVert > 0.0) && (y2NormVert < 0.0) ){
									//	angleVert2 = 2*3.1416 - angleVert2;}
									// calculate the angles between the lines centroid - Vert1
									// and centroid - Vert2
									angleVert1Vert2 = AngleIncluded(x1NormVert,y1NormVert,x2NormVert,y2NormVert);
									// calculate the angles between the lines centroid - corner
									// and centroid - Vert1
									angleCornerVert1 = AngleIncluded(x1NormVert,y1NormVert,cornerX,cornerY);							
									// calculate the angles between the lines centroid - corner
									// and centroid - Vert2
									angleCornerVert2 = AngleIncluded(x2NormVert,y2NormVert,cornerX,cornerY);
									
									//decide if corner is "in" the triangle

									// GMnSKY2008MLE
									//if ((cornAngle > angleVert1)||(cornAngle < angleVert2)) {
									if ((angleVert1Vert2 > angleCornerVert1) && (angleVert1Vert2 > angleCornerVert2)) {
										//get out of the while loop while the getting is good
										goodness++;
										break;
									}
									else {
										//move onto the next pair of edges
										x1NormVert = x2NormVert;
										
										// GMnSKY2008MLE
										//y1NormVert = x2NormVert;
										y1NormVert = y2NormVert;

										cnt++;
									}//else
									
								// GMnSKY2008MLE
								} // end else statement (cnt > (nv-1))

							}//while

							// GMnSKY2008MLE
							if (skip == 0) 
							{
						
								//-------------------------------------------------------------------------
								//find areas of different constructed triangles and decide if the corner is
								//  (i)	all the way out
								//  (ii)	all the way in		
						
								//Find the area of the polygon triangle
								//    This is the triangle formed by two vertices and centroid. 
								//    The two vertices are those that lay on either side of the
								//    line from centroid to grid corner.
							
								// GMnSKY2008MLE
								//baseMag = pow(pow(x2NormVert,2.0)+pow(cornerY,2.0),0.5);
								baseMag = pow(pow(x2NormVert,2.0)+pow(y2NormVert,2.0),0.5);

								sMag = pow(pow(x1NormVert,2.0)+pow(y1NormVert,2.0),0.5);
						
								// GMnSKY2008MLE
								//x2Unit = x2NormVert/baseMag;
								//y2Unit = y2NormVert/baseMag;
								//xProj = sMag*cos(angleVert2 - angleVert1)*x2Unit;
								//yProj = sMag*cos(angleVert2 - angleVert2)*y2Unit;
								//heightMag = pow(pow(xProj,2.0) + pow(yProj,2.0),0.5);
								//areaTriPoly = 0.5*baseMag*heightMag;
								areaTriPoly = sMag*baseMag*sin(angleVert1Vert2)*0.5;
						
								//Find the area of the rightward corner triangle
								//    This is the triangle formed by the centroid, the corner of
								//    interest and the vertex that is a LEFT (as opposed to Alex's right) rotation from
								//    the line from centroid to grid corner.
								baseMag = pow(pow(cornerX,2.0)+pow(cornerY,2.0),0.5);
								sMag = pow(pow(x1NormVert,2.0)+pow(y1NormVert,2.0),0.5);
						
								// GMnSKY2008MLE
								//x2Unit = cornerX/baseMag;
								//y2Unit = cornerY/baseMag;
								//xProj = sMag*cos(angleCorn - angleVert1)*x2Unit;
								//yProj = sMag*cos(angleCorn - angleVert1)*y2Unit;
								//heightMag = pow(pow(xProj,2.0) + pow(yProj,2.0),0.5);
								//areaTriCorn = 0.5*baseMag*heightMag;
								areaTriCorn = sMag*baseMag*sin(angleCornerVert1)*0.5;
							
								//Find the area of the leftward corner triangle
								//    This is the triangle formed by the centroid, the corner of
								//    interest and the vertex that is a LEFT rotation from
								//    the line b/t centroid and corner.
								baseMag = pow(pow(cornerX,2.0)+pow(cornerY,2.0),0.5);
							
								// GMnSKY2008MLE
								//sMag = pow(pow(x1NormVert,2.0)+pow(y1NormVert,2.0),0.5);
								sMag = pow(pow(x2NormVert,2.0)+pow(y2NormVert,2.0),0.5);
								//x2Unit = cornerX/baseMag;
								//y2Unit = cornerY/baseMag;						
								//xProj = sMag*cos(-angleCorn + angleVert1)*x2Unit;
								//yProj = sMag*cos(-angleCorn + angleVert1)*y2Unit;
								//heightMag = pow(pow(xProj,2.0) + pow(yProj,2.0),0.5);
								//Add left triangle area to right triangle area to get total area
								//areaTriCorn += 0.5*baseMag*heightMag;
								areaTriCorn += sMag*baseMag*sin(angleCornerVert2)*0.5;
							
								//compare the areas and decide if in or out
								//    If the triangle formed by the vertices and centroid is larger than the
								//    sum of the triangles formed by vertices, centroid AND corner, then the 
								//    corner is on the interior of the cell.
								//
	
								// GMnSKY2008MLE
								//if (areaTriCorn > areaTriPoly)
								if (areaTriCorn <= areaTriPoly + TOL)
									InOrOut[k] = 1;	//in
								else 
									InOrOut[k] = 0;	//out
							
							// GMnSKY2008MLE
							}  // end if (skip == 0)
							else	InOrOut[k] = 0;	//out  --- Case when cnt>nv-1

						}//k-for

						// GMnSKY2008MLE
						//loop through the polygon vertexes and find which are inside/outside of grid
						for (cnt = 0; cnt < nv-1; cnt++) 
						{ // till nv-1 because the last vertex is equal to the first
							// get the vertex
							x1NormVert = VoronX[cnt] - xCentroid;
							y1NormVert = VoronY[cnt] - yCentroid;

							if ((x1NormVert >= cornerXUL) && 
									(x1NormVert <= cornerXUL+ base->dR) && 
									(y1NormVert <= cornerYUL) && 
									(y1NormVert >= cornerYUL - base->dR) )
								// update Vert_InOrOut
								Vert_InOrOut[cnt] = 1;
						}
					
						//find sum of InOrOut
						sumInOrOut = 0;
						for (k = 0; k < 4; k++) {
							sumInOrOut += InOrOut[k];
						}

						// GMnSKY2008MLE
						//find sum of Vert_InOrOut
						Vert_sumInOrOut = 0;
						for (cnt = 0; cnt < nv-1; cnt++) {
							Vert_sumInOrOut += Vert_InOrOut[cnt];
						}

						// GMnSKY2008MLE: FOLLOWING LINES OF CODE REPLACE EARLIER LINES OF CODE COMMENTED OUT FURTHER BELOW
						//check the different cases.	    
						if ((sumInOrOut == 0) && (Vert_sumInOrOut == 0))
						{
							//all the way out
							fracOfGrid += 0.0;
							weightedAveSum += 0.0*base->gridIn[i_row][i_col]; 
						}
						else if (sumInOrOut == 4) 
						{
							//all the way in
							fracOfGrid += 1.0;
							weightedAveSum += 1.0*base->gridIn[i_row][i_col];						
						}
						else 
						// here we have different cases:
						// 1) some corner grid/s inside the polygon and no polygon vertex inside the grid
						// 2) some corner grid/s inside the polygon and some polygon vertexes inside the grid
						// 3) some polygon vertex/es inside the grid and no corner grid inside the polygon 
						//
						{
							// allocate memory for polyXY_1
							polyXY_1 = new double* [2];
							assert(polyXY_1 != 0);
							for (i_all=0; i_all < 2; i_all++) {
								polyXY_1[i_all] = 0;
								polyXY_1[i_all] = new double[nv+2*nv+4];
								assert(polyXY_1[i_all] != 0);
							}
							// initialize polyXY_1
							for (i_all=0; i_all < 2; i_all++) {
								for (k=0; k<(nv+2*nv+4) ; k++)
								{
									polyXY_1[i_all][k] = 0;
									polyXY_1[i_all][k] = 0;
								}
							}							
							// initialize the counter for polyXY_1
							ind_poly = 0;
							// insert in polyXY_1 the eventual corner grid/s that are inside the polygon
							if (sumInOrOut > 0)
							{
								//find which ONE corner is IN the polygon
								for (k = 0; k < 4; k++) {
									if (InOrOut[k] == 1)
									{
										// get the coordinates of the corner grid
										cornerX = cornerXUL;
										cornerY = cornerYUL;					
							        		if  (k == 1)
						      					cornerY = cornerY - base->dR;
							 			else if (k == 2)
										{
							      				cornerX = cornerX + base->dR;
											cornerY = cornerY - base->dR;
										}
							 			else if (k == 3)
							      				cornerX = cornerX + base->dR;
									 	// update poly_XY_1
										polyXY_1[0][ind_poly] = cornerX;
										polyXY_1[1][ind_poly] = cornerY;		
										ind_poly++;
									}
								}
							} // end if (sumInOrOut > 0)
							// insert in polyXY_1 the eventual polygon vertexes inside the grid
							if (Vert_sumInOrOut > 0)
							{
								//find which vertex is IN the grid
								for (cnt = 0; cnt < nv-1; cnt++) {
									if (Vert_InOrOut[cnt] == 1)
									{
										// get the coordinates of the polygon vertex
										polyXY_1[0][ind_poly] = VoronX[cnt] - xCentroid;
							 			polyXY_1[1][ind_poly] = VoronY[cnt] - yCentroid;
										ind_poly++;
									}
								}
							} // end if (Vert_sumInOrOut > 0)
							// find the intersection points between the grid edges and the polygon edges
							for (cnt = 0; cnt < nv-1; cnt++) 
							{
								for (k = 0; k < 4; k++) 
								{ 
									skip = 0;	
									if (k==0)
									{
										// check if the polygon edge is vertical
										if (fabs(VoronX[cnt]-VoronX[cnt+1]) > TOL)
										{
											// intersection with the vertical line XUL
											coord_int = cornerXUL + xCentroid;
											flag_int = 2;
										}				
										else
											skip = 1; // lines are assumed parallel
									}
									else if	(k==1)
									{
										// check if the polygon edge is horizontal
										if (fabs(VoronY[cnt]-VoronY[cnt+1]) > TOL)
										{
											// intersection with the horizontal line YUL
											coord_int = cornerYUL  + yCentroid;
											flag_int = 1;
										}
										else
											skip = 1;
									}
									else if	(k==2)
									{
										// check if the polygon edge is vertical
										if (fabs(VoronX[cnt]-VoronX[cnt+1]) > TOL)
										{
											// intersection with the vertical line XUL+dR		
											coord_int = cornerXUL + base->dR  + xCentroid;
											flag_int = 2;
										}
										else
											skip = 1;
									}
									else if	(k==3)
									{
										// check if the polygon edge is horizontal
										if (fabs(VoronY[cnt]-VoronY[cnt+1]) > TOL)
										{
											// intersection with the horizontal line YUL-dR
											coord_int = cornerYUL - base->dR  + yCentroid;
											flag_int = 1;
										}	
										else
											skip = 1;
									}
									// if skip = 0, calculate the intersection point and determine if 
									// it is included in the polygon edge					
									if (skip == 0)
									{
										// calculate the intersection point
										findIntersectionPoint(VoronX, VoronY, cnt+1, coord_int, 
													&X_inters, &Y_inters, flag_int);
										// traslate the reference system
										X_inters = X_inters - xCentroid;
										Y_inters = Y_inters - yCentroid;
										x1NormVert = VoronX[cnt] - xCentroid;
										y1NormVert = VoronY[cnt] - yCentroid;
										x2NormVert = VoronX[cnt+1] - xCentroid;
										y2NormVert = VoronY[cnt+1] - yCentroid;
										angleVert1Vert2 = AngleIncluded(x1NormVert,y1NormVert,x2NormVert,y2NormVert);
										angleCornerVert1 = AngleIncluded(x1NormVert,y1NormVert,X_inters,Y_inters);
										angleCornerVert2 = AngleIncluded(x2NormVert,y2NormVert,X_inters,Y_inters);
										// Area X1-X2					
										baseMag = pow(pow(x2NormVert,2.0)+pow(y2NormVert,2.0),0.5);
										sMag = pow(pow(x1NormVert,2.0)+pow(y1NormVert,2.0),0.5);
										areaTriPoly = sMag*baseMag*sin(angleVert1Vert2)*0.5;
										// Area X1-Inters
										baseMag = pow(pow(X_inters,2.0)+pow(Y_inters,2.0),0.5);
										sMag = pow(pow(x1NormVert,2.0)+pow(y1NormVert,2.0),0.5);
										areaTriCorn = sMag*baseMag*sin(angleCornerVert1)*0.5;
										// Area X2-Inters
										baseMag = pow(pow(X_inters,2.0)+pow(Y_inters,2.0),0.5);
										sMag = pow(pow(x2NormVert,2.0)+pow(y2NormVert,2.0),0.5);
										areaTriCorn += sMag*baseMag*sin(angleCornerVert2)*0.5;
										if ((areaTriCorn <= areaTriPoly+TOL) && 
										 		(X_inters <= cornerXUL+base->dR + TOL) &&
												(X_inters>= cornerXUL - TOL) &&
												(Y_inters <= cornerYUL + TOL) && 
												(Y_inters >= cornerYUL - base->dR - TOL))
										{
											polyXY_1[0][ind_poly] = X_inters;
											polyXY_1[1][ind_poly] = Y_inters;
											ind_poly++;
										}			
									} // end (skip == 0)
								} // end loop for (k = 0; k < 4; k++)
							} // end for (cnt = 0; cnt < nv-1; cnt++)
							// end istructions that search for intersection points between grid and polygon edges
							// check if we have more than 2 points in polyXY_1
							if (ind_poly < 3)
							{
								// we have less than 3 points - no contribution
								fracOfGrid += 0.0;
								weightedAveSum += 0.0*base->gridIn[i_row][i_col]; 
							}
							else
							{
								// sort the points in counter-clockwise order and eventually
								// eliminate from poly_XY_1 the points that are very close
								// NOTE: after the previous instructions, ind_poly has been 
								//       incremented, so that it represents the number of 
								//       points in poly_XY_1.
								Ord_CounterClock(polyXY_1, &ind_poly);
								// check if we have more than 2 points in polyXY_1
								if (ind_poly < 3)
								{
									// we have less than 3 points - no contribution
									fracOfGrid += 0.0;
									weightedAveSum += 0.0*base->gridIn[i_row][i_col];  
								}
								else
								{
									//calculate area								
									area1 = polygonArea(polyXY_1,ind_poly);
									if (area1 < 0)
										area1 = fabs(area1);
									//find fraction of area IN the polygon
									fracOfGrid += area1/(base->dR * base->dR);
									weightedAveSum += (area1/(base->dR * base->dR)) * base->gridIn[i_row][i_col];	
								}
							}//end else of the if (ind_poly < 3)							
							// deallocate memory							
							for (i_deall=0; i_deall < 2; i_deall++) 
							{
								delete [] polyXY_1[i_deall];						
								polyXY_1[i_deall] = 0;
							}
							delete [] polyXY_1;
							polyXY_1 = 0;
						} // end else of the different cases for sumInOrOut and Vert_sumInOrOut			
					} // end else of gridvalue == dummy
				}//i_row-for
			}//i_col-for			
            if (fracOfGrid > 1e-9) { // Use a small tolerance
                value = (weightedAveSum / fracOfGrid);
            }
            else {
                // If no overlap was found, the algorithm failed.
                // Fallback to the dummy value. The higher-level
                // FindNeighbValue() function will then handle it.
                if (simCtrl->Verbose_label == 'Y') {
                    cout << "\n\n\tWarning! tResample::convertToVoronoiFormat: CASE #General\n"
                         << "The complex intersection algorithm found zero overlapping area." << endl;
                    printCellData();
                }
                value = base->dummy;
            }
			// GMnSKY2008MLE: END OF ABOVE LINES OF CODE THAT REPLACE EARLIER LINES OF CODE COMMENTED OUT BELOW

		}//else if
		
		//END OF EDITS BY ALEX RINEHART JULY 17TH 2007 @ NMT
		// END OF EDITS BY GIUSEPPE MASCARO MAY 2008 @ NMT
		//--------------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------------
		
		
		// Finding centroid of polygon: VoronX & VoronY ARE CLOSED polygons  
		else { 
			ki = polyCentroid(VoronX, VoronY, nv, &xCentroid, &yCentroid, &areaT);
			if (ki > 0) { 
				cout<<"\nError in function 'polyCentroid'!"<<endl;
				exit(2);
			}
			if (areaT < 0)
				areaT = fabs(areaT);
			
			L = XmaxInd;
			M = YmaxInd;
			findCorrespInd(xCentroid, yCentroid, &L, &M);
			
			// If the center of mass is outside of the domain
			// use vertices to define if it crosses at all
			ii = -1;
			while (((base->gridIn[M][L]) == base->dummy) && (ii < (nv-1)) ) {
				L = XmaxInd;
				M = YmaxInd;
				ii++;   //loops through the vertices
				findCorrespInd(VoronX[ii], VoronY[ii], &L, &M);
			}
			
			if (base->gridIn[M][L] != base->dummy)
				value = base->gridIn[M][L];
			// Current Voronoi cell is not within the grid domain
			// A dummy value is accepted (will be modified later)
			else {
				printCellData();
				value = base->dummy; 
			}
		}
  } // STATEMENT (flag == 1)
  
	/*These lines are being commented out because they cause a segmentation fault on the cluster
	for (i=0; i < 2; i++) {
		delete [] polyXY_1[i];
		delete [] polyXY_2[i];
	}
	delete [] polyXY_1;
	delete [] polyXY_2;
	*/

	return value;
}

/*********************************************************************
**
**  vCell:: polyCentroid()
**
**  Calculates the centroid (xCentroid, yCentroid) and area
**  of a polygon, given its vertices (x[0], y[0]) ... (x[n-1], y[n-1]). It
**  is assumed that the contour is closed, i.e., that the vertex following
**  (x[n-1], y[n-1]) is (x[0], y[0]).  The algebraic sign of the area is
**  positive for counterclockwise ordering of vertices in x-y plane;
**  otherwise negative.
**
**  Returned values:  0 for normal execution;  1 if the polygon is
**  degenerate (number of vertices < 3);  and 2 if area = 0 (and the
**  centroid is undefined).
**
**********************************************************************/
int vCell::polyCentroid(double x[], double y[], int n,
						double *xCentroid, double *yCentroid, double *area) 
{
	int i, j;
	double ai, atmp = 0.0, xtmp = 0.0, ytmp = 0.0;
	double x_i, x_j, y_i, y_j; // Added by Giuseppe Mascaro - OCTOBER 2012	
	if (n < 3) return 1;
	for (i = n-1, j = 0; j < n; i = j, j++)
    {
		// Instructions modified by Giuseppe Mascaro to deal with potential 
		// ploblems arising when the coordinates are large and the polygons 
		// have small areas.
		x_i = x[i] - x[1];
		y_i = y[i] - y[1];
		x_j = x[j] - x[1];
		y_j = y[j] - y[1];
		// End of instructions added by Giuseppe Mascaro
		
		ai = x_i * y_j - x_j * y_i;
		atmp += ai;
		xtmp += (x_j + x_i) * ai;
		ytmp += (y_j + y_i) * ai;
		
		/* OLD INSTRUCTION- ERASED ON OCTOBER 2012
		ai = x[i] * y[j] - x[j] * y[i];
		atmp += ai;
		xtmp += (x[j] + x[i]) * ai;
		ytmp += (y[j] + y[i]) * ai;
		 */
    }
	*area = atmp / 2.0;
	if (atmp != 0.0)
    {
		*xCentroid =  xtmp / (3.0 * atmp) + x[1];
		*yCentroid =  ytmp / (3.0 * atmp) + y[1];
		return 0;
    }
	return 2;
}

/*********************************************************************
**
**  vCell::findCorrespInd()
**  
**  Function:  findCorrespInd
**  Arguments: - X & Y of a polygon's centroid 
**             - ix * iy indices of a box bounding current voronoi cell
**  Objective:   to set indices to appropriate values
**  Return value: void
**  Algorithm: - if centroid is not in the "corner" -> change indices  
**
**********************************************************************/
void vCell::findCorrespInd(double xCentroid, double yCentroid, 
                           int *ix, int *iy) 
{
	while (base->coorXG[*ix] > xCentroid) {
		*ix = *ix - 1;
	}
	
	while (base->coorYG[*iy] < yCentroid) {
		*iy = *iy - 1;
	}
	return;
}

/*********************************************************************
**
**  vCell::defineTwoSubArrays()
** 
** Function:  defineTwoSubArrays (INTERSECTION WITH OX OR OY)
** Arguments: - Arrays vorX & vorY of coordinates of vor. nodes either X-s&Y-s 
**            - Value defining EDGE (either on OX or OY)
** 	   - Number of nodes - nv
** 	   - flag defining either current cell was intersected by OX or OY
** 
**         -> Arrays VoronY & VoronX (coordinates of voronoi cell) MUST be 
** 	      defined and permitted for use
** 	   -> Arrays coorYG & coorXG (coordinates of grid nodes) MUST be 
** 	      defined and permitted for use
** 	   -> Memory must be allocated for polyXY_1 & polyXY_2
** Objective: to find intersection points and assign coordninates of 2 new
**            polygons 
** Return value: void
**               - filled arrays polyXY_1 & polyXY_2 and dimensions L & M
** Algorithm: - kind of Sutherland-Hodgeman Cliping Algorithm applied for 
**              particular case of a convex polygon crossed with ONE edge
** 
**********************************************************************/
void vCell::defineTwoSubArrays(double **polyXY_1, double **polyXY_2, 
                               int *L, int *M,
                               double *vorX, double *vorY, 
                               int nn, double coor, int flag) 
{
	int i, bumts = 0;
	int l = 0, m = 0;
	double X0=-999.0, Y0=-999.0;
	double *ptr;     //Used as run-time arrays
	
	if (flag == 1)      // Intersection with _Y_ axis 
		ptr = vorY;       // Sets the current one to appropriate array 
	else if (flag == 2) // Intersection with _X_ axis 
		ptr = vorX;       // Sets the current one to appropriate array 
	
	for  (i=0; i < nn; i++) { 
		if (ptr[i] > coor) {
			if (bumts == 2) {  // has come from the OUTside
				findIntersectionPoint(vorX, vorY, i, coor, &X0, &Y0, flag);
				polyXY_1[0][l] = polyXY_2[0][m] = X0;
				polyXY_1[1][l] = polyXY_2[1][m] = Y0;
				l++;
				m++;
			}
			polyXY_1[0][l] = vorX[i];
			polyXY_1[1][l] = vorY[i];
			l++;
			bumts = 1;
		}
		else if (ptr[i] < coor) {
			if (bumts == 1) {  // has come from the INside
				findIntersectionPoint(vorX, vorY, i, coor, &X0, &Y0, flag);
				polyXY_1[0][l] = polyXY_2[0][m] = X0;
				polyXY_1[1][l] = polyXY_2[1][m] = Y0;
				l++;
				m++;
			}
			polyXY_2[0][m] = vorX[i];
			polyXY_2[1][m] = vorY[i];
			m++;
			bumts = 2;
		}
		else if (ptr[i] == coor) {
			polyXY_1[0][l] = polyXY_2[0][m] = vorX[i];
			polyXY_1[1][l] = polyXY_2[1][m] = vorY[i];
			l++;
			m++;
		}
	}
	*L = l;  // Dimension of 'polyXY_1'
	*M = m;  // Dimension of 'polyXY_2'
	return;
}

/*********************************************************************
**
**  vCell::findIntersectionPoint()
**
**  Function:  findIntersectionPoint
**  Arguments: vorX - array of X coordinates of the polygon vertices
**             vorY - array of X coordinates of the polygon vertices
**  Objective: calculates intersection point of a polygon with a line
**  Return value: void (through pointers --> coordinates X & Y)
**
**********************************************************************/
void vCell::findIntersectionPoint(double *vorX, double *vorY, int i, 
								  double coor, double *X, double *Y, int flag) 
{
	// Define a small tolerance for floating point comparisons
    const double TOLERANCE = 1e-9;

    if (flag == 1) { // Intersecting a HORIZONTAL grid line (y = coor)
        double dy = vorY[i] - vorY[i-1];
        if (std::abs(dy) < TOLERANCE) {
            // The polygon edge is horizontal. An intersection is either the entire
            // edge (if it lies on the line) or non-existent.
            // This case should be handled by the calling logic, but to be safe,
            // we can return a value that indicates failure or just return one of the points.
            *X = vorX[i-1]; // Return a known point to avoid undefined behavior
            *Y = vorY[i-1];
            // Optionally, add a warning message here for debugging.
            std::cout << "Warning: Horizontal polygon edge detected in findIntersectionPoint." << std::endl;
            return;
        }
        *X = vorX[i-1] + (coor - vorY[i-1]) * (vorX[i] - vorX[i-1]) / dy;
        *Y = coor;
    }
    else if (flag == 2) { // Intersecting a VERTICAL grid line (x = coor)
        double dx = vorX[i] - vorX[i-1];
        if (std::abs(dx) < TOLERANCE) {
            // The polygon edge is vertical.
            *Y = vorY[i-1];
            *X = vorX[i-1];
            std::cout << "Warning: Vertical polygon edge detected in findIntersectionPoint." << std::endl;
            return;
        }
        *Y = vorY[i-1] + (vorY[i] - vorY[i-1]) * (coor - vorX[i-1]) / dx;
        *X = coor;
    }
	return;
}

// GMnSKY2008MLE
/*********************************************************************
**
**  vCell::Ord_CounterClock()
**
**  Function:  Ord_CounterClock
**  Arguments: points_coord - pointers to array of coordinates the points 
**	       to be ordered in counter-clockwise order
	       dim - number of points
**  Objective: order points in counter-clockwise order 
**  Return value: void (through pointers --> coordinates points_coord)
**
**********************************************************************/
void vCell::Ord_CounterClock(double **points_coord, int *dim) 
{

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

	double x0, y0;
	int ind, count_valid;
	std::vector<std::vector<double>> tempArr(2, std::vector<double>(*dim));
	double maxX, maxY, minX, minY;
	std::vector<double> angles(*dim);
	double TOL_angle = 1e-5;
	std::vector<int> differ(*dim - 1);

	// store the values of points_coord in tempArr
	for (ind = 0; ind < *dim; ind++)
	{
		tempArr[0][ind] = points_coord[0][ind];
		tempArr[1][ind] = points_coord[1][ind];
	}

	// find a point x0 and y0 in the middle of the points
	minX = points_coord[0][0];
	maxX = points_coord[0][0];
	minY = points_coord[1][0];
	maxY = points_coord[1][0];
	for (ind = 1; ind < *dim; ind++)
	{
		minX = min(minX,points_coord[0][ind]);
		maxX = max(maxX,points_coord[0][ind]);
		minY = min(minY,points_coord[1][ind]);
		maxY = max(maxY,points_coord[1][ind]);
	}

	x0 = minX + (maxX - minX)/2;
	y0 = minY + (maxY - minY)/2;

	// translate the coordinates with respect to x0 and y0
	for (ind = 0; ind < *dim; ind++)
	{
		tempArr[0][ind] = tempArr[0][ind] - x0;
		tempArr[1][ind] = tempArr[1][ind] - y0;
	}


	for (ind = 0; ind < *dim; ind++)
	{
	// determine the angle of the point with respect to the x-axis
		if (tempArr[0][ind] == 0.0)
    			angles[ind] = 3.1416/2;		
	    	else
		{		
			angles[ind] = abs(atan((tempArr[1][ind] / tempArr[0][ind])));
		}
	

 		// update the angle value so that it varies between 0 and 2pi   
		if (tempArr[0][ind]<=0.0 && tempArr[1][ind]>=0.0)
			angles[ind] = 3.1416 - angles[ind];
	    	else if (tempArr[0][ind]<=0.0 && tempArr[1][ind]<0.0)
			angles[ind] = 3.1416 + angles[ind];
	    	else if (tempArr[0][ind]>0.0 && tempArr[1][ind]<0.0)
			angles[ind] = 2*3.1416 - angles[ind];
	}
	


	// sort the angles in increasing order	
	// and the points coordintes accordingly
	for (ind = 0; ind < *dim-1; ind++)
	{
    		double temp;
		int pos_min = ind;
		for (int ind_bis = ind+1; ind_bis < *dim; ind_bis++)
			if (angles[pos_min] > angles[ind_bis])
				pos_min = ind_bis;
		if (pos_min != ind)
		{
			temp = angles[ind];
			angles[ind] = angles[pos_min];
			angles[pos_min] = temp;
			temp = tempArr[0][ind];
			tempArr[0][ind] = tempArr[0][pos_min];
			tempArr[0][pos_min] = temp;
			temp = tempArr[1][ind];
			tempArr[1][ind] = tempArr[1][pos_min];
			tempArr[1][pos_min] = temp;
		}
	}


	// find if there is some angles with close values
	for (ind = 0; ind < *dim-1; ind++)		
	{
		if (fabs(angles[ind]-angles[ind+1]) < TOL_angle)
			differ[ind] = 1;
		else
			differ[ind] = 0;
	}
	
	
	count_valid = 0;
	// find if there is some angles with close values
	for (ind = 0; ind <*dim-1; ind++)		
	{
		if (differ[ind] == 0) // keep the value
		{
			points_coord[0][count_valid] = tempArr[0][ind] + x0;
			points_coord[1][count_valid] = tempArr[1][ind] + y0;
			count_valid++;
		}		
	}

	points_coord[0][count_valid] = tempArr[0][*dim-1] + x0;
	points_coord[1][count_valid] = tempArr[1][*dim-1] + y0;
	
	*dim = count_valid + 1;

}

// GMnSKY2008MLE
//*********************************************************************
//**
//**  vCell::AngleIncluded()
//**
//**  Function:  AngleIncluded
//**  Arguments: x1,y1 - coordinates of point 1
//**             x2,y2 - coordinates of point 2
//**  Objective: calculates the angle between the vectors O1 and O2 
//**             (where O is the origin of the coordinate system)
//**  Return value: double angle in radiant
//**
//**********************************************************************/
double vCell::AngleIncluded(double X1, double Y1, double X2, double Y2) 
{
double InnProd; 
double Mod1, Mod2;
double angle;

// Inner product
InnProd = (X1*X2 + Y1*Y2);
// modules of vectors O1 and O2
Mod1 = pow(pow(X1,2)+pow(Y1,2),0.5);
Mod2 = pow(pow(X2,2)+pow(Y2,2),0.5);

if ((Mod1 > 0) && (Mod2 > 0))
	angle = acos(InnProd / (Mod1*Mod2));
else
	angle = 0.0;

return angle;
}

/*********************************************************************
**
**  vCell::polygonArea()
**  
**  Function:  polygonArea
**  Arguments: *poly[2] - coordinates of NOT closed polygon
**              - poly[0][] - X coordinates
**              - poly[1][] - Y coordinates
**  Objective: calculates polygon area: polygon is considered UNclosed
**  Return value: - area of a polygon
**  Algorithm: [O'Rourke], page 21 
**
**********************************************************************/
double vCell::polygonArea(double **poly, int L) 
{
	double sum = 0.0;
	for (int i=0; i < L; i++) {
		if (i == (L-1))
			sum += (poly[0][i]*poly[1][0] - poly[1][i]*poly[0][0]);
		else 
			sum += (poly[0][i]*poly[1][i+1] - poly[1][i]*poly[0][i+1]);
	}
	return(fabs(sum/2.0));
}

/*********************************************************************
**
**  vCell:printCellData()
**
**  Outputs coordinates of voronoi vertices, and indices of bounding box
**
**********************************************************************/
void vCell::printCellData() 
{
	int i;
	
	if (simCtrl->Verbose_label == 'Y') {
		cout<<"--------------------------------------------"<<endl;
		cout<<"XminInd = "<<XminInd<<"\tXmaxInd = "<<XmaxInd<<endl<<flush;
		cout<<"YminInd = "<<YminInd<<"\tYmaxInd = "<<YmaxInd<<endl<<flush;
		
		cout<<"inX[i] = [ ";
		for (i=0; i < nv-1; i++)
			cout<<inX[i]<<" ";
		cout<<"]"<<endl<<flush;
		
		cout<<"inY[i] = [ ";
		for (i=0; i < nv-1; i++)
			cout<<inY[i]<<" ";
		cout<<"]"<<endl<<flush;
		
		cout<<"VoronX[i] = [ ";
		for (i=0; i < nv; i++)
			cout<<VoronX[i]<<" ";
		cout<<"]"<<endl<<flush;
		
		cout<<"VoronY[i] = [ ";
		for (i=0; i < nv; i++)
			cout<<VoronY[i]<<" ";
		cout<<"]"<<endl<<flush;
		cout<<endl;
		cout<<"--------------------------------------------"<<endl;
	}
	return;
}

//=========================================================================
//
//
//                       End of tResample.cpp
//
//
//=========================================================================

