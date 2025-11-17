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
**  meshElements.cpp: Functions for mesh element classes tNode, tEdge,
**                    and tTriangle. 
**
***************************************************************************/

#include "src/tMeshElements/meshElements.h"

#include "src/Headers/TemplDefinitions.h"
#include "src/Headers/globalIO.h"

#ifdef PARALLEL_TRIBS
#include "src/tParallel/tParallel.h"
#endif

//====================================================================
//
//
//                  Section 1: Mesh Elements Global Functions
//
//
//====================================================================

int PointsCCW( tArray< double > &, tArray< double > &, tArray< double > & );

/*************************************************************************
**
**  FindIntersectionCoords
**
**  Finds and returns intersection of line segments
**  defined by endpoints given as arguments; 1st seg. endpts are xy1, xy2;
**  2nd seg. endpts are xy3, xy4. (SL) 
**
*************************************************************************/

tArray< double > FindIntersectionCoords( tArray< double > xy1, tArray< double > xy2,
                                         tArray< double > xy3, tArray< double > xy4 )
{
	double dxa, dxb, dya, dyb, a, b, c, f, g, h;
	
	tArray< double > intxy(2);
	
	dxa = xy2[0] - xy1[0];
	dxb = xy4[0] - xy3[0];
	dya = xy2[1] - xy1[1];
	dyb = xy4[1] - xy3[1];
	a = dya;
	b = -dxa;
	c = dxa * xy1[1] - dya * xy1[0];
	f = dyb;
	g = -dxb;
	h = dxb * xy3[1] - dyb * xy3[0];  
	
	if( fabs(dxa) > THRESH && fabs(dxb) > THRESH ){
		if( fabs(f - g * a / b) > 0 ){
			intxy[0] = (g * c / b - h) / (f - g * a / b);
			intxy[1] = (-c - a * intxy[0]) / b;
		}
	}
	else{
		if( fabs(dya) > THRESH && fabs(dyb) > THRESH ){
			if( fabs(g - f * b / a) > 0 ){
				intxy[1] = (f * c / a - h) / (g - f * b / a);
				intxy[0] = (-c - b * intxy[1]) / a;
			}
		}
		else {
			if( fabs(dya) <= THRESH ){
				intxy[0] = xy3[0];
				intxy[1] = xy1[1];
			}
			else{
				intxy[0] = xy1[0];
				intxy[1] = xy3[1];
			}
		}
	}
	return intxy;
}


//=========================================================================
//
//
//                  Section 2: tNode Class Functions
//
//
//=========================================================================            


/***********************************************************************
**
**  tNode::insertFrontSpokeList
**  tNode::insertBackSpokeList
**
**  Places eptr at the front or back of the spoke list (respectively)
**  and makes the list circular.
**
***********************************************************************/

void tNode::insertFrontSpokeList( tEdge *eptr ){
	spokeList.insertAtFront( eptr );
	assert( spokeList.getFirst() != 0 );
	makeWheel();
}

void tNode::insertBackSpokeList( tEdge *eptr ){
	spokeList.insertAtBack( eptr );
	assert( spokeList.getFirst() != 0 );
	makeWheel();
}

/***********************************************************************
**
**  tNode::NextSpoke()
**
***********************************************************************/          

const tEdge *tNode::NextSpoke( tPtrListNode< tEdge > * current ){
	if( current == 0 || current == spokeList.getLast() )
		current = spokeList.getFirstNC();
	else current = current->getNextNC();
	return current->getPtr();
}

/*****************************************************************************
**
**  tNode::AttachFirstSpoke
**
**  Attaches the first spoke to the node by pointing edg to that spoke,
**  and then telling the spoke to point to itself. thespoke is the edge
**  being added.
**
*****************************************************************************/

void tNode::AttachFirstSpoke( tEdge *thespoke ){
	assert( thespoke!=0 );
	assert( thespoke->getOriginPtr()==this );
	edg = thespoke;
	thespoke->setCCWEdg( thespoke );
}

/*****************************************************************************
**
**  tNode::Dist
**
**  Computes the distance of the node from the line formed
**  by points p0 p1 using x y.
**
*****************************************************************************/

double tNode::Dist( tNode * p0, tNode * p1 ){
	double a,b,c,res;
	
	a=(p1->y)-(p0->y);  
	b=-((p1->x)-(p0->x));
	c=-((a*(p0->x))+(b*(p0->y)));
	res=(a*x + b*y + c) / sqrt(a*a + b*b);
	if (res<0) res=-res;
	return(res);
}

/*****************************************************************************
**
**  tNode::EdgToNod
**
**  Finds and returns the spoke (edge) that connects the current node to _nod_,
**  or zero if no such spoke is found.
**
*****************************************************************************/

tEdge *tNode::EdgToNod( tNode * nod ){

	if (this->spokeList.getSize() != 0) {
		tPtrListIter< tEdge > spokIter( this->spokeList );
		tEdge * ce;
	
		for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() ){
			if( ce->getDestinationPtr()->getID() == nod->getID() ) return ce;
		}
	}

	// 
	// This version for mesh builder files which have thrown out the spokes
	// but have retained the connectivity through the ccw edge pointers PKF
	// 
	else {
		const tNode* destNode = (tNode*) edg->getDestinationPtr();
		if (destNode->getID() == nod->getID())
			return edg;
		tEdge* nextedge = edg->getCCWEdg();
		while (nextedge != edg) {
			destNode = nextedge->getDestinationPtr();
			if (destNode->getID() == nod->getID())
				return nextedge;
			nextedge = nextedge->getCCWEdg();
		}  
	}
	return 0;
}     

/*****************************************************************************
**
**  tNode::ComputeVoronoiArea
**
**  Computes the node's Voronoi area by summing the area of embedded
**  triangles, and also calls CalcSpokeVEdgLengths to compute the length of
**  the sides of the Voronoi cell.
**
**  The basic Voronoi polygon is described by the set of "right-hand
**  Voronoi vertices" associated with each spoke (edge). These vertices
**  are computed by setVoronoiVertices() as the intersection of the
**  perpendicular bisectors of consecutive spokes. However, in some cases
**  the basic polygon will be distorted, with consecutive vertices NOT
**  being counter-clockwise from one another. (This seems to be the result
**  of numerical errors (possibly in estimating the intersection of two
**  nearly parallel bisectors); according to Sugihara and Iri
**  [J. Comp. Geom. & Appl., 1994, v. 4, p. 179], each Delaunay triangle
**  should be associated with one Voronoi vertex). This is handled by
**  detecting "loops" in the Voronoi polygon and cutting them off by
**  taking the area of the closest (counterclockwise) intersection of 
**  perpendicular bisectors.
**
*****************************************************************************/

double tNode::ComputeVoronoiArea()
{
	int cw;
	double area = 0;
	double a, b, c, dx, dy, dx0, dx1, dy0, dy1, dx2, dy2;
	double vx, vy, x0, y0, x1, y1, x2, y2, m1, m2;
	double temparea; 
	tEdge *ce, *cen, *cenn, *cennn, *edgptr;
	tPtrList< tEdge > vedgList;
	tPtrListIter< tEdge > vtxIter( vedgList );
	tList< tArray< double > > vcL;  
	tListIter< tArray< double > > vcI( vcL ); 
	tArray< double > xy, xyn, xynn, xynnn, xy1, xy2, xy3, xy4;
	int i;
	
	// Create a duplicate list of edges; we will modify this list to obtain
	// the correct vertices. In some cases, we may need to delete an edge
	// to get the correct list of vertices; we don't want to delete the
	// spoke ptr, so we make a duplicate list.
	
	ce = edg;
	do{
		assert( ce != nullptr );
		vedgList.insertAtBack( ce );
		ce = ce->getCCWEdg();
	} while( ce != edg );
	
	vedgList.makeCircular();
	
	// Check boundary status: Voronoi area only defined for non-boundary nodes
	
	if(( boundary == 0 ) || (boundary == 3))
	{
		cw = TRUE;
		do
		{
			cw = FALSE;
			// For each edge, look at the next two to see if clockwise
			for( ce=vtxIter.FirstP(); !( vtxIter.AtEnd() ); ce=vtxIter.NextP() )
			{
				// First edge
				xy = ce->getRVtx();

				// Are any of the initial voronoi vertices equal to the centroid
				if (x == xy[0] && y == xy[1])
					cerr << "Cell " << id << ": voronoi vertex = centroid" << endl;

				// Second edge
				cen = vtxIter.NextP();
				xyn = cen->getRVtx();

				// Third edge
				cenn = vtxIter.NextP();
				xynn = cenn->getRVtx();

				// Are vertices 1,2,3 clockwise
				// Difference between voronoi vertex 2 and 3
				dx0 = xynn[0] - xyn[0];	
				dy0 = xynn[1] - xyn[1];

				// Difference between voronoi vertex 1 and 2
				dx1 = xy[0] - xyn[0];
				dy1 = xy[1] - xyn[1];

				// If the differences are very small set them to be 0
				if (fabs(dx0) < THRESH) dx0 = 0.0;
				if (fabs(dy0) < THRESH) dy0 = 0.0;
				if (fabs(dx1) < THRESH) dx1 = 0.0;
				if (fabs(dy1) < THRESH) dy1 = 0.0;

				// Voronoi vertices 1,2,3 clockwise
				if( dy0 * dx1 > dx0 * dy1 )
				{
					// Fourth edge
					cennn = vtxIter.NextP();
					xynnn = cennn->getRVtx();

					// Are vertices 2,3,4 clockwise
					dx0 = xynnn[0] - xynn[0];
					dy0 = xynnn[1] - xynn[1];
					dx1 = xyn[0] - xynn[0];
					dy1 = xyn[1] - xynn[1];
					
					// If the differences are very small set them to be 0
					if (fabs(dx0) < THRESH) dx0 = 0.0;
					if (fabs(dy0) < THRESH) dy0 = 0.0;
					if (fabs(dx1) < THRESH) dx1 = 0.0;
					if (fabs(dy1) < THRESH) dy1 = 0.0;

					// Voronoi vertices 2,3,4 clockwise also
					if( dy0 * dx1 > dx0 * dy1 )
					{
						//two consecutive clockwise vertices=>want intersection
						//of bisectors of edges 2 and 4
						cw = TRUE;

						// Centroid
						x0 = x;
						y0 = y;

						// Location of first edge destination node
						xy1 = ce->getDestinationPtr()->get2DCoords();

						// Location of third edge destination node
						xy2 = vtxIter.PrevP()->getDestinationPtr()->get2DCoords();

						// Midpoint between centroid and first voronoi vertex
						x1 = ( x0 + xy1[0] ) / 2;
						y1 = ( y0 + xy1[1] ) / 2;

						// Midpoint between centroid and second voronoi vertex
						x2 = ( x0 + xy2[0] ) / 2;
						y2 = ( y0 + xy2[1] ) / 2;

						// Differences between centroid and two midpoints
						dx1 = x1 - x0;
						dy1 = y1 - y0;
						dx2 = x2 - x0;
						dy2 = y2 - y0;

						if( fabs(dy1)>0 && fabs(dy2) > 0 ){
							m1 = -dx1/dy1;
							m2 = -dx2/dy2;
							
							if (m1 == m2){
								cout<<"\n\n1: Point Correction";
								cout<<"\nOrigin X = "<<ce->getOriginPtr()->getX();
								cout<<"\nOrigin Y = "<<ce->getOriginPtr()->getY();
								cout<<"\nOrigin Z = "<<ce->getOriginPtr()->getZ();
								cout<<"\nOrigin B = "<<ce->getOriginPtr()->getBoundaryFlag();
								cout<<"\n\nDestination X = "<<ce->getDestinationPtr()->getX();
								cout<<"\nDestination Y = "<<ce->getDestinationPtr()->getY();
								cout<<"\nDestination Z = "<<ce->getDestinationPtr()->getZ();
								cout<<"\nDestination B = "<<ce->getDestinationPtr()->getBoundaryFlag();
								cout<<"\n\n";
							}
							
							assert ( m1 != m2 );
							vx = (y2-m2*x2-y1+m1*x1) / (m1-m2);
							vy = m1*(vx-x1)+y1;
						}
						else{
							if( fabs(dx1) > 0 && fabs(dx2) > 0 ){
								m1 = dy1/dx1;
								m2 = dy2/dx2;
								
								if (m1 == m2){
									cout<<"\n\n2: Point Correction";
									cout<<"\nOrigin X = "<<ce->getOriginPtr()->getX();
									cout<<"\nOrigin Y = "<<ce->getOriginPtr()->getY();
									cout<<"\nOrigin Z = "<<ce->getOriginPtr()->getZ();
									cout<<"\nOrigin B = "<<ce->getOriginPtr()->getBoundaryFlag();
									cout<<"\n\nDestination X = "<<ce->getDestinationPtr()->getX();
									cout<<"\nDestination Y = "<<ce->getDestinationPtr()->getY();
									cout<<"\nDestination Z = "<<ce->getDestinationPtr()->getZ();
									cout<<"\nDestination B = "<<ce->getDestinationPtr()->getBoundaryFlag();
									cout<<"\n\n";
								}
								
								assert ( m1 != m2 );
								vy=(m1*y1+x1-m2*y2-x2)/(m1-m2);
								vx= -vy*m1+m1*y1+x1;
							}
							else if( fabs(dx1) > 0 ){
								vx = x1;
								vy = y2;
							}
							else{
								vx = x2;
								vy = y1;
							}
						}

						// Second edge
						edgptr = vtxIter.PrevP();
						xyn[0] = vx;
						xyn[1] = vy;
						dx = xy[0] - vx;
						dy = xy[1] - vy;
						
						// Is new voronoi vertex equal to the centroid
						if (x == vx && y == vy)
							cerr << "Cell " << id << ": new vertex = centroid" << endl;

						edgptr->setVEdgLen( sqrt( dx*dx + dy*dy ) );
						edgptr->setRVtx( xyn );

						// Third edge
						edgptr = vtxIter.ReportNextP();
						edgptr->setVEdgLen(0.0);
						edgptr->setRVtx( xynnn );
						edgptr = 0;
						
						// Remove third edge
						vedgList.removeNext( edgptr, vtxIter.NodePtr() );
					}
				}
				vtxIter.Get( ce->getID() );
			}
		} while( cw ); //while we're still finding loops in the polygon
		
		//Before the next step, make a list of V. vertex coord. arrays.
		//In doing so, check for parts of the V. area lying outside the
		//mesh domain and cut them off by substituting coordinates of
		//two intersections of V. edges with boundary edge for the V.
		//vertex lying outside the boundary. This should take care of any
		//outlying area as long as all boundaries are convex.
		// Go through spokes and put RVtx of ccw edge in coord list, but
		// first check that the vtx lies within the bndies
		
		tEdge *ne, *nne;
		tNode *bn0, *bn1;
		for( ce = vtxIter.FirstP(); !(vtxIter.AtEnd()); ce = vtxIter.NextP() )
		{
			ne = ce->getCCWEdg();
			xy1 = ne->getRVtx();
			//checking polygon edge is on boundary and ccw edge's RVtx is on
			//wrong side of bndy edge...
			// CJC 2025 debug if statement wouldn't execute for voronoi connected to outlet point causing area to be equal to last node
			if( ce->getDestinationPtrNC()->getBoundaryFlag() != kNonBoundary && ne->getDestinationPtrNC()->getBoundaryFlag() != kNonBoundary )
			{
				bn0 = ce->getDestinationPtrNC();
				bn1 = ne->getDestinationPtrNC();
				xy2 = bn0->get2DCoords();
				xy3 = bn1->get2DCoords();
				if( !PointsCCW( xy1, xy2, xy3 ) )
				{					
					xy = FindIntersectionCoords( ce->getRVtx(), xy1, xy2, xy3 );
					vcL.insertAtBack( xy );
					nne = ne->getCCWEdg();
					xy = FindIntersectionCoords( xy1, nne->getRVtx(), xy2, xy3 );
					vcL.insertAtBack( xy );
				}
				else vcL.insertAtBack( xy1 );
			}
			else vcL.insertAtBack( xy1 );
		}
		
		// Now that we've found the correct vertices, make triangles to
		// fill the polygon; the sum of the tri areas is the v. area.
		// For a convex polygon, we can compute the total area as the
		// sum of the area of triangles [P(1) P(i) P(i+1)] for i=2,3...N-1.
		
		// coords of first vertex:
		xy = *(vcI.FirstP()); //ce = vtxIter.FirstP();
		
		// Find out # of vertices in polygon:
		int nverts = vcL.getSize(); 
		for( i=2; i<=nverts-1; i++ )
		{
			xyn = *(vcI.NextP()); 
			xynn = *(vcI.NextP());
			
			dx = xyn[0] - xy[0];
			dy = xyn[1] - xy[1];
			
			a = sqrt( dx*dx + dy*dy );
			dx = xynn[0] - xyn[0];
			dy = xynn[1] - xyn[1];
			
			b = sqrt( dx*dx + dy*dy );
			dx = xynn[0] - xy[0];
			dy = xynn[1] - xy[1];
			
			c = sqrt( dx*dx + dy*dy );
			
			//Added by smr to take care of rounding error, leading to negative value
			
			temparea = 4*a*a*b*b - (c*c - (b*b + a*a))*(c*c - (b*b + a*a));
			if (temparea < 0) temparea = 0;
			
			area += 0.25*sqrt( temparea);
			
			vcI.Prev();
		}
	}
	varea = area;
	varea_rcp = 1.0/varea;
	
	return area;
}

/*******************************************************************
**
**  tNode::getVoronoiVertexList()
**
**  Creates and returns a list of (x,y) coordinates for the
**  Voronoi vertices associated with the node. The list is 
**  created by moving around the spokes and adding each spoke's
**  right-hand Voronoi vertex to the back of the list.
**  A pointer to the vertex list is passed as a parameter; any
**  prior contents are flushed before the list of points is created.
**
*******************************************************************/

void tNode::getVoronoiVertexList( tList<Point2D> * vertexList ){
	tArray<double> vtxarr;
	Point2D vtx;
	assert( !boundary );
	vertexList->Flush();
	
	// Loop around spokes, adding the right-hand Voronoi vertex of
	// each to the list, until we've gone all the way around
	tEdge *ce = edg;
	do{
		vtxarr = ce->getRVtx();
		vtx.x = vtxarr[0];
		vtx.y = vtxarr[1];
		vertexList->insertAtBack( vtx );
		ce = ce->getCCWEdg();
	}
	while( ce!=edg );
	
	assert( vertexList->getSize()!=0 );
}


/*******************************************************************
**
**  tNode::getVoronoiVertexXYZList()
**
**  Creates and returns a list of (x,y,z) coordinates for the
**  Voronoi vertices associated with the node. The list is 
**  created by moving around the spokes and adding each spoke's
**  right-hand Voronoi vertex to the back of the list. The z
**  coordinate is obtained by linear interpolation from the 3
**  points of the triangle of which the vertex is the circumcenter.
**  A pointer to the vertex list is passed as a parameter; any
**  prior contents are flushed before the list of points is created.
**
*******************************************************************/

void tNode::getVoronoiVertexXYZList( tList<Point3D> * vertexList )
{
	tArray<double> vtxarr, zvals(3), p0(2), p1(2), p2(2);
	Point3D vtx;
	tNode *n1, *n2;
	assert( !boundary );
	vertexList->Flush();
	
	// Loop around spokes, adding the right-hand Voronoi vertex of
	// each to the list, until we've gone all the way around
	tEdge *ce = edg;
	n2 = ce->getDestinationPtrNC();
	do{
		ce = ce->getCCWEdg();
		n1 = n2;
		n2 = ce->getDestinationPtrNC();
		vtxarr = ce->getRVtx();
		vtx.x = vtxarr[0];
		vtx.y = vtxarr[1];
		zvals[0] = this->z;
		zvals[1] = n1->getZ();
		zvals[2] = n2->getZ();
		vtx.z=PlaneFit(vtx.x,vtx.y,this->get2DCoords(),n1->get2DCoords(),n2->get2DCoords(),zvals);
		vertexList->insertAtBack( vtx );
		
	}
	while( ce!=edg );
	
	assert( vertexList->getSize()!=0 );
}


/*******************************************************************
**
**  tNode::makeCCWEdges
**
**  This function provides for compatibility between the CCW Edge
**  data structure and the Spoke List data structure. It sets up
**  CCW edge connectivity from the spoke list data (which is 
**  assumed to be up to date) by: (1) setting the node's edg 
**  pointer to the first spoke on the list, and (2) setting the
**  ccwedg pointer for each spoke.
**
*******************************************************************/

void tNode::makeCCWEdges()
{
	tEdge *ce, *ccwe;
	tPtrListIter< tEdge > spokIter( spokeList );
	
	ce = spokIter.FirstP();
	assert( ce != 0 );
	setEdg( ce );
	
	for( ; !(spokIter.AtEnd()); ce = spokIter.NextP() ){
		ccwe = spokIter.ReportNextP();
		assert( ccwe != 0 );
		ce->setCCWEdg( ccwe );
	}
}

/*******************************************************************
**
**  tNode::ConvertToClosedBoundary
**
**  Makes the node into a closed boundary by setting its boundary
**  status flag. The function also updates the boundary ("flow
**  allowed") status of the node's spokes and their complements.
**
******************************************************************/

void tNode::ConvertToClosedBoundary()
{
	tEdge *ce;   // an edge and its complement
	
	// Reset boundary flag
	boundary = kClosedBoundary;
	
	// Signal all connected edges and their complements to become no-flow
	ce = edg;
	do{
		assert( ce!=0 );
		if( ce->getBoundaryFlag()==kFlowAllowed ){
			ce->setFlowAllowed( 0 );
		}  
	} while( (ce=ce->getCCWEdg()) != edg );
}

//#ifdef PARALLEL_TRIBS
/**************************************************************************
**
**  get or set ResIndex
**
**  The resample index may not be the same as the node ID 
**  when running in parallel. So an additional index is
**  available for resampling.
**
**************************************************************************/
                                                                                 
void tNode::setResIndex(int rind) {
  resIndex = rind;
}  
                                                                                             
int tNode::getResIndex() {
    return resIndex;
}  
   
//#endif 
  
/**************************************************************************
**
**  TellAll
**
**  Debugging routine.
**
**************************************************************************/

#ifndef NDEBUG
void tNode::TellAll()
{
	cout << " NODE " << id << ":\n";
	cout << "  x=" << x << " y=" << y << " z=" << z;
	cout << "  boundary: " << boundary
        << "\n  varea: " << varea << endl;
	if( edg )
		cout << "  points to edg #" << edg->getID() << endl;
	else cout << "  edg is undefined!\n";
	
}
#endif

//=========================================================================
//
//
//                  Section 3: tEdge Class Functions
//                                                    
//                                                    
//
//=========================================================================

/**************************************************************************
**
**  tEdge::CalcLength
**
**  Computes the edge length and returns it. (Length is the projected
**  on the x,y plane). Assumes org and dest are valid.
**
**************************************************************************/

double tEdge::CalcLength(){
	assert( org!=0 );  
	assert( dest!=0 );
	
	double dx = org->getX() - dest->getX();
	double dy = org->getY() - dest->getY();
	len = sqrt( dx*dx + dy*dy );
	return len;
}

/**************************************************************************
**
**  tEdge::TellCoords
**
**  Debugging routine that reports edge ID and coords of endpoints.
**
**************************************************************************/
#ifndef NDEBUG
void tEdge::TellCoords(){
	cout << "EDGE " << id << ":\n";
	cout << "  " << org->getID() << " (" << org->getX() << ","
        << org->getY() << ") -> " << dest->getID() << " ("
        << dest->getX() << "," << dest->getY() << ")" << endl;
}
#endif

/**************************************************************************
**
**  tEdge::FindComplement
**
**  Finds and returns the edge's complement edge. Does this by checking
**  the spokes connected to its destination node and returning the one
**  that connects back to its origin.
**
**  Returns:  ptr to the complement edge
**  Assumes:  valid connectivity (one of destination's spokes connects
**            back to origin)
**
**************************************************************************/

tEdge * tEdge::FindComplement(){
	assert( org!=0 && dest!=0 && dest->getEdg()!=0 );
	
	tEdge * ce = dest->getEdg();
	while( ce->getDestinationPtrNC() != org ){
		ce = ce->getCCWEdg();
		assert( ce!=0 && ce!=dest->getEdg() );
	}
	return ce;
}

//=========================================================================
//
//
//                  Section 4: tTriangle Class Functions
//                                                    
//                                                    
//
//=========================================================================


/*****************************************************************************
**
**  tTriangle::FindCircumcenter
**
**  Finds the circumcenter of the triangle by finding the intersection of
**  the perpendicular bisectors of sides (p0,p1) and (p0,p2). Returns the
**  coordinates of the circumcenter as a 2-element array. Note that the
**  circumcenter is also the Voronoi cell vertex associated with the
**  triangle's three nodes (that's the point of computing it).
**
*****************************************************************************/

tArray< double > tTriangle::FindCircumcenter()
{
	double x1, y1, x2, y2, dx1, dy1, dx2, dy2, m1, m2;
	tArray< double > xyo, xyd1, xyd2, xy(2);
	
	assert( pPtr(0) && pPtr(1) && pPtr(2) );
	
	// Coordinates of triangle's nodes p0, p1, and p2
	xyo = pPtr(0)->get2DCoords();
	xyd1 = pPtr(1)->get2DCoords();
	xyd2 = pPtr(2)->get2DCoords();
	
	// Find the midpoints of the two sides (p0,p1) and (p0,p2) and store them 
	// in (x1,y1) & (x2,y2). Then get the distance between p0 and the 
	// midpoints of each side
	
	x1 = (xyo[0] + xyd1[0]) / 2;
	y1 = (xyo[1] + xyd1[1]) / 2;
	x2 = (xyo[0] + xyd2[0]) / 2;
	y2 = (xyo[1] + xyd2[1]) / 2;
	dx1 = x1-xyo[0];
	dy1 = y1-xyo[1];
	dx2 = x2-xyo[0];
	dy2 = y2-xyo[1];
	
	// Compute the intercept of the bisectors of the two sides:
	// Case: neither spoke is horizontal (ok to divide by dy1 and dy2)
	
	if( fabs(dy1)>0 && fabs(dy2)>0 ){
		assert( dy1!=0 && dy2!=0 );
		m1= -dx1/dy1;
		m2 = -dx2/dy2;
		if (m1 == m2){
			cout<<"\nPoint 0 X = "<<pPtr(0)->getX();
			cout<<"\nPoint 0 Y = "<<pPtr(0)->getY();
			cout<<"\nPoint 0 Z = "<<pPtr(0)->getZ();
			cout<<"\nPoint 0 B = "<<pPtr(0)->getBoundaryFlag();
			cout<<"\n\nPoint 1 X = "<<pPtr(1)->getX();
			cout<<"\nPoint 1 Y = "<<pPtr(1)->getY();
			cout<<"\nPoint 1 Z = "<<pPtr(1)->getZ();
			cout<<"\nPoint 1 B = "<<pPtr(0)->getBoundaryFlag();
			cout<<"\n\nPoint 2 X = "<<pPtr(2)->getX();
			cout<<"\nPoint 2 Y = "<<pPtr(2)->getY();
			cout<<"\nPoint 2 Z = "<<pPtr(2)->getZ();
			cout<<"\nPoint 2 B = "<<pPtr(0)->getBoundaryFlag();
			cout<<"\n\n";
		}
		assert( m1!=m2 ); // should never happen; means edges are parallel
		xy[0] = (y2 - m2 * x2 - y1 + m1 * x1) / (m1 - m2);
		xy[1] = m1 * (xy[0] - x1) + y1;
	}
	
	// Case: one spoke is horizontal, but neither are vertical
	
	else if( dx1!=0 && dx2!=0 ){
		assert( dx1!=0 && dx2!=0 );
		m1 = dy1/dx1;
		m2 = dy2/dx2;
		if (m1 == m2){
			cout<<"\nPoint 0 X = "<<pPtr(0)->getX();
			cout<<"\nPoint 0 Y = "<<pPtr(0)->getY();
			cout<<"\nPoint 0 Z = "<<pPtr(0)->getZ();
			cout<<"\nPoint 0 B = "<<pPtr(0)->getBoundaryFlag();
			cout<<"\n\nPoint 1 X = "<<pPtr(1)->getX();
			cout<<"\nPoint 1 Y = "<<pPtr(1)->getY();
			cout<<"\nPoint 1 Z = "<<pPtr(1)->getZ();
			cout<<"\nPoint 1 B = "<<pPtr(0)->getBoundaryFlag();
			cout<<"\n\nPoint 2 X = "<<pPtr(2)->getX();
			cout<<"\nPoint 2 Y = "<<pPtr(2)->getY();
			cout<<"\nPoint 2 Z = "<<pPtr(2)->getZ();
			cout<<"\nPoint 2 B = "<<pPtr(0)->getBoundaryFlag();
			cout<<"\n\n";
		}
		assert( m1!=m2 );
		xy[1] = (m1 * y1 + x1 - m2 * y2 - x2) / (m1 - m2);
		xy[0] = -xy[1] * m1 + m1 * y1 + x1;
	}
	
	// Special case: one is vertical, the other horizontal
	
	else{
		if( dx1!=0 ){
			xy[0] = xyo[0] + dx1;
			xy[1] = xyo[1] + dy2;
			assert( dx2==0 && dy1==0 );
		}
		else{
			xy[0] = xyo[0] + dx2;
			xy[1] = xyo[1] + dy1;
			assert( dx1==0 && dy2==0 );
		}
	}
	
	return xy;
}

/*****************************************************************************
**
**  tTriangle::TellAll()
** 
**  Debugging Routing
**
*****************************************************************************/

#ifndef NDEBUG
void tTriangle::TellAll()
{
	int i;
	
	//assert( this!=0 );//WR--09192023:'this' pointer cannot be null in well-defined C++ code; comparison may be assumed to always evaluate to true
	cout << "TRIANGLE #" << id << ":\n";
	for( i=0; i<3; i++ ){
		cout << "  P" << i << " ";
		if( p[i]!=0 ) cout << p[i]->getID() << " (" << p[i]->getX() << ","
			<< p[i]->getY() << ")";
		else cout << "(ndef)";
		cout << "  E" << i << " ";
		if( e[i]!=0 ) cout << e[i]->getID();
		else cout << "(ndef)";
		cout << "  T" << i << " ";
		if( t[i]!=0 ) cout << t[i]->getID();
		else cout << "(ndef)";
		cout << endl;
	}
}

#endif

//=========================================================================
//
//
//              End of MeshElements.cpp
//
//
//=========================================================================
