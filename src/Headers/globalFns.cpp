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
**
**  globalFns.cpp: Global tRIBS Class Functions
**  
**  Uses tNode objects in NewTriCCW, InNewTri, Intersect
**
***************************************************************************/

#include "src/Headers/globalFns.h"
#include "src/Headers/TemplDefinitions.h"

//=========================================================================
//
//
//                  Section 1: globalFns Functions
//
//
//=========================================================================


/***************************************************************************
**  
**  globalFns::UnitVector( tEdge* ePtr )
** 
***************************************************************************/
tArray< double > UnitVector( tEdge* ePtr )
{
	assert( ePtr != 0 );
	tArray< double > oxy( ePtr->getOriginPtr()->get2DCoords() );
	tArray< double > dxy( ePtr->getDestinationPtr()->get2DCoords() );
	double dx = dxy[0] - oxy[0];
	double dy = dxy[1] - oxy[1];
	double mag = sqrt( dx * dx + dy * dy );
	tArray< double > uvect(2);
	uvect[0] = dx / mag;
	uvect[1] = dy / mag;
	return uvect;
}

/***************************************************************************
**  
**  globalFns::FindCosineAngle0_2_1(  )
** 
***************************************************************************/
double FindCosineAngle0_2_1( tArray< double > &p0,
                             tArray< double > &p1,
                             tArray< double > &p2 )
{
	//assert( (&p0 != 0) && (&p1 != 0) && (&p1 != 0) );
	double dx0, dx1, dy0, dy1;
	double dotp, magp;
	dx0 = p0[0] - p2[0];
	dx1 = p1[0] - p2[0];
	dy0 = p0[1] - p2[1];
	dy1 = p1[1] - p2[1];
	dotp = dx0 * dx1 + dy0 * dy1;
	magp = sqrt( dx0 * dx0 + dy0 * dy0 ) * sqrt( dx1 * dx1 + dy1 * dy1 );
	return dotp / magp;
}

/***************************************************************************
**  
**  globalFns::TriPasses(  )
**
**  Check for Delauny Triangulation
** 
***************************************************************************/
int TriPasses( tArray< double > &ptest,
               tArray< double > &p0,
               tArray< double > &p1,
               tArray< double > &p2 )
{
	//assert( (&ptest != 0) && (&p0 != 0) && (&p1 != 0) && (&p1 != 0) );//WR--09192023
	double dx0, dx1, dy0, dy1;
	double crossp, dotp, angle0_2_1, angle0_test_1;
	
	// Find angle p0-p2-p1
	dx0 = p0[0] - p2[0];
	dx1 = p1[0] - p2[0];
	dy0 = p0[1] - p2[1];
	dy1 = p1[1] - p2[1];
	crossp = dx0 * dy1 - dx1 * dy0;
	dotp = dx0 * dx1 + dy0 * dy1;
	angle0_2_1 = atan2( crossp, dotp );
	
	// Find angle p0-ptest-p1
	dx0 = p0[0] - ptest[0];
	dx1 = p1[0] - ptest[0];
	dy0 = p0[1] - ptest[1];
	dy1 = p1[1] - ptest[1];
	crossp = dx0 * dy1 - dx1 * dy0;
	dotp = dx0 * dx1 + dy0 * dy1;
	angle0_test_1 = atan2( crossp, dotp );
	
	// Compare and return the result
	if( angle0_2_1 < angle0_test_1 ){
		return 0;
	}
	else{
		return 1;
	}
}

/***************************************************************************
**  
**  globalFns::PointsCCW(  )
** 
**  Determines whether points are counter-clockwise
** 
***************************************************************************/
int PointsCCW( tArray< double > &p0,
               tArray< double > &p1,
               tArray< double > &p2 )
{
	//assert( &p0 != 0 && &p1 != 0 && &p1 != 0 ); //WR--09192023
	double* a0 = p0.getArrayPtr();
	double* a1 = p1.getArrayPtr();
	double* a2 = p2.getArrayPtr();
	
	return ( predicate.orient2d( a0, a1, a2 ) > 0 );
}

/***************************************************************************
**  
**  globalFns::NewTriCCW(  tTriangle *ct )
** 
**  Determines whether triangle's "new" coords are CCW
**  Modified for use in tRIBS. Meandering removed
** 
***************************************************************************/
int NewTriCCW( tTriangle *ct )
{
	assert( ct != 0 );
	
	tNode *cn;
	cn = (tNode *) ct->pPtr(0);
	tArray< double > p0( cn->get2DCoords() );
	
	//Meandering off in tRIBS
	//-----------------------
	//if( cn->Meanders() ) p0 = cn->getNew2DCoords();
	
	cn = (tNode *) ct->pPtr(1);
	tArray< double > p1( cn->get2DCoords() );
	
	//Meandering off in tRIBS
	//-----------------------
	//if( cn->Meanders() ) p1 = cn->getNew2DCoords();
	
	cn = (tNode *) ct->pPtr(2);
	tArray< double > p2( cn->get2DCoords() );
	
	//Meandering off in tRIBS
	//-----------------------
	//if( cn->Meanders() ) p2 = cn->getNew2DCoords();
	
	if( PointsCCW( p0, p1, p2 ) ) return 1;
	else return 0;
}


/***************************************************************************
**  
**  globalFns::InNewTri( )
** 
**  Determines whether coords are in "new" triangle
**  Meandering turned off
** 
***************************************************************************/
int InNewTri( tArray< double > &xy, tTriangle *ct )
{
	int j;
	tNode *vtx;
	tArray< double > xy1, xy2;
	for( j=0; j<3; j++ ){
		vtx = (tNode *) ct->pPtr(j);
		
		//Meandering off in tRIBS
		//-----------------------
		//if( vtx->Meanders() ) xy1 = vtx->getNew2DCoords();
		//else xy1 = vtx->get2DCoords();
		
		xy1 = vtx->get2DCoords();
		
		vtx = (tNode *) ct->pPtr( (j+1)%3 );
		
		//Meandering off in tRIBS
		//-----------------------
		//if( vtx->Meanders() ) xy2 = vtx->getNew2DCoords();
		//else xy2 = vtx->get2DCoords();
		
		xy2 = vtx->get2DCoords();
		
		if ( ( (xy1[1] - xy[1]) * (xy2[0] - xy[0]) ) >
			 ( (xy1[0] - xy[0]) * (xy2[1] - xy[1])) )
			break;
	}
	if( j == 3) return 1;
	else return 0;
}


/***************************************************************************
**  
**  globalFns::Intersect( tEdge * ae, tEdge * be )
** 
**  Edge Intersection routine. Meandering turned off.
**
***************************************************************************/
int Intersect( tEdge * ae, tEdge * be )
{
	tNode * lnode;
	
	if( !ae || !be )
	{
		cout<<"Intersect: Warning: invalid edge(s)"<<endl<<flush;
		return( 0 );
	}
	
	if( !ae->getOriginPtr() || !ae->getDestinationPtr() ||
		!be->getOriginPtr() || !be->getOriginPtr() ){
		cout<<"Intersect: Warning: invalid org or dest"<<endl<<flush;
		return( 0 );
	}
	
	lnode = (tNode *) ae->getOriginPtrNC();
	tArray< double > A( lnode->get2DCoords() );
	
	//Meandering off in tRIBS
	//-----------------------
	//if( lnode->Meanders() ) A = lnode->getNew2DCoords();
	
	lnode = (tNode *) ae->getDestinationPtrNC();
	tArray< double > B( lnode->get2DCoords() );
	
	//Meandering off in tRIBS
	//-----------------------
	//if( lnode->Meanders() ) B = lnode->getNew2DCoords();
	
	lnode = (tNode *) be->getOriginPtrNC();
	tArray< double > C( lnode->get2DCoords() );
	
	//Meandering off in tRIBS
	//-----------------------
	//if( lnode->Meanders() ) C = lnode->getNew2DCoords();
	
	lnode = (tNode *) be->getDestinationPtrNC();
	tArray< double > D( lnode->get2DCoords() );
	
	//Meandering off in tRIBS
	//-----------------------
	//if( lnode->Meanders() ) D = lnode->getNew2DCoords();
	
	double rnumsign;
	double snumsign;
	
	double denomsign =
		predicate.DifferenceOfProductsOfDifferences( B[0], A[0],
													 D[1], C[1],
													 B[1], A[1],
													 D[0], C[0] );
	if( denomsign == 0.0 ) // segments parallel
	{
		rnumsign =
		predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
													 D[0], C[0],
													 A[0], C[0],
													 D[1], C[1] );
		
		
		// The following needs to be re-checked. SMR commented out max/min
		// statements since macros did not work on SGI compiler. This routine
		// is therefore not currently working.
		
		if( rnumsign != 0.0 ) return 0;
		else{
			if( A[0] != B[0] ) //segments not vertical
			{
				return 1; // segments overlap
			}
			else // vertical
			{
				return 1; // segments overlap
			}
		}
		
	}
	else if( denomsign > 0.0 ) // not parallel
	{
		if( A[0] == C[0] && A[1] == C[1] ) return 0; // segments
		if( A[0] == D[0] && A[1] == D[1] ) return 0; // intersect 
		if( B[0] == C[0] && B[1] == C[1] ) return 0; // at
		if( B[0] == D[0] && B[1] == D[1] ) return 0; // endpoints
		
		rnumsign =
			predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
                                                         D[0], C[0],
                                                         A[0], C[0],
                                                         D[1], C[1] );
		if( rnumsign < 0.0 ) return 0; // r < 0
		snumsign =
			predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
                                                         B[0], A[0],
                                                         A[0], C[0],
                                                         B[1], A[1] );
		if( snumsign < 0.0 ) return 0; // s < 0
	}
	
	else{
		if( A[0] == C[0] && A[1] == C[1] ) return 0; // segments
		if( A[0] == D[0] && A[1] == D[1] ) return 0; // intersect 
		if( B[0] == C[0] && B[1] == C[1] ) return 0; // at
		if( B[0] == D[0] && B[1] == D[1] ) return 0; // endpoints
		
		rnumsign =
			predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
                                                         D[0], C[0],
                                                         A[0], C[0],
                                                         D[1], C[1] );
		if( rnumsign > 0.0 ) return 0; // r < 0
		snumsign =
			predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
                                                         B[0], A[0],
                                                         A[0], C[0],
                                                         B[1], A[1] );
		if( snumsign > 0.0 ) return 0; // s < 0
	}
	
	// reverse directions of segments so we can still detect
	// intersection by the signs of the numerators and denomenators:
	// A -> B, B->A, C->D, D->C
	
	denomsign = 
		predicate.DifferenceOfProductsOfDifferences( A[0], B[0],
													 C[1], D[1],
													 A[1], B[1],
													 C[0], D[0] );
	if( denomsign > 0.0 )
	{
		rnumsign =
		predicate.DifferenceOfProductsOfDifferences( B[1], D[1],
													 C[0], D[0],
													 B[0], D[0],
													 C[1], D[1] );
		if( rnumsign < 0.0 ) return 0; // r > 1
		snumsign =
			predicate.DifferenceOfProductsOfDifferences( B[1], D[1],
                                                         A[0], B[0],
                                                         B[0], D[0],
                                                         A[1], B[1] );
		if( snumsign < 0.0 ) return 0; // s > 1
	}
	else if( denomsign < 0.0 )
	{
		rnumsign =
		predicate.DifferenceOfProductsOfDifferences( B[1], D[1],
													 C[0], D[0],
													 B[0], D[0],
													 C[1], D[1] );
		if( rnumsign > 0.0 ) return 0; // r > 1
		snumsign =
			predicate.DifferenceOfProductsOfDifferences( B[1], D[1],
                                                         A[0], B[0],
                                                         B[0], D[0],
                                                         A[1], B[1] );
		if( snumsign > 0.0 ) return 0; // s > 1
	}
	//else segments must intersect other than at endpoints:
	return 1;  
}



/***************************************************************************
**  
**  globalFns::IntersectsAnyEdgeInList
** 
**  Returns the first edge in the list (passed by pointer) which intersects 
**  "edge" or NULL if "edge" intersects no other edges.
**
***************************************************************************/
tEdge* IntersectsAnyEdgeInList( tEdge* edge, tPtrList< tEdge >& edglistRef )
{
	tEdge * ce;
	tPtrListIter< tEdge > edgIter( edglistRef );
	if( !edge ){
		cout<<"IntersectsAnyEdge: Warning: invalid edge"<<endl<<flush;
		return( NULL );
	}
	
	// Check all edges in list
	
	for( ce = edgIter.FirstP(); !(edgIter.AtEnd()); ce = edgIter.NextP() ){
		if( edge != ce ){
			if( Intersect( edge, ce ) ) return( ce );
		}     
	}
	assert( edgIter.AtEnd() );
	return( NULL );
}

/***************************************************************************
**  
**  globalFns::InterpSquareGrid( )
** 
**  Interpolation of Grid used in MakeRandomPointsFromArcGrid()
**  Use Tetzlaff & Harbaugh, '88, grid interpolation:
**  z(z1, z2, z3, z4, X, Y) = z1 * XY + z2 * (1-X)Y + z3 * X(1-Y) 
**  Assumes passed coordinates are normalized to grid spacing = 1.							
**
***************************************************************************/
double InterpSquareGrid( double xgen, double ygen, tMatrix< double >& elev,
                         int nodata )
{
	int nodatacount = 0;
	int nx = elev.getNumCols();
	int ny = elev.getNumRows();
	
	int ix = (int)xgen;                              // z3----z4
	int iy = (int)ygen;                              // |     |
	double xrem = xgen - (double)ix;                 // |     |
	double yrem = ygen - (double)iy;                 // z1----z2
	double yrows = elev.getNumRows();
	int i = ix;
	int j = (int)(yrows - 1.0 - ygen);
	
	double z3 = elev( j, i );
	double z4 = ( i < nx - 1 ) ? elev( j, i + 1 ) : nodata;
	double z1 = ( j < ny - 1 ) ? elev( j + 1, i ) : nodata;
	double z2 = ( i < nx - 1 && j < ny - 1 ) ? elev( j + 1, i + 1 ) : nodata;
	
	if( z1 == nodata ) nodatacount++;
	if( z2 == nodata ) nodatacount += 2;
	if( z3 == nodata ) nodatacount += 4;
	if( z4 == nodata ) nodatacount += 8;
	if( nodatacount == 15 ) return nodata;	     
	if( nodatacount == 0 )			      
		return z2 + (z4 - z2) * xrem + (z1 - z2) * yrem 
			- (z1 + z4 - z2 - z3) * xrem * yrem;
	if( nodatacount == 9 &&			     
		( ( xrem <= 0.5 && yrem >= 0.5 ) || ( xrem >= 0.5 && yrem <= 0.5 ) ) )
		return z2 + 0.5 * ( z3 - z2 ) * ( xrem + yrem ); 
	if( nodatacount == 6 &&			     
		( ( xrem <= 0.5 && yrem <= 0.5 ) || ( xrem >= 0.5 && yrem >= 0.5 ) ) )
		return 0.5 * ( z1 + z4 + ( z4 - z1 ) * ( xrem - yrem ) ); 
	if( nodatacount == 1 &&			     
		!( xrem < 0.5 && yrem < 0.5 ) )
		return z2 + ( z4 - z2 ) * xrem + 0.5 * ( z3 - z2 ) * yrem
			- ( z4 - 0.5 * ( z2 + z3 ) ) * xrem * yrem; 
	if( nodatacount == 2 &&			    
		!( xrem > 0.5 && yrem < 0.5 ) )
		return 0.5 * ( z1 + z4 + ( z4 - z1 ) * ( xrem - yrem ) )
			- ( 0.5 * ( z1 + z4 ) - z3 ) * xrem * yrem; 
	if( nodatacount == 4 &&			       
		!( xrem < 0.5 && yrem > 0.5 ) )
		return z2 + ( z4 - z2 ) * xrem + ( z1 - z2 ) * yrem
			- ( 0.5 * ( z1 + z4 ) - z2 ) * xrem * yrem; 
	if( nodatacount == 8 &&			       
		!( xrem > 0.5 && yrem > 0.5 ) )
		return z2 + 0.5 * ( z3 - z2 ) * xrem + ( z1 - z2 ) * yrem 
			- ( z1 - 0.5 * ( z2 + z3 ) ) * xrem * yrem;
	if( nodatacount == 3 && yrem >= 0.5 )       
		return z4 + (z3 - z4) * yrem;              
	if( nodatacount == 10 && xrem <= 0.5 )      
		return z1 + (z3 - z1) * xrem;             
	if( nodatacount == 12 && yrem <= 0.5 )       
		return z2 + (z1 - z2) * yrem;             
	if( nodatacount == 5 && xrem >= 0.5 )        
		return z2 + (z4 - z2) * xrem;             
	if( nodatacount == 14 && xrem <= 0.5 && yrem <= 0.5 ) 
		return z1;
	if( nodatacount == 13 && xrem >= 0.5 && yrem <= 0.5 ) 
		return z2;
	if( nodatacount == 11 && xrem <= 0.5 && yrem >= 0.5 ) 
		return z3;
	if( nodatacount == 7 && xrem >= 0.5 && yrem >= 0.5 )  
		return z4;
	return nodata;
}


/**********************************************************************
** 
**  globalFns::PlaneFit 
** 
**  A plane is fit given the x,y,z coordinates of three points.
**  Returns the z value on this plane at the x and y location
**  that is sent to it.
**
**  tx, ty are the location where you want to know z
**  p0, p1, p2 are arrays containing x and y values (p0[0]=x; p0[1]=y)
**  zs contain the z values at p0, p1, and p2, respectively in the array
**
**********************************************************************/
double PlaneFit(double x, double y, tArray<double> p0,
                tArray<double> p1, tArray<double> p2, tArray<double> zs)
{
	double a, b, c;
	double y0, y1, y2, x0, x1, x2, z0, z1, z2;
	y0=p0[1];
	y1=p1[1];
	y2=p2[1];
	x0=p0[0];
	x1=p1[0];
	x2=p2[0];
	z0=zs[0];
	z1=zs[1];
	z2=zs[2];
	
	a=(-y1*z2+z2*y0+z1*y2-y2*z0+z0*y1-y0*z1)/(y2*x1-x1*y0-x2*y1+y1*x0-x0*y2+y0*x2);
	b=-(x2*z1-z1*x0-z2*x1+x1*z0-z0*x2+x0*z2)/(y2*x1-x1*y0-x2*y1+y1*x0-x0*y2+y0*x2);
	c=(y2*x1*z0-z0*x2*y1+z2*y1*x0-y2*z1*x0+y0*x2*z1-z2*x1*y0)/(y2*x1-x1*y0-x2*y1+y1*x0-x0*y2+y0*x2);
	
	return(a*x+b*y+c);  
}

/**********************************************************************
** 
**  globalFns::LineFit 
**
**********************************************************************/
double LineFit(double x1, double y1, double x2, double y2, double nx)
{
	double slope = (y2-y1)/(x2-x1);
	return(slope*(nx-x1)+y1);
} 

/**********************************************************************
** 
**  globalFns::DistanceBW2Points
**
**********************************************************************/
double DistanceBW2Points(double x1, double y1, double x2, double y2)
{
	return(pow(pow(x1-x2,2)+pow(y1-y2,2),0.5));
}

//=========================================================================
//
//
//                        End of globalFns.cpp
//
//
//=========================================================================
