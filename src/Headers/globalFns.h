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
**  globalFns.h: Global tRIBS Class Header
**
***************************************************************************/

#ifndef GLOBALFNS_H
#define GLOBALFNS_H

//=========================================================================
//
//
//                  Section 1: globalFns Include Statements
//
//
//=========================================================================

#include "src/tArray/tArray.h"
#include "src/tMeshElements/meshElements.h"
#include "src/tPtrList/tPtrList.h"
#include "src/Mathutil/predicates.h"
#include "src/tArray/tMatrix.h"

#include <iostream>
#include <cmath>

using namespace std;

//=========================================================================
//
//
//                  Section 2: Function Declarations
//
//
//=========================================================================

extern Predicates predicate; 

double ran3( long * ); 

tArray< double > UnitVector( tEdge* );

double FindCosineAngle0_2_1( tArray< double > &, tArray< double > &,
                             tArray< double > & );

int TriPasses( tArray< double > &, tArray< double > &,
               tArray< double > &, tArray< double > & );

int PointsCCW( tArray< double > &, tArray< double > &, tArray< double > & );

int NewTriCCW( tTriangle * );

int InNewTri( tArray< double > &, tTriangle * );

int Intersect( tEdge *, tEdge * );

tEdge* IntersectsAnyEdgeInList( tEdge*, tPtrList< tEdge >& );

double InterpSquareGrid( double, double, tMatrix< double >&, int );

tArray< double > FindIntersectionCoords( tArray< double >,
                                         tArray< double >,
                                         tArray< double >,
                                         tArray< double > );

double PlaneFit(double x, double y, tArray<double> p0,
                tArray<double> p1, tArray<double> p2, tArray<double> zs);

double LineFit(double x1, double y1, double x2, double y2, double nx);

double DistanceBW2Points(double x1, double y1, double x2, double y2 );

//=========================================================================
//
//
//                  Section 3: Templated Function Declarations
//
//
//=========================================================================

template< class tSubNode >
int Next3Delaunay( tPtrList< tSubNode > &, tPtrListIter< tSubNode > & );

template< class tSubNode >
int PointAndNext2Delaunay( tSubNode &, tPtrList< tSubNode > &,
                           tPtrListIter< tSubNode > & );
template< class T > 
ostream &operator<<( ostream &, const tArray< T > & );

/**********************************************************************
**
**  Templated function BinaryWrite
**
**********************************************************************/

template< class outDataType >
inline void BinaryWrite(ostream& outStream, const outDataType& outData)
{
  outStream.write(
    reinterpret_cast<const char*>(&outData), sizeof(outDataType));
}

/**********************************************************************
**
**  Templated function BinaryRead
**
**********************************************************************/

template< class inHolderType >
inline istream& BinaryRead(istream& inStream, inHolderType& inHolder)
{
   return inStream.read(
      reinterpret_cast<char*>(&inHolder), sizeof(inHolderType));
}

#endif

//=========================================================================
//
//
//                          End of globalFns.h 
//
//
//=========================================================================
