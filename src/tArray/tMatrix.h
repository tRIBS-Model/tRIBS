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
**  tMatrix.h: Header file for class tMatrix.
**
**  Class tMatrix implements matrices (2D arrays) as 1D arrays of tArray
**  objects. 
**
***************************************************************************/

#ifndef TMATRIX_H
#define TMATRIX_H

#include "src/tArray/tArray.h"

//=========================================================================
//
//
//                  Section 1: tMatrix Class Definition
//
//
//=========================================================================

template < class T >
class tMatrix
{
public:
    tMatrix();
    tMatrix( int nr, int nc );
    ~tMatrix();
    T &operator()( int row, int col );
   int getNumRows() {return nrows;}
   int getNumCols() {return ncols;}
    
private:
    tArray<T> *ptr;
    int nrows;
    int ncols;  
};

#endif

//=========================================================================
//
//
//                           End of tMatrix.h
//
//
//=========================================================================
