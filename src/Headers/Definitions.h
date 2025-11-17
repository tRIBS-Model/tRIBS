/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) 2025. tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

//=========================================================================
//
//
//                  Section 1: Definitions and Macros
//
//
//=========================================================================

// Definitions

#define VERSION "tRIBS 2.0"
#define TRUE 1
#define FALSE 0
#define kMaxNameLength 200
#define kMaxNameSize 200
#define kName 200
#define kCount 200000
#define sCount 20
#define kMaxExt 20
#define kCommentMark '#'
#define kTimeLineMark ' '
#define kUniformGrid 0    
#define kPerturbedGrid 1
#define kRandomGrid 2
#define kCornerOutlet 0   
#define kOpenSide 1
#define kOppositeSidesOpen 2
#define kAllSidesOpen 3
#define kSpecifyOutlet 4
#define kClosedBoundary 1
#define kOpenBoundary 2
#define kNonBoundary 0
#define kStream 3             
#define kFlowAllowed 1
#define kFlowNotAllowed 0
#define kRepairMesh 1
#define kNoRepair 0
#define kNoUpdateMesh 0
#define kUniformMesh 0
#define kPerturbedMesh 1
#define kRandomMesh 2
#define PI 3.1415926
#define EPS 2.2204e-16
#define THRESH  1e-6
#define kFlooded 1  	
#define kNotFlooded 0  	
#define kCurrentLake 2  	
#define kSink 3  	
#define kOutletFlag 4  	
#define kOutletPreFlag 5  	
#define kVeryHigh 100000  	

// Macros
#define ROUND(x) (int)(x+0.5)
#define SIGN(x)  ( x>0 ? 1 : 0 )

#endif

//=========================================================================
//
//
//                         End of Definitions.h
//
//
//=========================================================================
