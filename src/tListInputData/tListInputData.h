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
**  tListInputData.h: 	Header file for class tListInputData
**
**  This class is used to read in lists of Delaunay-triangulated
**  mesh elements from three user-provided input files, which contain
**  the nodes, directed edges, and triangles in the mesh, respectively.
**
**************************************************************************/

#ifndef TLISTINPUTDATA_H
#define TLISTINPUTDATA_H

#include "src/Headers/Inclusions.h"
#include "src/Headers/Definitions.h"

#include <iostream>
#include <fstream>
using namespace std;

//=========================================================================
//
//
//                  Section 1: tListInputData Class Declaration
//
//
//=========================================================================

/**************************************************************************
**  
**  class tListInputData()
**
**  tListInputData reads an established triangulation from a set of four
**  files, and stores the data in a series of arrays. The files are:
**    <name>.nodes  --  node (point) data
**    <name>.edges  --  directed edge data
**    <name>.tri    --  triangle data
**    <name>.z      --  "z" value data (elevation or other)
**
**************************************************************************/

template< class tSubNode >
class tListInputData{
    friend class tMesh< tSubNode >; 

public:
    tListInputData( tInputFile & ); 
    ~tListInputData();

private:
    void GetKeyEntry();         // not currently supported
    void GetFileEntry();        // read data from files
    int nnodes, nedges, ntri;  	// # nodes, edges, & triangles
    ifstream nodeinfile;   	// node input file
    ifstream edgeinfile;   	// edge input file
    ifstream triinfile;    	// triangle input file
    ifstream zinfile;      	// "z" input file
    tArray< double > x;      	// node x coords
    tArray< double > y;      	// node y coords
    tArray< double > z;      	// node z values
    tArray< int > edgid;     	// node edge ID #s
    tArray< int > boundflag; 	// node boundary codes
    tArray< int > orgid;     	// directed edge origin node ID #s
    tArray< int > destid;    	// directed edge destination node ID #s
    tArray< int > nextid;    	// ID #s of next counter-clockwise edges
    tArray< int > p0;     	// IDs of triangle node 0
    tArray< int > p1;     	// IDs of triangle node 1
    tArray< int > p2;     	// IDs of triangle node 2
    tArray< int > e0;     	// IDs triangle clockwise-oriented edge 0
    tArray< int > e1;     	// IDs triangle clockwise-oriented edge 1
    tArray< int > e2;   	// IDs triangle clockwise-oriented edge 2
    tArray< int > t0;     	// IDs of neighboring tri's opposite node 0
    tArray< int > t1;     	// IDs of neighboring tri's opposite node 1
    tArray< int > t2;     	// IDs of neighboring tri's opposite node 2
};

#endif

//=========================================================================
//
//
//                      End of tListInputData.h
//
//
//=========================================================================

