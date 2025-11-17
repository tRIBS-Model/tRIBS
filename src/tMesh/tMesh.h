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
**  tMesh.h: Header file for class tMesh
**
**  tMesh is the master class that handles the implementation of the
**  triangulated mesh. The class includes lists of the mesh elements
**  (nodes, triangles, and edges; see meshElements.h/.cpp), and provides
**  functionality to:
**    - read in or create meshes, either from scratch, from a list of
**      points, from a pre-existing set of triangulation files (e.g., a
**      previous run), or an Arc/Info files
**    - move, add, and/or delete nodes
**    - update Delaunay and Voronoi geometry
**
\***************************************************************************/

#ifndef TMESH_H
#define TMESH_H

#include "src/Headers/Inclusions.h"
#include "src/tMesh/tTriangulator.h"

#ifdef ALPHA_64
  #include <stdlib.h>
  #include <strings.h>
#elif defined LINUX_32
  #include <stdlib.h>
  #include <strings.h>
  #include <iostream>
#elif defined WIN
  #ifdef _WIN32
  #define srand48(x) srand(x)
  #define drand48() double(rand())/RAND_MAX
  #define strcasecmp(x,y) _stricmp(x,y)
  #endif
#else 
  #include <stdlib.h>
  #include <strings.h>
#endif

using namespace std;

//=========================================================================
//
//
//                  Section 1: tMesh Class Declarations
//
//
//=========================================================================

template< class tSubNode >
class tMesh
{
  tMesh(const tMesh&);
  tMesh& operator=(const tMesh&);
public:
   tMesh();
   tMesh( SimulationControl* );
   tMesh( SimulationControl*, tInputFile & );
   tMesh( tMesh * );
   ~tMesh();
 
   SimulationControl *simCtrl;    // Pointer to simulation control 

   void MakeMeshFromScratch( tInputFile & );   		// creates a new mesh
   void MakeMeshFromInputData( tInputFile & ); 		// reads in existing mesh
   void MakeMeshFromPoints( tInputFile & );    		// creates mesh from pts
   void MakeRandomPointsFromArcGrid( tInputFile & ); 	// mesh from arc (rand)
   void MakeHexMeshFromArcGrid( tInputFile & );	        // mesh from arc (hex)
   void MakeLayersFromInputData( tInputFile & );
   void MakePointFromFileArcInfo( tInputFile & );       //creates pts from .net 
   void MakePointFromFileArcInfoGen( tInputFile & );    //creates pts from .pnt and .line
   void MakeMeshFromTriangulator( tInputFile & );    	//creates mesh using Tipper triang
   void MakeMeshFromMeshBuilder( tInputFile & );	// creates from .meshb
   void Print();

   void MakeCCWEdges();
   void setVoronoiVertices();
   void CalcVoronoiEdgeLengths();
   void CalcVAreas();
   tTriangle *LocateTriangle( double, double );
   tTriangle *LocateNewTriangle( double, double );
   tTriangle *TriWithEdgePtr( tEdge * );
   int DeleteNode( tListNode< tSubNode > *, int repairFlag=1 );
   int DeleteNode( tSubNode *, int repairFlag=1 );
   int ExtricateNode( tSubNode *, tPtrList< tSubNode > & );
   int DeleteEdge( tEdge * );
   int ExtricateEdge( tEdge * );
   int DeleteTriangle( tTriangle * );
   int ExtricateTriangle( tTriangle * );
   int RepairMesh( tPtrList< tSubNode > & );
   int AddEdgeAndMakeTriangle( tPtrList< tSubNode > &,
                               tPtrListIter< tSubNode > & );
   int MakeTriangle( tPtrList< tSubNode > &,
                     tPtrListIter< tSubNode > & );
   int AddEdge( tSubNode *, tSubNode *, tSubNode * );
   tSubNode *AddNode( tSubNode &, int updatemesh = 0, double time = 0.0 );
   tSubNode *AddNodeAt( tArray< double > &, double time = 0.0 );
   void AddNodesAround( tSubNode *, double time=0.0 ); 
   tMeshList<tEdge> * getEdgeList()
   { return &edgeList; }
   tMeshList<tSubNode> * getNodeList()
   { return &nodeList; }
   tMeshList<tSubNode> * getUnsortList();  
   tList< tTriangle > * getTriList();  
   tEdge *getEdgeComplement( tEdge * );

   // MeshBuilder Option 9 requires access to NodeTable
   tSubNode* getNodeFromID(int id) { return NodeTable[id]; }
   void setNodeFromID(int id, tSubNode* node) { NodeTable[id] = node; }
   void allocateNodeTable(int sz) { NodeTable = new tSubNode*[sz]; }
   void deleteNodeTable() { delete [] NodeTable; }
  
   void CheckMeshConsistency( int boundaryCheckFlag=1 );
   void CheckMeshConsistency(  tInputFile &, int boundaryCheckFlag=1 );
   int ChangePointOrder(tInputFile &, tList<double> XY );
   void UpdateMesh();
   int CheckForFlip( tTriangle *, int, int );
   void FlipEdge( tTriangle *, tTriangle *, int, int );
   tEdge * IntersectsAnyEdge( tEdge * );
   void CheckTriEdgeIntersect();
   void CheckLocallyDelaunay();
   void MoveNodes( double time = 0.0 );
   void TellAboutNode(tSubNode *);
   void writeRestart(fstream &);
   void readRestart(fstream &);
  
#ifndef NDEBUG
   void DumpEdges();
   void DumpSpokes( tSubNode * );
   void DumpTriangles();
   void DumpNodes();
#endif
   
protected:
   int nnodes, nedges, ntri;       	// # of nodes, edges, and tri's
   tMeshList< tSubNode > nodeList; 	// list of nodes
   tMeshList< tSubNode > unsortList;  	// list of unsorted nodes 
   tMeshList< tEdge > edgeList;       	// list of directed edges
   tList< tTriangle > triList;     	// list of triangles  (revert from tPtrList)
   int miNextNodeID;               	// next ID for added triangle
   int miNextEdgID;                	// next ID for added triangle
   int miNextTriID;                	// next ID for added triangle
   long seed;                      	// random seed
   tTriangle* mSearchOriginTriPtr; 	// ptr to tri. from which to start searches

   tSubNode** NodeTable;		// lookup table for node pointers
					// must be available to tFlowNet and
					// tGraph if using MeshBuilder files
   
};

#endif

//=========================================================================
//
//
//                          End of tMesh.h
//
//
//=========================================================================
