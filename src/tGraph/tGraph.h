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
**  tGraph.h: Header for tConnectivity class and objects
**
**  tGraph Class used in tRIBS for the parallel version, providing 
**  information about upstream and downstream basins
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: tGraph Include and Define Statements
//
//
//=========================================================================

#ifndef TGRAPH_H
#define TGRAPH_H

#include <iostream>
#include <vector>
#include <set>

#include "src/tSimulator/tSimul.h"
#include "src/tCNode/tCNode.h"
#include "src/tMesh/tMesh.h"
#include "src/tFlowNet/tKinemat.h"
#include "src/tFlowNet/tReservoir.h" // JECR2015
#include "src/tFlowNet/tResData.h" // JECR2015
#include "src/tGraph/tGraphNode.h"
#include "src/tInOut/tInputFile.h"

struct IDOrder {
  bool operator() ( const tCNode* n1, const tCNode* n2 )
  const {
    return n1->getID() < n2->getID();
  }
};

//=========================================================================
//
//
//                  Section 2: tGraph Class Definitions
//
//
//=========================================================================

class tGraph {

public:
  /// Constructor
  tGraph();
  /// Destructor
  ~tGraph();

  /// Initialize
  static void initialize(SimulationControl* s, tMesh<tCNode>* m, tKinemat* f,
    tInputFile& InputFile);

  /// Initialize for MeshBuilder input where reach partition happens first
  static void initialize(SimulationControl* s, tMesh<tCNode>* m,
    tInputFile& InputFile);

  /// Set the tKinemat separately for MeshBuilder because it has not been
  /// constructed at the point that initialize has been called
  static void setFlowNet(tKinemat* f)	{ flow = f; }

  /// Finalize
  static void finalize();

  /// Build graph, partition and mesh from MeshBuilder files (option 9)
  static void ReadDirectoryFromMeshBuilder();
  static void ReadFlowMesh();
  static void ReadFlowNode(fstream&, tCNode*);
  static void ReadFlowNode(fstream&, int* id, int* edge, int* flow, int* node);
  static void ReadFlowEdge(fstream&, tEdge*, int* orig, int* dest);
  static void ReadFlowEdge(fstream&, int* id, int* ccw);

  /// Determine stream reach connectivity
  static void connectivity();
  /// Write out stream reach connectivity
  static void outputConnectivity(tInputFile& InputFile);

  /// Partition graph based on # of processors
  static void partition(tInputFile& InputFile);
  /// Partition graph into n pieces
  static void partition(int n, tInputFile& InputFile);
  /// Read reach-based partitions from file
  static void readReachPartitionFromFile(char* pfile);
  /// Read inlet/outlet-based partitions from file
  static void readInletOutletPartitionFromFile(char* pfile);
  /// Create default partitions
  static void createDefaultPartition(int np);

  /// List ids of all active nodes
  static void listActiveNodes();

  /// Is stream reach in this partition?
  static bool inLocalPartition(int r);
  /// Get partition of selected stream reach
  static int getPartition(int r);
  /// Return list of partition members
  static std::vector<int> getPartitionMembers(int p);

  /// Update Mesh/Flow per partition
  /// Updates node list in Mesh
  /// Updates stream reach list in Flow
  static void update();

  /// Calculate overlapping nodes for flux exchange
  static void calculateOverlap();

  /// Calculate runoff/runon nodes
  static void calculateRunFlux();

  /// Display node counts per reach
  static void reachNodeCounts();
  /// Display stream nodes in each reach
  static void listStreamNodes();

  /// Does the local partition have the last reach
  static bool hasLastReach();
  /// Is this the last reach in the local partition
  /// with a common destination reach in a remote
  /// partition?
  static bool lastLocalReachWithCommonDest(int r, int rdest);
  /// Is this the last reach in a remote partition
  /// with a common destination reach in the local
  /// partition?
  static bool lastRemoteReachWithCommonDest(int r, int rdest);

  /// Is first reach upstream of second reach?
  static bool isUpstreamOf(int r1, int r2);
  /// Is first reach downstream of second reach?
  static bool isDownstreamOf(int r1, int r2);

  /// Return reach that node is in 
  static int inWhichReach(tCNode* cn);
  /// Return reach with given head and outlet ids
  static int whichReach(int iHead, int iOutlet);

  /// Is this node an overlap node?
  static bool isOverlapNode(tCNode* cn);

  /// Is this node a downstream node?
  static bool isDownstreamNode(tCNode* cn);
  /// Is this node an upstream node?
  static bool isUpstreamNode(tCNode* cn);

  /// Is this a remote flux node?
  static bool isRemotefluxNode(tCNode* cn);
  /// Is this a local flux node?
  static bool isLocalfluxNode(tCNode* cn);

  /// Are there upstream stream reach(s)?
  static bool hasUpstream(int r);
  /// Are there downstream stream reach(s)?
  static bool hasDownstream(int r);

  /// Send initial data to overlapping nodes
  static void sendInitial();
  /// Receive initial data from overlapping nodes
  static void receiveInitial();

  /// Send point data to downstream stream reach(s)
  static void sendDownstream(int rid, tCNode* snode, double value);
  /// Receive point data from upstream stream reach(s)
  static void receiveUpstream(int rid, tCNode* rnode);

  /// Send point data to downstream stream reach(s)
  static void sendQpin(int rid, tCNode* snode, double value);
  /// Receive point data from upstream stream reach(s)
  static void receiveQpin(int rid, tCNode* rnode);

  /// Send groundwater
  static void sendGroundWater();
  static void receiveGroundWater();

  /// Send and receive overlap after unsaturated zone calculation
  static void sendNwt();
  static void receiveNwt();

  /// Send data to overlapping nodes
  static void sendOverlap();
  /// Receive data from overlapping nodes
  static void receiveOverlap();
  /// Reset data in overlapping nodes
  static void resetOverlap();

  /// Send data to upstream overlapping flow nodes
  static void sendUpstreamFlow();
  /// Receive data from downstream overlapping flow nodes
  static void receiveDownstreamFlow();

  /// Send runon flux data to downstream nodes
  static void sendRunFlux(tCNode* cn);
  /// Receive runon flux data from upstream nodes
  static void receiveRunFlux(tCNode* cn);

private:

  static SimulationControl*      sim;             //!< SimulationControl
  static tMesh<tCNode>*          mesh;            //!< Mesh
  static tKinemat*               flow;            //!< Kinematic flow

  static int                     numGlobalReach;  //!< # of stream reaches
  static int                     numGlobalPart;   //!< # of partitions
  static int                     localPart;       //!< Local partition

  static std::vector<tGraphNode> conn;            //!< Connectivity of reaches

  static std::vector<int>        reach2partition; //!< Reach id to partition #
  static std::vector<int>        localReach;      //!< List of local reaches 
  static std::vector<int>        pointsPerReach;  //!< # of points per reach

  static int*                    hid;             //!< Reach head node IDs
  static int*                    oid;             //!< Reach outlet node IDs
  static int*                    aboveid;         //!< node above outlet IDs

  static std::vector<tCNode*>    nodeAboveOutlet; //!< Node above outlet
  static std::set<tCNode*,IDOrder>* upFlow;      //!< Upstream flow nodes
  static std::set<tCNode*,IDOrder>* downFlow;    //!< Downstream flow nodes

  static std::set<tCNode*,IDOrder>* localFlux;   //!< Send overlap flux nodes
  static std::set<tCNode*,IDOrder>* remoteFlux;  //!< Receive overlap flux nodes
  static bool                    lastReach;      //!< Contains last reach

  // MeshBuilder required variables
  static int  numGlobalNodes;           //!< # of nodes in problem
  static int  numGlobalEdges;           //!< # of edges in problem
  static int  nodeBytes;                //!< Size of MeshBuilder node
  static int  edgeBytes;                //!< Size of MeshBuilder edge

  static int* nodesPerReach;            //!< Nodes in each reach
  static int* internalEdgesPerReach;    //!< Edges internal to each reach
  static int* externalEdgesPerReach;    //!< Edges external to each reach

  static int* fluxNodesPerReach;        //!< Flux external nodes in each reach
  static int* fluxEdgesPerReach;	//!< Flux external edges to each reach

  static int* nodeOffset;               //!< Offset within nodes file
  static int* edgeOffset;               //!< Offset within edges file

  static int* fluxNodeOffset;           //!< Offset within flux node file
  static int* fluxEdgeOffset;		//!< Offset within flux edge file
};

#endif

//=========================================================================
//
//
//                          End of tGraph.h 
//
//
//=========================================================================
