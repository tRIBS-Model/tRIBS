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
**  tGraph.cpp: Functions for class tGraph (see tGraph.h)
**
***************************************************************************/

#include "src/tGraph/tGraph.h"
#include "src/Headers/globalIO.h"
#include "src/Headers/Definitions.h"
#include "src/tMeshList/tMeshList.h"

#ifdef PARALLEL_TRIBS
#include "src/tParallel/tParallel.h"
#endif

#include <cassert>
#include <map>

SimulationControl* tGraph::sim = 0;
tMesh<tCNode>* tGraph::mesh = 0;
tKinemat* tGraph::flow = 0;

int tGraph::numGlobalReach = 0;
int tGraph::numGlobalPart = 1;
int tGraph::localPart = 0;

std::vector<tGraphNode> tGraph::conn;

std::vector<int> tGraph::reach2partition;
std::vector<int> tGraph::localReach;
std::vector<int> tGraph::pointsPerReach;

int* tGraph::hid = 0;;
int* tGraph::oid = 0;
int* tGraph::aboveid = 0;

std::vector<tCNode*> tGraph::nodeAboveOutlet;
std::set<tCNode*,IDOrder>* tGraph::upFlow = 0;
std::set<tCNode*,IDOrder>* tGraph::downFlow = 0;
std::set<tCNode*,IDOrder>* tGraph::localFlux = 0;
std::set<tCNode*,IDOrder>* tGraph::remoteFlux = 0;

bool tGraph::lastReach = false;

// Message tags
const int INITIAL     = 1000;
const int DOWNSTREAM  = 3000;
const int UPSTREAM    = 4000;
const int OVERLAP     = 5000;
const int RUNON       = 6000;
const int QPIN        = 7000;
const int GROUNDWATER = 9000;
const int NWT         = 10000;

// MeshBuilder variables
int tGraph::nodeBytes = 0;
int tGraph::edgeBytes = 0;
int tGraph::numGlobalNodes = 0;
int tGraph::numGlobalEdges = 0;

int* tGraph::nodesPerReach = 0;
int* tGraph::internalEdgesPerReach = 0;
int* tGraph::externalEdgesPerReach = 0;
int* tGraph::fluxNodesPerReach = 0;
int* tGraph::fluxEdgesPerReach = 0;

int* tGraph::nodeOffset = 0;
int* tGraph::edgeOffset = 0;
int* tGraph::fluxNodeOffset = 0;
int* tGraph::fluxEdgeOffset = 0;

tGraph::tGraph() {}
tGraph::~tGraph() {}

/*************************************************************************
**
** Finalize
**
*************************************************************************/

void tGraph::finalize(){
  sim = nullptr;
  mesh = nullptr;
  flow = nullptr;
  numGlobalReach = 0;
  
  conn.erase(conn.begin(), conn.end());
  reach2partition.erase(reach2partition.begin(), reach2partition.end());

  localPart = -1;
  localReach.erase(localReach.begin(), localReach.end());

  pointsPerReach.erase(pointsPerReach.begin(), pointsPerReach.end());

  nodeAboveOutlet.erase(nodeAboveOutlet.begin(), nodeAboveOutlet.end());

  if (hid != NULL) delete [] hid;
  if (oid != NULL) delete [] oid;
  if (aboveid != NULL) delete [] aboveid;

  for (int i = 0; i < numGlobalPart; i++) {
    upFlow[i].erase(upFlow[i].begin(), upFlow[i].end());
    downFlow[i].erase(downFlow[i].begin(), downFlow[i].end());
    localFlux[i].erase(localFlux[i].begin(), localFlux[i].end());
    remoteFlux[i].erase(remoteFlux[i].begin(), remoteFlux[i].end());
  }
  delete [] upFlow;
  delete [] downFlow;
  delete [] localFlux;
  delete [] remoteFlux;
 
  numGlobalPart = 0;
  lastReach = false;
} 

/*************************************************************************
**
** Initialize
**
*************************************************************************/

void tGraph::initialize(SimulationControl* s, tMesh<tCNode>* m, tKinemat* f,
  tInputFile& InputFile) {

  sim = s;
  mesh = m;
  flow = f;

  // Create stream connectivity table
  Cout << "\nCreating stream reach connectivity table..." << endl;
  connectivity();

  // Partition stream reach graph
  Cout << "\nPartitioning stream reach graph..." << endl;
  partition(InputFile);

  // Update stream reach and node list information
  Cout << "\nUpdating stream reach and node list..." << endl;
  update();

  // Write connectivity to file
#ifdef PARALLEL_TRIBS
  if (tParallel::isMaster()) {
#endif
    Cout << "\nWriting connectivity file..." << endl;
    outputConnectivity(InputFile);
#ifdef PARALLEL_TRIBS
  }
#endif
}

/*************************************************************************
**
** Initialize for MeshBuilder input where no tMesh exists when this is called
** and tGraph partitions and reads only information needed for local reaches
**
*************************************************************************/

void tGraph::initialize(SimulationControl* s, tMesh<tCNode>* m,
  tInputFile& InputFile)
{
  sim = s;
  mesh = m;

  // Create stream connectivity table
  Cout << "\nRead reach directory..." << endl;
  ReadDirectoryFromMeshBuilder();

  // Partition stream reach graph
  Cout << "\nPartitioning stream reach graph..." << endl;
  partition(InputFile);

  // Update stream reach and node list information
  Cout << "\nRead partitioned reach and node list..." << endl;
  ReadFlowMesh();

  // Calculate overlapping nodes for flux exchange
  Cout << "\nCalculate overlap ..." << endl;
  calculateOverlap();

  // Set the node above outlet flux information (replaces calculateRunFlux())
  // Information came from meshBuilder directory
  Cout << "\nSet node above run flux ..." << endl;
  for (int i = 0; i < numGlobalReach; i++)
      nodeAboveOutlet.push_back(NULL); 

  for (int reach = 0; reach < numGlobalReach; reach++) {
    for (int i = 0; i < localReach.size(); i++) {

      if (hid[localReach[i]] == oid[reach] &&
          reach2partition[reach] != localPart) {
            nodeAboveOutlet[reach] = mesh->getNodeFromID(aboveid[reach]);
      }
    }
  }

  // Write connectivity to file
  Cout << "\nWrite connectivity..." << endl;
  outputConnectivity(InputFile);
}

/*************************************************************************
**
** Graph information calculated by MeshBuilder is read in
** Includes counts of nodes and edges for each reach and offsets within file
**
*************************************************************************/

void tGraph::ReadDirectoryFromMeshBuilder()
{
  fstream reachStr("reach.meshb", ios::in);
  std::string dummy;
  int reach;

  reachStr >> numGlobalReach            // Number of reaches 
           >> numGlobalNodes            // Total number of nodes
           >> nodeBytes                 // Size of node information in bytes
           >> numGlobalEdges            // Total number of edges
           >> edgeBytes;                // Size of edge information in bytes

  // Read reach node and edge counts (put boundary nodes in extra reach)
  nodesPerReach = new int[numGlobalReach + 1];
  nodeOffset = new int[numGlobalReach + 1];

  fluxNodesPerReach = new int[numGlobalReach + 1];
  fluxNodeOffset = new int[numGlobalReach + 1];

  internalEdgesPerReach = new int[numGlobalReach + 1];
  externalEdgesPerReach = new int[numGlobalReach + 1];
  edgeOffset = new int[numGlobalReach + 1];

  fluxEdgesPerReach = new int[numGlobalReach + 1];
  fluxEdgeOffset = new int[numGlobalReach + 1];

  for (int i = 0; i < numGlobalReach + 1; i++) {
    reachStr >> reach
             >> nodesPerReach[i]
             >> nodeOffset[i]

             >> fluxNodesPerReach[i]
             >> fluxNodeOffset[i]

             >> internalEdgesPerReach[i]
             >> externalEdgesPerReach[i]
             >> edgeOffset[i]

             >> fluxEdgesPerReach[i]
             >> fluxEdgeOffset[i];

    pointsPerReach.push_back(nodesPerReach[i]);
  }

  int size, upsize, downsize, value, id;
  tCNode* curnode;

  // Read reach head node IDs and reach outlet node IDs
  hid = new int[numGlobalReach];
  oid = new int[numGlobalReach];
  aboveid = new int[numGlobalReach];

  for (int i = 0; i < numGlobalReach; i++) {
    reachStr >> dummy >> id >> hid[i] >> oid[i] >> aboveid[i];
  }

  // Read connectivity of stream reaches
  for (int i = 0; i < numGlobalReach; i++) {
    reachStr >> dummy >> id;
    tGraphNode rnode(id);
    conn.push_back(rnode);

    reachStr >> dummy >> upsize;
    for (int j = 0; j < upsize; j++) {
      reachStr >> value;
      conn[i].addUpstream(value);
    }

    reachStr >> dummy >> downsize;
    for (int j = 0; j < downsize; j++) {
      reachStr >> value;
      conn[i].addDownstream(value);
    }
  }
}

/*************************************************************************
**
** Create stream reach connectivity table from tKinemat object.
**
*************************************************************************/

void tGraph::connectivity() {
  // Get reach heads and outlets
  tPtrList< tCNode >& hlist = flow->getReachHeadList();
  tPtrList< tCNode >& olist = flow->getReachOutletList();

  
  // Set number of stream reaches
  // Initialize point count per reach
  numGlobalReach = hlist.getSize();
  for (int i = 0; i < numGlobalReach; i++)
      pointsPerReach.push_back(0);

  // Figure out connectivity
  tCNode *chead, *coutlet;
  tPtrListIter< tCNode > HeadIter(hlist);
  tPtrListIter< tCNode > OutletIter(olist);

  int i;
  hid = new int[numGlobalReach];
  oid = new int[numGlobalReach];
  // First collect head and outlet IDs
  for (chead = HeadIter.FirstP(), coutlet = OutletIter.FirstP(), i=0;
       !(HeadIter.AtEnd());
    chead = HeadIter.NextP(), coutlet = OutletIter.NextP(), i++) {
    if (sim->debug == 'Y') {
      Cout << "Reach " << i << " head " 
           << chead->getID() << " " << endl;
      Cout << "Reach " << i << " outlet "          
           << coutlet->getID() << " " << endl;
    }
    tGraphNode rnode(i);
    conn.push_back(rnode);
    hid[i] = chead->getID();
    oid[i] = coutlet->getID();
  }                  	                                                              
  // Figure out which are connected
  // For each reach, find all outlets == head, upstream
  //                 find all heads == outlet, downstream
  for (i = 0; i < numGlobalReach; i++) { 
    for (int j = 0; j < numGlobalReach; j++) {
      if (j != i) {
        if (oid[j] == hid[i]) conn[i].addUpstream(j);
        if (hid[j] == oid[i]) conn[i].addDownstream(j); 
      }
    }
  }
  
  // For debugging, show reach nodes and stream nodes
  if (sim->debug == 'Y') {
    for (i = 0; i < numGlobalReach; i++) 
      cout << conn[i] << endl;
  }
}

/*************************************************************************
**
** Partition graph based on number of partitions previously set.
**
*************************************************************************/

void tGraph::outputConnectivity(tInputFile& InputFile) {

  // Write out reach connectivity to a file
  char connName[kMaxNameSize];
  InputFile.ReadItem(connName, "OUTFILENAME");
  strcat(connName, ".connectivity");
  ofstream connFile;
  connFile.open(connName);
  if (!connFile.good()) {
      cout << "connFile problem" << endl;
      exit(5);
  }
  connFile << "# ReachID PointCount HeadNodeID OutletNodeID Partition NumberDownstream [REACHID1 REACHID2 ...] NumberFlux [FLUXID1 FlUXID2 ...]" << endl;
  connFile << numGlobalReach << " " << numGlobalPart << endl;
  for (int i = 0; i < numGlobalReach; i++) {
    connFile << conn[i].getID() << " "
             << pointsPerReach[i] << " "
             << hid[i] << " "
             << oid[i] << " "
             << reach2partition[i] << " ";
    std::vector<int> downReach = conn[i].getDownstream();
    connFile << downReach.size();
    for (int j = 0; j < downReach.size(); j++) {
      connFile << " " << downReach[j];
    }
    std::vector<int> flux = conn[i].getFlux();
    connFile << " " << flux.size();
    for (int j = 0; j < flux.size(); j++) {
        connFile << " " << flux[j]; 
    }
    connFile << "\n";
  }
  connFile.close();
}

/*************************************************************************
**
** Partition graph based on number of partitions previously set.
**
*************************************************************************/

void tGraph::partition(tInputFile& InputFile) {

#ifdef PARALLEL_TRIBS
  // If running in parallel, get # of processors as # of partitions
  // and processor number as local partition
  numGlobalPart = tParallel::getNumProcs();
  localPart = tParallel::getMyProc();
#endif

  assert(numGlobalPart > 0);
  partition(numGlobalPart, InputFile);

  // Create array of sets for upstream/downstream overlapping flow nodes
  upFlow = new std::set<tCNode*,IDOrder>[numGlobalPart];
  downFlow = new std::set<tCNode*,IDOrder>[numGlobalPart];

  // Create array of sets for flux overlapping nodes
  localFlux = new std::set<tCNode*,IDOrder>[numGlobalPart];
  remoteFlux = new std::set<tCNode*,IDOrder>[numGlobalPart];
}

/*************************************************************************
**
** Partition graph based on given number of partitions.
**
*************************************************************************/

void tGraph::partition(int np, tInputFile& InputFile) {

  assert(np > 0);

  if (sim->debug == 'Y') {
    Cout << "np = " << np << endl;
  }

  for (int i = 0; i < numGlobalReach; i++) 
      reach2partition.push_back(0);

  // Special case: only 1 partition
  // All reaches in partition 0
  if (np == 1) {
    localPart = 0;
    for (int i = 0; i < numGlobalReach; i++) 
      localReach.push_back(i);
  }
  
  // More than 1 partition
  else {
    // Check for partitioning file option
    int optgfile = 0;
    optgfile = InputFile.ReadItem( optgfile, "GRAPHOPTION"); 

    // Check for partitioning file
    if (optgfile == 1 || optgfile == 2) {
      char pfile[256];
      strcpy(pfile,"");
      InputFile.ReadItem( pfile, "GRAPHFILE" );

      // If a file was given, read in partitioning based
      // on file format type
      if (optgfile == 1)
        readReachPartitionFromFile(pfile);
      else
        readInletOutletPartitionFromFile(pfile);
    }

    // Else no partition file was given,
    // create a simple default partitioning
    else {
        createDefaultPartition(np);
    }

    // Collect local reaches together
    assert(localPart >= 0);
    for (int i = 0; i < numGlobalReach; i++) {
      if (reach2partition[i] == localPart) {
        localReach.push_back(i);
      }
    }
  }

  if (sim->debug == 'Y') {
    // Print out local reaches
    cout << "\nLocal reaches:" << endl;
    for (int i = 0; i < localReach.size(); i++)
      cout << localReach[i] << endl;

    // Print out reach to partition relationships
    for (int i = 0; i < numGlobalReach; i++)
      Cout << "Reach " << i << " to Partition " << reach2partition[i] 
           << endl;
  }
}

/*************************************************************************
**
** Read graph partitioning from file.
** The inlet/outlet format is as follows:
**
** A line for each reach should contain:
**   <partition number> <inlet ID> <outlet ID>
**
** Example (7 reaches on 4 processors):
**
**  0 1082 1145
**  0 1083 1145
**  1 1085 1156
**  1 1145 1156
**  2 1084 1190
**  2 1156 1190
**  3 1190 1198
**
*************************************************************************/

void tGraph::readInletOutletPartitionFromFile(char* pfile) {

    ifstream partFile;
    partFile.open(pfile);
    int whichPart, headID, outletID;
    int nr = 0;
    while (nr < numGlobalReach && !partFile.eof()) {
        partFile >> whichPart >> headID >> outletID;
        int preach = whichReach(headID, outletID);
        if (sim->debug == 'Y') {
            Cout << "Partition = " << whichPart
                 << " head = " << headID
                 << " outlet = " << outletID
                 << " reach = " << preach
                 << endl;
        }
        nr++;
        assert(preach >= 0 && preach < numGlobalReach);
        reach2partition[preach] = whichPart;
    }
    partFile.close();
    Cout << "\nPartitioning read from inlet/output file " << pfile << endl;
}


/*************************************************************************
**
** Read graph partitioning from file.
** The reach format is as follows:
**
** A line for each reach should contain:
**   <partition number> <reach ID>
**
** Example (7 reaches on 4 processors):
**
**  0 0 
**  0 1
**  1 3
**  1 4
**  2 2
**  2 5
**  3 6
**
*************************************************************************/

void tGraph::readReachPartitionFromFile(char* pfile) {

    ifstream partFile;
    partFile.open(pfile);
    int whichPart, preach;
    int nr = 0;
    while (nr < numGlobalReach && !partFile.eof()) {
        partFile >> whichPart >> preach;
        if (sim->debug == 'Y') {
            Cout << "Partition = " << whichPart
                 << " reach = " << preach
                 << endl;
        }
        nr++;
        assert(preach >= 0 && preach < numGlobalReach);
        reach2partition[preach] = whichPart;
    }
    partFile.close();
    Cout << "\nPartitioning read from reach file " << pfile << endl;
}


/*************************************************************************
**
** Create default partitions.
**
*************************************************************************/

void tGraph::createDefaultPartition(int np) {

    // Split reaches as evenly as possible between partitions
    // Last partition may have less than others
    assert(numGlobalReach >= np);
    int rstart = 0;
    int rend = -1;
    int rp = np;
    int nreach = numGlobalReach;
    for (int i = 0; i < np; i++) {
      rstart = rend + 1;
      rend = (i == (np - 1)) ? numGlobalReach - 1 :
                              rstart + (nreach + rp - 1) / rp - 1;
      nreach = nreach - (rend - rstart + 1);
      rp--;
      // Assign partition to reaches
      for (int j = rstart; j <= rend; j++)
        reach2partition[j] = i;
    }

}

/*************************************************************************
**
** Check if stream reach in the local partition.
**
*************************************************************************/

bool tGraph::inLocalPartition(int r) {
  if (r >= 0 && r < numGlobalReach)
    if (reach2partition[r] == localPart) return true;
  return false;
}

/*************************************************************************
**
** Check if last stream reach in the local partition.
**
*************************************************************************/

bool tGraph::hasLastReach() {
  return lastReach;
}

/*************************************************************************
**
** Is this the last reach in the local partition with a common 
** destination reach in a remote partition?
**
*************************************************************************/

bool tGraph::lastLocalReachWithCommonDest(int r, int rdest) {
  assert(r >= 0 && r < numGlobalReach);        // Legal reach id
  assert(reach2partition[r] == localPart);     // In local partition
  assert(reach2partition[rdest] != localPart); // In remote partition

  // Is destination actually downstream of reach r
  if ( isDownstreamOf(rdest, r) ) {
    std::vector<int> ureach = conn[ rdest ].getUpstream();
    for (int j = 0; j < ureach.size(); j++) {

      if ( ( hid[ureach[j]] < hid[r] ) && inLocalPartition( ureach[j] ) ) 
        return false;
    }
    // Last one
    return true;
  }

  // Reaches are not connected
  return false;
}

/*************************************************************************
**
** Is this the last reach in a remote partition with a common
** destination reach in the local partition?
**
*************************************************************************/

bool tGraph::lastRemoteReachWithCommonDest(int r, int rdest) {
  assert(r >= 0 && r < numGlobalReach);        // Legal reach id
  assert(reach2partition[r] != localPart);     // In remote partition
  assert(reach2partition[rdest] == localPart); // In local partition
  // Is destination actually downstream of reach r
  if ( isDownstreamOf(rdest, r) ) {
    std::vector<int> ureach = conn[ rdest ].getUpstream();
    for (int j = 0; j < ureach.size(); j++) {
      if ( ( hid[ureach[j]] < hid[r] ) && 
           ( reach2partition[ ureach[j] ] == reach2partition[r]) )
        return false;
    }
    // Last one
    return true;
  }

  // Reaches not connected
  return false;
}

/*************************************************************************
**
** Is first reach downstream of second reach?
**
*************************************************************************/

bool tGraph::isDownstreamOf(int r1, int r2) {
  assert(r1 >= 0 && r1 < numGlobalReach); // Legal reach id
  assert(r2 >= 0 && r2 < numGlobalReach); // Legal reach id
  std::vector<int> dreach = conn[ r2 ].getDownstream();
  for (int i = 0; i < dreach.size(); i++)
    if ( dreach[i] == r1 ) return true; // yes

  return false;                         // no
}

/*************************************************************************
**
** Is first reach upstream of second reach?
**
*************************************************************************/

bool tGraph::isUpstreamOf(int r1, int r2) {
  assert(r1 >= 0 && r1 < numGlobalReach); // Legal reach id
  assert(r2 >= 0 && r2 < numGlobalReach); // Legal reach id
  std::vector<int> ureach = conn[r2].getUpstream();
  for (int i = 0; i < ureach.size(); i++)
    if ( ureach[i] == r1 ) return true; // yes
                                                                                
  return false;                         // no
}

/*************************************************************************
**
** Return partition number of selected stream reach.
**
*************************************************************************/

int tGraph::getPartition(int r) {
  assert(r >= 0 && r < numGlobalReach);
  assert(reach2partition.size() >= r);
  return reach2partition[r];
}

/*************************************************************************
**
** Return list of partition members for selected partition.
**
*************************************************************************/

std::vector<int> tGraph::getPartitionMembers(int p) {
  assert(p >= 0 && p < numGlobalPart);
  std::vector<int> pmember;
  for (int i = 0; i < reach2partition.size(); i++) {
    if (reach2partition[i] == p)
      pmember.push_back(i);
  }
  return pmember; 
}

/*************************************************************************
**
** Print out node IDs for stream nodes in each reach.
**
*************************************************************************/

void tGraph::listStreamNodes() {

  // Get reach heads and outlets
  tPtrList< tCNode >& hlist = flow->getReachHeadList();
  tPtrList< tCNode >& olist = flow->getReachOutletList();
  tPtrListIter< tCNode > HeadIter(hlist);
  tPtrListIter< tCNode > OutletIter(olist);
  int i;
  tCNode *chead, *coutlet, *creach;

  for (chead = HeadIter.FirstP(), coutlet = OutletIter.FirstP(), i=0;
     !(HeadIter.AtEnd());
     chead = HeadIter.NextP(), coutlet = OutletIter.NextP(), i++) {
     creach = chead;
     while (creach != coutlet) {
        Cout << "Reach " << i << " Stream node " << creach->getID() << endl;
        creach = creach->getDownstrmNbr();
     }
     Cout << endl;
   }
}

/*************************************************************************
**
** Update Mesh and Flow so only data in the local partition is active.
** For stream reaches in the local partition, make all the internal points
** that flow to them active in the Mesh node list. Make the rest inactive.
** Remove or inactivate non-local stream reaches in tKinemat.
**
*************************************************************************/

void tGraph::update() {

  // Determine if last reach is on this processor
  lastReach = inLocalPartition(numGlobalReach-1);

  // If only running on 1 processor with 1 partition,
  // there is no need to do the rest of this
  if (numGlobalPart == 1) return;

  // Get reach heads and outlets
  tPtrList< tCNode >& hlist = flow->getReachHeadList();
  tPtrList< tCNode >& olist = flow->getReachOutletList();

  // Create iterator for all nodes, heads, and outlets
  tMeshList<tCNode> *nlist = mesh->getNodeList();
  tMeshListIter<tCNode> niter(nlist);
  tListNode<tCNode> *prevNode = 0;
  tMeshList<tEdge> *elist = mesh->getEdgeList();
  tMeshListIter<tEdge> eiter(elist);
  tListNode<tEdge> *prevEdge = 0;
  int i, ia, ea;
  tCNode *cn, *co;
  tCNode cnode;
  tEdge *ce;
  tEdge cedge;

  // Go through all the active nodes
  int ncount = 0;
  int nlimit = 10000;
  int ninc = 10000;
  int creach;
  ia = 0;
  int rc = 0;
  cn = niter.FirstP();
  while (niter.IsActive() && cn != 0) {
    // Count up points per reach
    creach = cn->getReach();
    if (creach >= 0 && creach < numGlobalReach)
      pointsPerReach[creach]++;
    else {
        cout << "Illegal reach number." << endl;
        exit(3);
    }

    // Check if node is associated with a reach not in this partition
    // Make nodes in reaches not in the local partition inactive
    // The previous node stays the same
    if (!inLocalPartition(creach)) {
      int xid = cn->getID();
      int xreach = cn->getReach();
      cn = niter.NextP();
      if(prevNode != 0) {
        nlist->nextToBack(prevNode);
      }
      // Case of first node
      else {
        nlist->frontToBack();
      }
      ia++;
    }
    // Node not made inactive, change previous node
    else {
      prevNode = niter.NodePtr();
      cn = niter.NextP();
    }

    ncount++;
    if (ncount == nlimit) {
      Cout << "Processed " << ncount << " nodes " << endl;
      nlimit += ninc;
    } 
  }

  Cout << "\nUpdated node list based on local reaches..." << endl;

  // Make any edge with an inactive origin node, inactive
  std::map<int,std::map<int,int> > reachFlux; 
  int ecount = 0;
  int elimit = 50000;
  int einc = 50000;
  ea = 0; 
  ce = eiter.FirstP();
  while (eiter.IsActive() && ce != 0) {

    co = (tCNode *)ce->getOriginPtrNC();
    cn = (tCNode *)ce->getDestinationPtrNC();

    int coReach = co->getReach();
    int cnReach = cn->getReach();
    if (coReach != cnReach) reachFlux[coReach][cnReach] = 1;

    bool coActive = inLocalPartition(coReach);
    bool cnActive = inLocalPartition(cnReach);
    
    // If origin node is inactive, make edge inactive
    // The previous edge stays the same
    if (!coActive) { 
      ce = eiter.NextP();
      if (prevEdge != 0) {
        elist->nextToBack(prevEdge);
      }
      // Case of first edge
      else {
        elist->frontToBack();
      }
      ea++;

      // Check for flow nodes
      // Determine which reaches the nodes are in
      if ((cn->getID() == hid[cnReach])
                   && cnActive 
                   && !coActive
                   && (cnReach > coReach)) {
        if (sim->debug == 'Y') {
          cout << "upFlow " << cn->getID() << " " << co->getID()
               << " -> " << cn->getID() << " "
               << coReach << " -> " << cnReach << endl;
        }
        upFlow[getPartition(coReach)].insert(cn);
      }
    }
    
    // Check for flow nodes
    else {

      // Determine which reaches the nodes are in
      // If destination node is inactive and
      // destination node is in a downstream reach
      if ((cn->getID() == hid[cnReach])
                    && !cnActive 
                    && coActive
                    && (cnReach > coReach)) {
        if (sim->debug == 'Y') {
            cout << "downFlow " << cn->getID() << " " << co->getID()
                 << " -> " << cn->getID() << " "
                 << coReach << " -> " << cnReach << endl;
        }
        downFlow[getPartition(cnReach)].insert(cn);
      }

      // Edge not made inactive, set new previous edge.
      prevEdge = eiter.NodePtr();
      ce = eiter.NextP();
    }

    ecount++;
    if (ecount == elimit) {
      Cout << "Processed " << ecount << " edges " << endl;
      elimit += einc;
    }
  }

  Cout << "\nUpdated edge list based on local reaches..." << endl;

  // Set reaches that exchange flux
  for (int i = 0; i < numGlobalReach; i++) {
    for (int j = 0; j < numGlobalReach; j++) {
      if ((i != j) && (reachFlux[i][j] == 1))
        conn[i].addFlux(j);
    }
  }
  reachFlux.clear();

  // Display after stats
  if (sim->debug == 'Y') {
    cout << "Partition " << localPart << " # nodes made inactive = " 
         << ia << endl;
    cout << "Partition " << localPart << " # edges made inactive = " 
         << ea << endl;
    reachNodeCounts();
  }

  // Calculate overlapping nodes for flux exchange
  calculateOverlap();

  // Calculate runoff/runon flux nodes
  calculateRunFlux();
}


/*************************************************************************
**
** Load the nodes and edges belonging to reaches in this processors
** partition, and load all boundary nodes and edges onto every processor
**
*************************************************************************/

void tGraph::ReadFlowMesh()
{
   int nnodes, nedges;
   int id, firstID, flowID, streamID, ccwID;
   int origID, destID, origReach, destReach, origBoundary, destBoundary;
   int origPartition, destPartition;
   int reach;

   tCNode curnode, origNode, destNode;
   tEdge curedge;
   tCNode* cn;

   fstream nodeStr, edgeStr, fluxNodeStr, fluxEdgeStr;

#ifdef PARALLEL_TRIBS
   tParallel::barrier();
#endif

   // Files with one copy of every node and edge ordered by reach
   nodeStr.open("nodes.meshb", ios::in | ios::binary);
   edgeStr.open("edges.meshb", ios::in | ios::binary);
   BinaryRead(edgeStr, nedges);
   BinaryRead(nodeStr, nnodes);

   // Files with duplicate nodes and edges per reach for partitioned mesh
   fluxNodeStr.open("fluxnodes.meshb", ios::in | ios::binary);
   fluxEdgeStr.open("fluxedges.meshb", ios::in | ios::binary);

   // Boundary directory information is in the last slot of reach information
   int boundaryReach = numGlobalReach;

   // Map of reaches to connecting reaches
   std::map<int,std::map<int,int> > reachFlux;

   // Get node and edge list from the tMesh we are going to fill in
   tMeshList<tCNode>* nodeList = mesh->getNodeList();
   tMeshList<tEdge>* edgeList = mesh->getEdgeList();

   // Determine if last reach is on this processor
   lastReach = inLocalPartition(numGlobalReach - 1);

   //////////////////////////////////////////////////////////////////////////
   //
   // Pass 1 read of the node file
   // Creates all nodes and inserts in nodelist but doesn't have edge pointers
   // to match up with edge ids for first edge and flow edge
   //
   cout << "Read reach nodes: Pass 1" << endl;

   // Local reach nodes added to mesh
   std::vector<int>::iterator riter;
   for (riter = localReach.begin(); riter != localReach.end(); riter++) {
      reach = (*riter);
      nodeStr.seekg(nodeOffset[reach], ios::beg);

      for (int node = 0; node < nodesPerReach[reach]; node++) {
         ReadFlowNode(nodeStr, &curnode);
         nodeList->insertAtActiveBack(curnode);
      }
   }

   // Boundary (inactive) reach nodes added to mesh
   nodeStr.seekg(nodeOffset[boundaryReach], ios::beg);
   for (int node = 0; node < nodesPerReach[boundaryReach]; node++) {
      ReadFlowNode(nodeStr, &curnode);
      nodeList->insertAtBack(curnode);
   }

   // Create the lookup table for currently loaded nodes which are the local
   // nodes for the reaches on this partition and the boundary nodes
   mesh->allocateNodeTable(nnodes);
   tMeshListIter< tCNode > nodIter( nodeList );
   for (cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP())
      mesh->setNodeFromID(cn->getID(), cn);

   // Local flux nodes added to mesh may be duplicates because both reaches
   // are local to this partition, so check before adding
   for (riter = localReach.begin(); riter != localReach.end(); riter++) {
      origReach = (*riter);
      fluxNodeStr.seekg(fluxNodeOffset[origReach], ios::beg);

      // Read the flux nodes attached to this reach
      for (int node = 0; node < fluxNodesPerReach[origReach]; node++) {
         ReadFlowNode(fluxNodeStr, &destNode);

         // Add flux node if not in the mesh already
         if (mesh->getNodeFromID(destNode.getID()) == 0) {
            nodeList->insertAtBack(destNode);
            destReach = destNode.getReach();

            // Note that this is not good because we don't have the pointer
            // to the destNode coming out of the insertAtBack
            mesh->setNodeFromID(destNode.getID(), &destNode);
         }
      }
   }

   // Flux nodes for boundary nodes
   fluxNodeStr.seekg(fluxNodeOffset[boundaryReach], ios::beg);

   for (int node = 0; node < fluxNodesPerReach[boundaryReach]; node++) {
      ReadFlowNode(fluxNodeStr, &destNode);

      // Add flux node if not in the mesh because of another reach
      if (mesh->getNodeFromID(destNode.getID()) == 0) {
         nodeList->insertAtBack(destNode);
         destReach = destNode.getReach();

         // Note that this is not good because we don't have the pointer
         // to the destNode coming out of the insertAtBack
         mesh->setNodeFromID(destNode.getID(), &destNode);
      }
   }

   // Create node table with addition of flux nodes
   mesh->deleteNodeTable();
   mesh->allocateNodeTable(nnodes);
   for (cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP())
      mesh->setNodeFromID(cn->getID(), cn);

   //////////////////////////////////////////////////////////////////////////
   //
   // Pass 1 read of the edge file
   //
   cout << "Read reach edges: Pass 1" << endl;
   // Read every edge in local reaches
   for (riter = localReach.begin(); riter != localReach.end(); riter++) {
      origReach = (*riter);
      edgeStr.seekg(edgeOffset[origReach], ios::beg);

      // Internal edges where origin and destination are in same reach
      for (int edge = 0; edge < internalEdgesPerReach[origReach]; edge++) {
         ReadFlowEdge(edgeStr, &curedge, &origID, &destID);

         tCNode *origNode = mesh->getNodeFromID(origID);
         curedge.setOriginPtr(origNode);

         tCNode *destNode = mesh->getNodeFromID(destID);
         curedge.setDestinationPtr(destNode);

         curedge.setFlowAllowed(1);
         edgeList->insertAtActiveBack(curedge);
      }

      // External edges where origin and destination in different reach
      // If the reach = -1 it is the boundary nodes and must always be loaded
      for (int edge = 0; edge < externalEdgesPerReach[origReach]; edge++) {
         ReadFlowEdge(edgeStr, &curedge, &origID, &destID);

         tCNode *origNode = mesh->getNodeFromID(origID);
         curedge.setOriginPtr(origNode);
         origBoundary = origNode->getBoundaryFlag();
         origPartition = -1;
         if (origReach != -1)
            origPartition = getPartition(origReach);

         tCNode *destNode = mesh->getNodeFromID(destID);
         curedge.setDestinationPtr(destNode);
         destBoundary = destNode->getBoundaryFlag();
         destReach = destNode->getReach();
         destPartition = -1;
         if (destReach != -1)
            destPartition = getPartition(destReach);

         // set the "flowallowed" status and put in edge list appropriately
         if (origBoundary == kClosedBoundary ||
             destBoundary == kClosedBoundary ||
             (origBoundary == kOpenBoundary && destBoundary == kOpenBoundary)) {
            curedge.setFlowAllowed(0);
            edgeList->insertAtBack(curedge);
         } else {
            curedge.setFlowAllowed(1);
            edgeList->insertAtActiveBack(curedge);

            // Record the existence of flux sharing
            reachFlux[origReach][destReach] = 1;
         }

         // Set upflow nodes which are always attached to flux nodes
         if (origNode->getID() == hid[origReach] && 
             origPartition == localPart &&
             destPartition != localPart &&
             destPartition != -1 &&
             origReach > destReach)
            upFlow[destPartition].insert(origNode);

         // Set downflow nodes which are always attached to flux nodes
         if (destNode->getID() == hid[destReach] && 
             origPartition == localPart &&
             destPartition != localPart &&
             destPartition != -1 &&
             destReach > origReach)
            downFlow[destPartition].insert(destNode);
      }

      // Flux edges are completely in another reach and connect a flux node
      // to its flow edge which is required by groundwater
      fluxEdgeStr.seekg(fluxEdgeOffset[origReach], ios::beg);
      for (int edge = 0; edge < fluxEdgesPerReach[origReach]; edge++) {
         ReadFlowEdge(fluxEdgeStr, &curedge, &origID, &destID);

         tCNode *origNode = mesh->getNodeFromID(origID);
         curedge.setOriginPtr(origNode);

         tCNode *destNode = mesh->getNodeFromID(destID);
         curedge.setDestinationPtr(destNode);

         curedge.setFlowAllowed(0);
         edgeList->insertAtBack(curedge);
      }
   }

   // Read boundary (inactive) edges (internal and external treated the same)
   edgeStr.seekg(edgeOffset[boundaryReach], ios::beg);
   int totalBoundaryEdges = internalEdgesPerReach[boundaryReach] +
                            externalEdgesPerReach[boundaryReach];

   for (int edge = 0; edge < totalBoundaryEdges; edge++) {
      ReadFlowEdge(edgeStr, &curedge, &origID, &destID);

      tCNode *origNode = mesh->getNodeFromID(origID);
      curedge.setOriginPtr(origNode);
      origBoundary = origNode->getBoundaryFlag();

      tCNode *destNode = mesh->getNodeFromID(destID);
      curedge.setDestinationPtr(destNode);
      destBoundary = destNode->getBoundaryFlag();

      // set the "flowallowed" status and put in edge list appropriately
      curedge.setFlowAllowed(0);
      edgeList->insertAtBack(curedge);
   }

   // Create the lookup table of edges indexed by edge id
   tEdge** EdgeTable = new tEdge*[nedges];
   tEdge* ce; 
   tMeshListIter< tEdge > edgIter( edgeList );
   for (ce = edgIter.FirstP(); !(edgIter.AtEnd()); ce = edgIter.NextP())
      EdgeTable[ce->getID()] = ce;
      
   // Read flux edges which might already be in the edge list
   fluxEdgeStr.seekg(fluxEdgeOffset[boundaryReach], ios::beg);
   for (int edge = 0; edge < fluxEdgesPerReach[boundaryReach]; edge++) {
      ReadFlowEdge(fluxEdgeStr, &curedge, &origID, &destID);

      if (EdgeTable[curedge.getID()] == 0) {
         EdgeTable[curedge.getID()] = &curedge;
         tCNode *origNode = mesh->getNodeFromID(origID);
         curedge.setOriginPtr(origNode);
         origBoundary = origNode->getBoundaryFlag();

         tCNode *destNode = mesh->getNodeFromID(destID);
         curedge.setDestinationPtr(destNode);
         destBoundary = destNode->getBoundaryFlag();

         curedge.setFlowAllowed(0);
         edgeList->insertAtBack(curedge);
      }
   }

   // Create the lookup table of edges indexed by edge id
   delete [] EdgeTable;
   EdgeTable = new tEdge*[nedges];
   for (ce = edgIter.FirstP(); !(edgIter.AtEnd()); ce = edgIter.NextP())
      EdgeTable[ce->getID()] = ce;
      
   //////////////////////////////////////////////////////////////////////////
   //
   // Pass 2 read of the node file to fill in node pointers
   //
   cout << "Read reach nodes: Pass 2" << endl;

   // Local reach nodes read to retrieve ids which can be looked up in tables
   for (riter = localReach.begin(); riter != localReach.end(); riter++) {
      reach = (*riter);
      nodeStr.seekg(nodeOffset[reach], ios::beg);

      for (int node = 0; node < nodesPerReach[reach]; node++) {
         ReadFlowNode(nodeStr, &id, &firstID, &flowID, &streamID);

         // Get the pointer to the existing node in the mesh
         cn = mesh->getNodeFromID(id);
         if (cn != 0) {
            cn->setEdg(EdgeTable[firstID]);
            if (flowID >= 0)
               cn->setFlowEdg(EdgeTable[flowID]);
            if (streamID >= 0)
               cn->setStreamNode(mesh->getNodeFromID(streamID));
         }
      }
   }

   // Boundary reach nodes read to retrieve ids which can be looked up in tables
   nodeStr.seekg(nodeOffset[boundaryReach], ios::beg);

   for (int node = 0; node < nodesPerReach[boundaryReach]; node++) {
      ReadFlowNode(nodeStr, &id, &firstID, &flowID, &streamID);

      cn = mesh->getNodeFromID(id);
      if (cn != 0) {
         cn->setEdg(EdgeTable[firstID]);
         if (flowID >= 0)
            cn->setFlowEdg(EdgeTable[flowID]);
         if (streamID >= 0)
            cn->setStreamNode(mesh->getNodeFromID(streamID));
      }
   }

   // Flux boundary nodes read to retrieve ids which can be looked up in tables
   fluxNodeStr.seekg(fluxNodeOffset[boundaryReach], ios::beg);

   for (int node = 0; node < fluxNodesPerReach[boundaryReach]; node++) {
      ReadFlowNode(fluxNodeStr, &id, &firstID, &flowID, &streamID);

      cn = mesh->getNodeFromID(id);
      if (cn != 0) {

         // When a boundary node goes to a node within a reach that node is
         // called a flux node and it is stored along with its ccw edges
         // Additionally the destination node for the ccw edges are stored
         // as flux nodes BUT not their ccw edges, so we must test before
         // assigning setEdg()  Note that if something else fails this
         // might be the problem

         if (EdgeTable[firstID] != 0)
            cn->setEdg(EdgeTable[firstID]);
         if (flowID >= 0)
            cn->setFlowEdg(EdgeTable[flowID]);
         if (streamID >= 0)
            cn->setStreamNode(mesh->getNodeFromID(streamID));
      }
   }

   // Flux nodes read to retrieve ids which can be looked up in tables
   for (riter = localReach.begin(); riter != localReach.end(); riter++) {
      reach = (*riter);
      fluxNodeStr.seekg(fluxNodeOffset[reach], ios::beg);

      for (int node = 0; node < fluxNodesPerReach[reach]; node++) {
         ReadFlowNode(fluxNodeStr, &id, &firstID, &flowID, &streamID);
         cn = mesh->getNodeFromID(id);
         if (cn != 0) {
            if (EdgeTable[firstID] != 0)
               cn->setEdg(EdgeTable[firstID]);
            if (flowID >= 0)
               cn->setFlowEdg(EdgeTable[flowID]);
            if (streamID >= 0)
               cn->setStreamNode(mesh->getNodeFromID(streamID));
         }
      }
   }

   //////////////////////////////////////////////////////////////////////////
   //
   // Pass 2 read of the edge file to fill in edge pointers
   //
   cout << "Read reach edges: Pass 2" << endl;
   tEdge *curedg, *ccwedg;

   // Local reach edges
   for (riter = localReach.begin(); riter != localReach.end(); riter++) {
      reach = (*riter);
      edgeStr.seekg(edgeOffset[reach], ios::beg);

      // Local reach internal edges
      for (int edge = 0; edge < internalEdgesPerReach[reach]; edge++) {
         ReadFlowEdge(edgeStr, &id, &ccwID);
         if (ccwID >= 0 && ccwID < nedges)
            EdgeTable[id]->setCCWEdg(EdgeTable[ccwID]);
      }

      // Local reach external edges
      for (int edge = 0; edge < externalEdgesPerReach[reach]; edge++) {
         ReadFlowEdge(edgeStr, &id, &ccwID);
         if (ccwID >= 0 && ccwID < nedges)
            EdgeTable[id]->setCCWEdg(EdgeTable[ccwID]);
      }

      // Boundary flux edges
      fluxEdgeStr.seekg(fluxEdgeOffset[reach], ios::beg);
      for (int edge = 0; edge < fluxEdgesPerReach[boundaryReach]; edge++) {
         ReadFlowEdge(fluxEdgeStr, &id, &ccwID);
         if (ccwID >= 0 && ccwID < nedges && EdgeTable[ccwID] != 0)
            EdgeTable[id]->setCCWEdg(EdgeTable[ccwID]);
      }
   }

   // Boundary edges
   edgeStr.seekg(edgeOffset[boundaryReach], ios::beg);
   for (int edge = 0; edge < totalBoundaryEdges; edge++) {
      ReadFlowEdge(edgeStr, &id, &ccwID);
      if (ccwID >= 0 && ccwID < nedges && EdgeTable[ccwID] != 0)
         EdgeTable[id]->setCCWEdg(EdgeTable[ccwID]);
   }

   // Boundary flux edges
   fluxEdgeStr.seekg(fluxEdgeOffset[boundaryReach], ios::beg);
   for (int edge = 0; edge < fluxEdgesPerReach[boundaryReach]; edge++) {
      ReadFlowEdge(fluxEdgeStr, &id, &ccwID);
      if (ccwID >= 0 && ccwID < nedges && EdgeTable[ccwID] != 0)
         EdgeTable[id]->setCCWEdg(EdgeTable[ccwID]);
   }

   // Set reaches that exchange flux
   for (int i = 0; i < numGlobalReach; i++) {
      for (int j = 0; j < numGlobalReach; j++) {
         if ((i != j) && (reachFlux[i][j] == 1))
            conn[i].addFlux(j);
      }
   }
   reachFlux.clear();

   edgeStr.close();
   nodeStr.close();
   fluxNodeStr.close();
   fluxEdgeStr.close();

   delete [] EdgeTable;

#ifdef PARALLEL_TRIBS
   tParallel::barrier();
#endif
}

/***************************************************************************
**      
** Read the flow node information
** Filled tCNode is returned along with id variables which must be looked
** up in EdgeTable or NodeTable to locate actual pointers
**      
***************************************************************************/
          
void tGraph::ReadFlowNode(fstream& nodeStr, tCNode* curnode)
{       
   int id, firstEdgeID, flowEdgeID, streamNodeID;
   int boundary, flood, tracer, reach;
   double hillpath, traveltime, srf, hsrf, psrf, satsrf, sbsrf;
   double SoilMoistureSC, SoilMoistureUNSC, RootMoistureSC;
   double EvapoTranspiration, ContrArea, NwtNew, Rain, Curvature, streampath;
   double x, y, z, varea;

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

   curnode->setID(id);
   curnode->setBoundaryFlag(boundary);
   curnode->setX(x);
   curnode->setY(y);
   curnode->setZ(z);
   curnode->setVArea(varea);
   curnode->setVArea_Rcp(1.0 / varea);
   curnode->setFloodStatus(flood);
   curnode->setTracer(tracer);
   curnode->setHillPath(hillpath);
   curnode->setTTime(traveltime);
   curnode->setsrf(srf);
   curnode->sethsrf(hsrf);
   curnode->setpsrf(psrf);
   curnode->setsatsrf(satsrf);
   curnode->setsbsrf(sbsrf);
   curnode->setEvapoTrans(EvapoTranspiration);
   curnode->setSoilMoistureSC(SoilMoistureSC);
   curnode->setSoilMoistureUNSC(SoilMoistureUNSC);
   curnode->setRootMoistureSC(RootMoistureSC);
   curnode->setContrArea(ContrArea);
   curnode->setNwtNew(NwtNew);
   curnode->setRain(Rain);
   curnode->setCurvature(Curvature);
   curnode->setStreamPath(streampath);
   curnode->setReach(reach);
}   

void tGraph::ReadFlowNode(fstream& nodeStr,
                          int* id, int* firstEdgeID, 
                          int* flowEdgeID, int* streamNodeID)
{       
   int boundary, flood, tracer, reach;
   double hillpath, traveltime, srf, hsrf, psrf, satsrf, sbsrf;
   double SoilMoistureSC, SoilMoistureUNSC, RootMoistureSC;
   double EvapoTranspiration, ContrArea, NwtNew, Rain, Curvature, streampath;
   double x, y, z, varea;

   BinaryRead(nodeStr, *id);
   BinaryRead(nodeStr, boundary);
   BinaryRead(nodeStr, x);
   BinaryRead(nodeStr, y);
   BinaryRead(nodeStr, z);
   BinaryRead(nodeStr, varea);
   BinaryRead(nodeStr, *firstEdgeID);
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
   BinaryRead(nodeStr, *flowEdgeID);
   BinaryRead(nodeStr, *streamNodeID);
}

    
/***************************************************************************
**  
** Read the flow edge information
**    
***************************************************************************/
    
void tGraph::ReadFlowEdge(fstream& edgeStr, tEdge* curedge,
                          int* origID, int* destID)
{   
   int id, ccwID;
   double x, y, z, varea, rvtx[2], length, slope, vedglen;
   tArray<double> RVtx(2);

   BinaryRead(edgeStr, id);
   BinaryRead(edgeStr, rvtx[0]);
   BinaryRead(edgeStr, rvtx[1]);
   BinaryRead(edgeStr, length);
   BinaryRead(edgeStr, slope);
   BinaryRead(edgeStr, vedglen);
   BinaryRead(edgeStr, *origID);
   BinaryRead(edgeStr, *destID);
   BinaryRead(edgeStr, ccwID);

   curedge->setID(id);
   RVtx[0] = rvtx[0];
   RVtx[1] = rvtx[1];
   curedge->setRVtx(RVtx);
   curedge->setLength(length);
   curedge->setSlope(slope);
   curedge->setVEdgLen(vedglen);
} 

void tGraph::ReadFlowEdge(fstream& edgeStr, int* id, int* ccwID)
{   
   int origID, destID;
   double x, y, z, varea, rvtx[2], length, slope, vedglen;
   tArray<double> RVtx(2);

   BinaryRead(edgeStr, *id);
   BinaryRead(edgeStr, rvtx[0]);
   BinaryRead(edgeStr, rvtx[1]);
   BinaryRead(edgeStr, length);
   BinaryRead(edgeStr, slope);
   BinaryRead(edgeStr, vedglen);
   BinaryRead(edgeStr, origID);
   BinaryRead(edgeStr, destID);
   BinaryRead(edgeStr, *ccwID);
}

/*************************************************************************
**
** List nodes in partition.
**
*************************************************************************/

void tGraph::listActiveNodes() {

  // Open file and write out node ids
  ofstream PartNodeOut;
  char fname[100];
  snprintf(fname,sizeof(fname),"partition%d.nodes", localPart);//WR--09192023: 'sprintf' is deprecated: This function is provided for compatibility reasons only.
  PartNodeOut.open(fname);

  // Create iterator for all nodes
  tMeshList<tCNode> *nlist = mesh->getNodeList();
  tMeshListIter<tCNode> niter(nlist);
  tCNode *cn, *cdest;
  tEdge *ce;
                                                                               
  // Go through all the active nodes
  // Print out ids
  cn = niter.FirstP();
  while (niter.IsActive()) {
      ce = cn->getFlowEdg();
      cdest = (tCNode*)ce->getDestinationPtrNC();
      int fromreach = cn->getReach();
      int toreach = cdest->getReach();
      PartNodeOut << cn->getID() << ":" << fromreach 
         << " " << cn->getBoundaryFlag() << " -> " 
         << cdest->getID() << " " << toreach
         << " " << cdest->getBoundaryFlag() 
         << ((fromreach != toreach) ? "*****" : " ") << endl;
      cn = niter.NextP();
  }

  PartNodeOut.close();
}

/*************************************************************************
**
** Determine runon/runoff flux nodes.
**
*************************************************************************/

void tGraph::calculateRunFlux() {

  // Collect local nodes above outlets
  for (int i = 0; i < numGlobalReach; i++)
      nodeAboveOutlet.push_back(NULL);

  // Get reach heads and outlets
  tPtrList< tCNode >& hlist = flow->getReachHeadList();
  tPtrList< tCNode >& olist = flow->getReachOutletList();
  tCNode *chead, *coutlet, *cnext, *crflux;
  tPtrListIter< tCNode > HeadIter(hlist);
  tPtrListIter< tCNode > OutletIter(olist);

  int i;
  // Loop through reach heads and outlets, looking for 
  // a head in another partition and outlet in the local partition
  // The node upstream from the outlet is the flux node (node above outlet)

  chead = HeadIter.FirstP();
  coutlet = OutletIter.FirstP();

  for (i=0;!(HeadIter.AtEnd()) && i < numGlobalReach-1;i++) { //WR--09192023: modified for loop to address this warning: left operand of comma operator has no effect [-Wunused-value]

    int hpart = getPartition(chead->getReach());
    int opart = getPartition(coutlet->getReach());
    // Outlet is local
    if (hpart != localPart && opart == localPart) {    
      cnext = chead;
      crflux = chead;
      while (cnext != coutlet) {
        crflux = cnext;
        cnext = cnext->getDownstrmNbr();
      }
      nodeAboveOutlet[crflux->getReach()] = crflux;

      if (sim->debug == 'Y') {
        cout << " RunFlux " << crflux->getID() 
             << " " << hpart
             << "   -> " << coutlet->getID() 
             << " " << localPart << endl;
      }
    }
    chead = HeadIter.NextP();
    coutlet = OutletIter.NextP();
  }
}

/*************************************************************************
**
** Determine overlapping nodes for flux exchange.
**
*************************************************************************/

void tGraph::calculateOverlap() {
  if (numGlobalPart == 1) return;
 
  tCNode *cnorg;
  tCNode *cndest;
  tEdge *ce;
  tMeshListIter<tEdge> edgIter(mesh->getEdgeList());

  for (ce = edgIter.FirstP(); edgIter.IsActive(); ce = edgIter.NextP() ) {
    //Destination and Origin Nodes
    cnorg = (tCNode*)ce->getOriginPtrNC();
    cndest = (tCNode*)ce->getDestinationPtrNC();

    //Excluding calculation of flux to the outlet point
    if ( (cnorg->getBoundaryFlag() != kOpenBoundary) &&
         (cndest->getBoundaryFlag() != kOpenBoundary) &&
         (cnorg->getBoundaryFlag() != kClosedBoundary) &&
         (cndest->getBoundaryFlag() != kClosedBoundary) ) {

       int cnReach = cndest->getReach();
       bool cnActive = inLocalPartition(cnReach);

       // If the destination node is inactive, save destination as a remote 
       // flux node. Save origin as a local flux node for another partition.
       if (!cnActive) {
         int cnPart = getPartition(cnReach);
         remoteFlux[cnPart].insert(cndest);
         localFlux[cnPart].insert(cnorg);

         if (sim->debug == 'Y') {
             cout << "local " << cnorg->getID() << " " << cnorg->getBoundaryFlag()
                  << " remote " << cndest->getID() << " " << cndest->getBoundaryFlag()
                  << endl;
         }
       }
    }
  }

  // Display after stats
  if (sim->debug == 'Y') {

    for (int i = 0; i < numGlobalPart; i++) {
      if (i != localPart) {
        cout << "Partition " << localPart << " number remoteFlux = " 
             << i << " " << remoteFlux[i].size() << endl;
        cout << "Partition " << localPart << " number localFlux = " 
             << i << " " << localFlux[i].size() << endl;
      }
    }

    std::set<tCNode*>::iterator ir;
    for (int i = 0; i < numGlobalPart; i++) {
      for (ir = remoteFlux[i].begin(); ir != remoteFlux[i].end(); ++ir) {
        cout << "Local Partition " << localPart << " remoteFlux from partition "
             << i << " " << (*ir)->getID() << " " << endl;
      }
    }

    std::set<tCNode*>::iterator il;
    for (int i = 0; i < numGlobalPart; i++) {
      for (il = localFlux[i].begin(); il != localFlux[i].end(); ++il) {
        cout << "Local Partition " << localPart << " localFlux to partition " 
             << i << " " << (*il)->getID() << " " << endl;
      }
    }

  }
}

/*************************************************************************
**
** Print node counts for reaches.
**
*************************************************************************/

void tGraph::reachNodeCounts() {
  // Get reach heads and outlets
  tPtrList< tCNode >& hlist = flow->getReachHeadList();
  tPtrList< tCNode >& olist = flow->getReachOutletList();

  // Create iterator for all nodes, heads, and outlets
  tMeshList<tCNode> *nlist = mesh->getNodeList();
  tMeshListIter<tCNode> niter(nlist);
  int i;

  // Collect information on node list
  int nsize = nlist->getSize();
  cout << "# nodes in list = " << nsize << endl;
  int asize = nlist->getActiveSize();
  cout << "# active nodes in list = " << asize << endl;

  // Count # of nodes associated with each reach
  // The outlet is only associated with the last reach
  int *ncount = new int[numGlobalReach];
  for (i = 0; i < numGlobalReach; i++)
    ncount[i] = 0;

  tCNode *cn, *chead, *coutlet, *creach;
  int rnum;
  // Figure out how many nodes are associated with each reach
  for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP()) {
    rnum = cn->getReach();
    if (rnum >= 0) ncount[rnum]++;
  }
  for (i = 0; i < numGlobalReach; i++)
    cout << "Reach " << i << " # nodes = " << ncount[i] << endl;
}

/*************************************************************************
**
** Return which reach a node is in.
**
*************************************************************************/

int tGraph::inWhichReach(tCNode* cn) {
  // Get reach heads and outlets, iterators
  tPtrList< tCNode >& hlist = flow->getReachHeadList();
  tPtrList< tCNode >& olist = flow->getReachOutletList();
  tPtrListIter< tCNode > HeadIter(hlist);
  tPtrListIter< tCNode > OutletIter(olist);
  int i;
  tCNode *creach, *coutlet, *chead;

  // Figure out which reach a node is in
  for (chead = HeadIter.FirstP(), coutlet = OutletIter.FirstP(), i=0;
       !(HeadIter.AtEnd());
       chead = HeadIter.NextP(), coutlet = OutletIter.NextP(), i++) {
    creach = chead;
    // Final outlet associated with last reach
    if ((i == numGlobalReach-1) &&
       ((cn->getStreamNode() == coutlet) || (cn == coutlet))) {
      return i;
    }
    // Check if in a reach
    while (creach != coutlet) {
      if ((cn->getStreamNode() == creach) || (cn == creach)) {
        return i;
      }
      creach = creach->getDownstrmNbr();
    }
  }
  return -1;
}

/*************************************************************************
**
** Return which reach is defined by the given head and outlet nodes.
**
*************************************************************************/

int tGraph::whichReach(int iHead, int iOutlet) {

  // Get reach heads and outlets
  tPtrList< tCNode >& hlist = flow->getReachHeadList();
  tPtrList< tCNode >& olist = flow->getReachOutletList();

  // Create iterator for all heads, and outlets
  tPtrListIter< tCNode > HeadIter(hlist);
  tPtrListIter< tCNode > OutletIter(olist);
  int i; // reach counter

  tCNode *coutlet, *chead;
  for (chead = HeadIter.FirstP(), coutlet = OutletIter.FirstP(), i=0;
     !(HeadIter.AtEnd()); 
     chead = HeadIter.NextP(), coutlet = OutletIter.NextP(), i++) {
    if (iHead == chead->getID() && iOutlet == coutlet->getID())
        return i;
  }
  
  // Not found
  return -1;  
}


/*************************************************************************
**
** Check if a stream reach has upstream stream reaches.
**
*************************************************************************/

bool tGraph::hasUpstream(int r) {
  assert(r >= 0 && r < numGlobalReach);
  assert(conn.size() >= r);
  return conn[r].hasUpstream();
}

/*************************************************************************
**
** Check if a stream reach has downstream stream reaches.
**
*************************************************************************/

bool tGraph::hasDownstream(int r) {
  assert(r >= 0 && r < numGlobalReach);
  assert(conn.size() >= r);
  return conn[r].hasDownstream();
}

/*************************************************************************
**
** Check if a node is in one of the overlap lists.
**
*************************************************************************/

bool tGraph::isOverlapNode(tCNode* cn) {

  for (int i = 0; i < numGlobalPart; i++) {
    set<tCNode*>::iterator ilocal = localFlux[i].find(cn);
    set<tCNode*>::iterator iremote = remoteFlux[i].find(cn);
    // If found, return true
    if ((ilocal != localFlux[i].end()) || (iremote != remoteFlux[i].end()))
      return true;
  }

  // Not found
  return false;
}

/*************************************************************************
**
** Check if a node is in the downstream list.
**
*************************************************************************/

bool tGraph::isDownstreamNode(tCNode* cn) {

  for (int i = 0; i < numGlobalPart; i++) {
    set<tCNode*>::iterator idwn = downFlow[i].find(cn);
    // If found, return true
    if (idwn != downFlow[i].end())
      return true;
  }

  // Not found
  return false;
}

/*************************************************************************
**
** Check if a node is in the upstream list.
**
*************************************************************************/

bool tGraph::isUpstreamNode(tCNode* cn) {

  for (int i = 0; i < numGlobalPart; i++) {
    set<tCNode*>::iterator iup = upFlow[i].find(cn);
    // If found, return true
    if (iup != upFlow[i].end())
      return true;
  }

  // Not found
  return false;
}

/*************************************************************************
**
** Check if a node is in the local flux list.
**
*************************************************************************/

bool tGraph::isLocalfluxNode(tCNode* cn) {

  for (int i = 0; i < numGlobalPart; i++) {
    set<tCNode*>::iterator iflux = localFlux[i].find(cn);
    // If found, return true
    if (iflux != localFlux[i].end())
      return true;
  }

  // Not found
  return false;
}

/*************************************************************************
**
** Check if a node is in the remote flux list.
**
*************************************************************************/

bool tGraph::isRemotefluxNode(tCNode* cn) {

  for (int i = 0; i < numGlobalPart; i++) {
    set<tCNode*>::iterator iflux = remoteFlux[i].find(cn);
    // If found, return true
    if (iflux != remoteFlux[i].end())
      return true;
  }

  // Not found
  return false;
}

/*************************************************************************
**
** Send data to downstream stream reach(s).
**
*************************************************************************/

void tGraph::sendDownstream(int rid, tCNode* snode, double value) {
  assert(rid >= 0 && rid < numGlobalReach);
#ifdef PARALLEL_TRIBS
  // double* ndata = new double[1]; //WR debug moved to inside scope of if statment to prevent memory leak
  // Get list of downstream reaches
  std::vector<int> dreach = conn[rid].getDownstream();
  // For each on another processor, send data
  for (int i = 0; i < dreach.size(); i++) {

    // Send if last reach
    if ( !inLocalPartition(dreach[i])) {
        double* ndata = new double[1]; // WR debug added
      int to_proc = reach2partition[ dreach[i] ]; // To processor
      // Collect node data and send
      ndata[0] = value;
      tParallel::send(to_proc, DOWNSTREAM+snode->getReach(), ndata, 1);
    }
  }
#endif
}

/*************************************************************************
**
** Receive data from upstream stream reach(s).
**
*************************************************************************/

void tGraph::receiveUpstream(int rid, tCNode* rnode) {
  assert(rid >= 0 && rid < numGlobalReach);
#ifdef PARALLEL_TRIBS
  double* ndata = new double[1];
  // Get list of upstream reaches
  std::vector<int> ureach = conn[rid].getUpstream();
  // For each on another processor, receive data
  for (int i = 0; i < ureach.size(); i++) {

    // Check for last remote reach coming to this local reach
    if ( !inLocalPartition(ureach[i])) {
    
      int from_proc = reach2partition[ ureach[i] ]; // From processor
      // Receive data and update node
      tParallel::receive(from_proc, DOWNSTREAM+rnode->getReach(), ndata, 1);
      rnode->addQstrm(ndata[0]);
    }
  }
  delete [] ndata;
  tParallel::freeBuffers();// WR debug: put this at end of each receive call, it checks to see which previously assinged pointer  arrays can be safely deleted
#endif
}

/*************************************************************************
**
** Reset overlap node values to zero.
**
*************************************************************************/

void tGraph::resetOverlap() {

#ifdef PARALLEL_TRIBS
  if (numGlobalPart == 1) return;

  for (int i = 0; i < numGlobalPart; i++) {
    std::set<tCNode*>::iterator ir;
    for (ir = remoteFlux[i].begin(); ir != remoteFlux[i].end(); ++ir) {
      (*ir)->setQpin(0.0);
      (*ir)->setGwaterChng(0.0);
    }
  }
  
#endif
}

/*************************************************************************
**
** Send data to overlapping flow nodes.
**
*************************************************************************/

void tGraph::sendOverlap() {

#ifdef PARALLEL_TRIBS
  int dsizeN = 0;
  for (int i = 0; i < numGlobalPart; i++) {
    dsizeN = 3 * upFlow[i].size(); 
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];
      int d = 0;

      // Pack data for upstream reach outlet nodes
      std::set<tCNode*>::iterator iup;
      for (iup = upFlow[i].begin(); iup != upFlow[i].end(); ++iup) {
        ndata[d++] = (*iup)->getQstrm();
        ndata[d++] = (*iup)->getNwtOld();
        ndata[d++] = (*iup)->getNfOld();
      }

      tParallel::send(i, OVERLAP, ndata, dsizeN);
    }
  }

#endif
}

/*************************************************************************
**
** Receive data from overlapping flow nodes.
**
*************************************************************************/

void tGraph::receiveOverlap() {

#ifdef PARALLEL_TRIBS

  int dsizeN = 0;
  for (int i = 0; i < numGlobalPart; i++) {
   dsizeN = 3 * downFlow[i].size();
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];

      tParallel::receive(i, OVERLAP, ndata, dsizeN);
      int d = 0;

      // Unpack data from downstream reach head nodes
      std::set<tCNode*>::iterator idw;
      for (idw = downFlow[i].begin(); idw != downFlow[i].end(); ++idw) {

        (*idw)->setQstrm(ndata[d++]);
        (*idw)->setNwtOld(ndata[d++]);
        (*idw)->setNfOld(ndata[d++]);

      }
      delete [] ndata;
    }
  }
  tParallel::freeBuffers();// WR debug: put this at end of each receive call, it checks to see which previously assinged pointer  arrays can be safely deleted
#endif
}

/*************************************************************************
**
** Send data to downstream stream reach(s).
**
*************************************************************************/

void tGraph::sendQpin(int rid, tCNode* snode, double value) {
  assert(rid >= 0 && rid < numGlobalReach);

#ifdef PARALLEL_TRIBS
  //double* ndata = new double[1]; //WR debug added to scope of if statment to prevent memory leak
  // Get list of downstream reaches
  std::vector<int> dreach = conn[rid].getDownstream();

  // For each on another processor, send data
  for (int i = 0; i < dreach.size(); i++) {

    // Send from last reach
    if ( !inLocalPartition(dreach[i])) {
        double* ndata = new double[1];
        int to_proc = reach2partition[ dreach[i] ]; // To processor

      // Collect node data and send
      ndata[0] = value;
      tParallel::send(to_proc, QPIN+snode->getReach(), ndata, 1);
    }
  }
#endif
}

/*************************************************************************
**
** Receive data from upstream stream reach(s).
**
*************************************************************************/

void tGraph::receiveQpin(int rid, tCNode* rnode) {
  assert(rid >= 0 && rid < numGlobalReach);

#ifdef PARALLEL_TRIBS
  double* ndata = new double[1];
  // Get list of upstream reaches
  std::vector<int> ureach = conn[rid].getUpstream();

  // For each on another processor, receive data
  for (int i = 0; i < ureach.size(); i++) {

    // Check for last remote reach coming to this local reach
    if ( !inLocalPartition(ureach[i])) {

      int from_proc = reach2partition[ ureach[i] ]; // From processor

      // Receive data and update node
      tParallel::receive(from_proc, QPIN+rnode->getReach(), ndata, 1);
      rnode->addQpin(ndata[0]);
    }
  }
  delete [] ndata;
  tParallel::freeBuffers();// WR debug: put this at end of each receive call, it checks to see which previously assinged pointer  arrays can be safely deleted
#endif
}

/*************************************************************************
**
** Send initial data to overlapping flux nodes.
**
*************************************************************************/

void tGraph::sendInitial() {
#ifdef PARALLEL_TRIBS
  if (numGlobalPart == 1) return;

  int dsizeN = 0;
  for (int i = 0; i < numGlobalPart; i++) {
    dsizeN = 2 * localFlux[i].size() + 2 * upFlow[i].size();
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];
      std::set<tCNode*>::iterator il;
      int d = 0;

      for (il = localFlux[i].begin(); il != localFlux[i].end(); ++il) {
        int soilID = (*il)->getSoilID();
        double vArea = (*il)->getVArea();
        ndata[d++] = static_cast<double>(soilID);//*(reinterpret_cast<double*>(&soilID)); //WR debug converts int pointer to double pointer and derefs as double--possible source undefined behavior
        ndata[d++] = vArea;
      }

      for (il = upFlow[i].begin(); il != upFlow[i].end(); ++il) {
        int soilID = (*il)->getSoilID();
        double vArea = (*il)->getVArea();
        ndata[d++] = static_cast<double>(soilID);//*(reinterpret_cast<double*>(&soilID));//WR debug converts int pointer to double pointer and derefs as double--possible source undefined behavior
        ndata[d++] = vArea;
      }
      tParallel::send(i, INITIAL, ndata, dsizeN);
    }
  }
#endif
}

/*************************************************************************
**
** Receive initial data from overlapping flux nodes.
**
*************************************************************************/

void tGraph::receiveInitial() {
#ifdef PARALLEL_TRIBS
  if (numGlobalPart == 1) return;

  int dsizeN = 0;
  for (int i = 0; i < numGlobalPart; i++) {
    dsizeN = 2 * remoteFlux[i].size() + 2 * downFlow[i].size();
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];

      tParallel::receive(i, INITIAL, ndata, dsizeN);

      std::set<tCNode*>::iterator ir;
      int d = 0;
      for (ir = remoteFlux[i].begin(); ir != remoteFlux[i].end(); ++ir) {
          int soilID = static_cast<int>(ndata[d++]);//int soilID = *(reinterpret_cast<int*>(&ndata[d++]));//WR debug converts double pointer to int pointer and derefs as int--possible source undefined behavior
        (*ir)->setSoilID(soilID);
        (*ir)->setVArea(ndata[d++]);
      }

      for (ir = downFlow[i].begin(); ir != downFlow[i].end(); ++ir) {
          int soilID = static_cast<int>(ndata[d++]);// int soilID = *(reinterpret_cast<int*>(&ndata[d++])); //WR debug converts int pointer to double pointer and derefs as int--possible source undefined behavior
        (*ir)->setSoilID(soilID);
        (*ir)->setVArea(ndata[d++]);
      }
      delete [] ndata;
    }
  }
  tParallel::freeBuffers();// WR debug: put this at end of each receive call, it checks to see which previously assinged pointer  arrays can be safely deleted
#endif
}

/*************************************************************************
**
** Send runon flux data to downstream nodes.
**
*************************************************************************/

void tGraph::sendRunFlux(tCNode* cn) {
#ifdef PARALLEL_TRIBS
  int dsizeN = 2;
  int d;
  double* ndata = new double[dsizeN];

  // Get list of downstream reaches
  std::vector<int> dreach = conn[cn->getReach()].getDownstream();

  // For each on another processor, send data
  for (int i = 0; i < dreach.size(); i++) {

    // Send from last reach
    if ( !inLocalPartition(dreach[i])) {

      int to_proc = reach2partition[ dreach[i] ]; // To processor

       d = 0;
       ndata[d++] = cn->getSrf();
       ndata[d++] = cn->getVArea();
       tParallel::send(to_proc, RUNON+cn->getReach(), ndata, dsizeN);

   }
 }

  delete [] ndata;
#endif
}

/*************************************************************************
**
** Receive runon flux data from upstream nodes.
**
*************************************************************************/

void tGraph::receiveRunFlux(tCNode* cn) {
#ifdef PARALLEL_TRIBS
  int dsizeN = 2;
  int d;
  double* ndata = new double[dsizeN];
  int creach = cn->getReach();

  // Get list of upstream reaches
  std::vector<int> ureach = conn[cn->getReach()].getUpstream();

  // For each on another processor, receive data
  for (int i = 0; i < ureach.size(); i++) {

    // Check for last remote reach coming to this local reach
    if ( !inLocalPartition(ureach[i])) {

      int from_proc = reach2partition[ ureach[i] ]; // From processor

      // Receive data and update node
      d = 0;
      tParallel::receive(from_proc, RUNON+nodeAboveOutlet[ureach[i]]->getReach(), ndata, dsizeN);
      nodeAboveOutlet[ureach[i]]->setsrf(ndata[d++]);
      nodeAboveOutlet[ureach[i]]->setVArea(ndata[d++]);

    }
  }
  delete [] ndata;
  tParallel::freeBuffers();// WR debug: put this at end of each receive call, it checks to see which previously assinged pointer  arrays can be safely deleted
#endif
}

/*************************************************************************
**
** Send data to upstream overlapping nodes.
**
*************************************************************************/

void tGraph::sendUpstreamFlow() {
#ifdef PARALLEL_TRIBS
  int dsizeN = 0;
  int d;
  for (int i = 0; i < numGlobalPart; i++) {
    dsizeN = upFlow[i].size();
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];
      d = 0;
      // Pack up flow data
      std::set<tCNode*>::iterator iup;
      for (iup = upFlow[i].begin(); iup != upFlow[i].end(); ++iup) {
        ndata[d++] = (*iup)->getQstrm();
      }
      tParallel::send(i, UPSTREAM, ndata, dsizeN);
    }
  }
#endif
}

/*************************************************************************
**
** Receive data from downstream overlapping nodes.
**
*************************************************************************/

void tGraph::receiveDownstreamFlow() {
#ifdef PARALLEL_TRIBS
  int dsizeN = 0;
  int d;
  for (int i = 0; i < numGlobalPart; i++) {
    dsizeN = downFlow[i].size();
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];
      tParallel::receive(i, UPSTREAM, ndata, dsizeN);
      d = 0;
      // Unpack flow data
      std::set<tCNode*>::iterator idw;
      for (idw = downFlow[i].begin(); idw != downFlow[i].end(); ++idw) {
        (*idw)->setQstrm(ndata[d++]);
      }
      delete [] ndata;
    }
  }
  tParallel::freeBuffers();// WR debug: put this at end of each receive call, it checks to see which previously assinged pointer  arrays can be safely deleted
#endif
}

/*************************************************************************
**
** Send groundwater data from remote to local flux nodes.
**
*************************************************************************/

void tGraph::sendGroundWater() {

#ifdef PARALLEL_TRIBS
  int dsizeN = 0;
  for (int i = 0; i < numGlobalPart; i++) {
    dsizeN = 11 * remoteFlux[i].size();
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];
      int c = 0;

      // Pack data for remote saturated flux nodes
      std::set<tCNode*>::iterator iflux;
      for (iflux = remoteFlux[i].begin(); iflux != remoteFlux[i].end();
          ++iflux) {
        list<double>& rflist = (*iflux)->getGwaterChngList();
        int count = rflist.size();
        ndata[c++] = static_cast<double>(count);//*(reinterpret_cast<double*>(&count));//WR debug converts int pointer to double pointer and derefs as int--possible source undefined behavior
        list<double>::iterator iter;
        for (iter = rflist.begin(); iter != rflist.end(); ++iter) {
          ndata[c++] = (*iter);
        }
      }
      tParallel::send(i, GROUNDWATER, ndata, dsizeN);
    }
  }
#endif 
}

/*************************************************************************
**
** Receive groundwater data from remote to local flux nodes.
**
*************************************************************************/

void tGraph::receiveGroundWater() {

#ifdef PARALLEL_TRIBS
  int dsizeN = 0;
  for (int i = 0; i < numGlobalPart; i++) {
    dsizeN = 11 * localFlux[i].size();
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];
      tParallel::receive(i, GROUNDWATER, ndata, dsizeN);
      int c = 0;

      // Unpack data for local saturated flux nodes
      std::set<tCNode*>::iterator iflux;
      for (iflux = localFlux[i].begin(); iflux != localFlux[i].end();
          ++iflux) {
        int count = static_cast<int>(ndata[c++]);//*(reinterpret_cast<int*>(&ndata[c++]));//WR debug: static_cast is safer
        for (int j = 0; j < count; j++) {
          (*iflux)->addGwaterChng(ndata[c++]);
        }
      }
      delete [] ndata; //WR debug: ndata can be deleted here as tParallel:receive uses MPI_Recv which is a blocking process
    }
  }
    tParallel::freeBuffers();// WR debug: put this at end of each receive call, it checks to see which previously assinged pointer  arrays can be safely deleted
#endif
}

/*************************************************************************
**
** Send Nwt data to overlapping flow nodes.
**
*************************************************************************/

void tGraph::sendNwt() {

#ifdef PARALLEL_TRIBS
  int dsizeN = 0;
  for (int i = 0; i < numGlobalPart; i++) {
    dsizeN = 1 * localFlux[i].size();
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];
      int d = 0;
      std::set<tCNode*>::iterator iflux;
      for (iflux = localFlux[i].begin(); iflux != localFlux[i].end(); 
          ++iflux) {
        ndata[d++] = (*iflux)->getNwtOld();
      }
      tParallel::send(i, NWT, ndata, dsizeN);
    }
  }

#endif
}

/*************************************************************************
**
** Receive data from overlapping flow nodes.
**
*************************************************************************/

void tGraph::receiveNwt() {

#ifdef PARALLEL_TRIBS
  int dsizeN = 0;
  for (int i = 0; i < numGlobalPart; i++) {
   dsizeN = 1 * remoteFlux[i].size();
    if (dsizeN > 0) {
      double* ndata = new double[dsizeN];
      tParallel::receive(i, NWT, ndata, dsizeN);
      std::set<tCNode*>::iterator iflux;
      int d = 0;
      // Unpack flux data from downstream
      for (iflux = remoteFlux[i].begin(); iflux != remoteFlux[i].end(); 
          ++iflux) {
        (*iflux)->setNwtOld(ndata[d++]);
      }
      delete [] ndata;
    }
  }
    tParallel::freeBuffers();// WR debug: put this at end of each receive call, it checks to see which previously assinged pointer  arrays can be safely deleted
#endif
}

//=========================================================================
//
//
//                        End of tGraph.cpp
//
//
//=========================================================================
