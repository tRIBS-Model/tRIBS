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
**  tFlowNet.h: Header for class tFlowNet (see tFlowNet.cpp) based on CHILD
**             and RIBS routines for Stream Net Routing
**
***************************************************************************/

#ifndef TFLOWNET_H
#define TFLOWNET_H

//=========================================================================
//
//
//                  Section 1: tFlowNet Include Statements
//
//
//=========================================================================

#include "src/Headers/Definitions.h"
#include "src/Headers/Classes.h"
#include "src/tMesh/tMesh.h"
#include "src/tMeshElements/meshElements.h"
#include "src/tInOut/tInputFile.h"
#include "src/tCNode/tCNode.h"
#include "src/tSimulator/tRunTimer.h"
#include "src/tFlowNet/tFlowResults.h"

//=========================================================================
//
//
//                  Section 2: tFlowNet Class Definitions
//
//
//=========================================================================

class tFlowResults;

class tFlowNet
{
public:
  tFlowNet();
  tFlowNet(SimulationControl*, tMesh< tCNode> *, tInputFile &, tRunTimer *);
  ~tFlowNet();

  int MaxTravel();
  int FindLakeNodeOutlet(tCNode *);
  int IsConnected(tCNode*, tCNode*);  
  int IsConfluence(tCNode*, tCNode*);           
  int IsStreamHead(tCNode*);                
  int IsStreamDisjoint(tEdge*,tCNode*,tPtrList<tEdge> &,tPtrList<tEdge> &);
  int IsToEliminate(tCNode*);

  int CheckNeighbor(tCNode*, tEdge*,tPtrList<tCNode> &,tPtrList<tEdge> &);
  int FindStreamDisjoints(tCNode*, int, tPtrList<tCNode> &);
  int FindConfluenceDisjoints(tPtrList<tCNode> &);

  void SetFlowVariables(tInputFile &);
  void SetBasinOutlet();
  void CalcSlopes();
  void InitFlowDirs();
  void FlowDirs();
  void SortNodesByNetOrder();
  void FillLakes();
  void initializeTravelTime();
  void initializeTravelTimeOnly();     
  void setTravelTime();
  void setMaxTravelTime();
  void setTravelVelocity(double);
  void SurfaceFlow();
  void DrainAreaVoronoi();
  void RouteFlowArea(tCNode *, double);
  void DeriveStreamReaches(tInputFile &); 
  void ComputeDistanceToStream();         
  void SortStreamNodes();                 
  void TellAboutNode(tCNode *cn);       
  void PrintArcInfoLinks(tInputFile &);  
  void WeightedShortestPath(); 
  void AddUnsettledNeighbors(tCNode*,tPtrList<tCNode> &,tPtrList<tEdge> &);
  void UpdatePathVariable(tCNode*);
  void CheckVDrainageWidths();
  void FixVoronoiEdgeWidth(tCNode *);
  void DeriveCurvature();

  void ReadFlowNetFromMeshBuilder();

  int IsBetweenEndPnts(tArray<double> &,tArray<double> &,
		       tArray<double> &,tArray<double> &, 
		       double, double);
  int AreSegmentsParallel(tArray<double> &,tArray<double> &,
			  tArray<double> &,tArray<double> &);
  int IsInTriangle(tArray< double > &,tArray< double > &, 
		   tArray< double > &, double, double);
  double ComputeEdgeWeight(tEdge*, double);
  double FindAngle(tCNode*, tCNode*, tCNode*);     
  double FindDistance(double, double, double, double); 
  double CompElevationTendency(tCNode*, tCNode*, int*, int);   
  double polygonArea(double **, int);
  double getCurrDischarge(int);

  void SetReachInformation();
  void writeRestart(fstream &) const;
  void readRestart(fstream &);

  tPtrList< tCNode >& getReachHeadList() { return NodesLstH; }
  tPtrList< tCNode >& getReachOutletList() { return NodesLstO; }
  tList< int >& getReachSizeList() { return NNodes; }


  tFlowResults* getResultsPtr() { return res; }
  tCNode*       getOutletPtr()  { return OutletNode;}

  SimulationControl *simCtrl;   

protected:
  tMesh<tCNode>  *gridPtr;
  tRunTimer      *timer;
  tFlowResults   *res;
  tCNode         *OutletNode;    // Ptr to the outlet node, 1 so far...
  tPtrList< tCNode > HeadsLst;   // List of basin stream heads
  tPtrList< tCNode > NodesLstH;  // NODELIST #1 Heads of stream reaches 
  tPtrList< tCNode > NodesLstO;  // NODELIST #2 Outlets of stream reaches 
  tList< int >       NNodes;     // NODELIST #3 # of nodes in each reach 

  int flowboxes; 	        // Size of current discharge array
  int maxttimeInitial;		// added by Ara Ko in 2017 to keep the minitial maxttime 
  double hillvel;   	        // Hillslope velocity, [m/sec]
  double streamvel; 	        // Stream velocity, [m/sec]
  double velratio;  	        // Ratio Stream/Hillslope
  double velcoef;   	        // Coeffcient velocity-discharge
  double flowexp;   	        // Power in velocity-discharge relationship
  double baseflow;  	        // Baseflow: leftover from event-based
  double dOtp;      	        // Time interval in hydrograph output, [hour]
  double timespan;  	        // Total runtime during simulation, [hour]
  double flowout;   	        // Discharge at the outlet, [m^3/sec]
  double maxttime;  	        // MAX travel time defined for watershed, [hour]
  double dist_hill_max;  	// MAX distance on hillslope, [m] 
  double dist_stream_max;	// MAX distance in stream, [m]
  double BasArea;               // Total Basin Area, [m^2]
  int percolationOption;	// ASM percolation option
};

#endif
	
//=========================================================================
//
//
//                          End tFlowNet.h
//
//
//=========================================================================
