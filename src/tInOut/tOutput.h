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
**  tOutput.h: Header for the tOutput class 
**
***************************************************************************/

#ifndef TOUTPUT_H
#define TOUTPUT_H

//=========================================================================
//
//
//                  Section 1: tOutput Include and Define Statement
//
//
//========================================================================= 

#include "src/Headers/Definitions.h"
#include "src/tMesh/tMesh.h"
#include "src/tMeshList/tMeshList.h"
#include "src/tMeshElements/meshElements.h"
#include "src/tInOut/tInputFile.h"
#include "src/tSimulator/tRunTimer.h"
#include "src/tRasTin/tResample.h"

using namespace std;

class tResample;

//=========================================================================
//
//
//                  Section 2: tOutput Class Definition
//
//
//=========================================================================

template< class tSubNode > 
class tOutput 
{
public:
  tOutput(SimulationControl*, tMesh<tSubNode>*, tInputFile&, tResample*);
  tOutput(SimulationControl*, tMesh<tSubNode>*, tInputFile&, 
	  tResample*, tRunTimer*);
  virtual ~tOutput();
  SimulationControl *simCtrl;    

  void WriteOutput(double);
  void CreateAndOpenFile(ofstream*, char*);
  void CreateAndOpenVizFile(ofstream*, char*);
  void ReadNodeOutputList();
  void CreateAndOpenPixel();
  void CreateAndOpenDynVar();
  void end_simulation();
  void SetInteriorNode();
 
  virtual void WriteDynamicVars(double); 
  virtual void WriteDynamicVarsBinary(double); 
  virtual void WriteDynamicVar(double);
  virtual void WriteNodeData(double);    
  virtual void WriteNodeData(double, tResample*);
  virtual void WriteGeometry(tResample*);
  virtual void WritePixelInfo(double);


  int numNodes;
  int *nodeList;

protected:
  tMesh<tSubNode> *g;        
  tRunTimer *timer;    
  tSubNode **uzel;     
  tResample *respPtr;

  char baseName[kMaxNameSize]; 
  char nodeFile[kMaxNameSize];
  char vizName[kMaxNameSize]; 

  int vizOption;
   
  ofstream nodeofs;
  ofstream edgofs;
  ofstream triofs;
  ofstream zofs;
  ofstream *pixinfo;
  ofstream *dynvars;              
};

//=========================================================================
//
//
//                  Section 3: tCOutput Class Definition
//
//
//=========================================================================

template< class tSubNode > 
class tCOutput : public tOutput<tSubNode>
{
public:
  tCOutput(SimulationControl*, tMesh<tSubNode>*, tInputFile&, tResample*);
  tCOutput(SimulationControl*, tMesh<tSubNode>*, tInputFile&,
	      tResample*, tRunTimer*);
  ~tCOutput();
  SimulationControl *simCtrl;    

  int  numOutlets;
  int *OutletList;

  char outletName[kMaxNameSize];

  void WriteDynamicVars(double);
  void WriteDynamicVarsBinary(double);
  void WriteDynamicVar(double);
  void WriteIntegrVars(double);
  void WritePixelInfo(double);
  void WriteNodeData(double);
  void WriteNodeData(double, tResample*); 
  void WriteGeometry(tResample*);
  void UpdateForNewRun(tInputFile &); 

  void WriteOutletInfo(double);
  void ReadOutletNodeList(char *);
  void CreateAndOpenOutlet();
  void SetInteriorOutlet();

private:
  tSubNode **Outlets;    //Pointer to an array of tCNode objects
  ofstream *outletinfo;     
  ofstream arcofs;
  ofstream vorofs;
  ofstream intofs;
  ofstream drareaofs;
  ofstream widthsofs;
};

#endif

//=========================================================================
//
//
//                          End of tOutput.h 
//
//
//=========================================================================
