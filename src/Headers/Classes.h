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
**  Classes.h: General tRIBS Class Declarations
**
***************************************************************************/

#ifndef CLASSES_H
#define CLASSES_H

//=========================================================================
//
//
//                  Section 1: Classes 
//
//
//=========================================================================

class tTriangle;
class tNode;
class tEdge;
class tInputFile;
class tRunTimer;
class tCNode;
class tHydroModel;
class tStorm;
class tRainGauge;
class tRainfall;
class tEvapoTrans;
class tIntercept;
class tHydroMet;
class tHydroMetStoch;
class tHydroMetConvert;
class tWaterBalance;
class tSnowPack; // SKY2008Snow from AJR2007
class tSnowIntercept; // SKY2008Snow from AJR2007
class tShelter; // SKY2008Snow from AJR2007
class GenericSoilData;
class SoilType;
class GenericLandData;
class LandType;
class tResample;
class vCell; 
class tVariant;
class Simulator;
class SimulationControl; 
class tFlowNet;
class tFlowResults;
class Predicates;
class Point2D;
class Point3D;
class tPreProcess;

//=========================================================================
//
//
//                  Section 2: Templated Classes
//
//
//=========================================================================

template< class tSubNode > class tMesh;
template< class T > class tArray;
template< class T > class tMatrix;
template< class NodeType > class tListNode;
template< class NodeType > class tList;
template< class NodeType > class tMeshList;
template< class NodeType > class tPtrListNode;
template< class NodeType > class tPtrList;
template< class NodeType > class tListIter;
template< class NodeType > class tMeshListIter;
template< class NodeType > class tPtrListIter;
template< class tSubNode > class tListInputData;
template< class tSubNode > class tOutput;
template< class tSubNode > class tCOutput;

#endif

//=========================================================================
//
//
//                        End of Classes.h
//
//
//=========================================================================
