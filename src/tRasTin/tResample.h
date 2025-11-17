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
**  tResample.h:   Header File for tResample class
**
***************************************************************************/

#ifndef  TRESAMPLE_H
#define  TRESAMPLE_H

#include "src/Headers/Inclusions.h"

class vCell;

//=========================================================================
//
//
//                  Section 1: tResample Class Declarations
//
//
//=========================================================================

class tResample 
{
 public:
  tResample(SimulationControl*, tMesh<tCNode> *);
  ~tResample();

  tMesh<tCNode> *mew; // Pointer to tMesh 
  vCell *eta;         // Pointer to the base level unit object
  SimulationControl *simCtrl;  // Pointer to simulation control

  int NVor;          //  Total number of Voronoi cells in the basin
  int      *nPoints; //  Contains number of points for each Voronoi cell
  double   **vXs;    //  Ragged array containing voronoi Xs off ALL cells
  double   **vYs;    //  Ragged array containing voronoi Ys off ALL cells
                   
  double   *varFromGrid;  // Assigned Voronoi cell values from GRID source
  int      *varFromPoint; // Assigned Voronoi cell values from POINT source 

  int      NR;       //  Rows in the INPUT grid
  int      MR;       //  Columns in the INPUT grid
  double   **gridIn; //  Values of the grid

  double   dR;       //  Spatial resolution of the GRID 
  double   dummy;    //  NODATA value in the GRID
  double   xllcR;    //  X -LOWER left corner of the GRID 
  double   yllcR;    //  Y -LOWER left corner of the GRID 
  double   xulcR;    //  X -UPPER left corner of the GRID 
  double   yulcR;    //  Y -UPPER left corner of the GRID   
  double   *coorXG;  //  X -Coordinates of all points of the GRID 
  double   *coorYG;  //  Y -Coordinates of all points of the GRID 

  char   fileIn[30];
  char   fileOut[30];

  double* doIt(char *, int);  // Returns array 'varFromGrid'
  int*    doIt(int *, double *, double *, int); // Returns array 'varFromGrid'

  void    VerticesNoAccBndEff();     
  void    MakeBoundaryPolygons();    
  void    MakeInteriorPolygons();    
  void    FindNeighbValue(tCNode *);
  void    PrintEdgeInfo(tEdge *); 

  void    readInputGrid(char *);
  void    DestrtResample();
  void    allocMemory(double **, int, int);
  void    In_Mrain_Name (char *, char *, int, int, int, int);
  void    Out_Mrain_Name(char *, char *, int);
  void    printToFile(ofstream&, double **, int, int);
  int     IsInTriangle(tArray<double> &,tArray<double> &,tArray<double> &, double, double);
  void    FixVoronoiPolygon(tCNode *, tCNode *, tArray< double > &, tArray< double > &);
  tArray<double> FindNormal(double, double, double, double, double);

};

//=========================================================================
//
//
//                  Section 2: vCell Class Declarations
//
//
//=========================================================================

class vCell 
{
 public:
  vCell();
  vCell(SimulationControl*, tResample *, double *, double *, int);
  ~vCell();

  int    nv;         // Number of vertices of the Voronoi cell: CLOSED
  double *inX;       // X-s of the Voronoi vertices: NOT CLOSED polygon
  double *inY;       // Y-s of the Voronoi vertices: NOT CLOSED polygon
  double *VoronX;    // X-s of the Voronoi vertices: CLOSED polygon
  double *VoronY;    // Y-s of the Voronoi vertices: CLOSED polygon

  double xMax;      // X max-min & Y max-min of current voronoi cell
  double yMax;
  double xMin;
  double yMin;

  int    YminInd;   // Indices of *coorXY bounding box with Voronoi cell 
  int    YmaxInd;
  int    XminInd;
  int    XmaxInd;

  tResample *base;
  SimulationControl *simCtrl;  // Pointer to simulation control

  void initializeVCell(SimulationControl*, tResample *, double *, 
		       double *, int);
  void initializeVCell(double *, double *, int);

  void DestrtvCell();
  void printCellData();

  void getClosedPolygon(double *, double *, int);
  void findMaxMin();
  void findCorrespInd(double, double, int *, int *);
  void defineTwoSubArrays(double **, double **, int *, int *,
                     double *, double *, int, double, int);
  void findIntersectionPoint(double *, double *, int, 
                     double, double *, double *, int);

  // GMnSKY2008MLE to fix memory errors
  //int findDirectionIntersectionPoint(int, double, double, double *, double *, double *, double *, int);//AJR 2007
  void Ord_CounterClock(double **, int *); // new function - Giuseppe
  double AngleIncluded(double , double , double , double); // new function - Giuseppe
  
  int  getBoxXs(double);
  int  getBoxYs(double);
  int  polyCentroid(double *, double *, int, double *, double *, double *);

  double polygonArea(double **, int);
  double findDistance(double, double, double, double);
  double convertToVoronoiFormat(int flag);
  int    convertPointData(int *, double *, double *, int);

};

#endif

//=========================================================================
//
//
//                         End of tResample.h
//
//
//=========================================================================
