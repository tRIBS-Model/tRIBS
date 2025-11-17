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
**  tShelter.h:   Header File for tInvariant classes
**
**  tShelter class is used for radiation sheltering and wind sheltering,
**  The radiation sheltering will be implemented for all energy balance
**  calculations involving net radiation if so desired. 
**
**  Currently, only the radiation sheltering algorithm is implemented. It is
**  entirely contained in the object constructor, as we only need to do this 
**  process once at the initialization of the model.
**
**  Currently, there are only 16 directions used to find horizon angles, which
**  will introduce error in larger regions. The algorithm is fairly straight
**  forward.
**
**    Initialize DEM array
**    Choose a direction
**    Pixelate and set up navigation variables
**    For each grid cell and in each direction, find steepest neigbor.
**    Resample finished grid to Voronoi cells.
**    Calculate sky view and land view factors.
**
**    RINEHART 2007 @ NEW MEXICO TECH
**  
**
***************************************************************************/

#ifndef  TSHELTER_H
#define  TSHELTER_H

#include "src/tRasTin/tResample.h"
#include "src/Headers/Inclusions.h"

//=========================================================================
//
//
//                  Section 1: tShelter Class Declarations
//
//
//=========================================================================

class tShelter : public tResample {
public:

  tShelter();
  tShelter( SimulationControl *, tMesh<tCNode> *, tInputFile &);
  ~tShelter();
  

  char GridInPath[50]; //file pathname of DEM
  double angleDiv; //angle by which we divide 360 degrees
  double horAngle;
  double maxTan, tempTan;
  double maxRow, maxCol;
  double maxZ, tempZ, initZ;
  double sv, lv;
  double slope, aspect, elevation;
  double **horAngleX; //pointer to temporary horizon angle grid
  double **tempGrid;
  int radSheltOpt, windSheltOpt;

};

#endif
