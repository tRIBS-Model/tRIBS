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
**  tVariant.h:   Header file for class tVariant.cpp 
**
**  tVariant class is an generic class for time-varying meteorological input
**  grids based on the tRainfall class.
** 
***************************************************************************/

#ifndef  TVARIANT_H
#define  TVARIANT_H

#include "src/Headers/Inclusions.h"

using namespace std;

//=========================================================================
//
//
//                  Section 1: tVariant Class Declaration
//
//
//=========================================================================

class tVariant{
public:
  tVariant();
  tVariant(tMesh<tCNode> *, tResample *);

  void newVariable(char *);
  void updateVariable(char *);
  void setFileNames(char *, char*);
  void noData(char *);
  int  composeFileName(tRunTimer *);
  char fileIn[kName];

  // SKYnGM2008LU
  ifstream infile2;
  char *getInputName();
  char *getExtension();
  void updateLUVarOfPrevGrid(const char *, char *);
  void updateLUVarOfBothGrids(const char *, char *);

protected:
  tMesh<tCNode> *gridPtr;
  tResample *respPtr;
  char inputName[kName];
  char extension[20];
  ifstream infile; 
};

#endif

//=========================================================================
//
//
//                         End of tVariant.h
//
//
//=========================================================================
