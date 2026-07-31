/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/

/***************************************************************************
**
**  tInvariant.h:   Header File for tInvariant classes
**
**  tInvariant class is used for soil and land use attributes by defining
**  two Generic_Data and Data_Type classes.
**
***************************************************************************/

#ifndef  TINVARIANT_H
#define  TINVARIANT_H

#include "src/Headers/Inclusions.h"

using namespace std;

class SoilType;
class LandType;   

//=========================================================================
//
//
//                  Section 1: tInvariant Generic Class Declarations
//
//
//=========================================================================

class GenericSoilData 
{
public:
    GenericSoilData(tMesh<tCNode> *, tInputFile *, tResample *);
    ~GenericSoilData();
    SoilType **SoilClass;
    
    int numClass;    			// Number of soil classes
    int currClass;   			// Indicates current class
    char soilGrid[kMaxNameSize];  	// Pathname to soil map
    char soilTable[kMaxNameSize]; 	// Pathname to soil reference table
    char scfile[kMaxNameSize]; 	// Pathname to soil grid file // Giuseppe 2016
    
    void setSoilPtr(int soilID);
    void printSoilPars();
    void SetSoilParameters(tMesh<tCNode>*, tResample*, tInputFile&, int);
    double getSoilProp(int);

    // WR debug convert read grid flow to use smart pointes
    // Define a typedef for convenience
    using CharArrayPtr = std::unique_ptr<char[]>;

    // Create vectors of unique_ptr to manage your character arrays
    vector<CharArrayPtr> SCgridParamNames;
    vector<CharArrayPtr> SCgridBaseNames;
    vector<CharArrayPtr> SCgridExtNames;
    vector<CharArrayPtr> SCgridName;

    // Changes added by Giuseppe in August 2016 to allow reading soil grids WR debug converted to smart pointers to minimize memory leak
    // char **SCgridParamNames, **SCgridBaseNames, **SCgridExtNames, **SCgridName;
	
	
};

class GenericLandData 
{      
public:
  GenericLandData(tMesh<tCNode> *, tInputFile *, tResample *);
  ~GenericLandData();

  LandType **LandClass;   
  int numClass;      
  int currClass;     
  char landGrid[kMaxNameSize];  
  char landTable[kMaxNameSize]; 

  void setLandPtr(int landID);
  void SetLtypeParameters(tMesh<tCNode>*, tResample*, tInputFile&, int); 
  double getLandProp(int);
};


//=========================================================================
//
//
//                  Section 2: tInvariant Type Class Declarations
//
//
//=========================================================================

class SoilType 
{
public:
  SoilType();
  SoilType(double *, int);
  ~SoilType();

  int numProps;       	// Number of soil properties to store
  double *sProperty;    // Array containing soil properties

  void   setProperty(int, double);
  double getProperty(int);
};


class LandType 
{         
public:     
  LandType();
  LandType(double *, int);
  ~LandType();

  int numProps;         // Number of land properties to store
  double *lProperty;    // Array containing land properties
  
  void   setProperty(int, double);
  double getProperty(int);
};

#endif

//=========================================================================
//
//
//                         End of tInvariant.h 
//
//
//=========================================================================
