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
**  tStorm.h: Header for class tStorm
**
**  A tStorm object generates random storms assuming an exponential
**  distribution of rainfall intensity, storm duration, and time to the
**  next storm (Eagleson, 1978b). Its services include reading the 
**  necessary parameters from a tInputFile, generating a new      
**  storm, and reporting its various values. Inherited from CHILD with
**  minor modifications. 
**
**  If you want to provide an option for NOT having storms vary
**  randomly, you can do so by setting optVariable to zero on initialization.
**
**  The GammaDev() function is provided for future reference; it is not
**  actually used in the current version.
**    
**  The storm parameters can also be varied sinusoidally with time to 
**  simulate long-term climatic fluctuations (GenerateStorm)
**
**************************************************************************/

#ifndef TSTORM_H
#define TSTORM_H

//=========================================================================
//
//
//                  Section 1: tStorm Include Statements
//
//
//=========================================================================

#include "src/Headers/Classes.h"
#include "src/Mathutil/mathutil.h"
#include "src/tInOut/tInputFile.h"

//=========================================================================
//
//
//                  Section 2: tStorm Class Definitions
//
//
//=========================================================================

class tStorm
{
public:
  tStorm();
  tStorm( SimulationControl*, tInputFile & );
  ~tStorm();
  SimulationControl *simCtrl;    

  void   SetStormVariables(tInputFile &);  
  void   GenerateStorm( double tm, double minp=0.0, double mind=0.0);
  int    getoptStorm() const;  
  int    DefineSeason(int);
  int    getSeasonID();
  double getStormDuration() const;
  double interstormDur() const;
  double getRainrate() const;
  double getMeanStormDur() const;
  double getMeanInterstormDur() const;
  double getMeanPrecip() const;

  void  updateSeasonVars();
  void  setRainrate(double );
  void  setSeasonID(int);
  void  setSeasonMonth(int, int, int); 
  void  setSeasonPMean(int, double);
  void  setSeasonStmDurMean(int, double); 
  void  setSeasonSplDurMean(int, double); 
  void  allocSeasonMemory(int);
  void  writeRestart(fstream &) const;
  void  readRestart(fstream &);

  int    RealRn;
  int    rid;
  int    *stbeg;    // Only when REAL rain drives all other vars
  int    *stend;   
  double *rainIN;   
   
private:
  double ExpDev( long * );
   
  int optStoch;       // Flag for stochastic mode
 
  double stdurMean;   // Mean duration
  double istdurMean;  // Mean time between storms
  double pMean;       // Mean rainfall intensity
  double p;           // Actual rainfall intensity for the current storm
  double stdur;       // Actual storm duration
  double istdur;      // Actual time between storms
  double p0;          // Climatological mean: the "weather" means can 
  double stdur0;      //  themselves vary over geologic time; p0, etc, are
  double istdur0;     //  the means of the means so to speak.
  double pdev;        // Absolute magnitude of deviation from p0, etc, under
  double stdurdev;    //  sinusoidal variation.
  double istdurdev;
  double twoPiLam;    // Parameter for sinusoidal variation: 2pi / period
  long   seed;        // Random seed
  double endtm;       // The end time of the run
                 
  int currSeasID;
  int NumSeas;
  int **SeasMo;
  double *MStmDr;
  double *MSplDr;
  double *MRate;

  ofstream stormfile; // File containing history of storm events
};


#endif

//=========================================================================
//
//
//			End of tStorm.h                   
//
//
//=========================================================================
