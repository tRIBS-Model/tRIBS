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
**  tEvapoTrans.h:   Header file for tEvapoTrans Class
**
**  This class encapsulates the evapotranspiration routines necessary for
**  interstorm computations. Three ET estimation methods are currently
**  implemented based on the Penman-Monteith equation,the Deardorff equation
**  and the Priestly-Taylor equation. Alternatively, specified from pan 
**  evaporation data.
**
***************************************************************************/

#ifndef EVAPOTRANS_H
#define EVAPOTRANS_H

#include "src/Headers/Inclusions.h"
#include "src/tRasTin/tRainfall.h"

class tRainfall;

//=========================================================================
//
//
//            Section 1: tEvapoTrans Class Declaration
//
//
//=========================================================================

class tEvapoTrans
{
 public:
  tEvapoTrans();
  tEvapoTrans(SimulationControl*,tMesh<tCNode> *, tInputFile &, tRunTimer *, 
	tResample *, tHydroModel *, tRainfall *);
  ~tEvapoTrans();
  SimulationControl *simCtrl;    // Pointer to simulation control

  void Debug(int, int);
  void DeleteEvapoTrans(); 
  void deleteLUGrids();
  void SetEvapTVariables(tInputFile &, tHydroModel *);
  void CreateHydroMetAndLU(tInputFile &); // SKYnGM2008LU
  void SetSunVariables();
  void SetEnvironment();
  void callEvapoTrans(tIntercept *, int);
  void callEvapoPotential();
  void initializeVariables();
  void assignStationToNode();
  void resampleGrids(tRunTimer *);
  void setCoeffs(tCNode *);
  void setTime(int);
  void readHydroMetData(int);
  void readHydroMetStat(char*);
  void readHydroMetGrid(char*);
  void readLUGrid(char*); // SKYnGM2008LU: added by AJR 2007
  void newHydroMetData(int);
  void newHydroMetStochData(int);
  void newHydroMetGridData(tCNode *);
  void newLUGridData(tCNode *); // SKYnGM2008LU: added by AJR 2007
  void createVariant();
  void createVariantLU(); // SKYnGM2008LU: added by AJR 2007
  void EvapPenmanMonteith(tCNode *);
  void EvapDeardorff(tCNode *);
  void EvapPriestlyTaylor(tCNode *);
  void EvapPan();
  void robustNess(double *, int);
  void setToNode(tCNode *);
  void betaFunc(tCNode *);
  void betaFuncT(tCNode *);
  void ComputeETComponents(tIntercept *, tCNode *, int, int);
  void DeriveAspect();
  void FunctionAndDerivative(tCNode *, double, double &, double &,  double);
  void HeatTransferProperties(tCNode *);
  void initialLUGridAssignment();
  void LUGridAssignment();
  void interpolateLUGrids(tCNode *);
  void constantLUGrids(tCNode* );
  void integratedLUVars(tCNode *, double);

  int  getEToption();
  int  julianDay();
  void DirectDiffuse(double);
  double latentHeat();
  double clausClap();
  double satVaporPress();
  double satVaporPress(double);
  double vaporPress();
  double psychoMetric();
  double totalPress();
  double densityMoist();
  double aeroResist();
  double stomResist();
  double inLongWave(tCNode *);
  double inShortWave(tCNode *);
  double energyBalance(tCNode*);

  // SKY2008Snow from AJR2007
  double compSkyCover();//find sky cover if marker is there -- RINEHART 2007 @ NMT
  double aboveHorizon(int);//check if sun is above horizon -- RINEHART 2007 @ NMT
  double ForceRestore(double, int); 
  double ComputeHourAngle(double, double);
  double ApproximateEP();
  double getinShortWave() const;
  double getinLongWave() const;
  double getoutLongWave() const;
  double getNetRad() const;
  double getDeltaT() const;
  double getSunRiseHour() const;
  double getSunSetHour() const;
  double getDayLength() const;
  double getDeltaAngle() const;
  double getPhiAngle() const;
  double getTauAngle() const;
  double rtsafe_mod_energy(tCNode*, double, double, double, 
			   double, double, double *, int *);

  void writeRestart(fstream &) const;
  void readRestart(fstream &);

  tHydroMetStoch *weatherSimul;

  // SKYnGM2008LU
  void SetGridTimeInfoVariables(tVariant *, char *);

 protected:
  tMesh<tCNode> *gridPtr;
  tResample *respPtr;
  tRunTimer *timer;
  tRainfall *rainPtr;
  tHydroMet *weatherStations;
  tHydroModel *hydrPtr;
  GenericLandData *landPtr;
  GenericSoilData *soilPtr;
  char stationFile[kName]{};

  char luFile[kName]{}; //SKYnGM2008LU: added by AJR 2007

  int *currentTime, *assignedStation;
  int VerbID{};
  int evapotransOption{}, metdataOption{}, Ioption{};
  int snowOption{}; // NEW FOR SNOW.... SKY2008Snow from AJR2007
  int shelterOption{}; // NEW FOR SHELTERING... SKY2008Snow from AJR2007
  int gFluxOption{}, dewHumFlag{}, ID{};
  int gmt{}, nodeHour{}, thisStation{}, oldTimeStep{};
  int numStations{}, arraySize{}, hourlyTimeStep{}, nParm{}, gridgmt{};
  int LUgridgmt{}; //SKYnGM2008LU: added by AJR 2007
  int vapOption{}, tsOption{}, nrOption{};
  int luOption{}, nParmLU;  //SKYnGM2008LU: added by AJR 2007
  int luInterpOption{};  //SKYnGM2008LU
  int timeCount;

  double timeStep{}, gridlat{}, gridlong{};
  double IfNotFirstTStepLU{}; //SKYnGM2008LU: added by AJR 2007
  double metHour{}, etHour{}, rainInt{}; // SKY2008Snow from AJR2007
  double LUgridlat{}, LUgridlong{}; //SKYnGM2008LU: added by AJR 2007

  double coeffH{}, coeffKt{}, coeffAl{}, coeffRs{}, coeffV{};
  double coeffKs{}, coeffCs{}, coeffPan{};
  // SKY2008Snow from AJR2007
  double coeffLAI{};
  double Rah{}, Rstm{};
  double SoilHeatCondTh{}, SoilHeatCpctTh{}, SoilHeatDiffTh{};
  double potEvap{}, actEvap{}, panEvap{}, betaS{}, betaT{};
  double airTemp{}, dewTemp{}, surfTemp{}, Tso{}, Tlo{};
  double rHumidity{}, atmPress{}, windSpeed{}, skyCover{};
  double netRad{}, latitude{}, longitude{}, vPress{};
  double inLongR{}, inShortR{}, outLongR{};
  double Tlinke{};
  double Is{}, Ic{}, Ics{}, Id{}, Ids{};
  double elevation{}, slope{}, aspect{};
  double atmPressC{}, surfTempC{}, skyCoverC{}, rHumidityC{}, dewTempC{};
  double windSpeedC{}, netRadC{}, vPressC{}, gFlux{}, hFlux{}, lFlux{}, Epot{}, rain{}, Gso{};
  double Io{}, alphaD{}, sinAlpha{}, del{}, phi{}, tau{}, circ{}, sunaz{};
  double SunRisHrLoc{}, SunSetHrLoc{}, DayLength{}, deltaT{};
  double RadDirObs{}, RadDifObs{};
  // CJC2025: New parameters
  double coeffSE{}, coeffST{};
  // SKY2008Snow from AJR2007
  //new for sheltering algorithm
  //	RMK: THE HA* VARIABLES SHOULD ACTUALLY BE HANDLED IN AN ARRAY
  //	BUT I COULDN'T GET THE INTERACTION W/ A POINTER IN TCNODE TO
  //	WORK CORRECTLY
  //
  //	RINEHART 2007 @ NEW MEXICO TECH
  double shelterFactorGlobal{}, landRefGlobal{};
  double horizonAngle{};
  double ha0000{}, ha0225{}, ha0450{}, ha0675{}, ha0900{}, ha1125{};
  double ha1350{}, ha1575{}, ha1800{}, ha2025{}, ha2250{}, ha2475{};
  double ha2700{}, ha2925{}, ha3150{}, ha3375{};
  //information for lapse rates
  //  RINEHART 2007 @ NEW MEXICO TECH
  double tempLapseRate{}; //K/m -- make sure that time steps are consistent
  //for output of cumulative number of hours of sunlight
  //  RINEHART 2007 @ NEW MEXICO TECH
  double SunHour{};

  char **gridParamNames, **gridBaseNames, **gridExtNames;
  char **LUgridParamNames, **LUgridBaseNames, **LUgridExtNames; // SKYnGM2008LU: added by AJR 2007

  tVariant *airpressure, *dewtemperature, *skycover, *windspeed;
  tVariant *airtemperature, *surftemperature, *netradiation, *incomingsolar; //E.R.V 3/6/2012
  tVariant *evapotranspiration, *relhumidity, *vaporpressure;

  // SKYnGM2008LU: added by AJR 2007
  tVariant *LandUseAlbGrid, *ThroughFallGrid, *VegHeightGrid; 
  tVariant *StomResGrid, *VegFractGrid; 

  // SKYnGM2008LU
  tVariant *CanStorParamGrid; 
  tVariant *IntercepCoeffGrid, *CanFieldCapGrid, *DrainCoeffGrid;
  tVariant *DrainExpParGrid, *OptTransmCoeffGrid, *LeafAIGrid;
  tVariant *EvapThreshGrid, *TransThreshGrid; // CJC2025

  // SKYnGM2008LU
  int numALfiles{}, numTFfiles{}, numVHfiles{}, numSRfiles{};
  int numVFfiles{}, numCSfiles{}, numICfiles{}, numCCfiles{};
  int numDCfiles{}, numDEfiles{}, numOTfiles{}, numLAfiles{};
  int numSEfiles{}, numSTfiles{}; // CJC2025 
  int *ALgridhours, *TFgridhours, *VHgridhours, *SRgridhours;
  int *VFgridhours, *CSgridhours, *ICgridhours, *CCgridhours;
  int *DCgridhours, *DEgridhours, *OTgridhours, *LAgridhours;  
  int *SEgridhours, *STgridhours; // CJC2025
  int NowTillWhichALgrid{}, NowTillWhichTFgrid{}, NowTillWhichVHgrid{}, NowTillWhichSRgrid{};
  int NowTillWhichVFgrid{}, NowTillWhichCSgrid{}, NowTillWhichICgrid{}, NowTillWhichCCgrid{};
  int NowTillWhichDCgrid{}, NowTillWhichDEgrid{}, NowTillWhichOTgrid{}, NowTillWhichLAgrid{};
  int NowTillWhichSEgrid{}, NowTillWhichSTgrid{}; // CJC2025
  char **ALgridFileNames, **TFgridFileNames, **VHgridFileNames, **SRgridFileNames;
  char **VFgridFileNames, **CSgridFileNames, **ICgridFileNames, **CCgridFileNames;
  char **DCgridFileNames, **DEgridFileNames, **OTgridFileNames, **LAgridFileNames;
  char **SEgridFileNames, **STgridFileNames; // CJC2025
  int AtFirstTimeStepLUFlag{};

  int skycover_flag; // intended for when nodata is set for XC gridded data so that skycover is estimated.

};

inline double tEvapoTrans::getDeltaT()      const {return deltaT;}
inline double tEvapoTrans::getSunRiseHour() const {return SunRisHrLoc;}
inline double tEvapoTrans::getSunSetHour()  const {return SunSetHrLoc;}
inline double tEvapoTrans::getDayLength()   const {return DayLength;}
inline double tEvapoTrans::getDeltaAngle()  const {return del;}
inline double tEvapoTrans::getPhiAngle()    const {return phi;}
inline double tEvapoTrans::getTauAngle()    const {return tau;}
inline double tEvapoTrans::getinShortWave() const {return inShortR;}
inline double tEvapoTrans::getinLongWave()  const {return inLongR;}
inline double tEvapoTrans::getoutLongWave() const {return outLongR;}
inline double tEvapoTrans::getNetRad()      const {return netRad;}


#endif

//=========================================================================
//
//
//                    End of tEvapoTrans.h
//
//
//=========================================================================
