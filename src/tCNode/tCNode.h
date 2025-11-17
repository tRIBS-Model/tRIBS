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
**  tCNode.h: Header for tCNode.cpp 
**
**  tCNode class is derived class from tNode used to store dynamic 
**  hydrologic variables and data. Updating data members through get/set
**  functions is performed at model time steps.
**
***************************************************************************/

#ifndef TCNODE_H
#define TCNODE_H

//=========================================================================
//
//
//                  Section 1: tCNode Include and Define Statements
//
//
//=========================================================================

#include "src/tMeshElements/meshElements.h"
#include "src/tInOut/tInputFile.h"
#include "src/Headers/globalFns.h"
#include <memory> // WR - added 09192023 :)
#include <list> //SMM - added 09232008

#ifdef ALPHA_64
  #include <iostream>
  #include <fstream>
  #include <assert.h>
  #include <math.h>
#elif defined LINUX_32
  #include <iostream>
  #include <fstream>
  #include <cassert>
  #include <cmath>

#elif defined MAC
  #include <iostream>
  #include <fstream>
  #include <cassert>
  #include <cmath>

#elif defined WIN
  #include <iostream>
  #include <fstream>
  #include <assert.h>
  #include <math.h>
#else 
  #include <iostream>
  #include <fstream>
  #include <assert.h>
  #include <math.h>
#endif

//=========================================================================
//
//
//                  Section 2: tCNode Class Definition
//
//
//=========================================================================

class tCNode : public tNode
{
public:
  tCNode();
  tCNode( tInputFile &infile );
  ~tCNode();

  double getNwtOld();                  //Unsat and Sat Zone Members
  double getNwtNew();
  double getMuOld();
  double getMuNew();
  double getMiOld();
  double getMiNew();
  double getNtOld();
  double getNtNew();
  double getNfOld();
  double getNfNew();
  double getRuOld();
  double getRuNew();
  double getRiOld();
  double getRiNew();
  double getQin();                     //Runoff Members
  double getQout();
  double getQpin();
  double getQpout();
  double getRain();
  double getSrf_Hr();
  double getCumSrf(); // added CJC2021
  double getSrf();
  double getHsrf();
  double getPsrf();
  double getSatsrf();
  double getSbsrf();
  double getEsrf();
  double getRsrf();
  double getRunOn();                   //Runon
  double getGwaterChng();
  std::list<double>& getGwaterChngList();  //SMM added 10132008
  double getTTime();
  double getHillPath();                //Routing Members
  double getStreamPath(); 
  double getInterceptLoss();           //Interception Members
  double getNetPrecipitation(); 
  double getCanStorage();    
  double getPotEvap();        
  double getActEvap();        
  double getStormLength();      
  double getCumIntercept();       
  double getEvapoTrans();              //Evapotranspiration Members
  double getEvapWetCanopy();   
  double getEvapDryCanopy();    
  double getEvapSoil();         
  double getSoilMoisture();     
  double getSoilMoistureSC();     
  double getSoilMoistureUNSC();     
  double getRootMoisture();
  double getRootMoistureSC();
  double getTransmiss();               //Meteorological Members
  double getAirPressure();      
  double getDewTemp();	   
  double getRelHumid();   
  double getVapPressure();     
  double getSkyCover();         
  double getWindSpeed();       
  double getAirTemp();         
  double getSurfTemp();         
  double getSoilTemp();
  double getNetRad();           
  double getGridET(); 

  double getShortRadIn();
  double getShortRadIn_dir();
  double getShortRadIn_dif();
  double getShortAbsbVeg();
  double getShortAbsbSoi();
  double getShortRadSlope();    /// Get back the raw slope‐corrected irradiance, JB2025

  double getLongRadIn();
  double getLongRadOut();

  double getGFlux();
  double getHFlux();
  double getLFlux();
  double getGnod();

  double getQgwOut();
  double getQgwIn();

  // ASM 2/10/2017
  double getChannelPerc();
  double getFt(); //ASM
  double getPercOccur(); //ASM
  double getavPerc(); //ASM

  // SKY2008Snow from AJR2007
  //snowpack
  double getLiqWE();//state
  double getIceWE();//state
  double getDU();//state
  double getSnTempC();//state
  double getLiqRouted();//flux
  double getCrustAge();//state
  double getDensityAge();//state
  double getEvapoTransAge();//state
  double getSnLHF();//flux
  double getSnSHF();//flux
  double getSnSub();//flux // CJC2020
  double getSnEvap();//flux // CJC2020
  double getSnGHF();//flux
  double getSnPHF();//flux
  double getSnRLout();//flux
  double getSnRLin();//flux
  double getSnRSin();//flux
  double getUnode();//flux
  double getUerror();//flux
  double getCumLHF();//integrated output
  double getCumMelt();//integrated output
  double getCumSHF();//integrated output
  double getCumSnSub();//integrated output // CJC2020
  double getCumSnEvap();//integrated output // CJC2020
  double getCumTotEvap();//integrated output // CJC2020
  double getCumBarEvap();//integrated output // CJC2020
  double getCumPHF();//integrated output
  double getCumRLin();//integrated output
  double getCumRLout();//integrated output
  double getCumRSin();//integrated output
  double getCumGHF();//integrated output
  double getCumUerror();//integrated output
  double getPersTimeMax();//integrated output
  double getPersTime();//integratedoutput
  double getPeakSWE();//integratedoutput
  double getPeakSWETemp();//integratedoutput
  double getInitPackTime();//integratedoutput
  double getInitPackTimeTemp();
  double getPeakPackTime();//integratedoutput

  //snowintercept
  double getIntSWE();//state
  double getIntPrec();//flux
  double getIntSnUnload();//flux
  double getIntSub();//flux
  double getCumIntSub();//integrated output
  double getCumIntUnl();//integrated output

  //tShelter members
  double getHorAngle0000();
  double getHorAngle0225();
  double getHorAngle0450();
  double getHorAngle0675();
  double getHorAngle0900();
  double getHorAngle1125();
  double getHorAngle1350();
  double getHorAngle1575();
  double getHorAngle1800();
  double getHorAngle2025();
  double getHorAngle2250();
  double getHorAngle2475();
  double getHorAngle2700();
  double getHorAngle2925();
  double getHorAngle3150();
  double getHorAngle3375();
  double getSheltFact();
  double getLandFact();
  double getCumHrsSun();
  double getCumHrsSnow();

  double getIntStormVar();
  double getBedrockDepth();
  double getHlevel();
  double getQstrm();
  double getChannelWidth();
  double getRoughness();
  double getFlowVelocity();

  //double getCanopyStorage();           //tWaterBalance Members
  // SKYnGM2008LU: getCanopyStorage() now replaced by getCanopyStorVol() below. 
  //tWaterBalance Members
  double getCanopyStorVol(); // SKYnGM2008LU
  double getUnSaturatedStorage(); 
  double getSaturatedStorage();     
  double getRecharge(); 
  double getUnSatFlowOut();
  double getUnSatFlowIn();

  double getContrArea();
  double getCurvature();
  double getBasinArea();
  double getAspect();
  double getVegFraction();

  // SKYnGM2008LU
  // additions for landuse grids AJR 2007
  double getLandUseAlb();
  double getThroughFall();
  double getVegHeight();
  double getStomRes();

  // SKYnGM2008LU
  double getIntercepCoeff();
  double getCanFieldCap();
  double getDrainCoeff();
  double getDrainExpPar();
  double getOptTransmCoeff();
  double getLeafAI();
  double getCanStorParam();
  // CJC2025: New Parameters
  double getEvapThresh();
  double getTransThresh();

  // SKYnGM2008LU
  double getLandUseAlbInPrevGrid();
  double getLandUseAlbInUntilGrid();
  double getThroughFallInPrevGrid();
  double getThroughFallInUntilGrid();
  double getVegHeightInPrevGrid();
  double getVegHeightInUntilGrid();
  double getStomResInPrevGrid();
  double getStomResInUntilGrid();
  double getVegFractionInPrevGrid();
  double getVegFractionInUntilGrid();
  double getCanStorParamInPrevGrid();
  double getCanStorParamInUntilGrid();
  double getIntercepCoeffInPrevGrid();
  double getIntercepCoeffInUntilGrid();
  double getCanFieldCapInPrevGrid();
  double getCanFieldCapInUntilGrid();
  double getDrainCoeffInPrevGrid();
  double getDrainCoeffInUntilGrid();
  double getDrainExpParInPrevGrid();
  double getDrainExpParInUntilGrid();
  double getOptTransmCoeffInPrevGrid();
  double getOptTransmCoeffInUntilGrid();
  double getLeafAIInPrevGrid();
  double getLeafAIInUntilGrid();
  // CJC2025: New Parameters
  double getEvapThreshInPrevGrid();
  double getEvapThreshInUntilGrid();
  double getTransThreshInPrevGrid();
  double getTransThreshInUntilGrid();



  // SKYnGM2008LU
  double getAvCanStorParam();
  double getAvIntercepCoeff();
  double getAvThroughFall();
  double getAvCanFieldCap();
  double getAvDrainCoeff();
  double getAvDrainExpPar();
  double getAvLandUseAlb();
  double getAvVegHeight();
  double getAvOptTransmCoeff();
  double getAvStomRes();
  double getAvVegFraction();
  double getAvLeafAI();
  // CJC2025: New Parameters
  double getAvEvapThresh();
  double getAvTransThresh();

  double getAvSoilMoisture();         // Integral characteristics
  double getAvEvapFract();
  double getAvET();

  tEdge * getFlowEdg();

  int getFloodStatus();
  int getSoilID();                     //Invariant Members
  int getLandUse();
  int NoMoreTracers();
  int getReach();
  // Added by Giuseppe Mascaro in 2016 to allow ingestion of soil grids
  double getKs();
  double getThetaS();
  double getThetaR();
  double getPoreSize();
  double getAirEBubPres();
  double getDecayF();
  double getSatAnRatio();
  double getUnsatAnRatio();
  double getPorosity();
  double getVolHeatCond();
  double getSoilHeatCap();
  
  double getCentroidX();				// Geometric Methods
  double getCentroidY();
  int  polyCentroid(double *, double *, int, double *, double *, double *);
  
  void setMuOld(double);               //Unsat and Sat Zone Members
  void setMuNew(double);
  void setRiOld(double);
  void setRiNew(double);
  void setRuOld(double);
  void setRuNew(double);
  void setNfOld(double);
  void setNfNew(double);
  void setNtOld(double);
  void setNtNew(double);
  void setNwtOld(double);
  void setNwtNew(double);
  void setMiOld(double);
  void setMiNew(double);
  void setQpout(double);               //Runoff Members
  void setQpin(double);
  void setRain(double);
  void setSrf_Hr(double);
  void setsrf(double);
  void sethsrf(double);
  void setesrf(double);
  void setpsrf(double);
  void setsatsrf(double);
  void setrsrf(double);
  void setsbsrf(double);
  void setRunOn(double);               //Runon
  void setGwaterChng(double);
  void setFlowEdg(tEdge *);
  void setFloodStatus( int status );
  void setIntStormVar(double);
  void setSoilID(int);                 //Invariant Members
  void setLandUse(int);               
  void setReach(int);
  
  // Added by Giuseppe Mascaro in 2016 to allow ingestion of soil grids
  void setKs(double);
  void setThetaS(double);
  void setThetaR(double);
  void setPoreSize(double);
  void setAirEBubPres(double);
  void setDecayF(double);
  void setSatAnRatio(double);
  void setUnsatAnRatio(double);
  void setPorosity(double);
  void setVolHeatCond(double);
  void setSoilHeatCap(double);
  
  void setTTime(double);               //Routing Members
  void setHillPath(double);   
  void setStreamPath(double); 
  void setInterceptLoss(double);       //Interception Members
  void setNetPrecipitation(double);
  void setCanStorage(double);      
  void setPotEvap(double);        
  void setActEvap(double);
  void setStormLength(int,double);        
  void setCumIntercept(double);   
  void setEvapoTrans(double);          //Evapotranspiration Members
  void setEvapWetCanopy(double);   
  void setEvapDryCanopy(double);   
  void setEvapSoil(double);        
  void setSoilMoisture(double);    
  void setSoilMoistureSC(double);    
  void setSoilMoistureUNSC(double);    
  void setRootMoisture(double);
  void setRootMoistureSC(double);
  void setTransmiss(double);      
  void setAirPressure(double);         //Meteorological Members
  void setDewTemp(double);        
  void setRelHumid(double); 
  void setVapPressure(double);       
  void setSkyCover(double);        
  void setWindSpeed(double);       
  void setAirTemp(double);        
  void setSurfTemp(double);        
  void setSoilTemp(double);
  void setNetRad(double);          
  void setGridET(double);

  void setShortRadIn(double);
  void setShortRadIn_dir(double);
  void setShortRadIn_dif(double);
  void setShortAbsbVeg(double);
  void setShortAbsbSoi(double);
  void setShortRadSlope(double); /// Set the slope‐corrected incoming shortwave (before canopy/albedo), JB2025

  void setLongRadIn(double);
  void setLongRadOut(double);

  void setGFlux(double);
  void setHFlux(double);
  void setLFlux(double);
  void setGnod(double);
 
  // ASM 2/10/2017
  void setChannelPerc(double);
  void setFt(double); //ASM

  void setContrArea(double);
  void setCurvature(double); 
  void setBasinArea(double);
  void setBedrockDepth(double);
  void setAspect(double);
  void setVegFraction(double);

  // SKYnGM2008LU
  // additions for landuse grids AJR 2007
  void setLandUseAlb(double);
  void setThroughFall(double);
  void setVegHeight(double);
  void setStomRes(double);

  // SKYnGM2008LU
  void setIntercepCoeff(double);
  void setCanFieldCap(double);
  void setDrainCoeff(double);
  void setDrainExpPar(double);
  void setOptTransmCoeff(double);
  void setLeafAI(double);
  void setCanStorParam(double);
  // CJC2025: New Parameters
  void setEvapThresh(double);
  void setTransThresh(double);

  void setLandUseAlbInPrevGrid(double);
  void setLandUseAlbInUntilGrid(double);
  void setThroughFallInPrevGrid(double);
  void setThroughFallInUntilGrid(double);
  void setVegHeightInPrevGrid(double);
  void setVegHeightInUntilGrid(double);
  void setStomResInPrevGrid(double);
  void setStomResInUntilGrid(double);
  void setVegFractionInPrevGrid(double);
  void setVegFractionInUntilGrid(double);
  void setCanStorParamInPrevGrid(double);
  void setCanStorParamInUntilGrid(double);
  void setIntercepCoeffInPrevGrid(double);
  void setIntercepCoeffInUntilGrid(double);
  void setCanFieldCapInPrevGrid(double);
  void setCanFieldCapInUntilGrid(double);
  void setDrainCoeffInPrevGrid(double);
  void setDrainCoeffInUntilGrid(double);
  void setDrainExpParInPrevGrid(double);
  void setDrainExpParInUntilGrid(double);
  void setOptTransmCoeffInPrevGrid(double);
  void setOptTransmCoeffInUntilGrid(double);
  void setLeafAIInPrevGrid(double);
  void setLeafAIInUntilGrid(double);
  // CJC2025: New Parameters
  void setEvapThreshInPrevGrid(double);
  void setEvapThreshInUntilGrid(double);
  void setTransThreshInPrevGrid(double);
  void setTransThreshInUntilGrid(double);

  // SKYnGM2008LU
  void setAvCanStorParam(double);
  void setAvIntercepCoeff(double);
  void setAvThroughFall(double);
  void setAvCanFieldCap(double);
  void setAvDrainCoeff(double);
  void setAvDrainExpPar(double);
  void setAvLandUseAlb(double);
  void setAvVegHeight(double);
  void setAvOptTransmCoeff(double);
  void setAvStomRes(double);
  void setAvVegFraction(double);
  void setAvLeafAI(double);
  // CJC2025: New Parameters
  void setAvEvapThresh(double);
  void setAvTransThresh(double);

  void setAvSoilMoisture(double);     // Integral characteristics
  void setAvEvapFract(double);
  void setAvET(double);

  void setHlevel(double);              //tKinemat Members
  void setQstrm(double);
  void setChannelWidth(double);
  void setRoughness(double);
  void setFlowVelocity(double);

  // SKY2008Snow from AJR2007
  //snow pack
  void setLiqWE(double);
  void setIceWE(double);
  void setDU(double);
  void setLiqRouted(double);
  void setSnTempC(double);
  void setCrustAge(double);
  void setDensityAge(double);
  void setEvapoTransAge(double);
  void setSnLHF(double);
  void setSnSHF(double);
  void setSnSub(double); // Snowpack sublimation CJC2020
  void setSnEvap(double); // Snowpack evaporation CJC2020
  void setSnGHF(double);
  void setSnPHF(double);
  void setSnRLout(double);
  void setSnRLin(double);
  void setSnRSin(double);
  void setUnode(double);
  void setUerror(double);
  void setPersTimeMax(double);//integrated output
  void setPersTime(double);//integratedoutput
  void setPeakSWE(double);//integratedoutput
  void setPeakSWETemp(double);//integratedoutput
  void setInitPackTime(double);//integratedoutput
  void setInitPackTimeTemp(double);  
  void setPeakPackTime(double);//integratedoutput
  //snowintercept
  void setIntSWE(double);
  void setIntPrec(double);
  void setIntSnUnload(double);
  void setIntSub(double);
  //tShelter members
  void setHorAngle0000(double);
  void setHorAngle0225(double);
  void setHorAngle0450(double);
  void setHorAngle0675(double);
  void setHorAngle0900(double);
  void setHorAngle1125(double);
  void setHorAngle1350(double);
  void setHorAngle1575(double);
  void setHorAngle1800(double);
  void setHorAngle2025(double);
  void setHorAngle2250(double);
  void setHorAngle2475(double);
  void setHorAngle2700(double);
  void setHorAngle2925(double);
  void setHorAngle3150(double);
  void setHorAngle3375(double);
  void setSheltFact(double);
  void setLandFact(double);

  // void setCanopyStorage(double);       //tWaterBalance Members
  // SKYnGM2008LU: setCanopyStorage(double) now replaced by setCanopyStorVol() below.
  // tWaterBalance Members
  void setCanopyStorVol(double); // SKYnGM2008LU
  void setUnSaturatedStorage(double);  
  void setSaturatedStorage(double);   
  void setRecharge(double);  
  void setUnSatFlowOut(double);
  void setUnSatFlowIn(double);

  void addGwaterChng(double);
  void addQpin(double);
  void addTTime(double);
  void addIntStormVar(double);
  void addContrArea(double);
  void addQgwOut(double);
  void addQgwIn(double);
  void addQstrm(double);
  void addAvSoilMoisture(double);
  void addCumSrf(double); //added CJC2021
  void addSrf_Hr(double);

  // SKY2008Snow from AJR2007
  //snow
  void addLatHF(double);
  void addSHF(double);
  void addSnSub(double); // Snowpack sublimation CJC2020
  void addSnEvap(double); // Snowpack evaporation CJC2020
  void addTotEvap(double); // Snowpack evaporation CJC2020
  void addBarEvap(double); // Snowpack evaporation CJC2020
  void addPHF(double);
  void addRLin(double);
  void addRLout(double);
  void addRSin(double);
  void addGHF(double);
  void addMelt(double);
  void addCumHrsSun(double);
  void addCumHrsSnow(double);
  void addCumUerror(double);
  //snInt
  void addIntSub(double);
  void addIntUnl(double);

  int  getTracer();                    //Tracer routines
  void setTracer(int);
  void ActivateSortTracer();
  void DeactivateTracer();
  void MoveSortTracerDownstream();
  void AddTracer();
  void setStreamNode(tCNode *);

  tCNode * getDownstrmNbr();
  tCNode * getStreamNode();

  tList< int >    * getTimeIndList();
  tList< double > * getQeffList();
  void   allocDataStack();
  void   deleteDataStack();

  int    nVerts;                         //tResample Members
  double *VertsX;   
  double *VertsY;  
  tEdge  *bndEdge1; 
  tEdge  *bndEdge2; 
  void   deleteVertArrays();   
  void   allocVertArrays(int); 
  void   writeRestart(fstream&) const;
  void   readRestart(fstream&);
  void   printVariables();

  int    satOccur;              // Surface saturation occurence 
  double hsrfOccur;             // Infiltration excess runoff Occurence
  double psrfOccur; 		// Perched Saturation Runoff Occurence 
  double satsrfOccur; 	        // Groundwater Saturation Occurence
  double sbsrfOccur; 		// Saturation from Below runoff Occurence
  double RechDisch;
  double percOccur;		// ASM Channel Percolation Occurence
  double avPerc;		// ASM Channel Percolation Occurence

  // for beta cutoff statements
  void setSoilCutoff(double);
  void setRootCutoff(double);

  double getSoilCutoff();
  double getRootCutoff();


protected:
  double ContrArea;             // Surface contributing area for the node 
  double Curvature;             // Topographic curvature of the element 
  double BasinArea;
  double alpha; 		// Slope angle
  double NwtOld, NwtNew; 	// Water table depth in mm
  double MuOld, MuNew; 		// Moisture Content above WT in mm
  double MiOld, MiNew; 		// Initialization Moisture in mm
  double NtOld, NtNew; 		// Top Front in mm
  double NfOld, NfNew; 		// Wetting Front in mm
  double RuOld, RuNew; 		// Recharge rate in mm
  double RiOld, RiNew; 		// Recharge rate in mm 
  double Qin, Qpin, Qout, Qpout; 
  double gwchange;
  list<double> gwc;  //SMM - added 09232008 for groundwater flux contributions 
  double srf_hr; 		// Runoff produced in one hour in mm
  double srf; 			// Total Runoff Generation in mm
  double cumsrf; 			// added CJC2021
  double hsrf; 			// Hortonian Runoff in mm
  double esrf; 			// Exfiltration in mm
  double psrf; 			// Perched Saturation Runoff in mm
  double satsrf; 	        // Groundwater Saturation in mm
  double rsrf; 			// Return Flow in mm
  double sbsrf; 		// Saturation from Below runoff in mm
  double RunOn;
  double Rain;                  // mm/hr
  double intstorm;              // hour
  double traveltime;            // hour
  double hillpath;              // m
  double streampath;            // m
  double Interception;    
  double NetPrecipitation; 
  double ActEvaporation;
  double CanStorage;    
  double PotEvaporation;   
  double StormLength;     
  double CumIntercept;     
  double EvapWetCanopy;   
  double EvapDryCanopy;   
  double EvapSoil;            
  double EvapoTranspiration;   
  double SoilMoisture;
  double SoilMoistureSC;       
  double SoilMoistureUNSC;       
  double RootMoisture;
  double RootMoistureSC;
  double Transmissivity;   
  double AirTemp;
  double DewTemp;
  double RelHumid;
  double VapPress;
  double SurfTemp;
  double SoilTemp;
  double WindSpeed;
  double SkyCover;
  double NetRad;
  double AirPressure;
  double GridET;
  double BedrockDepth;
  double ShortRadIn, ShortRadIn_dir, ShortRadIn_dif;
  double ShortAbsbVeg, ShortAbsbSoi;
  double ShortRadSlope;   // raw incoming shortwave on the slope (pre‐canopy/albedo), JB2025 @ ASU
  double LongRadIn;
  double LongRadOut;
  double gFlux;
  double hFlux;
  double lFlux;
  double Gnod;
  double Aspect;

  //ASM 2/10/2017
  double ChannelPerc;
  double Ft; //ASM

  // SKY2008Snow from AJR2007
  //snowpack
  double liqWEq;//state var
  double iceWEq;//state var
  double dU;//output
  double Unode;//state var
  double Uerror;//output
  double cumUError;//output
  double liqRoute;//output
  double snTemperC;//state var
  double crAge;//state var
  double densAge;//state var
  double ETage;//state var
  double snLHF;//output
  double snSHF;//output
  double snSub;//output // CJC2020
  double snEvap;//output // CJC2020
  double snGHF; //output
  double snPHF; //output
  double snRLin; //output
  double snRLout; //output
  double snRSin; //output
  double persTime; //output
  double persTimeTemp; //output
  double peakSWE; //output
  double peakSWEtemp; //output
  double initPackTime; //output
  double initPackTimeTemp;
  double peakPackTime; //output
  //snowintercept
  double intSWEq;//state var
  double intSnUnload;//output var
  double intSub;//output var
  double intPrec;//output var
  //shelter
  double horizonAngle0000;
  double horizonAngle0225;
  double horizonAngle0450;
  double horizonAngle0675;
  double horizonAngle0900;
  double horizonAngle1125;
  double horizonAngle1350;
  double horizonAngle1575;
  double horizonAngle1800;
  double horizonAngle2025;
  double horizonAngle2250;
  double horizonAngle2475;
  double horizonAngle2700;
  double horizonAngle2925;
  double horizonAngle3150;
  double horizonAngle3375;
  double sfact;
  double lfact;

  double VegFraction;

  // SKYnGM2008LU
  // additions for landuse grids AJR 2007
  double LandUseAlb;
  double ThroughFall;
  double VegHeight;
  double StomRes;

  // SKYnGM2008LU
  double IntercepCoeff;
  double CanFieldCap;
  double DrainCoeff;
  double DrainExpPar;
  double OptTransmCoeff;
  double LeafAI;
  double CanStorParam;
  // CJC2025: New Parameters
  double EvapThresh;
  double TransThresh;

  // SKYnGM2008LU
  double LandUseAlbInPrevGrid, LandUseAlbInUntilGrid, ThroughFallInPrevGrid, ThroughFallInUntilGrid;
  double VegHeightInPrevGrid, VegHeightInUntilGrid, StomResInPrevGrid, StomResInUntilGrid;
  double VegFractionInPrevGrid, VegFractionInUntilGrid, CanStorParamInPrevGrid, CanStorParamInUntilGrid;
  double IntercepCoeffInPrevGrid,IntercepCoeffInUntilGrid, CanFieldCapInPrevGrid, CanFieldCapInUntilGrid;
  double DrainCoeffInPrevGrid, DrainCoeffInUntilGrid, DrainExpParInPrevGrid, DrainExpParInUntilGrid;
  double OptTransmCoeffInPrevGrid, OptTransmCoeffInUntilGrid, LeafAIInPrevGrid, LeafAIInUntilGrid;
  double EvapThreshInPrevGrid, EvapThreshInUntilGrid, TransThreshInPrevGrid, TransThreshInUntilGrid; // CJC2025: New Parameters

  // SKYnGM2008LU
  double AvCanStorParam, AvIntercepCoeff, AvThroughFall, AvCanFieldCap;
  double AvDrainCoeff, AvDrainExpPar, AvLandUseAlb, AvVegHeight;
  double AvOptTransmCoeff, AvStomRes, AvVegFraction, AvLeafAI;
  double AvEvapThresh, AvTransThresh; // CJC2025: New Parameters

  double AvSoilMoisture;        // Integral characteristics
  double AvEvapFract;
  double AvET;

  // SKY2008Snow from AJR2007
  //snow
  double cumHrsSun;
  double cumLHF;
  double cumSHF;
  double cumSnSub; // Define snowpack sublimation CJC2020
  double cumSnEvap; // Define snowpack evaporation CJC2020
  double cumTotEvap; // Define total ET CJC2020
  double cumBarEvap; // Define bare soil evaporation CJC2020
  double cumPHF;
  double cumRLin;
  double cumRLout;
  double cumRSin;
  double cumGHF;
  double cumMelt;
  double cumIntSub;
  double cumIntUnl;
  double cumHrsSnow;

  double QgwIn;
  double QgwOut;

  double Hlevel;                //tKinemat Members
  double Qstrm;
  double Width;
  double Roughness;
  double FlowVelocity;
    
  // Added by Giuseppe Mascaro in 2016 to allow ingestion of soil grids
  double Ks;
  double ThetaS;
  double ThetaR;
  double PoreSize;
  double AirEBubPres;
  double DecayF;
  double SatAnRatio;
  double UnsatAnRatio;
  double Porosity;
  double VolHeatCond;
  double SoilHeatCap;

  // double CanopyStorage;         //tWaterBalance Members
  // SKYnGM2008LU: CanopyStorage now replaced by CanopyStorVol below
  // tWaterBalance Members
  double CanopyStorVol; // SKYnGM2008LU
  double UnSaturatedStorage;
  double SaturatedStorage;
  double Recharge;
  double UnSatFlowIn;
  double UnSatFlowOut;

  tEdge * flowedge; 
  tCNode * StreamPtr;

  int tracer; 
  int flood; 
  int soiID;     
  int LandUse;  
  int Reach;

  std::shared_ptr<tList<double> > Qeff; //WR debug convert to smart shared pointer to prevent memory leak
  std::shared_ptr<tList<int> > TimeInd;

  double xC;
  double yC;


  // update to beta functions, variable store cutoff of soil and dry canopy evap once Nwt = Bedrock
  double soil_cutoff;
  double root_cutoff;


};

#endif

//=========================================================================
//
//
//                           End of tCNode.h
//
//
//=========================================================================
