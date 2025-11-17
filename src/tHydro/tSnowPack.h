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
**  tSnowPack.h:   Header file for tSnowPack Class
**
**  This class creates a class derived from tEvapoTrans that basically 
**  operates as the snowpack in the model. There are a few different kinds of
**  functions. The first are the constructors/destructors. Second come the 
**  initialization/finalization schemes; this is where tSnowPack primarily
**  comes into contact with tCNode. Next, is the calling function, where 
**  the primary physics are outlined. A large amount of code from 
**  tEvapoTrans is recreated in the calling function. Then, there are
**  ancillary functions that define the material properties. Lastly come the
**  physical algorithms, which define or constrain the energy and mass fluxes
**  into the pack.
**
**  Variables similarly fall into a few classifications. These include 
**  material properties, state variables, output/verification variables
**  and energy/mass fluxes. I also construct unit conversions. This is
**  something that probably would be more appropriate to do one a global 
**  scope, but I have choosen to do here for simplicity.
**
**  The rest of the code only interacts with the snow module in three ways. 
**  The first in the initialization of the object. Secondly, it calls the
**  module. Lastly, the tCNode objects are adjusted.
**
**  When the module is called, it checks to see if (1) snow is on the ground
**  (2) snow is falling, or (3) snow is unloading from the canopy. If none of
**  these are the case, the tEvapoTrans::callEvapoPotential() and 
**  tEvapoTrans::callEvapoTrans() are recreated. Before this, it has already
**  dealt with the canopy interception. If there is snow, then the precipitation
**  and latent mass fluxes are calculated and the pack is adjusted appropriately.
**  Then, the current energy balance is calculated, and the change in energy is
**  recorded and the total energy of the pack is updated. If the total energy
**  is less than 0J/m^2, then the new temperature of the pack is calculated.
**  Otherwise, the liquid content is calculated. If the liquid content is 
**  greater than 40% of the total pack, the excess is routed and the total
**  SWE and the liquid WE are adjusted appropriately.
**
**  Also, the mean age of the pack and the age of the surface are tracked 
**  throughout in order to have an estimate of pack density and surface
**  albedo.
**
**  09 July 2007 -- AJR @ New Mexico Tech
**
**************************************************************************/

#ifndef TSNOWPACK_H
#define TSNOWPACK_H

//=========================================================================
//
//		Section 1: tSnowPack Include Statements
//
//=========================================================================

#include "src/Headers/Inclusions.h"
#include "src/tHydro/tEvapoTrans.h"



//=========================================================================
//
//		Section 2: tSnowPack Class Definitions
//
//=========================================================================

class tSnowPack : public tEvapoTrans
{

	
public:

	
  //constructors
  tSnowPack();
  tSnowPack(SimulationControl *,tMesh<tCNode> *, tInputFile &, tRunTimer *,
		  tResample *, tHydroModel *, tRainfall *);

  //destructor
  ~tSnowPack();

  //initialization and update routine
  void SetSnowPackVariables(tInputFile &);
  void SetSnowVariables(tInputFile &);
  void SetSnowInterceptVariables();
  void checkShelter(tCNode *cNode);

  //calling functions
  void callSnowPack(tIntercept *, int);
  void callSnowIntercept(tCNode *, tIntercept *,int count);

  //initialization, interact w/ tCNode
  void getFrNodeSnP(tCNode *);
  void setToNodeSnP(tCNode *);
  
  //physical routines
  double densityFromAge();
  void computeSub();
  void computeUnload();
  void updateRipeSnowPack(double);
  void updateSolidSnowPack(double);

  //EB functions

  //basic calculations
  double latentHFCalc(double);
  double sensibleHFCalc(double);
  double snowFracCalc();
  double precipitationHFCalc();
  double agingAlbedo();
  double resFactCalc();
  double inShortWaveSn(tCNode *);
  double inShortWaveCan();
  double emmisSn();

  
  //EB function
  void snowEB(int, tCNode *); // AJR2008, SKY2008Snow

  //conversion functions
  double CtoK(double);
  double KtoC(double);
  
  //communication functions
  int getSnowOpt();
  
  // Restart functions
  void writeRestart(fstream &) const;
  void readRestart(fstream &);

protected:
  int hillAlbedoOption;
  double densityAge; //hr
  double rainTemp;
  double ETAge; //min

  //discretization
  double timeSteph, timeSteps, timeStepm;
  double minutelyTimeStep;  
  
  //state variables (intrinsic)
  double liqWE, iceWE, snWE; //cm
  double snSub, snEvap; //cm added by CJC2020
  double canWE; //cm
  double liqRoute; //cm
  double liqWEm, iceWEm, snWEm; //m
  double Utot, Usn, Uwat, Utotold; //internal energy (kJ/m^2), set to 0 at T=0 C
  double liqWatCont; // degree of saturation
  double liqTempC, iceTempC, snTempC; //Celsius
  double liqTempK, iceTempK, snTempK; //Kelvin
  double crustAge; //hrs  

  //fluxes, changes in energy/mass
  double H,L,G,Prec,Rn; //sensible, latent, ground, precipitation, net radiative heat fluxes (kJ/(m^2)(s))
  double dUint, RLin, RLout, RSin, Uerr; //change in energy, radiative fluxes, EB error (kJ/m^2s)
  double snPrec, liqPrec; //cm
  double snPrecm, liqPrecm; //m
  double snPrecmm, liqPrecmm; //mm  
  double snUnload, snCanWE; //cm
  double vapPressSmb, vapPresskSPa; //vapor pressure (mb and Pa)

  //intercept variables (from Liston and Elder 2006, section 3)
  double Qcs, Ce, I, Iold, psiS; //
  double Imax, prec, LAI; //
  double kc, iceRad, dmdt; //
  double Omega, Sp, RH, D, rhoVap; //
  double Sh, Nu, Re; //
  double KtAtm, Ta, Mwater, R; //
  double RdryAir, esatIce, nu, beta; //
  double acoefficient; //
  double Lm; // unloading
  double airTempK; //
  double effPrecip; //


    //density parameters
  double rholiqcgs, rhoicecgs, rhosncgs; //g/cm^3
  double rholiqkg, rhoicekg, rhosnkg; //kg/m^3
  double rhoAir; //kg/m^3
  double phfOnOff;
  
  //thermal propoerties
  double cpsnowkJ, cpicekJ, cpwaterkJ, cpairkJ; //heat capacity kJ/Kg
  double latFreezekJ, latVapkJ, latSubkJ; //latent heats kJ/g

  //EB/Surface properties
  double resFact;
  double albedo;
  double hillalbedo;
  double compactParam, rhoSnFreshkg;
  double minSnTemp;
  double snliqfrac; // Added by CJC2020

  //output variables
  double snDepth,snDepthm; //snow depths (cm,m)
  double snOnOff;
  double peakSnWE, peakSnWEtemp; //maximum SWE (cm)
  double persMax, persMaxtemp; //time of peristence of pack (hours)
  double inittime, inittimeTemp, peaktime; //time of initial bulk of snow pack and time of peak

  //conversion factors
  double naughttokilo, kilotonaught, cgsRHOtomks, mksRHOtocgs;
  double naughttocm, cmtonaught, ctom, mtoc;


};

#endif

//==========================================================================
//
//		tSnowPack -- End of tSnowPack.h
//
//==========================================================================
