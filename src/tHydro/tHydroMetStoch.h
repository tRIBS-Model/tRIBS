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
**  tHydroMetStoch.h:   Header file for class tHydroMetStoch
**
**  This class encapsulates the hydrometeorological forcing generated
**  using functions of weather generator implemented within this class
**  Data members include current day weather data.
**
***************************************************************************/

#ifndef THYDROMETSTOCH_H
#define THYDROMETSTOCH_H

#include "src/Headers/Inclusions.h"

#define  PI12  0.261799387799149

//=========================================================================
//
//
//                  Section 1: tHydroMetStoch Class Declaration
//
//
//=========================================================================

class tHydroMetStoch
{
 public:
  tHydroMetStoch();
  tHydroMetStoch(tMesh<tCNode> *, tRunTimer *, tInputFile &, 
		 tEvapoTrans *, tRainfall *);
  ~tHydroMetStoch();
  void   setLat(double,int);
  void   setLong(double,int);
  void   setGmt(int);
  void   setSunH(double);
  void   setSinH(double);
  void   setIo(double);
  void   setIdir(double);
  void   setIdir_vis(double);
  void   setIdir_nir(double);
  void   setIdif(double);
  void   setIdif_vis(double);
  void   setIdif_nir(double);

  void   SetStochasticHydroMet(tInputFile &);
  void   SetSkyVars();
  void   SetWindVars();
  void   SetTemperatureVars();
  void   SetCloudinessVars();
  void   SetCondBetaPars(double);
  void   ReadWeatherParameters(char *);
  void   OutputHydrometVars();
  void   PrintOutArray(double *, int);

  void   SimulateHydrometVars();

  // SKY2008Snow from AJR2007
  void   SimulateAirDPTemperature();

  double SimulateAirDPTemperatureCurtis();
  double SimulateDewPointTemperature(double, double);
  double SimulateSkyCover(int);
  double SimulateWindCurtis();

  void   ComputeDailyTempCoeffs();
  void   ComputeDailyAvailEnergy(double);
  void   ComputeDailyEpCld(double, double);
  void   ComputeEF(double, double,  double);
  double ComputeDailyTdewKimball(double, double);

  void   Function0(int, double, double, double);
  double Function1(double, double);
  double Function2(double);
  double Function3(double);
  double ExpectedCloudiness(int, double, double);
  double GetCloudTransitValue();
  double inLongWave(double, double);

  void   writeRestart(fstream &) const;
  void   readRestart(fstream &);

  // SKY2008Snow from AJR2007
  void   EstimateResidualsVector();
  void   ExpectedTmaxTmin();
  void   SimulateTmaxTmin(double, double, double, double);
  void   SetDisplaceCorrects(double *, double *);
  void   SetDisplaceVector(double *, int, double, double);
  void   AddDisplaceVector(double *, int, double *);
  void   DiurnalCycleTemperature(double *, double *);
  void   DiurnalCycleDewPointTemperature(double *);
  void   SetTemperatAmplitude(double, double);
  void   SetTemperatParameters();
  void   SetPeakTemperatShift();
  void   StoreCurrValues();
  void   GetVectRandomNorm(double *, int, double, double);

  int    getGmt();
  int    getYear();
  int    getMonth(); 
  int    getDay();
  int    getHour(); 
  double getLat(int);
  double getLong(int);
  double getSunH();
  double getSinH();
  double getIo();
  double getIdir();
  double getIdir_vis();
  double getIdir_nir();
  double getIdif();
  double getIdif_vis();
  double getIdif_nir();

  double getAirTemp();
  double getAirTemp_1();
  double getDewTemp();
  double getSurfTemp();
  double getAtmPress();
  double getSkyCover();
  double getSkyCover_1();
  double getRHumidity();
  double getWindSpeed();
  double getNetRad();
  double getPanEvap();
  double getVaporPress();
  double getNlm();
  double getNci();
  double getF0();

 protected:
  tMesh<tCNode> *gridPtr;
  tRunTimer     *timer;
  tEvapoTrans   *etPtr; 
  tRainfall     *rainPtr;

  ofstream hout;

  int gmt;
  double latitude, longitude;
  double basinLat, basinLong;

  // Previous temperature model -- NOT COUPLED
  double Tmon[12];    // Mean monthly temperatures
  double Stdmon[12];  // Temperatures STDs
  double **A;         // Parameter matrices of multivariate AR(1) model    
  double **B;
  double T_ro1;       // Temperature autocorrelation coefficient lag 1 day
  double tmone;       // Simulated from N(Tseas[i],Sseas[i]^2) mean monthly toC
  double sgma;        // Daily temperature variation
  double tmaxe, tmine;// Expected values of max and min temperatures 
  double td;          // Mean daily temperature
  double tmax;        // Max daily temperature
  double tmin;        // Min daily temperature
  double tdnext;      // Next day (t+1) mean daily temperature
  double tmaxnext;    // Next day (t+1) max daily temperature
  double tminnext;    // Next day (t+1) min daily temperature
  double deltt;
  double deltt_1;
  double tShift;
  double tShift_1;
  double hi[3];       // Standardized residuals for tmax and tmin
  double hi_1[3];
  double DiurnT[24];  // Diurnal temperature cycle (hour)
  double DiurnDP[24]; // Diurnal dew point temperature cycle (hour)


  // Variables for the radiation model (implemented in tEvapotrans.cpp)
  double SunH, SinH;
  double Io;
  double Idir;
  double Idir_vis;
  double Idir_nir;
  double Idif;
  double Idif_vis;
  double Idif_nir;

  // Variables for atmospheric pressure model (NOT COUPLED) 
  double atmPress;
  double Patm_mean;
  double Patm_std;

  // Variables for the cloudiness model 
  double skyCover, skyCover_1;
  double M0;
  double StdSky;
  double GammaSky;
  double m_ro1;
  double mt_1;
  double tAftrRn;
  double tNextRn;
  double Nlm, Nci, f0;

  // Wind model variables
  double windSpeed, windSpeed_1;
  double WindSp_mean, WindSp_std, WindSh_fact;
  double w_ro1, w_skew, w_skew_gam;

  // Current and (t-1) values of air and dew temperature
  double airTemp, airTemp_1;
  double dewTemp, dewTemp_1;

  // Curtis and Eagleson [1982] temperature model
  int TempOpt;     // Option what temperature model to use
  int CldParOpt, TempParOpt;
  int RainOpt;
  double **Cldi;
  double **Tdevi;
  double **Bi;
  double EtA[12][11];
  double EtB[12][11];
  double AnnR[12];
  double betaA, betaB;
  double b0, b1, b2, b3, b4;
  double b5, b6, b7; // currently NOT used 
  double K0, K1, K2, K3, K4, K5, K6;
  double I1, I2, I3, I4, I5, I6, I7;
  double I1_1, I2_1, I3_1, I4_1, I5_1, I6_1, I7_1;
  double temp11pm, qt_1;
  double MeanTDev, StdTDev, AutoCorTDev;
  double Tdev, Tdev_1; 

  // Auxiliary variables
  double AnnRain, EF;
  double AvEnergy, AvEnergy_1;
  double AvEp, AvEp_1;
  double AvCld, AvCld_1;
  double Tmax, Tmin, Tmax_1, Tmin_1;
  double cloudfactor;

  double vaporPress, rHumidity;
  int year, month, day, hour, minute;
  int HrCount;
};

#endif

//=========================================================================
//
//
//                        End of tHydroMetStoch.h
//
//
//=========================================================================
