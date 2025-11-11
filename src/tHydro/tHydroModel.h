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
**  tHydroModel.h:   Header file for tHydroModel Class
**
**  Core class for hydrologic computations modified as a combination of
**  RIBS hydrologic processes and CHILD software environment
**
***************************************************************************/

#ifndef THYDROMODEL_H 
#define THYDROMODEL_H

//=========================================================================
//
//
//                  Section 1: tHydroModel Include Statements
//
//
//=========================================================================

#include "src/Headers/Inclusions.h"

#define LAMBEPS 2.2204E-16

//=========================================================================
//
//
//                  Section 2: tHydroModel Class Definitions
//
//
//=========================================================================

class tHydroModel
{
public:     
  tHydroModel(SimulationControl*, tMesh < tCNode > *,
	      tInputFile &, tResample *, tWaterBalance *, tRunTimer *);
  ~tHydroModel();

  int    HydroNodesExist();

  void   InitSet(tResample *); 
  void   SetHydroMVariables(tInputFile &, tResample *, int);
  void   SetHydroNodes(char *); 
  void   InitIntegralVars(); 
  void   SetupNodeUSZ(tCNode *);
  void   SetupNodeSZ(tCNode *);
  void   CheckMoistureContent(tCNode *);

  void   UnSaturatedZone(double);
  void   SaturatedZone(double);
  void   Reset();
  void   ResetGW();
  void   ComputeFluxesNodes1D(); 
  void   ComputeFluxesEdgesND();
  void   PrintOldVars(tCNode *, tEdge *, double, int);
  void   PrintNewVars(tCNode *, double); 
  void   PrintNewGWVars(tCNode *, int); 
  
  double get_Total_Moist(double);         
  double get_Upper_Moist(double, double); 
  double get_Lower_Moist(double, double) const;
  
  //double get_Z1Z2_Moist(double, double, double);
  // SKY2008Snow from AJR2007
  double  get_Z1Z2_Moist(double, double, double, int);
  
  double get_InitMoist_depthZ(double) const;
  double get_EdgeMoist_depthZ(double) const;

  double get_RechargeRate(double, double) const;
  double get_UnSat_LateralFlow(double, double, double ) const;
  double get_Sat_LateralFlow(double , double, double, double) const;
  double getTransmissivityFinD(double) const;
  double getTransmissivityInfD(double) const;
  double GetCellRunon(tCNode *, double);
  double ComputeSurfSoilMoist(double);

  void    writeRestart(fstream &) const;
  void    readRestart(fstream &);

  void   set_Suction_Term(double);    
  void   SetCellRunon(tCNode *, double, double, double, int);

  void   polyn(double,double&,double&,double,double,double) const;
  void   polyn(double,double,double&,double&,double,double,double) const;
  double LambertW(double);
  double Newton(double, double);  
  double Newton(double, double, double, int) const;
  double rtsafe_mod(double, double, double, double, double, double, double);

  char   gwatfile[kMaxNameSize]{};
  char   bedrockfile[kMaxNameSize]{};

  GenericSoilData *soilPtr;   
  GenericLandData *landPtr;     
  tWaterBalance *balPtr;         // Pointer to water balance 
  SimulationControl *simCtrl;    // Pointer to simulation control

#ifdef PARALLEL_TRIBS
  double getDM100() { return dM100; }
  double getDMRt()  { return dMRt; }
  double getMTh100() { return mTh100; }
  double getMThRt()  { return mThRt; }
  double getFSoi100() { return fSoi100; }
  double getFTop100() { return fTop100; }
  double getFClm100() { return fClm100; }
  double getFGW100() { return fGW100; }
#endif

protected: 
  tArray<double> gwaterval;

private: 
  int ID{}, numNodes{};
  int *nodeList;
  int EToption{}, Ioption{};                // Options for hydrologic processs
  
  // SKY2008Snow from AJR2007
  int SnOpt{};
  
  int gFluxOption{}, BRoption{}, RdstrOption{}, GWoption{};
  int RunOnoption{};

  int percolationOption{}; //ASM 2/14/2017

  double NwtOld, NwtNew;   		// Water table depth in mm
  double MuOld,  MuNew;    		// Moisture Content above WT in mm
  double MiOld,  MiNew;    		// Initialization Moist above WT in mm
  double NfOld,  NfNew;    		// Wetting Front in mm
  double NtOld,  NtNew;    		// Top Front in mm
  double RuOld,  RuNew;    		// Recharge rate above wet front in mm
  double RiOld,  RiNew;    		// Recharge rate in mm 
  double QpIn{}, QpOut{}, QIn{}, QOut{};
  double IntStormVar{};
  double IntStormMAX{};

  double R{}, R1{}, Rain{};      		// Rainfall in mm/h;
  double BasArea{};                       // Total basin area in m^2
  double alpha{};  			// Slope angle of the node in radians
  double Cos{}, Sin{};
  double gwchange{};

  double srf{};    			// Total Runoff Generation
  double hsrf{};   			// Hortonian Runoff
  double esrf{};   			// Exfiltration
  double psrf{};   			// Perched Saturation Runoff
  double satsrf{}; 			// Groundwater Saturation
  double sbsrf{};  			// Saturation from Below runoff
  double qrunon{};                        // Runon amount [mm/hr]

  double G{};           			// Capillary drive across the wet front
  double SeIn{}, Se0{};   			// Effective saturation in the power
  double ThRiNf{}, ThReNf{};                // (3 + 1/lambda)

  double Ksat{};
  double F{};
  double Ths{};
  double Thr{};
  double Ar{};
  double UAr{};
  double PoreInd{};
  double Eps{};
  double Psib{};
  double porosity{};
  double Stok{};    			// Cumulative runoff value M^3
  double TotRain{}; 			// Cumulative rainfall value M^3

  // SKYnGM2008LU: Land Use Parameters
  double a_LU{};
  double b1_LU{};
  double P_LU{};
  double S_LU{};
  double K_LU{};
  double b2_LU{};
  double Al_LU{};
  double h_LU{};
  double Kt_LU{};
  double Rs_LU{};
  double V_LU{};
  double LAI_LU{};
  // CJC2025 Stress Thresholds
  double ST_LU{};
  double SE_LU{};

  // SKY2008Snow from AJR2007
  double snowMeltEx{};
  double swe{};

  double TotGWchange{}; 			// Cumulative GW storage change M^3
  double TotMoist{}; 			// Cumulative change in moisture storage
  double DtoBedrock{}; 			// Depth to bedrock
  double surfaceSoilDepth; // Depth for surface soil moisture [mm]
  double rootZoneDepth;    // Depth for root zone moisture [mm]
  
  ofstream fctout;
  double fSoi100{}, fTop100{}, fClm100{}, fGW100{}, dM100{}, dMRt{}, mTh100{}, mThRt{};

  tMesh<tCNode>   *gridPtr;      // Pointer to mesh
  tRunTimer *timer;              // Pointer to timer
   
};

#endif

//=========================================================================
//
//
//                       End of tHydroModel.h
//
//
//=========================================================================
