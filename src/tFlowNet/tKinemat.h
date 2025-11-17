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
**  tKinemat.h: Header for class tKinemat (see tKinemat.cpp) which inherets
**              tFlowNet and implements hydraulic channel routing.
**
***************************************************************************/

#ifndef TKINEMAT_H
#define TKINEMAT_H

//=========================================================================
//
//
//                  Section 1: tKinemat Include Statements
//
//
//=========================================================================

#include "src/tFlowNet/tFlowNet.h"
#include "src/tFlowNet/tReservoir.h"
#include "src/tFlowNet/tResData.h"
#include "src/Headers/Inclusions.h"

using namespace std;

//=========================================================================
//
//
//                  Section 2: tKinemat Define Statements
//
//
//=========================================================================

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
                  (maxarg1) : (maxarg2))
#define FREERETURN {delete [] X; delete [] X1; delete [] XR; delete [] F;\
       delete [] gradf; delete [] aa; delete [] bb; delete [] cc; return; }


//=========================================================================
//
//
//                  Section 3: tKinemat Class Definitions
//
//
//=========================================================================

class tKinemat : public tFlowNet
{
 public:
  tKinemat();
  tKinemat(char **);
  tKinemat(SimulationControl*, tMesh<tCNode> *, tInputFile &, tRunTimer *);
  ~tKinemat();

  void KinematWave(double *, double *, double *,  double *, 
                   double *, double *, double, double, int *);

  void SolveForTwoNodeReach(double *, double *, double *, double *,
                            double *, double *, double, double);

  void ComputeFunction(double *, double *, double *, double *, 
                       double *, double *, double *,  double *,
                       double, double, int);

  void ComputeJacobian(double *, double *, double *, double *,
                       double *, double *, double *, double *, int);
  void SolveLinearSystem(double *, double *, double *, double *, int);

  void GAUSS(double *, double **, int);

  void tridag(double *, double *, double *, double *, double *,
	      unsigned long);
  void bandec(double **, unsigned long, int, int, double **,
              unsigned long *, double *);
  void banbks(double **, unsigned long, int, int, double **,
              unsigned long *, double *);
 
  void lnsrch(int, double *, double, double *, double *, double *,
	      double *,  double, int *, double *,  double *,
              double *, double *, double *, double *, 
              double *, double, double);
  
  void ControlPrint(ofstream &, double *, int);
  void UpdateHsShifted(double *, double *, double, int);
  void AllocateMemory(int);
  void FreeMemory();
  void InitializeStreamReach(int, int);
  void AssignLateralInflux();
  void PrintFlowStacks(ofstream &, tCNode *);
  void ComputeCoefficientArrays();
  void AssignQin();
  void ComputeQout();
  void UpdateStreamVars();
  void initialize_values(tInputFile &, double); // JECR 2015
  void Reservoir_Routing(int); // JECR 2015
  void RunRoutingModel(int, int *, double);
  void RunHydrologicRouting();
  void SurfaceFlow();
  void setTravelVelocityKin(double, double);
  void UpdateForNewRun(tInputFile &, int);
  void AssignChannelWidths(tInputFile &);

  double fmin(double *, int);
  double getQout() { return Qout; }
  double getHn()  { return his[m]; }
  double ComputeNodeQstrm(int);
  double ComputeNodeFlowVel(int);
  double RetrieveQeff(tCNode *);

  void writeRestart(fstream &) const;
  void readRestart(fstream &);

  ifstream GeomtFile;       // Channel geometry input file
  ofstream theOFStream;     // Output file to store all the info
  ofstream ControlOut;      // Output file to control simulations

  SimulationControl *simCtrl;   

#ifdef PARALLEL_TRIBS
  int getId() { return id; }
  int getN() { return n; }
  double* getBis() { return bis; }
  double* getAis() { return ais; }
  double* getRifis() { return rifis; }
  double* getSiis() { return siis; }
  double* getC() { return C; }
  double* getY1() { return Y1; }
  double* getY2() { return Y2; }
  double* getY3() { return Y3; }

  // Open Outlet file on the processor that it resides
  void openOutletFile(tInputFile &);
#endif

protected:
  int     id;         // ID of the stream link
  int     n, m, m1;   // # nodes, # nodes-1, # nodes-2

  double  dt, dtReff; // Computational dt and time interval for lateral influx 
  double  qit, Qin;   // Qin is influx (m3/s) in the upper node (qit=0.5Qin)
  double  H0;         // H0 is the water level corresponding to Qin, m
  double  Qout;       // Discharge in the n-th node of the stream reach, (m3/s)
  double  maxH;       // max water level from the preceding time step
  double  maxReff;    // max lateral influx at the current time step
  double  Roughness;  // Uniform roughness
  double  Width;      // Uniform width
  double  kincoef;    // Uniform coefficient in power law relationship 
   
  double  *ais, *bis, *his, *reis, *siis, *rifis, *sumis;
  double  *C, *Y1, *Y2, *Y3;

  double  *OutletHlev; // Used for storage of the outlet H values

  tCNode  *cHead;      // Ptr to a current stream head node
  tCNode  *cOutlet;    // Ptr to a current outlet node
  
  int     TimeSteps;

/*** Start edits by JECR 2015 ***/
  int optres, checkID, checkNode;
  double  ChannelConduc, ReachLoss, poro, reis1, Pchannel, Preach, ChanWidth, TotWidth, ChanLength, channelPorosity; //ASM
  double  TransientConduc, TransientTime, TotChanLength; //ASM
  double  ExtraPerc, PotentialPerc, TotalPerc, Qin1; //ASM 2/14/2017 Variables for calculating channel loss
  int percolationOption, CountLimit, Count; //ASM perolcation option
  double  *clis, *NodeLoss; // ASM not sure if double is correct
  double  ParallelPerc, IndividualPerc; //ASM test to make parallel percolation
  double  Ft_init,Ft,ThetaResid,PoreInd,PsiB,PsiF,DeltTh,test,Ft_prime,Ft_previous,IntStormVar,IntStormMax,rate; //ASM Green Ampt variables
  double tempVariable;
  double resTimeStep;
  double resRunTime;
  tReservoir LevelPool;
  tPreProcess ResReadItem;
/**** End edits by JECR 2015 ****/

};

#endif

//=========================================================================
//
//
//                          End tKinemat.h
//
//
//=========================================================================
