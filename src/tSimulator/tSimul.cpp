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
**  tSimul.cpp: Functions for class Simulator and SimulationControl 
**              (see tSimul.h)
**
***************************************************************************/

#include <cmath>
#include <cstdint>
#include <sstream>
#include <vector>
#include "src/tSimulator/tSimul.h"
#include "src/Headers/TemplDefinitions.h"
#include "src/Headers/globalIO.h"

#ifdef PARALLEL_TRIBS
#include "src/tGraph/tGraph.h"
#include "src/tParallel/tParallel.h"
#include <mpi.h>
#endif

//=========================================================================
//
//
//                  Section 1: Simulator Constructors/Destructors
//
//
//=========================================================================

Simulator::Simulator(SimulationControl *simctrlptr, tRainfall *rainptr, 
					 tRunTimer *tmrptr, tCOutput<tCNode> *otpptr,
                tRestart<tCNode> *restartptr)
{  
	simCtrl = simctrlptr;
	rainIn  = rainptr;
	timer   = tmrptr;
	outp    = otpptr;
    restart = restartptr;

	// Time tag of initial time, hour
	begin_hour = timer->getCurrentTime(); 
	
	// Counter of time, used for GW model
	GW_label = 0.;

	// Output data on Mesh and Voronoi elements
	outp->WriteOutput( 0 );
	
	// Get rainsearch if rainfall used
	if (rainIn->rainfallType == 1) {
		searchRain = rainIn->searchRain;     // Rainfall search threshold
	}
	count = 0;
}

Simulator::~Simulator() 
{  
	simCtrl = nullptr;
	rainIn = nullptr;
	timer = nullptr;
	outp = nullptr;
    
	Cout<<"Simulator Object has been destroyed..."<<endl<<flush;
}


//=========================================================================
//
//
//                  Section 2: Simulator Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  Simulator::initialize_simulation()
**  
**  Carry out initial activities before simulation begins 
**
*****************************************************************************/
void Simulator::initialize_simulation(tEvapoTrans *EvapoTrans, tSnowPack *SnowPack,
                                      tInputFile& InFl) 
{
	simCtrl->first_time = 'Y';

    //Read in previous command line arguments that are now specified in the input file WR 08282023
    /*  removed command line arguments that should be specified in input file
    "OPTGROUNDWATER" -G    Run groundwater model: GW_model_label
    "OPTSPATIAL" -R    Write intermediate states (spatial output): inter_results
    */

    if (InFl.IsItemIn( "OPTGROUNDWATER" ))
        simCtrl->GW_model_label = InFl.ReadItem(simCtrl->GW_model_label, "OPTGROUNDWATER");
    else
        simCtrl->GW_model_label = true; //Default option

    if (InFl.IsItemIn( "OPTSPATIAL" ))
        simCtrl->inter_results = InFl.ReadItem(simCtrl->inter_results, "OPTSPATIAL");
    else
        simCtrl->inter_results = false; //Default option

	// Ouput pre-processing
	if (simCtrl->inter_results)
		outp->CreateAndOpenDynVar();

    //WR debug 01032024: this was setting the met forcing values to 0 at time 0 in the Dynamic and Pixel files
	// Output initial conditions
	//outp->WriteDynamicVars( timer->getCurrentTime() );
	//outp->WritePixelInfo(   timer->getCurrentTime() );

	// Prepare rainfall input
	if (rainIn->rainfallType == 1) {
		
		// Compose rainfall file name
		while ( !(rainIn->Compose_In_Mrain_Name(timer)) ) { 
			if (count == 0) {
				Cout<<"\nWarning: Next rainfall file "<<rainIn->mrainfileIn
				<<" is missing..."<<endl;
			}
			Cout<<"File "<<rainIn->mrainfileIn<<" was not found..."<<endl;
			
			timer->addRainTime();
			
			if ( timer->getRainTime()-timer->getEndTime() > searchRain ) {
				Cout<<"\nRainfall search threshold exceeded... "<<endl;
				Cout<<count+1<<" rainfall input files are missing..."<<endl; 
				Cout<<"Exiting Program..."<<endl<<endl<<endl;
				exit(2);
			}
			count++;
		}
		
		lmr_hour = timer->getRainTime();     //Time tag of LAST measured rain hour
		dt_rain  = timer->getRainDT();   
		
		//rainIn->NewRain(timer); //WR debug 01032024: this was effectively truncating the rainfall vector, shifting all values by one hour toward the initial runtime
		
		if (simCtrl->Verbose_label == 'Y') {
			Cout<<"\nNext rainfall input: "<<lmr_hour<<" hours in simulation."<<endl;
			Cout<<"Unsaturated zone time steps for interval: ";
			Cout<<timer->getElapsedSteps(lmr_hour)<<endl;
		}
		
		if ( dt_rain < timer->getTimeStep() ) {
			Cout <<"\nComputation DT for unsaturated zone must be ";
			Cout <<"less/equal to DT of Rainfall input data"<<endl;
			Cout <<"Exiting Program..."<<endl;
			exit(2);
		}
		count = 0;
	} else {
		// rainIn->callRainGauge(); 
		// rainIn->callRainGauge(timer); // SKY2008Snow//WR debug 01032024: this was effectively truncating the rainfall vector, shifting all values by one hour toward the initial runtime
	}
	
	met_hour = timer->getMetTime(1);
	eti_hour = timer->getMetTime(2);
	
   // SKY2008Snow
   // Read in weather station data
   //SMM - 09252008 moved from simulation_loop, needs to be done before restart
   if (SnowPack->getSnowOpt() == 0) {
      EvapoTrans->CreateHydroMetAndLU(InFl);
   }
   else { //snow active
      SnowPack->CreateHydroMetAndLU(InFl);
   }

}

/*****************************************************************************
**  
**  Simulator::simulation_loop()
**  
**  Rainfallloops for the whole basin for one cycle. It computes basin 
**  evolution with measured rain, write results if needed.
**
**  Algorithm:
**   get timer information and define number of steps
**   call function to read measured rain
**   for all steps:
**     call f_n to compute model evolution (computation loop)
**     if step corresponds to end of measured rain
**        if writing results is active
**            call function to basin state and results
**
*****************************************************************************/
void Simulator::simulation_loop(tHydroModel *Moisture, tKinemat *Flow,
								tEvapoTrans *EvapoTrans, tIntercept *Intercept, 
								tWaterBalance *Balance, tSnowPack *SnowPack, // SKY2008Snow from AJR2007
								tInputFile &InFl) // SKY2008Snow
{
	Cout<<"\nHydrologic Simulation begins...\n"<<endl;
	
#ifdef PARALLEL_TRIBS
   // Open Outlet file on the processor that it resides
   Flow->openOutletFile(InFl);

   // Exchange static data for ghost nodes
   tGraph::sendInitial();
   tGraph::receiveInitial();
#endif

   // Get the restart information
   double restartIntrvl = 0.0;
   double nextRestartDump = 0.0;
   char restartDir[kName];
   int optrestart = InFl.ReadItem(optrestart, "RESTARTMODE");

   if (optrestart == 1 || optrestart == 3) {
     restartIntrvl = InFl.ReadItem(restartIntrvl, "RESTARTINTRVL");
     InFl.ReadItem(restartDir, "RESTARTDIR");
     nextRestartDump = timer->getCurrentTime() + restartIntrvl;
   }

	while( !timer->IsFinished() ) {
		
		// Output current time info depending on I/O options
		timer->Advance(timer->getTimeStep());
        if (simCtrl->disp_time == 'Y') {
            PrintRunTimeVars(Moisture, 0);
        }
	
		// Check if precipitation variables have to be updated
		UpdatePrecipitationInput();

		// Simulate Interception, ET processes
		SurfaceHydroProcesses( EvapoTrans, Intercept, SnowPack); // SKY2008Snow from AJR2007
		
		// Simulate Infiltration, Groundwater processes
		SubSurfaceHydroProcesses( Moisture );

		// Simulate Hydraulic/ Hydrologic routing
		Flow->SurfaceFlow();

		// Output various simulated variables
		OutputSimulatedVars( Flow );

		// Update water balance variables
		UpdateWaterBalance( Balance );

		// Update the system 
		Moisture->Reset(); 

#ifdef PARALLEL_TRIBS
      // Reset overlap nodes
      tGraph::resetOverlap();
#endif

      // Write restart files
      if ( (optrestart == 1 || optrestart == 3) && timer->getCurrentTime() >= nextRestartDump) {
          writeRestart(restartDir);
          nextRestartDump += restartIntrvl;
      }
	}
	return;
}

/*****************************************************************************
**  
**  Simulator::end_simulation()
**  
**  Carry out final activities after simulation ends
**
*****************************************************************************/
void Simulator::end_simulation(tKinemat *Flow) 
{ 
	Flow->getResultsPtr()->
			writeAndUpdate( timer->getCurrentTime() );
	
	Flow->getResultsPtr()->
		whenTimeIsOver( timer->getCurrentTime() );
	
	double tend  = timer->getEndTime();
	double spout = timer->getSpatialOutputInterval();
	
	if (!simCtrl->inter_results ||
		(simCtrl->inter_results && spout > tend) ||
		(simCtrl->inter_results && (tend/spout-floor(tend/spout)) > 0))
		
		outp->WriteDynamicVars( timer->getCurrentTime() );
	
	outp->end_simulation();
	
	Cout<<"\nSimulation completed...\n"<<endl;
	
	return;
}

/*****************************************************************************
**  
**  Simulator::PrintRunTimeVars()
**  
**  Prints out the time variables as the simulation progresses
**  
*****************************************************************************/
void Simulator::PrintRunTimeVars(tHydroModel *Moisture, int opt)
{
	if (opt) 
		Cout<<"  "<<timer->year<<"\t"<<timer->month<<"\t"
			<<timer->day<<"\t"<<timer->hour<<"\t"<<timer->minute<<endl;


        if (Moisture->HydroNodesExist()) {
          Cout<<"\n\tx=x=x Current time: "<<timer->getCurrentTime()
            <<" hour x=x=x"<<endl;
        }
        else if (!fmod(timer->getCurrentTime(), timer->getGWTimeStep())) {
          Cout<<"\tx=x=x Current time: "<<timer->getCurrentTime()
                <<" hour x=x=x"<<endl;
        }


}

/*****************************************************************************
**  
**  Simulator::UpdatePrecipitationInput()
**  
**  Handles precipitation input provided to the model either as measured
**  or stochastically created values
**
*****************************************************************************/
void Simulator::UpdatePrecipitationInput()
{
	// Options for radar or rain gauges
	if (rainIn->rainfallType == 1)
		get_next_mrain(simCtrl->mode);
	
	else if (rainIn->rainfallType == 2) {
		if ( timer->isGaugeTime(timer->getRainDT()) ) {
			get_next_gaugerain();  // updates rain to nodes from station data
		}
	}
	return;
}

/*****************************************************************************
**  
**  Simulator::SurfaceHydroProcesses()
**  
**  Calls functions to simulate evapotranspiration, and interception.
**  Handles the time variables for the function calls
**
*****************************************************************************/
void Simulator::SurfaceHydroProcesses(tEvapoTrans *EvapoTrans, 
									  tIntercept  *Intercept, tSnowPack *SnowPack) // SKY2008Snow from AJR2007
{
	// Update meteorological and ET/I time
	get_next_met();

    // SKY2008Snow from AJR2007
	if (SnowPack->getSnowOpt() == 0) {

		// Possible combinations of Evapotrans and Intercept on/off
		// 1) Both ON
		if (EvapoTrans->getEToption() !=0 && Intercept->getIoption() != 0) {
			if ( timer->getCurrentTime() == met_hour ) {
				EvapoTrans->callEvapoPotential();
			}
			if ( timer->getCurrentTime() == eti_hour ) {
				EvapoTrans->callEvapoTrans( Intercept, 1);
			}
		}
		// 2) Interception ON
		if (EvapoTrans->getEToption() == 0 && Intercept->getIoption() != 0) {
			Cout<<"\nInterception Option "<<Intercept->getIoption()
				<<" not valid if "<<endl;
			Cout<<"Evaporation scheme turned off. \n"<<endl;
			Cout<<"Exiting Program...\n\n"<<endl;
			exit(1);
		}
		// 3) ET ON
		if (EvapoTrans->getEToption() !=0 && Intercept->getIoption() == 0) {
			if ( timer->getCurrentTime() == met_hour ) {
				EvapoTrans->callEvapoPotential();
			}
			if ( timer->getCurrentTime() == eti_hour )
				EvapoTrans->callEvapoTrans( Intercept, 0);
		}

	// SKY2008Snow from AJR2007 starts here
	} //end if (no snow)

	else { //snow active

		// ADDED BY RINEHART 2007 @ NMT
		//
		// Possible combinations of Evapotrans and Intercept on/off
		// 1) BOTH ON
		if (SnowPack->getEToption() !=0 && Intercept->getIoption() != 0) {
			if ( timer->getCurrentTime() == eti_hour ) {

				SnowPack->callSnowPack(Intercept,1);
            }
		}
		// 2) INTERCEPTION ON
		if (SnowPack->getEToption() == 0 && Intercept->getIoption() != 0) {

			Cout<<"\nInterception Option "<<Intercept->getIoption()
				<<" not valid if "<<endl;
			Cout<<"Snow scheme turned off." <<endl;
			Cout<<"\nExiting Program...\n\n"<<endl;
			exit(1);
		}
		// 3) ET ON
		if (SnowPack->getEToption() !=0 && Intercept->getIoption() == 0) {
			if ( timer->getCurrentTime() == eti_hour ) {
				SnowPack->callSnowPack(Intercept,0);
			}

		} //evapotrans options
	} //snow option
	// SKY2008Snow from AJR2007 ends here
	 
	return;
}

/*****************************************************************************
**  
**  Simulator::SubSurfaceHydroProcesses()
**  
**  Makes function calls to simulate infiltration and groundwater dynamics
**  
*****************************************************************************/
void Simulator::SubSurfaceHydroProcesses(tHydroModel *Moisture)
{
	// Call Unsaturated Zone in tHydroModel
	Moisture->UnSaturatedZone( timer->getTimeStep() );
    
	GW_label = fmod(timer->getCurrentTime(), timer->getGWTimeStep());
	
	// Call Saturated Zone in tHydroModel 
	if (simCtrl->GW_model_label) {
		if ( !GW_label ) {
			Moisture->ResetGW(); 
			Moisture->SaturatedZone( timer->getGWTimeStep() );
		}
	}
}

/*****************************************************************************
**  
**  Simulator::OutputSimulatedVars()
**  
**  Handles calls to tOutput for writing output files with simulated
**  variables: both pixel and catchment scale
**
*****************************************************************************/
void Simulator::OutputSimulatedVars(tKinemat *Flow)
{
	// If it's necessary -> Output PixelInfo
	if ( ! (fmod(timer->getCurrentTime(), timer->getEtIStep())) ) {
		if ( outp->nodeList )
			outp->WritePixelInfo( timer->getCurrentTime() );
	}

	// Write streamflow for interior outlets
	// TODO: Need to change this later to get an average flow, i.e., 1-hr step
	outp->WriteOutletInfo( timer->getCurrentTime() );
	
	// Write spatial output
	if ( timer->CheckSpatialOutputTime() ) {
		// If it's time -> Output DynVars     
		if ( simCtrl->inter_results )
			outp->WriteDynamicVars( timer->getCurrentTime() );
	}
	return;
}

/*****************************************************************************
**  
**  Simulator::UpdateWaterBalance()
**  
**  Assigns various water balance variables
**  
*****************************************************************************/
void Simulator::UpdateWaterBalance(tWaterBalance *Balance)
{ 
	Balance->UnSaturatedBalance();
	if (!GW_label)
		Balance->SaturatedBalance();
	if (timer->getCurrentTime() == met_hour)
		Balance->CanopyBalance();
	Balance->BasinStorage( timer->getCurrentTime() );
	return;
}

/*****************************************************************************
**  
**  Simulator::get_next_mrain(mode)
**  
**  Get the next measured file name and evaluate duration of rainfall loop 
**
**  Return value: int: error code
**                0: no error
**                -1: Time tag of measured rain smaller than beginning
**                1: Time tag of greater than end
**                10: there is no next file
**  Algorithm:
**   get next measured rainfall name from rain data structure
**
*****************************************************************************/
void Simulator::get_next_mrain(int mode) 
{  
	begin_hour = timer->getCurrentTime(); 
	
	// NODE: Need to redefine lmr_hour -->
	// In this implementation, it searches for the next rainfall file
	// incrementing each time by dtRain. It is assumed that rainfall 
	// for the next found file can be applied to ALL simulation periods 
	// preceding the end of the interval of found rainfall input
	if (lmr_hour < begin_hour && mode==AUTO_INPUT) { 
		
		timer->addRainTime();
		// Unless a file is detected - go through possible list
		while ( !(rainIn->Compose_In_Mrain_Name(timer)) ) { 
			if (count == 0) {
				Cout<<"\nWarning: Next rainfall file "<<rainIn->mrainfileIn
				<<" is missing..."<<endl;
			}
			Cout<<"File "<<rainIn->mrainfileIn<<" was not found..."<<endl;
			
			timer->addRainTime();
			
			if ( timer->getRainTime()-timer->getEndTime() > searchRain ) {
				Cout<<"\nRainfall search threshold exceeded... "<<endl;
				Cout<<count+1<<" rainfall input files are missing..."<<endl; 
				Cout<<"Exiting Program..."<<endl<<endl<<endl;
				exit(2);
			}
			count++;
		}
		
		lmr_hour = timer->getRainTime();
		rainIn->NewRain(timer);
		
		if (simCtrl->Verbose_label == 'Y') {
			Cout<<"Next rainfall input: "<<lmr_hour<<" hours in simulation.\n";
			Cout<<"Unsaturated zone time steps for interval: ";
			Cout<<timer->getElapsedSteps(lmr_hour)<<endl;
		}
	}
	
	else if (lmr_hour < begin_hour && mode==STD_INPUT) { 
		timer->addRainTime();
		
		if ( !(rainIn->Compose_In_Mrain_Name(timer)) ) { 
			Cout<<"\nFile "<<rainIn->mrainfileIn<<" was not found...";
			Cout<<"Exiting Program..."<<endl;
			exit(2);
		}
		
		lmr_hour = timer->getRainTime();
		rainIn->NewRain(timer);  
	}
	// else just use the same intensity values in tCNode
	return;
}

/*****************************************************************************
**  
**  Simulator::get_next_met()
**  
**  Update the meteorological time
**
*****************************************************************************/
void Simulator::get_next_met() 
{    
	if (met_hour < timer->getCurrentTime() ) { 
		timer->addMetTime(1);
		met_hour = timer->getMetTime(1);
	}
	
	if (eti_hour < timer->getCurrentTime() ) { 
		timer->addMetTime(2);
		eti_hour = timer->getMetTime(2);
	}
	return;
}

/*****************************************************************************
**  
**  Simulator::get_next_gaugerain()
**  
**  Call the tRainfall function that gets a new rain gauge value
**
*****************************************************************************/
void Simulator::get_next_gaugerain() 
{
	// rainIn->callRainGauge();
	rainIn->callRainGauge(timer); // SKY2008Snow 
	return;
}

/***************************************************************************
**
** Simulator::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/
void Simulator::writeRestart(char* directory) const
{
  Cout << "WRITE RESTART at time " << timer->getCurrentTime() << endl << endl;

  double   ct     = timer->getCurrentTime();
  int32_t  yr     = static_cast<int32_t>(timer->getYear());
  int32_t  mo     = static_cast<int32_t>(timer->getMonth());
  int32_t  dy     = static_cast<int32_t>(timer->getDay());
  int32_t  hr     = static_cast<int32_t>(timer->getHour());
  int32_t  snOpt  = static_cast<int32_t>(restart->getSnowOpt());

  stringstream sFile;
  sFile << directory << "/tRIBS_Rstrt_";
  sFile << setw(5) << setfill('0') << (int) ct;

#ifdef PARALLEL_TRIBS
  // ---- PARALLEL PATH ----
  // Each rank serializes its outlet records (nodeID+Hlev pairs) and node
  // records to separate in-memory buffers, then rank 0 gathers both and
  // writes a single flat file: header + all outlet records + all node records.
  // On read, every rank opens the same file and filters by nodeID, so the
  // file is completely rank-count-agnostic.

  ostringstream outletBuf(ios::binary);
  restart->writeRestartOutlets(outletBuf);
  string outletData = outletBuf.str();
  int outletLocalSize = static_cast<int>(outletData.size());

  ostringstream nodeBuf(ios::binary);
  restart->writeRestartNodes(nodeBuf);
  string nodeData = nodeBuf.str();
  int nodeLocalSize = static_cast<int>(nodeData.size());

  int numProcs = tParallel::getNumProcs();
  int32_t totalOutlets = static_cast<int32_t>(tParallel::sumBroadcast(restart->getNumOutlets()));
  int32_t totalNodes   = static_cast<int32_t>(tParallel::sumBroadcast(restart->getNumNodes()));

  // Gather outlet buffer sizes and bytes to rank 0.
  vector<int> outletSizes(numProcs, 0);
  MPI_Gather(&outletLocalSize, 1, MPI_INT, outletSizes.data(), 1, MPI_INT,
             MASTER_PROC, MPI_COMM_WORLD);

  vector<int> outletDispls(numProcs, 0);
  int totalOutletBytes = 0;
  if (tParallel::isMaster()) {
    for (int i = 0; i < numProcs; i++) { outletDispls[i] = totalOutletBytes; totalOutletBytes += outletSizes[i]; }
  }
  vector<char> allOutletData;
  if (tParallel::isMaster()) allOutletData.resize(totalOutletBytes);
  MPI_Gatherv(outletData.data(), outletLocalSize, MPI_CHAR,
              allOutletData.empty() ? nullptr : allOutletData.data(),
              outletSizes.data(), outletDispls.data(), MPI_CHAR,
              MASTER_PROC, MPI_COMM_WORLD);

  // Gather node buffer sizes and bytes to rank 0.
  vector<int> nodeSizes(numProcs, 0);
  MPI_Gather(&nodeLocalSize, 1, MPI_INT, nodeSizes.data(), 1, MPI_INT,
             MASTER_PROC, MPI_COMM_WORLD);

  vector<int> nodeDispls(numProcs, 0);
  int totalNodeBytes = 0;
  if (tParallel::isMaster()) {
    for (int i = 0; i < numProcs; i++) { nodeDispls[i] = totalNodeBytes; totalNodeBytes += nodeSizes[i]; }
  }
  vector<char> allNodeData;
  if (tParallel::isMaster()) allNodeData.resize(totalNodeBytes);
  MPI_Gatherv(nodeData.data(), nodeLocalSize, MPI_CHAR,
              allNodeData.empty() ? nullptr : allNodeData.data(),
              nodeSizes.data(), nodeDispls.data(), MPI_CHAR,
              MASTER_PROC, MPI_COMM_WORLD);

  if (tParallel::isMaster()) {
    fstream rStr;
    rStr.open(sFile.str().c_str(), ios::out|ios::binary);
    uint32_t magic = 0x54524853u, schema = 2u, verMaj = 6u, verMin = 4u;
    BinaryWrite(rStr, magic);   BinaryWrite(rStr, schema);
    BinaryWrite(rStr, verMaj);  BinaryWrite(rStr, verMin);
    BinaryWrite(rStr, ct);      BinaryWrite(rStr, yr);
    BinaryWrite(rStr, mo);      BinaryWrite(rStr, dy);
    BinaryWrite(rStr, hr);      BinaryWrite(rStr, snOpt);
    BinaryWrite(rStr, totalNodes);
    BinaryWrite(rStr, totalOutlets);
    rStr.write(allOutletData.data(), totalOutletBytes);
    rStr.write(allNodeData.data(), totalNodeBytes);
    rStr.close();
  }

#else
  // ---- SERIAL PATH ----
  int32_t numN   = static_cast<int32_t>(restart->getNumNodes());
  int32_t numOut = static_cast<int32_t>(restart->getNumOutlets());

  fstream rStr;
  rStr.open(sFile.str().c_str(), ios::out|ios::binary);
  uint32_t magic = 0x54524853u, schema = 1u, verMaj = 6u, verMin = 4u;
  BinaryWrite(rStr, magic);   BinaryWrite(rStr, schema);
  BinaryWrite(rStr, verMaj);  BinaryWrite(rStr, verMin);
  BinaryWrite(rStr, ct);      BinaryWrite(rStr, yr);
  BinaryWrite(rStr, mo);      BinaryWrite(rStr, dy);
  BinaryWrite(rStr, hr);      BinaryWrite(rStr, snOpt);
  BinaryWrite(rStr, numN);    BinaryWrite(rStr, numOut);
  restart->writeRestart(rStr);
  rStr.close();
#endif
}

/***************************************************************************
**
** Simulator::readRestart() Function
**
***************************************************************************/
void Simulator::readRestart(tInputFile &InFl)
{
  char restartFile[kName];
  InFl.ReadItem(restartFile, "RESTARTFILE");

  // All ranks (serial or parallel) open the same file and filter by nodeID.
  // No rank suffix — the file is rank-count-agnostic.
  fstream rStr;
  rStr.open(restartFile, ios::binary|ios::in);
  if (!rStr.is_open()) {
    cerr << "ERROR: Cannot open restart file: " << restartFile << endl;
    exit(1);
  }

  uint32_t magic, schema, verMaj, verMin;
  BinaryRead(rStr, magic);
  if (magic != 0x54524853u) {
    cerr << "ERROR: Restart file is not a valid v6 format (bad magic number).\n"
         << "  v5.x restart files are not compatible with v6.0.0." << endl;
    exit(1);
  }
  BinaryRead(rStr, schema);

#ifdef PARALLEL_TRIBS
  if (schema != 2u) {
    cerr << "ERROR: Restart file schema " << schema
         << " is not a parallel v6 file (expected 2)." << endl;
    exit(1);
  }
#else
  if (schema != 1u) {
    cerr << "ERROR: Restart file schema " << schema
         << " is not supported (expected 1)." << endl;
    exit(1);
  }
#endif

  BinaryRead(rStr, verMaj);
  BinaryRead(rStr, verMin);

  double   fileCt;
  int32_t  yr, mo, dy, hr, snOpt, totalN, totalOut;
  BinaryRead(rStr, fileCt);
  BinaryRead(rStr, yr);
  BinaryRead(rStr, mo);
  BinaryRead(rStr, dy);
  BinaryRead(rStr, hr);
  BinaryRead(rStr, snOpt);
  BinaryRead(rStr, totalN);
  BinaryRead(rStr, totalOut);

  if (snOpt != restart->getSnowOpt()) {
    cerr << "WARNING: Restart file snowOpt=" << snOpt
         << " differs from current config snowOpt=" << restart->getSnowOpt()
         << ". Snow state may be inconsistent." << endl;
  }

#ifdef PARALLEL_TRIBS
  int32_t globalNodes = static_cast<int32_t>(tParallel::sumBroadcast(restart->getNumNodes()));
  if (totalN != globalNodes) {
    cerr << "ERROR: Restart file node count (" << totalN
         << ") does not match mesh (" << globalNodes << ")." << endl;
    exit(1);
  }
#else
  if (totalN != restart->getNumNodes()) {
    cerr << "ERROR: Restart file node count (" << totalN
         << ") does not match mesh (" << restart->getNumNodes() << ")." << endl;
    exit(1);
  }
#endif

  Cout << "Restart loaded: "
       << yr << "-" << setfill('0') << setw(2) << mo << "-" << setw(2) << dy
       << " " << setw(2) << hr << ":00"
       << " (" << fixed << setprecision(3) << fileCt << " simulation hrs elapsed)"
       << endl << endl;

  // Each rank reads the full outlet block and node block, applying only the
  // records whose nodeIDs belong to this rank.
  restart->readRestartGlobal(rStr, totalOut, totalN);

  rStr.close();
}

//=========================================================================
//
//
//                          End of tSimul.cpp
//
//
//=========================================================================
