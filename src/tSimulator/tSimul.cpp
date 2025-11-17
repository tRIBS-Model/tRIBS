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

#include <sstream>
#include "src/tSimulator/tSimul.h"
#include "src/Headers/TemplDefinitions.h"
#include "src/Headers/globalIO.h"

#ifdef PARALLEL_TRIBS
#include "src/tGraph/tGraph.h"
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
	
	// For forecasted rainfall, turned off
	lfr_hour = 0; 
	
	// Counter of time, used for GW model
	GW_label = 0.;

	// Output data on Mesh and Voronoi elements
	outp->WriteOutput( 0 );
	
	// Get rainsearch if rainfall used
	if (rainIn->getoptStorm() == 0) {
		if (rainIn->rainfallType == 1 || rainIn->rainfallType == 2) {
			searchRain = rainIn->searchRain;     // Rainfall search threshold
		}
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
    "OPTINTERHYDRO")-H    Write intermediate hydrographs (.mrf): hydrog_results
    "OPTHEADER"); -M    Do NOT Write headers in pixel/hydrograph/voronoi output files: : Header_label
    */

    if (InFl.IsItemIn( "OPTGROUNDWATER" ))
        simCtrl->GW_model_label = InFl.ReadItem(simCtrl->GW_model_label, "OPTGROUNDWATER");
    else
        simCtrl->GW_model_label = true; //Default option

    if (InFl.IsItemIn( "OPTSPATIAL" ))
        simCtrl->inter_results = InFl.ReadItem(simCtrl->inter_results, "OPTSPATIAL");
    else
        simCtrl->inter_results = false; //Default option

    if (InFl.IsItemIn( "OPTINTERHYDRO" ))
        simCtrl->hydrog_results = InFl.ReadItem(simCtrl->hydrog_results, "OPTINTERHYDRO");
    else
        simCtrl->hydrog_results = false; //Default option

	// Ouput pre-processing
	if (simCtrl->inter_results)
		outp->CreateAndOpenDynVar();

    //WR debug 01032024: this was setting the met forcing values to 0 at time 0 in the Dynamic and Pixel files
	// Output initial conditions
	//outp->WriteDynamicVars( timer->getCurrentTime() );
	//outp->WritePixelInfo(   timer->getCurrentTime() );

	// Prepare rainfall input if stochastic rainfall Off
	if ( !rainIn->getoptStorm()) {
		if (rainIn->rainfallType == 1 || rainIn->rainfallType == 2) {
			
			// Check if time for rainfall forecast 
			if (fmod(timer->getCurrentTime(), timer->getRainDT())==0 &&
				timer->getoptForecast()!=0)
				fState = checkForecast();
			
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
	}
	
	// Stochastic Rainfall Initialization
	else {
		// Call storm generator, update time, and set rainfall
		rainIn->GenerateStorm( timer->getCurrentTime() );
		timer->UpdateStorm( rainIn->getStormDuration() + rainIn->interstormDur() );
		rainIn->NewRain(rainIn->getRainrate());
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
		UpdatePrecipitationInput( rainIn->getoptStorm() );

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

		// End simulation beyond forecast interval
		if ( !rainIn->getoptStorm() ) {
			if (timer->getoptForecast() != 0 && fState == 3)
				break; 
		}

      // Write restart files
      if ( optrestart > 0 && timer->getCurrentTime() == nextRestartDump) {
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
	if ( !simCtrl->hydrog_results )
		Flow->getResultsPtr()->
			writeAndUpdate( timer->getCurrentTime(), 0 );
	
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
void Simulator::UpdatePrecipitationInput(int opt)
{
	// ==============================================
	// For measured radar or raingauge rainfall input
	if ( !opt ) { 
		
		// Check if time for rainfall forecast
		if (fmod(timer->getCurrentTime(), timer->getRainDT())==0 && 
			timer->getoptForecast()!=0)
			fState = checkForecast();
		
		// Options for radar or rain gauges 
		if (rainIn->rainfallType == 1 || rainIn->rainfallType == 2)
			get_next_mrain(simCtrl->mode);
		
		else if (rainIn->rainfallType == 3) {
			if ( timer->isGaugeTime(timer->getRainDT()) ) {
				get_next_gaugerain();  // updates rain to nodes from station data
			}
		}
	}
	// ==============================================
	// For stochastic rainfall input
	else {  
		
		// Within Storm period
		if (timer->getCurrentTime() <= 
			(timer->getStormTime()-rainIn->interstormDur()))
			rainIn->NewRain(rainIn->getRainrate());        //Set Rainfall to value
		
		// Within Intestorm period
		else if (timer->getCurrentTime() > 
				 (timer->getStormTime()-rainIn->interstormDur()) && 
				 timer->getCurrentTime() <= timer->getStormTime()) {
			rainIn->NewRain();                            //Set Rainfall to zero
			rainIn->setRainrate(0.0);
		}
		
		// After interstorm period: call storm generator, update time, set rainfall
		else if (timer->getCurrentTime() >= timer->getStormTime()) {
			rainIn->GenerateStorm( timer->getCurrentTime() );	
			timer->UpdateStorm( rainIn->getStormDuration() + rainIn->interstormDur() );
			rainIn->NewRain(rainIn->getRainrate());
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
			if (Intercept->getIoption() == 1) {
				if ( timer->getCurrentTime() == eti_hour )
					EvapoTrans->callEvapoTrans( Intercept, 1 );
			}
			else {
				Cout<<"\nInterception Option "<<Intercept->getIoption()
				   <<" not valid if "<<endl;
				Cout<<"Evaporation scheme turned off. \n\tPlease use:"<<endl;
				Cout<<"\t\t(1) for Gray (1970) Method: Two Parameter Model"<<endl;
				Cout<<"Exiting Program...\n\n"<<endl;
				exit(1);
			}
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
			if ( timer->getCurrentTime() == met_hour ) {

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
			if ( timer->getCurrentTime() == met_hour ) {
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
	int forenum;

	// If it's necessary -> Output PixelInfo
	if ( ! (fmod(timer->getCurrentTime(), timer->getEtIStep())) ) {
		if ( outp->nodeList )
			outp->WritePixelInfo( timer->getCurrentTime() );
	}

	// Write streamflow for interior outlets
	// TODO: Need to change this later to get an average flow, i.e., 1-hr step
	outp->WriteOutletInfo( timer->getCurrentTime() );
	
	// If it's time -> Output Hydrograph  
	if ( timer->CheckOutputTime() ) {
		if (simCtrl->fore_rain_label == 'N')
			forenum=0;
		else 
			forenum=1;
        if ((simCtrl->hydrog_results) && (timer->getCurrentTime())) {
            Flow->getResultsPtr()->
                    writeAndUpdate( timer->getCurrentTime(), forenum );
        }
		
		// Write selected dynamic variables
		// if ( simCtrl->inter_results == 'Y' )
		//   outp->WriteDynamicVar( timer->getCurrentTime() );
	}
	
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

/*****************************************************************************
**  
**  Simulator::checkForecast()
**  
**  Check the forecast state. Returns integer representing state:
**  
**  0 = Before and up to forecast time, Use QPE
**  1 = In Forecast Period and up to lead time, Use QPF
**  2 = In Forecast Period and after lead time, Use Average Rainfall
**  3 = After Forecast Period, End simulation
**
*****************************************************************************/
int Simulator::checkForecast() 
{
	int state;
	
	if (timer->getCurrentTime() < timer->getfTime())
		state = 0;
	else if (timer->getCurrentTime() < (timer->getfTime() + timer->getfLead()) &&
			 timer->getCurrentTime() >= timer->getfTime())
		state = 1;
	else if (timer->getCurrentTime() < (timer->getfTime() + timer->getfLength()) &&
			 timer->getCurrentTime() >= timer->getfLead())
		state = 2;
	else if (timer->getCurrentTime() >= (timer->getfTime() + timer->getfLength()))
		state = 3;
	
	rainIn->setfState(state);
	
	return state;
}

/*****************************************************************************
**  
**  Simulator::check_mod_status()
**  
**  Checks if the model must stay on after a run by checking one of the 
**  SimulationControl flags: 0 - 'NO', 1 - 'YES'
**
*****************************************************************************/
int Simulator::check_mod_status() 
{ 
	if (simCtrl->mod_is_on == 'Y')
		return 1;
	else
		return 0;
}

/*****************************************************************************
**  
**  Simulator:: RunItAgain()
**
**  To run the model using previously constructed mesh and assigned to it 
**  various properties e.g. soils, landuse, etc.
**  
**  Algorithm: Re-initialize all the objects whithout deleting these
**           objects, accessing their data members through the functions
**           used in their constructors
**
*****************************************************************************/
void Simulator::RunItAgain( tInputFile &InFl, tHydroModel *Moisture, 
							tKinemat *Flow, tEvapoTrans *EvapoTrans, 
							tIntercept *Intercept, tWaterBalance *Balance,
							tPreProcess *PreProcessor, tSnowPack *SnowPack) // SKY2008Snow from AJR2007
{ 
	char wish = 'Z';
	char keep = 'Z';
	char filein[80];
	char yesno[20];
	
	simCtrl->num_simul++;
	
	cerr<<"\n\n----------------------------------------------------"<<endl
		<<"\tMODEL RUN #"<<simCtrl->num_simul <<" COMPLETED\n"<<endl
		<<"\n\tDo you want to continue (type 'y' or 'n')?\n\tEnter Option: "<<flush;
	
	while ( wish == 'Z' ) {
		cin>>yesno;
		
		if ((yesno[0] == 'n' || yesno[0] == 'N') && (yesno[1] == '\0')) {
			cerr<<"\nProgram Finishing..."<<endl<<flush;
			cerr<<"----------------------------------------------------"<<endl;
			simCtrl->mod_is_on = 'N';
			return;
		}
		else if ((yesno[0] == 'y' || yesno[0] == 'Y') && (yesno[1] == '\0'))
			wish = 'y';
		else {
			wish = 'Z';
			cerr<<"\n\tCommand not understood... \n\tEnter Option: "<<flush;
		}
	}
	cerr<<"\n\tEnter input data filename: "; 
	cin>>filein;
	
	ifstream source( filein );
	while ( !source && simCtrl->mod_is_on == 'Y') {
		cerr<<"\n\tFile does not exist... Check spelling...\n"
		<<"\t (To Exit, please type 'n') \n"
		<<"\n\tEnter input data filename: "; 
		cin >> filein;
		if ((filein[0] == 'n' || filein[0] == 'N') && (filein[1] == '\0')) {
			simCtrl->mod_is_on = 'N';
			cerr<<"\nProgram Finishing..."<<endl<<flush;     
			cerr<<"----------------------------------------------------"<<endl;
			return;
		}
		source.open( filein );
	}
	source.close();
	simCtrl->infile = filein;
	
	cerr<<endl<<flush;
	cerr<<"\tName of *.in file: '"<<simCtrl->infile<<"'"<<endl<<flush;
	cerr<<endl<<flush;
	cerr<<"\tPlease indicate if soil and landuse maps are changed\n"
		<<"\t(if so, new resampling will need to be carried out)\n"
		<<"\n\tPlease type ('y' or 'n'): "<<flush;
	
	wish = 'Z';
	while ( wish == 'Z' ) {
		cin>>yesno;
		
		if ((yesno[0] == 'n' || yesno[0] == 'N') && (yesno[1] == '\0'))
			wish = 0;
		else if ((yesno[0] == 'y' || yesno[0] == 'Y') && (yesno[1] == '\0'))
			wish = 1;
		else {
			wish = 'Z';
			cerr<<"\n\tCommand not understood. \n\tEnter Option: "<<flush;
		}
	}
	
	cerr<<endl<<flush;
	cerr<<"\tPlease indicate if state of the system should be used\n"
		<<"\tas initial condtion for the next run\n"
		<<"\tPlease type ('y' or 'n'): "<<flush;
	
	keep = 'Z';
	while ( keep == 'Z' ) {
		cin>>yesno;
		
		if ((yesno[0] == 'n' || yesno[0] == 'N') && (yesno[1] == '\0'))
			keep = 0;
		else if ((yesno[0] == 'y' || yesno[0] == 'Y') && (yesno[1] == '\0'))
			keep = 1;
		else {
			keep = 'Z';
			cerr<<"\n\tCommand not understood. \n\tType 'y' or 'n': "<<flush;
		}
	}
	cerr<<endl<<flush;
	cerr<<"----------------------------------------------------"<<endl;
	
	if ( !wish ) {
		cerr<<"\n\tNOTE: Previous soil and landuse maps are used. Continuing..."
		<<endl<<flush;
	}
	
	cout<<"\n---------------------------------------------------\n"
		<<"\tA new tRIBS run has been initiated..."<<endl<<flush;
	cout<<"-----------------------------------------------------\n";
	
	simCtrl->infile = filein;
	cerr<<"\n\tName of *.in file: '"<<simCtrl->infile<<"'"<<endl<<flush;
	
	// Re-initializing tInputFile
	InFl.CloseOldAndOpenNew( simCtrl->infile ); // Close previous IN file and startnew
	
	// Check validity of the input file 
	PreProcessor->CheckInputFile( InFl );
	
	// Re-initializing tRunTimer 
	timer->InitializeTimer( InFl );
	
	// Re-initializing tOutput 
	outp->UpdateForNewRun( InFl );
	
	// Re-initialize tWaterBalance
	Balance->DeleteWaterBalance();
	Balance->SetWaterBalance( InFl );
	
	// Re-initializing tFlowNet 
	Flow->SetFlowVariables( InFl );
	Flow->setTravelVelocity( 0.0 );
	Flow->initializeTravelTimeOnly();
	Flow->UpdateForNewRun( InFl , keep);
	
	// Re-initializing tFlowResults 
	Flow->getResultsPtr()->free_results();
	Flow->getResultsPtr()->SetFlowResVariables( InFl, Flow->MaxTravel() );
	
	// Re-initializing tInvariant & tHydroModel 
	Moisture->soilPtr->SetSoilParameters(rainIn->getMeshPtr(),rainIn->getRsmplPtr(),
										 InFl, wish );
	Moisture->landPtr->SetLtypeParameters(rainIn->getMeshPtr(),rainIn->getRsmplPtr(), 
										  InFl, wish );
	Moisture->SetHydroMVariables( InFl, rainIn->getRsmplPtr(), keep );
	
	// Re-initializing tRainfall 
	rainIn->SetStormVariables( InFl );
	if ( !rainIn->getoptStorm() )
		rainIn->SetRainVariables( InFl );
	
	// Re-initializing tEvapoTrans 
	EvapoTrans->DeleteEvapoTrans();
	EvapoTrans->SetEvapTVariables( InFl, Moisture );
	
	// Re-initializing tIntercept 
	Intercept->SetIntercpVariables( InFl, Moisture );
	
	// Re-initializing Simulator 
	begin_hour = timer->getCurrentTime();
	GW_label = 0.;
	
	// Initialize simulation
	initialize_simulation(EvapoTrans, SnowPack, InFl ); 
       //SMM 09252008 added parameters
	
	// Start simulation
	simulation_loop( Moisture, Flow, EvapoTrans, Intercept, Balance, SnowPack, InFl); // SKY2008Snow from AJR2007
	
	// Finish simulation 
	end_simulation( Flow );
	
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

  fstream rStr;
  stringstream sFile;
  sFile << directory << "/tRIBS_Rstrt_";
  sFile << setw(5) << setfill('0') << (int) timer->getCurrentTime();

#ifdef PARALLEL_TRIBS
  sFile << "_" << tParallel::getMyProc();
#endif 

  rStr.open(sFile.str().c_str(), ios::out|ios::binary);

  // Dump local simulator information
  BinaryWrite(rStr, count);
  BinaryWrite(rStr, fState);
  BinaryWrite(rStr, dt_rain);
  BinaryWrite(rStr, lfr_hour);
  BinaryWrite(rStr, lmr_hour);
  BinaryWrite(rStr, begin_hour);
  BinaryWrite(rStr, met_hour);
  BinaryWrite(rStr, eti_hour);
  BinaryWrite(rStr, GW_label);
  BinaryWrite(rStr, searchRain);

  // Dump information from objects controlled by tRestart
  restart->writeRestart(rStr);

  rStr.close();
}

/***************************************************************************
**
** Simulator::readRestart() Function
**
***************************************************************************/
void Simulator::readRestart(tInputFile &InFl)
{
  Cout << "READ RESTART at time " << timer->getCurrentTime() << endl << endl;

  char restartFile[kName];
  InFl.ReadItem(restartFile, "RESTARTFILE");

  fstream rStr;
  stringstream sFile;
  sFile << restartFile;

#ifdef PARALLEL_TRIBS
  sFile << "_" << tParallel::getMyProc();
#endif

  rStr.open(sFile.str().c_str(), ios::binary|ios::in);

  // Read local simulator information
  BinaryRead(rStr, count);
  BinaryRead(rStr, fState);
  BinaryRead(rStr, dt_rain);
  BinaryRead(rStr, lfr_hour);
  BinaryRead(rStr, lmr_hour);
  BinaryRead(rStr, begin_hour);
  BinaryRead(rStr, met_hour);
  BinaryRead(rStr, eti_hour);
  BinaryRead(rStr, GW_label);
  BinaryRead(rStr, searchRain);

  // Read information from objects controlled by tRestart
  restart->readRestart(rStr);

  rStr.close();
}

//=========================================================================
//
//
//                          End of tSimul.cpp
//
//
//=========================================================================
