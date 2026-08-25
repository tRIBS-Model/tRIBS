/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/


#include "src/Headers/globalIO.h"
#include "Headers/Inclusions.h"
#include "src/Mathutil/predicates.h"
#include "src/tFlowNet/tKinemat.h"
#include "tFlowNet/tReservoir.h" // JECR2014
#include "src/tRasTin/tRainfall.h"
#include "src/tRasTin/tShelter.h" // SKY2008Snow from AJR2007
#include "src/tSimulator/tSimul.h"
#include "src/Headers/TemplDefinitions.h"
#include "src/tHydro/tSnowPack.h" // SKY2008Snow from AJR2007


#ifdef PARALLEL_TRIBS
#include "tGraph/tGraph.h"
#include "tParallel/tParallel.h"
#endif

using namespace std;

Predicates predicate;

// Create global output streams for serial/parallel
tOstream Cout(cout);
tOstream Cerr(cerr);

int serialSimulation(int, char**);
int parallelSimulation(int, char**);

//=========================================================================
//
//
//                  Section 1: Main Program
//
//
//=========================================================================

int main( int argc, char **argv )
{  
#ifdef PARALLEL_TRIBS
	parallelSimulation(argc, argv);
#else
	serialSimulation(argc, argv);
#endif
}

//=========================================================================
//
// Serial simulation will process all nodes on a single processor so
// tGraph is not called to partition nodes and only read parts of the
// graph belonging to a processor
//
//=========================================================================

int serialSimulation( int argc, char **argv )
{
	// Check command-line arguments
	SimulationControl SimCtrl(argc, argv);
	
	Cout<<"\n\nPart 1: Read Input Parameters from "<<"'"<<argv[1]<<"'"<<endl;
	Cout<<"---------------------------------------------------"<<endl;
	tInputFile InputFile( SimCtrl.infile );
	
	// Preprocessing meteorological data
	tPreProcess PreProcessor( InputFile );

	// Partition-only mode needs the MPI process count and tGraph, neither of
	// which exist in the serial binary
	int optpar;
	optpar = InputFile.ReadItem( optpar, "PARALLELMODE" );
	if (optpar == 2) {
		Cerr<<"\nPARALLELMODE = 2 (partition-only) requires the parallel binary:"
			<<"\n   mpirun -np <N> tRIBSpar <input file>"
			<<"\nExiting."<<endl;
		exit(1);
	}

	// Timer, checks environmental variables
	tRunTimer Timer( InputFile );

	// Creating Mesh
	// Option 9 meshbuilder files can be handled like others in serial
	// so no special order has to be imposed as it is with the parallel simulation
	tMesh<tCNode> BasinMesh( &SimCtrl, InputFile );
	
	Cout<<"\n\nPart 3: Creating Stream Network from Mesh "<<endl;
	Cout<<"---------------------------------------------"<<endl;
	tKinemat Flow( &SimCtrl, &BasinMesh, InputFile, &Timer );

	Cout<<"\n\nPart 4: Creating Resampling Object"<<endl;
	Cout<<"--------------------------------------"<<endl;
	tResample RsmplMaster( &SimCtrl, &BasinMesh ); 

	Cout<<"\n\nPart 4a: Creating Sheltering Object for DEM Input" << endl;
	Cout<<"----------------------------------------------------"<<endl;
	tShelter Shelter( &SimCtrl, &BasinMesh, InputFile );// SKY2008Snow from AJR2007

	Cout<<"\n\nPart 5: Creating Output Files from "<<"'"<<argv[1]<<"'"<<endl;
	Cout<<"------------------------------------------"<<endl;
	tCOutput<tCNode> Output( &SimCtrl, &BasinMesh, InputFile, 
							 &RsmplMaster, &Timer );
	
	// Creating WaterBalance class
	tWaterBalance Balance( &SimCtrl, &BasinMesh, InputFile);
	
	Cout<<"\n\nPart 6: Creating Hydrologic System"<<endl;
	Cout<<"--------------------------------------"<<endl;
	tHydroModel Moisture( &SimCtrl, &BasinMesh, InputFile, 
						  &RsmplMaster, &Balance, &Timer );
	
	Cout<<"\nCreating Rainfall Setup...\n";
	tRainfall Rainfall( &SimCtrl, &BasinMesh, InputFile, &RsmplMaster, &Timer );
	
	Cout<<"\nCreating EvapoTranspiration Setup...\n";
	tEvapoTrans EvapoTrans( &SimCtrl, &BasinMesh, InputFile, &Timer, 
							&RsmplMaster, &Moisture, &Rainfall);
	
	Cout<<"\nCreating Interception Setup...\n";
	tIntercept Intercept( &SimCtrl, &BasinMesh, InputFile, &Timer, 
						  &RsmplMaster, &Moisture ); 

	Cout<<"\nInitializing SnowPack setup... \n";
	tSnowPack SnowPack( &SimCtrl, &BasinMesh, InputFile, &Timer, &RsmplMaster, &Moisture, &Rainfall); // SKY2008Snow from AJR2007

	Cout<<"\nCreating Restart setup...\n";
	tRestart<tCNode> Restart( &Timer, &BasinMesh, &Flow, &SnowPack);

	Cout<<"\n\nPart 7: Creating and Initializing Simulation"<<endl;
	Cout<<"------------------------------------------------"<<endl<<endl;
	Simulator Simulant( &SimCtrl, &Rainfall, &Timer, &Output, &Restart );
	Simulant.initialize_simulation(&EvapoTrans, &SnowPack, InputFile); 
	
	// Simulation starts new or from restart file
	int optrestart;
	optrestart = InputFile.ReadItem( optrestart, "RESTARTMODE");
	if (optrestart == 2 || optrestart == 3 ) {
		Simulant.readRestart(InputFile);
		EvapoTrans.flagRestartStart();
		SnowPack.flagRestartStart();
	}

	Cout<<"\n\nPart 8: Hydrologic Simulation Loop"<<endl;
	Cout<<"--------------------------------------"<<endl;
	Simulant.simulation_loop( &Moisture, &Flow, &EvapoTrans,
							  &Intercept, &Balance, &SnowPack, // SKY2008Snow from AJR2007
							  InputFile); // SKY2008Snow
	Simulant.end_simulation( &Flow );

	Cout<<"\n\nPart 9: Deleting Objects and Exiting Program"<<endl;
	Cout<<"------------------------------------------------"<<endl<<endl;
	return 0; // 04/07/2020 Added this to eliminate warning Clizarraga
}

//=========================================================================
//
// Parallel simulation first partitions reaches and then reads the portion
// of the files with nodes belonging to the partition
//
//=========================================================================

int parallelSimulation(int argc, char **argv)
{
#ifdef PARALLEL_TRIBS
	// Startup parallel communications
	tParallel::initialize(argc, argv);

	// Check command-line arguments
	SimulationControl SimCtrl(argc, argv);

	Cout<<"\n\nPart 1: Read Input Parameters from "<<"'"<<argv[1]<<"'"<<endl;
	Cout<<"---------------------------------------------------"<<endl;
	tInputFile InputFile( SimCtrl.infile );
        
	// Preprocessing meteorological data
	tPreProcess PreProcessor( InputFile );

	// Timer, checks environmental variables 
	tRunTimer Timer( InputFile );

	// Build the full mesh, then the stream network, then partition the reach
	// graph. tGraph reads the partition from GRAPHFILE when one is provided;
	// otherwise it generates the partition in-process via METIS and writes
	// it next to the model outputs (OUTFILENAME).
	int numberOfProcessors = tParallel::getNumProcs();

	// PARALLELMODE 2 = partition-only: build the mesh and stream network,
	// generate and write the .reach graphfile, print partition statistics,
	// and exit without running the hydrologic simulation
	int optpar;
	optpar = InputFile.ReadItem( optpar, "PARALLELMODE" );
	bool partitionOnly = (optpar == 2);

	{

	tMesh<tCNode> BasinMesh( &SimCtrl, InputFile );

	Cout<<"\n\nPart 3: Creating Stream Network from Mesh "<<endl;
	Cout<<"---------------------------------------------"<<endl;
	tKinemat Flow( &SimCtrl, &BasinMesh, InputFile, &Timer );

	Cout<<"\n\nPart 3b: Creating Stream Reach partitioning "<<endl;
	Cout<<"-----------------------------------------------"<<endl;
	tGraph::initialize( &SimCtrl, &BasinMesh, &Flow, InputFile, partitionOnly);

	if (partitionOnly) {
		if (numberOfProcessors == 1)
			Cout<<"\nNote: run with 1 MPI process, so all reaches are in "
				<<"partition 0 and no graphfile is written."
				<<"\nRun 'mpirun -np <N> tRIBSpar <input file>' to build a "
				<<"graphfile for N partitions."<<endl;
		Cout<<"\nPartition-only run (PARALLELMODE = 2): "
			<<"skipping the hydrologic simulation."<<endl;
	}
	else {

	Cout<<"\n\nPart 4: Creating Resampling Object"<<endl;
	Cout<<"--------------------------------------"<<endl;
	tResample RsmplMaster( &SimCtrl, &BasinMesh );  

	Cout<<"\n\nPart 4a: Creating Sheltering Object for DEM Input" << endl;
	Cout<<"----------------------------------------------------"<<endl;
	tShelter Shelter( &SimCtrl, &BasinMesh, InputFile );// SKY2008Snow from AJR2007

	Cout<<"\n\nPart 5: Creating Output Files from "<<"'"<<argv[1]<<"'"<<endl;
	Cout<<"------------------------------------------"<<endl;
	tCOutput<tCNode> Output( &SimCtrl, &BasinMesh, InputFile, 
                             &RsmplMaster, &Timer );     

	// Creating WaterBalance class
	tWaterBalance Balance( &SimCtrl, &BasinMesh, InputFile);

	Cout<<"\n\nPart 6: Creating Hydrologic System"<<endl;
	Cout<<"--------------------------------------"<<endl;
	tHydroModel Moisture( &SimCtrl, &BasinMesh, InputFile,
                              &RsmplMaster, &Balance, &Timer );

	Cout<<"\nCreating Rainfall Setup...\n";
	tRainfall Rainfall( &SimCtrl, &BasinMesh, InputFile, &RsmplMaster, &Timer );

	Cout<<"\nCreating EvapoTranspiration Setup...\n";
	tEvapoTrans EvapoTrans( &SimCtrl, &BasinMesh, InputFile, &Timer,
                                &RsmplMaster, &Moisture, &Rainfall);

	Cout<<"\nCreating Interception Setup...\n";
	tIntercept Intercept( &SimCtrl, &BasinMesh, InputFile, &Timer, 
                          &RsmplMaster, &Moisture ); 

	Cout<<"\nInitializing SnowPack setup... \n";
	tSnowPack SnowPack( &SimCtrl, &BasinMesh, InputFile, &Timer, &RsmplMaster, 
                        &Moisture, &Rainfall); // SKY2008Snow from AJR2007

	Cout<<"\nCreating Restart setup...\n";
	tRestart<tCNode> Restart( &Timer, &BasinMesh, &Flow, &SnowPack);

	Cout<<"\n\nPart 7: Creating and Initializing Simulation"<<endl;
	Cout<<"------------------------------------------------"<<endl<<endl;
	Simulator Simulant( &SimCtrl, &Rainfall, &Timer, &Output, &Restart );
	Simulant.initialize_simulation(&EvapoTrans, &SnowPack, InputFile); 

	// Simulation starts new or from restart file
	int optrestart;
	optrestart = InputFile.ReadItem( optrestart, "RESTARTMODE");
	if (optrestart == 2 || optrestart == 3 ) {
		Simulant.readRestart(InputFile);
		EvapoTrans.flagRestartStart();
		SnowPack.flagRestartStart();
	}

	Cout<<"\n\nPart 8: Hydrologic Simulation Loop"<<endl;
	Cout<<"--------------------------------------"<<endl;
	Simulant.simulation_loop( &Moisture, &Flow, &EvapoTrans,
                              &Intercept, &Balance, &SnowPack, // SKY2008Snow from AJR2007
                              InputFile); // SKY2008Snow
	Simulant.end_simulation( &Flow );

	}

	}

	Cout<<"\n\nPart 9: Deleting Objects and Exiting Program"<<endl;
	Cout<<"------------------------------------------------"<<endl<<endl;

  // Finalize graph partitioning
  tGraph::finalize();

  // Finalize parallel communications
  tParallel::finalize();

  return(1);
#endif
  return 0; // Added this to eliminate warning 04/07/2020 Clizarraga
}

//=========================================================================
//
//
//                         End of main.cpp
//
//
//=========================================================================
