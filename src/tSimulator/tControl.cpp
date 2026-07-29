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
**  tControl.cpp: Functions for class SimulationControl 
**              (see tControl.h)
**
***************************************************************************/

#include "src/tSimulator/tControl.h"
#include "src/Headers/globalIO.h"

#ifdef PARALLEL_TRIBS
#include "src/tParallel/tParallel.h"
#endif

//=========================================================================
//
//
//                  Section 1: SimulationControl Constructors/Destructors
//
//
//=========================================================================

SimulationControl::SimulationControl(int argc, char **argv)
{  
	Cout<<"\n-----------------------------------------------------------------"
	<<"-------";
	Cout<<"\n\n\t\t tRIBS Distributed Hydrologic Model";
	Cout<<"\n\t\t TIN-based Real-time Integrated Basin Simulator \n\n";
	Cout<<"-----------------------------------------------------------------"
		<<"-------"<<endl;
	
	char VERSION[50] = "6.0.0, Summer 2026";
	
	Cout <<"\n\ntRIBS Version "<< VERSION <<endl<<endl;
	
	static char usage[]=
		"Usage : %s [-A] [-V NodeID] [-T]\n";

	
	mode = STD_INPUT;       //Default: If file doesn't exist, assume zero rainfall
	Verbose_label    = 'N';
    debug            = 'N';
    disp_time = 'N';
	VerbID = -999;
	
#ifdef PARALLEL_TRIBS
  // Command line args are broadcast to all nodes
  tParallel::inputArgs(argc, argv);
#endif

	if( argc < 2 ){
		Cout<<"\n\nUsage: " << argv[0] <<" <input file>  [options]"<<endl<<endl;
		Cout<<"Options: "<<endl;
		Cout<<"\t-A    Automatic listing of rainfall files (zero if missing)"<<endl;
		Cout<<"\t-V [NodeID] Verbose mode (output run-time information)"<<endl<<endl;
		Cout<<"Provide name of an input file. Exiting program...\n"<<endl;
		exit(1);
	}
	
	infile = argv[1];  		//Name of input file 
	
    /* WR 08282023 removed command line arguments and specified in input file
    "OPTGROUNDWATER" -G    Run groundwater model: GW_model_label
    "OPTSPATIAL" -R    Write intermediate states (spatial output): inter_results
    */

	for (int i=2; i < argc; i++)  //Start from the 3rd argument
    { 
		if (*argv[i] != '-'){ 
			if (i > 2) {
				if ((*(argv[i-1]+1)) == 'V') {
					sscanf(argv[i],"%d",&VerbID);
					cout<<"\nVerbose Node ID = "<<VerbID<<endl;
				}
				else {
					printf("Simulation control: Syntax incorrect %s\n", argv[i]);
					printf(usage, argv[0]);
					exit(1);
				}
			}
			else {
				printf("Simulation control: Syntax incorrect %s\n", argv[i]);
				printf(usage, argv[0]);
				exit(1);
			}
		}
		
		switch(*(argv[i]+1)){ 
			case 'A': 		//Automatic listing
			{ 
				mode=AUTO_INPUT;      // if file does not exist -> assumes '0'-s!
				break; 
			}
			case 'V':                 //Verbose Output
			{
				Verbose_label = 'Y';
				break;
			}
            case 'T':                        //For tGraph debugging output
            {
                disp_time = 'Y';
                break;
            }
            case 'D':                        //For tGraph debugging output
            {
                debug = 'Y';
                break;
            }
      }
   }
}

SimulationControl::~SimulationControl() 
{  
	Cout<<"SimulationControl Object has been destroyed..."<<endl<<flush;
	
	Cout<<"\n-----------------------------------------------------------------"
		<<"-------";
	Cout<<"\n\n\t\t tRIBS Distributed Hydrologic Model";
	Cout<<"\n\t\t TIN-based Real-time Integrated Basin Simulator";
	Cout<<"\n\t\t Ralph M. Parsons Laboratory";
	Cout<<"\n\t\t Massachusetts Institute of Technology \n\n";
	Cout<<"-----------------------------------------------------------------"
		<<"-------"<<endl;
}

//=========================================================================
//
//
//                          End of tControl.cpp
//
//
//=========================================================================
