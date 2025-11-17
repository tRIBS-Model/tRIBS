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
**  tFlowResults.cpp: Functions for class tFlowResults
**			  (see tFlowResults.h)
**
***************************************************************************/

#include "src/tFlowNet/tFlowResults.h"
#include "src/Headers/globalFns.h"
#include "src/Headers/globalIO.h"

#ifdef PARALLEL_TRIBS
#include "src/tParallel/tParallel.h"
#endif

//=========================================================================
//
//
//                  Section 1: tFlowResults Constructors/Destructors
//
//
//=========================================================================

/***************************************************************************
** 
**  tFlowResults Constructor
**
**  Allocate memory and initialize the tFlowResults Class
**
**  Algorithm:
**	set timer *
**	results step of end simulation + add_time
**	Allocate memory for the hydrographs
**	Initialize hydrographs to zero
**
**  ribsOutput used as flag for compatibility with RIBS user interphase
**  visualization of the previous and current hydrograph. 
**              0: current hydrograph only
**              1: previous and current hydrograph
** 
***************************************************************************/
tFlowResults::tFlowResults(SimulationControl *simControl, tInputFile &infile, 
						   tRunTimer *timptr, double add_time)
{
	timer   = timptr;
	simCtrl = simControl;
	SetFlowResVariables(infile, add_time);
}

tFlowResults::~tFlowResults()
{
	timer = NULL;
	simCtrl = NULL;
	// Commented out temporarily to investigate destructor problems
	// free_results();
	free_results(); // GMnSKY2008MLE: uncommented to fix memory leaks 
	Cout<<"tFlowResults Object has been destroyed..."<<endl<<flush;
}

/***************************************************************************
**
**  tFlowResults::SetFlowResVariables()
**
** To initialize basic variables of tFlowResults and allocate 
** memory for streamflow arrays
**
***************************************************************************/
void tFlowResults::SetFlowResVariables(tInputFile &infile, double add_time)
{
	int i;
	int flag = 1;
	char tmp[kMaxNameSize];
	writeFlag = 0;
	count = 1;
	
	// Initialize char arrays
	for (i = 0; i < kMaxNameSize; i++)  {
		baseHydroName[i] = ' ';
		tmp[i] = ' ';
	}
	
	infile.ReadItem( baseHydroName,"OUTHYDROFILENAME"); 
	ribsOutput = infile.ReadItem( ribsOutput, "RIBSHYDOUTPUT");
	infile.ReadItem( Extension, "OUTHYDROEXTENSION");
	
	i=kMaxNameSize-1;
	while (flag) {
		if ( (baseHydroName[i] == '/') || (i == 0) )
			flag = 0;
		else
			tmp[i] = baseHydroName[i];
		i--;
	}
	strcpy( outlet, tmp );
	
	Cout<<"\nFlow Output Files:"<<endl<<endl;
	Cout<<"Hydrograph File Name: \t\t'"<<baseHydroName<<"'"<<endl;
	Cout<<"Hydrograph Outlet Name: \t\t "<<outlet<<endl;
	Cout<<"Hydrograph File Extension: \t"<<Extension<<endl;
	
	limit = (int)add_time; //TODO: this limit should really just be the length of runtime, rigth?
	
	if ((phydro = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: phydro failed..."<<endl;
	
	if ((mhydro = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: mhydro failed..."<<endl;
	
	if ((HsrfRout = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: HsrfRout failed..."<<endl;
	
	if ((SbsrfRout = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: SbsrfRout failed..."<<endl;
	
	if ((PsrfRout = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: PsrfRout failed..."<<endl;
	
	if ((SatsrfRout = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: SatsrfRout failed..."<<endl;
	
	if ((prr = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: prr failed..."<<endl;
	
	if ((crr = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: crr failed..."<<endl;
	
	if ((max = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: max failed..."<<endl;
    
	if ((min = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: min failed..."<<endl;
	
	if ((sat = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: sat failed..."<<endl;
	
	if ((msm = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: msm failed..."<<endl;
	
	if ((msmRt = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: msmRt failed..."<<endl;
	
	if ((msmU = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: msmU failed..."<<endl;
	
	if ((mgw = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: mgw failed..."<<endl;
	
	if ((met = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: met failed..."<<endl;

	// SKY2008Snow from AJR2007
	if ((swe = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: swe failed..."<<endl;
	if ((melt = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: melt failed..."<<endl;
	if ((snsub = (double*)calloc(limit,sizeof(double)))==NULL) // CJC2020
		cout <<"\ntFlowResults: snsub failed..."<<endl; // CJC2020
	if ((snevap = (double*)calloc(limit,sizeof(double)))==NULL) // CJC2020
		cout <<"\ntFlowResults: snevap failed..."<<endl; // CJC2020
	if ((stC = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: stC failed..."<<endl;
	if ((DUint = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: DUint failed..."<<endl;
	if ((slhf = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: slhf failed..."<<endl;
	if ((sshf = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: sshf failed..."<<endl;
	if ((sghf = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: sghf failed..."<<endl;
	if ((sphf = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: sphf failed..."<<endl;
	if ((srli = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: srli failed..."<<endl;
	if ((srlo = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: srlo failed..."<<endl;
	if ((srsi = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: srsi failed..."<<endl;
	if ((intsn = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: intsn failed..."<<endl;
	if ((intsub = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\ntFlowResults: intsub failed..."<<endl;
	if ((intunl = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\nFlowResults: intunl failed..."<<endl;
	if ((sca = (double*)calloc(limit,sizeof(double)))==NULL)
		cout <<"\nFlowResults: sca failed..."<<endl;

	if ((frac = (double*)calloc(limit,sizeof(double)))==NULL)
		cout<<"\ntFlowResults: frac failed..."<<endl;
	
	if ((fState = (int*)calloc(limit,sizeof(int)))==NULL)
		cout<<"\ntFlowResults: fState failed..."<<endl;
	
	if ((qunsat = (double*)calloc(limit,sizeof(double)))==NULL) // CJC2025
		cout<<"\ntFlowResults: qunsat failed..."<<endl;

    if ((Perc = (double*)calloc(limit,sizeof(double)))==NULL)
        cout<<"\ntFlowResults: Percolation failed..."<<endl; //ASM percolation option

    for (i=0; i<limit; i++) {
		mhydro[i] = phydro[i] = 0.0;
		prr[i] = crr[i] = 0.0;
		HsrfRout[i] = SbsrfRout[i] = 0.0;
		PsrfRout[i] = SatsrfRout[i] = 0.0;
		fState[i] = 0;
		max[i]=sat[i]=msm[i]=mgw[i]=msmU[i]=msmRt[i]=met[i]=frac[i]=0.0;
		      
		// SKY2008Snow from AJR2007
		swe[i] = melt[i] = intsn[i] = intsub[i] = intunl[i] = sca[i] = 0.0;
		stC[i] = DUint[i] = srsi[i] = srlo[i] = srli[i] = sghf[i] = sphf[i] = 0.0;
		sshf[i] = slhf[i] = snsub[i] = snevap[i] = 0.0; // Added snsub, snevap CJC2020
		qunsat[i] = 0.0; // CJC2025

                Perc[i]=0.0; //ASM Percolation opt

		min[i]=9999.99;
	}
	iimax=0;
	return;
}

/***************************************************************************
** 
**  tFlowResults: free_results()
**
**  Free memory for tFlowResults data structure. 
**  Points set to NULL
**
***************************************************************************/
void tFlowResults::free_results() 
{ 
	limit= 0;
	
	free(phydro);
	free(mhydro);
	free(HsrfRout);
	free(SbsrfRout);
	free(PsrfRout);
	free(SatsrfRout);
	free(prr);
	free(crr);
	free(max);
	free(min);
	free(sat);
	free(msm);
	free(msmU);
	free(msmRt);
	free(mgw);
	free(met);
	free(frac);
	free(fState);
	// SKY2008Snow from AJR2007
	free(swe);
	free(melt);
	free(snsub); // CJC2020
	free(snevap); // CJC2020
	free(stC);
	free(DUint);
	free(slhf);
	free(sshf);
	free(sghf);
	free(sphf);
	free(srli);
	free(srlo);
	free(srsi);
	free(intsn);
	free(intsub);
	free(intunl);
	free(sca);
	free(qunsat); // CJC2025

    //ASM 5/5/2016
    free(Perc); //ASM percolation option

	prr=nullptr;
	crr=NULL;
	phydro=NULL;
	mhydro=NULL;
	HsrfRout=NULL;
	SbsrfRout=NULL;
	PsrfRout=NULL;
	SatsrfRout=NULL;
	fState=NULL;
	max = NULL;
	min = NULL;
	msm = NULL;
	msmU = NULL;
	msmRt = NULL;
	mgw = NULL;
	met = NULL;
	sat = nullptr;
	frac = NULL;

	// SKY2008Snow from AJR2007
	swe = NULL;
	melt = NULL;
	snsub = NULL; // CJC2020
	snevap = NULL; // CJC2020
	stC = NULL;
	DUint = NULL;
	slhf = NULL;
	sshf = NULL;
	sghf = NULL;
	sphf = NULL;
	srli = NULL;
	srlo = NULL;
	srsi = NULL;
	intsn = NULL;
	intsub = NULL;
	intunl = NULL;
	sca = NULL;
    //ASM 5/5/2016
    Perc = NULL; //ASM percolation option
	qunsat = NULL; // CJC2025

	return;
}

//=========================================================================
//
//
//                  Section 2: tFlowResults: Write Hydrograph
//
//
//=========================================================================

/***************************************************************************
** 
**  tFlowResults: read_prev_hyd(char *filename, int index)
**
**  Read hydrographs
**
**  Algorithm:
**    open hydro file
**    read headings
**    until end of file:
**     if type tag is index
**        read streamflow corresponding to previous step
**     else
**        discard data
**    close file
**
***************************************************************************/
void tFlowResults::read_prev_hyd(char *filename, int index) 
{ 
	int j,ii,jl;
	FILE *ifile;
	char str[80];
	float flt1, flt2;
	
	if ((ifile=fopen(filename,"r")) == NULL) {  
		cout<<"\nError: Unable to open *.mrf file: "<<filename<<endl;
		cout<<"Exiting Program..."<<endl;
		exit(2);
	}
	
	fscanf(ifile,"%d",&jl);
	for (j=0;j<jl;j++) { 
		fscanf(ifile,"%s",str);
	}
	
	ii=0;
	while (!feof(ifile)) { 
		fscanf(ifile,"%d",&jl);
		if (jl==index) { 
			fscanf(ifile,"%s %f %f",str,&flt1,&flt2);
			phydro[ii]=(double)flt1;
			crr[ii]=prr[ii]=(double)flt2;
			ii++;
		}
		else fscanf(ifile,"%s %f %f",str,&flt1,&flt2);
	}
	
	iimax=ii;
	fclose(ifile);
	return;
}

/***************************************************************************
** 
**  tFlowResults: add_fore_hyd(filename, index)
**
**  Write hydrograph corresponding to forecasted rain
**
**  Algorithm:
**    open hydro file
**    get maximum range in array
**    for index in array range:
**      write streamflow corresponding to current step
**
***************************************************************************/
void tFlowResults::add_fore_hyd(char *filename, int index) 
{
	FILE *ifile;
	int ii;              	// Loop counter 
	int it_hour,it_min;  	// Hours and minutes to print results
	
	if ((ifile=fopen(filename,"a")) == NULL) {
		cout<<"\nError: Unable to open *.mrf file: "<<filename<<endl;
		cout<<"Exiting Program..."<<endl;
		exit(2);
	}
	for (ii=0; ii<iimax; ii++) { 
		timer->res_time_mid(ii, &it_hour, &it_min);
		fprintf(ifile," %d %d_%d %f %f\n",
				index+2, it_hour, it_min,
				phydro[ii]+mhydro[ii],
				crr[ii]);
	}
	fclose(ifile);     
	return;         
}

/***************************************************************************
** 
**  tFlowResults: writeAndUpdate(inter_hour, dt_rain)
**
**  Write basin state variables, output and result hydrographs
**
**  Algorithm:
**	compose time tag of current state
**	call basin model to compose path
**	call basin model to write state
**	call network model to compose path
**	call network model to write state and outputAlgorithm:
**
***************************************************************************/
void tFlowResults::writeAndUpdate( double time, int forenum )
{ 
	int  hour, minute;
	char timetag[20];
	writeFlag = 0;
	
	if (simCtrl->Verbose_label == 'Y')
		cout<<"\n\ttFlowResults: Time to write hydrograph; time = "
			<<time<<endl<<flush;
	
	hour   = (int)floor(time);
	minute = (int)floor((time-hour)*60);

	snprintf(timetag, sizeof(timetag),"%04d_%02d.", hour, minute);
	strcpy( currHydroName, baseHydroName );
	strcat( currHydroName, timetag );
	strcat( currHydroName, Extension);
	
	write_inter_hyd(currHydroName, outlet, forenum);
	
	update_prev_hyd();
	
	reset_meas_hyd();
}

/***************************************************************************
** 
**  tFlowResults: write_inter_hyd(char *filename, char *identification,
**                               int foreNum)
**
**  Write hydrographs corresponding to measured rain
** 
**  Algorithm:
**	open hydro file
**	get maximum range in array
**	if forecasted rain is active
**	   write headings for three variables
**	else
**	   write headings for two variables
**	for index in array range:
**	   write type tag 0
**         write streamflow corresponding to previous step
**      for index in array range:
**	   write type tag 1
**         write streamflowcorresponding to current step
**      copy hydro file name as last hydro file
**
***************************************************************************/
void tFlowResults::write_inter_hyd(char *filename, char *identification,
								   int foreNum) 
{ 
	FILE *ifile;
	int ii;              //Loop counter 
	int it_hour,it_min;  //Hours and minutes to print results 
	
	if ((ifile=fopen(filename,"w")) == nullptr) {
		cout<<"\nError: Unable to open *.mrf file: "<<filename<<endl;
		cout<<"Exiting Program..."<<endl;
		exit(2);
	}
	
	// Assign Forecast State if Option != 0
	if (timer->getoptForecast()!=0)
		fState[count] = checkForecast();
	
	// RIBS Compatibility
	
	if (ribsOutput == 1) {
		
		// Previous hydrograph
		for (ii=1; ii < iimax; ii++) {    
			timer->res_time_begin(ii, &it_hour, &it_min);
			fprintf(ifile,"0 %04d.%02d %f %f\n",
					it_hour, it_min, phydro[ii], prr[ii-1]);
		}
		// Current hydrograph
		for (ii=1; ii < iimax; ii++) {    
			timer->res_time_begin(ii, &it_hour, &it_min);
			fprintf(ifile," 1 %04d.%02d %f %f\n",
					it_hour, it_min, phydro[ii]+mhydro[ii], crr[ii-1]);
		}
	}
	
	else if (ribsOutput == 0) {
		
#ifdef PARALLEL_TRIBS

   // Variables for sums, mins, and maxs
   double *pPhydro, *pMhydro, *pCrr, *pMax, *pMin, *pMsm, *pMsmRt,
          *pMsmU, *pMgw, *pMet, *pSat, *pFrac,
          *pSwe, *pMelt, *pSnSub, *pSnEvap, *pStC, *pDUint, *pSlhf, *pSshf, *pSphf, // Added *pSnSub, *pSnEvap CJC2020
          *pSghf, *pSrli, *pSrlo, *pSrsi, *pIntsn, *pIntsub,
          *pIntunl, *pSca, *pPerc, *pQunsat; //ASM 5/5/2016

   // If running in parallel, collect sums, mins, maxs
   pPhydro = tParallel::sum(phydro, iimax);
   pMhydro = tParallel::sum(mhydro, iimax);
   pCrr = tParallel::sum(crr, iimax);
   pMax = tParallel::max(max, iimax);
   pMin = tParallel::min(min, iimax);
   pMsm = tParallel::sum(msm, iimax);
   pMsmRt = tParallel::sum(msmRt, iimax);
   pMsmU = tParallel::sum(msmU, iimax);
   pMgw = tParallel::sum(mgw, iimax);
   pMet = tParallel::sum(met, iimax);
   pSat = tParallel::sum(sat, iimax);
   pFrac = tParallel::sum(frac, iimax);
   pSwe = tParallel::sum(swe, iimax);
   pMelt = tParallel::sum(melt, iimax);
   pSnSub = tParallel::sum(snsub, iimax); // CJC2020
   pSnEvap = tParallel::sum(snevap, iimax); // CJC2020
   pStC = tParallel::sum(stC, iimax);
   pDUint = tParallel::sum(DUint, iimax);
   pSlhf = tParallel::sum(slhf, iimax);
   pSshf = tParallel::sum(sshf, iimax);
   pSphf = tParallel::sum(sphf, iimax);
   pSghf = tParallel::sum(sghf, iimax);
   pSrli = tParallel::sum(srli, iimax);
   pSrlo = tParallel::sum(srlo, iimax);
   pSrsi = tParallel::sum(srsi, iimax);
   pIntsn = tParallel::sum(intsn, iimax);
   pIntsub = tParallel::sum(intsub, iimax);
   pIntunl = tParallel::sum(intunl, iimax);
   pSca = tParallel::sum(sca, iimax);
   pPerc = tParallel::sum(Perc, iimax); //ASM percolation option
   pQunsat = tParallel::sum(qunsat, iimax);
   
   // Master processor writes file
   if (tParallel::isMaster()) {
#endif

		// Print out header information
		if (simCtrl->Header_label=='Y' && writeFlag == 0) {
			
			// fprintf(ifile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			// SKY2008Snow from AJR2007
			fprintf(ifile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",

					"Time", "Srf","MAP","RainMax", "RainMin","FState", "MSM100", "MSMRt", 
					// "MSMU", "MGW","MET", "%Sat", "%Rain");
					// SKY2008Snow from AJR2007
					"MSMU", "MDGW","MET", "SatPercent", "RainPercent",
					"AvSWE" , "AvMelt" , "AvSnSub" , "AvSnEvap" , "AvSTC" , "AvDUInt" , "AvSLHF" , "AvSSHF" , "AvSPHF" , "AvSGHF" , //added by AJR 2007 @ NMT // Added "AvSnSub" , "AvSnEvap" CJC2020
					"AvSRLI" , "AvSRLO" , "AvSRSI" , "AvInSn" , "AvInSu" , "AvInUn" , "SCA", "ChannelPercolation", "Qunsat");//added by AJR 2007 @ NMT
			fprintf(ifile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					"hr" , "m3/s" , "mm/hr" , "mm/hr" , "mm/hr" , "[]" , "[]" , "[]" , "[]" , "mm" , 
	    				"mm" , "[]" , "[]" , "cm" , "cm" , "cm" , "cm" , "C" , "kJ/m2" , "kJ/m2" , "kJ/m2" , //added by AJR 2007 @ NMT // added "cm" , "cm" CJC2020
	    				"kJ/m2" , "kJ/m2" , "kJ/m2" , "kJ/m2" , "kJ/m2" , "cm" , "cm" , "cm", "[]", "m3", "mm/hr" );//added by AJR 2007 @ NMT // removed extra header for kJ/m2 CJC2020 //akram: Need to confirm unit for Percolation

			// fprintf(ifile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			//		"hr","m3/s","mm/hr","mm/hr","mm/hr","[]", "[]", "[]", "[]",
			//		"mm","mm", "[]", "[]");
			writeFlag = 1;
		}
		
		for (ii=0; ii < iimax; ii++) {    
			timer->res_time_begin(ii+1, &it_hour, &it_min);
			//fprintf(ifile,"%d.%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",

#ifdef PARALLEL_TRIBS

         // Print min, max, and summ variables from all processors
         fprintf(ifile,"%d.%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
               it_hour, it_min, pPhydro[ii]+pMhydro[ii], pCrr[ii],
               pMax[ii], pMin[ii], fState[ii], pMsm[ii], pMsmRt[ii],pMsmU[ii], 
               pMgw[ii], pMet[ii], pSat[ii], pFrac[ii],
               pSwe[ii], pMelt[ii], pSnSub[ii], pSnEvap[ii], pStC[ii], pDUint[ii], pSlhf[ii], pSshf[ii], // Added pSnSub[ii], pSnEvap[ii] CJC2020
               pSphf[ii], pSghf[ii], pSrli[ii], pSrlo[ii], pSrsi[ii], pIntsn[ii],
               pIntsub[ii], pIntunl[ii], pSca[ii], pPerc[ii], pQunsat[ii]);

#else
			// SKY2008Snow from AJR2007
			fprintf(ifile,"%d.%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
					it_hour, it_min, phydro[ii]+mhydro[ii], crr[ii], 
					max[ii], min[ii], fState[ii], msm[ii], msmRt[ii],msmU[ii], mgw[ii], met[ii], sat[ii], frac[ii],
					swe[ii], melt[ii], snsub[ii], snevap[ii], stC[ii], DUint[ii], slhf[ii], sshf[ii], sphf[ii], sghf[ii],//added by AJR 2007 @ NMT // Added snsub[ii], snevap[ii] CJC2020
					srli[ii], srlo[ii], srsi[ii], intsn[ii], intsub[ii], intunl[ii], sca[ii], Perc[ii], qunsat[ii]);//added by AJR 2007 @ NMT


					//it_hour, it_min, phydro[ii]+mhydro[ii], crr[ii], 
					//max[ii], min[ii], fState[ii], msm[ii], msmRt[ii],
					//msmU[ii], mgw[ii], met[ii], sat[ii], frac[ii]);
#endif
		}

#ifdef PARALLEL_TRIBS
   }

     delete [] pPhydro;
   delete [] pMhydro;
   delete [] pCrr;
   delete [] pMax;
   delete [] pMin;
   delete [] pMsm;
   delete [] pMsmRt;
   delete [] pMsmU;
   delete [] pMgw;
   delete [] pMet;
   delete [] pSat;
   delete [] pFrac;
   delete [] pSwe;
   delete [] pMelt;
   delete [] pSnSub; // CJC2020
   delete [] pSnEvap; // CJC2020
   delete [] pStC;
   delete [] pDUint;
   delete [] pSlhf;
   delete [] pSshf;
   delete [] pSphf;
   delete [] pSghf;
   delete [] pSrli;
   delete [] pSrlo;
   delete [] pSrsi;
   delete [] pIntsn;
   delete [] pIntsub;
   delete [] pIntunl;
   delete [] pSca;
   delete [] pPerc; //ASM percolation option
   delete [] pQunsat; // CJC2025

#endif

	}
	Cout<<"\nCreating Hydrograph Output: '"<<filename<<"'"<<endl;

#ifdef PARALLEL_TRIBS
  // If running parallel, Master closes file
  if (tParallel::isMaster())
#endif
	fclose(ifile);                 
	
	// Copy current rain hyetograph to previous
	for (ii=0; ii < limit; ii++)
		prr[ii] = crr[ii];
	
	count++;
	return;
}

/***************************************************************************
** 
**  tFlowResults: write_extra_hyd(name, identification)
**
**  Write hydrographs for the extra gauge. The extra gauge is used to store
**  virtual variables for the user interface. Values stored are not neccessarily
**  hydrographs
**
**  Algorithm:
**     open hydro file
**     get maximum range in array
**     write headings for one variable
**     for index in array range:
**       write streamflow corresponding to current step
**
***************************************************************************/
void tFlowResults::write_extra_hyd(char *name, char *identification)
{
	FILE *ifile;
	int ii;               // Loop counter 
	int it_hour, it_min;  // Hours and minutes to print results 
	
	if ((ifile=fopen(name,"w")) == NULL) {
		cout<<"\nError: Unable to open extra hydro file: "<<name<<endl;
		cout<<"Exiting Program..."<<endl;
		exit(2);
	}
	
	fprintf(ifile,"  1  \n%s\n", identification);
	for (ii=0; ii < iimax; ii++) {
		timer->res_time_mid(ii, &it_hour, &it_min);
		fprintf(ifile," 0 %04d.%02d %f %f\n",
				it_hour, it_min, phydro[ii]+mhydro[ii], crr[ii]);
	}
	fclose(ifile);   
	return;
}

/***************************************************************************
** 
**  tFlowResults: write_Runoff_Types(char *filename, char *identification)
**
**  Write hydrographs for each runoff type to single file.
**  Option to write to separete files commented out.
**
***************************************************************************/
void tFlowResults::write_Runoff_Types(char *filename, char *)
{
	FILE *ifile;
	int ii;               
	int it_hour, it_min;  
	
#ifdef PARALLEL_TRIBS

  // Sum hydrographs for each runoff type across processors
  double* hr = tParallel::sum(HsrfRout, iimax);
  double* sbr = tParallel::sum(SbsrfRout, iimax);
  double* pr = tParallel::sum(PsrfRout, iimax);
  double* satr = tParallel::sum(SatsrfRout, iimax);


  if (tParallel::isMaster()) {
    for (int i = 0; i < iimax; i++) {
      HsrfRout[i] = hr[i];
      SbsrfRout[i] = sbr[i];
      PsrfRout[i] = pr[i];
      SatsrfRout[i] = satr[i];
    }

    delete [] hr;
    delete [] sbr;
    delete [] pr;
    delete [] satr;
  }

  // If running in parallel, the Master processor writes this file
  if (tParallel::isMaster()) {
#endif

	if ((ifile=fopen(filename,"w")) == NULL) {
		cout<<"\nError: Unable to open *.rft file: "<<filename<<endl;
		cout<<"Exiting Program..."<<endl;
		exit(2);
	}
	
	if (simCtrl->Header_label=='Y') {
		fprintf(ifile,"%s\t","Time");
		fprintf(ifile,"%s\t","Hsrf");
		fprintf(ifile,"%s\t","Sbsrf");
		fprintf(ifile,"%s\t","Psrf");
		fprintf(ifile,"%s\n","Satsrf");
		fprintf(ifile,"%s\t%s\t%s\t%s\t%s\n","hr","m3/s","m3/s","m3/s","m3/s");
	}
	
	// Current Hydrographs
	for (ii=0; ii<iimax; ii++)  {      // Current
		timer->res_time_begin(ii+1, &it_hour, &it_min);
		fprintf(ifile,"%04d.%02d ", it_hour, it_min);
		fprintf(ifile,"\t%f\t", HsrfRout[ii]);
		fprintf(ifile,"%f\t", SbsrfRout[ii]);
		fprintf(ifile,"%f\t", PsrfRout[ii]);
		fprintf(ifile,"%f\n", SatsrfRout[ii]);
	}   
	Cout<<"Creating Runoff Type Output: '"<<filename<<"'"<<endl;
	fclose(ifile);   

#ifdef PARALLEL_TRIBS
  }
#endif

	return;
}

/***************************************************************************
** 
**  tFlowResults::whenTimeIsOver( double time )
**
**  Call routine for write_Runoff_Types
**
***************************************************************************/
void tFlowResults::whenTimeIsOver( double time )
{
	int  hour, minute;
	char fullName[kMaxNameSize+kMaxExt];
	char timetag[20];
	
	if (simCtrl->Verbose_label == 'Y')
		cout<<"\n\ttFlowResults: Time to write runoff types; time = "
			<<time<<endl<<flush;
	
	hour   = (int)floor(time);
	minute = (int)floor((time-hour)*60);
	
    snprintf(timetag,sizeof(timetag),"%04d_%02d", hour, minute);
	strcpy( fullName, baseHydroName );
	strcat( fullName, timetag );
	strcat( fullName, ".rft" );
	
	write_Runoff_Types(fullName, outlet);
	return;
}

//=========================================================================
//
//
//                  Section 3: tFlowResults: Store Volume
//
//
//=========================================================================

/***************************************************************************
** 
**  tFlowResults: store_volume(double time, double value)
**
**  Store value of discharge (expressed in volume per cell) corresponding to
**  a lag equal to "time" from the beginning of the time step. The value is 
**  distributed equally during the simulation time step and stored 
**  proportionally in every results time that over laps with simulation time
**
**  Algorithm:
**      get first results step that overlaps with the simulation step
**	get last results step that overlaps with the simulation step
**   	if first=last
**     	   store value multiplied by dt_calc(simulation)/dt_res(results)
**  	else if difference last-first = 1
**         get fraction of first: end of first results step -
**                                begin of simulation time step
**         store fraction in first
**         get fraction of last:  simulation time step -
**                                  fraction of first
**         store fraction in last
**      else
**         get fraction of first: end of first results step -
**                                begin of simulation time step
**         store fraction in first
**         for all intermediate results steps
**           store fraction
**         get fraction of last:  end of simulation step -
**                                begin of last results step
**         store fraction in last
**
***************************************************************************/
void tFlowResults::store_volume(double time, double value) 
{       
	int init, end, ii;
	double dCalc, dRes, dInt;
	double Sum = 0;
	
	dCalc = timer->getTimeStep();
	dRes  = timer->getOutputInterval();
	
	init = timer->getResStep(time+.01*dCalc); 
	end  = timer->getResStep(time+.99*dCalc);
	
	if (init==end) {
		add_m_volume(value, init); 
	}
	else {
		if (end-init == 1) {  
			dInt = timer->res_hour_end(init) - timer->get_abs_hour(time);
			add_m_volume(value*dInt/dCalc, init); 
			add_m_volume(value*(dCalc-dInt)/dCalc, end); 
			if (dInt/dCalc > 1.) 
				cout<<"\nWarning: Improper time relationship in tFlowResults.cpp"<<endl;
		}
		else { // difference between init & end > 1
			dInt = timer->res_hour_end(init) - timer->get_abs_hour(time);
			if (dInt/dCalc > 1.) 
				cout<<"\nWarning: Improper time relationship in tFlowResults.cpp"<<endl;
			add_m_volume(value*dInt/dCalc, init); 
			Sum = dInt/dCalc;
			for (ii=init+1; ii < end; ii++) {
				add_m_volume(value*dRes/dCalc, ii);
				Sum += dRes/dCalc;
			}
			dInt=timer->get_abs_hour(time+dCalc) - timer->res_hour_begin(end); 
			add_m_volume(value*dInt/dCalc, end);
			Sum += dInt/dCalc;
			if (Sum > 1.05) {
				cout<<"\nWarning: Improper time relationship in tFlowResults.cpp"<<endl;
				cout<<"Sum = "<<Sum<<"\n\n";
			}
		}
	}
	return;
}

/***************************************************************************
** 
**  tFlowResults: store_volume_Type(double time, double value, int Type)
**
**  Store value of discharge (expressed in volume per cell) corresponding to
**  a lag equal to "time" from the beginning of the time step for each
**  runoff type
**
***************************************************************************/
void tFlowResults::store_volume_Type(double time, double value, int Type)
{
	int init, end, ii;
	double dCalc, dRes, dInt;
	double Sum = 0;
	
	dCalc = timer->getTimeStep();
	dRes  = timer->getOutputInterval();
	
	// To prevent synchronous stuff 
	init = timer->getResStep(time+.01*dCalc); 
	end  = timer->getResStep(time+.99*dCalc);
	
	if (init == end) {
		add_m_volume_Type(value, init, Type); 
	}
	else {
		if ((end-init) == 1) { 
			dInt = timer->res_hour_end(init) - timer->get_abs_hour(time);
			add_m_volume_Type(value*dInt/dCalc, init, Type); 
			add_m_volume_Type(value*(dCalc-dInt)/dCalc, end, Type); 
			if (dInt/dCalc > 1.) 
				cout<<"\nWarning: Improper time relationship in tFlowResults.cpp"<<endl;
		}
		else { // difference between init & end > 1 
			dInt = timer->res_hour_end(init) - timer->get_abs_hour(time);
			if (dInt/dCalc > 1.) 
				cout<<"\nWarning: Improper time relationship in tFlowResults.cpp"<<endl;
			add_m_volume_Type(value*dInt/dCalc, init, Type); 
			Sum = dInt/dCalc;
			
			for (ii=init+1; ii<end; ii++)  {
				add_m_volume_Type(value*dRes/dCalc, ii, Type);
				Sum += dRes/dCalc;
			}
			dInt=timer->get_abs_hour(time+dCalc) - timer->res_hour_begin(end); 
			add_m_volume_Type(value*dInt/dCalc, end, Type);
			Sum += dInt/dCalc;
			if (Sum > 1.05) {
				cout<<"\nWarning: Improper time relationship in tFlowResults.cpp"<<endl;
				cout<<"Sum = "<<Sum<<"\n\n";
			}
		}
	}
	return;
}

/***************************************************************************
** 
**  tFlowResults: add_m_volume_Type(double Value, int iStep, int Type)
**
**  Add value of discharge to volume for each runoff type
**
***************************************************************************/
void tFlowResults::add_m_volume_Type(double value, int iStep, int Type)
{
	// Hortonian runoff
	if (Type == 1)      
		HsrfRout[iStep] += value/timer->getOutputIntervalSec();
	
	// Saturation from below
	else if (Type == 2) 
		SbsrfRout[iStep] += value/timer->getOutputIntervalSec();
	
	// Perched saturation runoff
	else if (Type == 3) 
		PsrfRout[iStep] += value/timer->getOutputIntervalSec();
	
	// Ground water return flow
	else if (Type == 4) 
		SatsrfRout[iStep] += value/timer->getOutputIntervalSec();
	return;
}

//=========================================================================
//
//
//                  Section 4: tFlowResults: Rainfall Function
//
//
//=========================================================================

/***************************************************************************
** 
**  tFlowResults: store_rain(double time, double value)
**
**  Store value of rainfall corresponding to a lag equal to "time"
**  from the beginning of the time step. The value is distributed equally
**  during the simulation time step and stored proportionally in every
**  results time step that over laps with the simulation time step
**
** Algorithm:
**   get first results step that overlaps with the simulation step
**   get last results step that overlaps with the simulation step
**   if first=last
**     store value multiplied by dt_calc(simulation)/dt_res(results)
**   else if difference last-first = 1
**           get fraction of first: end of first results step -
**                                  begin of simulation time step
**           store fraction in first
**           get fraction of last:  simulation time step -
**                                    fraction of first
**           store fraction in last
**        else
**           get fraction of first: end of first results step -
**                                  begin of simulation time step
**           store fraction in first
**           for all intermediate results steps
**             store fraction
**           get fraction of last:  end of simulation step -
**                                  begin of last results step
**           store fraction in last
**
***************************************************************************/
void tFlowResults::store_rain(double time, double value) 
{    
	int init,end,ii;
	double dcalc,dres,dint;
	
	dcalc = timer->getTimeStep();
	dres  = timer->getOutputInterval();
	
	init = timer->getResStep(time-.01*dcalc); 
	end  = timer->getResStep(time-.99*dcalc);
	
	if (init==end) {
		crr[init] += value*dcalc/dres;
	}
	else {
		if (end-init == 1) { 
			dint = timer->res_hour_end(init) - timer->get_abs_hour(time);
			crr[init] += value*dint/dres; 
			crr[end]  += value*(dcalc-dint)/dres; 
		}
		else { 
			dint = timer->res_hour_end(init) - timer->get_abs_hour(time);
			crr[init] += value*dint/dres; 
			for (ii=init+1; ii < end; ii++) {
				crr[ii] +=value;
			}
			dint=timer->get_abs_hour(time+dcalc) - timer->res_hour_begin(end); 
			crr[end] += value*dint/dres;
		}
	}
	return;
}

/***************************************************************************
** 
**  tFlowResults: store_maxminrain(double time, double value)
**
***************************************************************************/
void tFlowResults::store_maxminrain(double time, double value, int flag) 
{   
	int init,end;
	double dcalc,dres;
	
	dcalc = timer->getTimeStep();
	dres  = timer->getOutputInterval();
	init = timer->getResStep(time-.01*dcalc); 
	end  = timer->getResStep(time-.99*dcalc);
	
	if (init==end) {
		if (flag == 0) {
			if (value > max[init])
				max[init] = value;
			if (value <= min[init])
				min[init] = value;    
		}
		else if (flag == 1) {
			frac[init] += value*dcalc/dres;
		}
	}
	return;
}

/***************************************************************************
** 
**  tFlowResults: store_saturation(double time, double value)
** TODO is it possible to optimize this part of the code?
***************************************************************************/
void tFlowResults::store_saturation(double time, double value, int flag) 
{   
	int init,end;
	double dcalc, dres;
	
	dcalc = timer->getTimeStep();
	dres  = timer->getOutputInterval();
	init  = timer->getResStep(time-.01*dcalc); 
	end   = timer->getResStep(time-.99*dcalc);
	
	if (init==end) {
		if (flag == 0)
			msm[init] += value*dcalc/dres;
		else if (flag == 1)
			msmRt[init] += value*dcalc/dres;
		else if (flag == 2)
			msmU[init] += value*dcalc/dres;
		else if (flag == 3)
			sat[init] += value*dcalc/dres;
		else if (flag == 4)
			mgw[init] += value*dcalc/dres;
		else if (flag == 5)
			met[init] += value*dcalc/dres;
		
		// SKY2008Snow from AJR2007
		else if (flag == 6)
			swe[init] += value*dcalc/dres;
		else if (flag == 7)
			melt[init] += value*dcalc/dres;
		else if (flag == 8)
			stC[init] += value*dcalc/dres;
		else if (flag == 9)
			DUint[init] += value*dcalc/dres;
		else if (flag == 10)
			slhf[init] += value*dcalc/dres;
		else if (flag == 11)
			sshf[init] += value*dcalc/dres;
		else if (flag == 12)
			sghf[init] += value*dcalc/dres;
		else if (flag == 13)
			sphf[init] += value*dcalc/dres;
		else if (flag == 14)
			srli[init] += value*dcalc/dres;
		else if (flag == 15)
			srlo[init] += value*dcalc/dres;
		else if (flag == 16)
			srsi[init] += value*dcalc/dres;
		else if (flag == 17)
			intsn[init] += value*dcalc/dres;
		else if (flag == 18)
			intsub[init] += value*dcalc/dres;
		else if (flag == 19)
			intunl[init] += value*dcalc/dres;
		else if (flag == 20)
			sca[init] += value*dcalc/dres;
		else if (flag == 21)
			snsub[init] += value*dcalc/dres; // CJC2020
		else if (flag == 22)
			snevap[init] += value*dcalc/dres; // CJC2020
		else if (flag == 24)
			qunsat[init] += value*dcalc/dres; // CJC2025

            //ASM percolation option
        else if (flag ==23)
            Perc[init] += value*225; //ASM the 225 converts to m3


    }
	return;
}

/*****************************************************************************
**  
**  tFlowResults::checkForecast()
**  
**  Check the forecast state. Returns integer representing state:
**  
**  0 = Before and up to forecast time, Use QPE
**  1 = In Forecast Period and up to lead time, Use QPF
**  2 = In Forecast Period and after lead time, Use Average Rainfall
**  3 = After Forecast Period, End simulation
**
*****************************************************************************/
int tFlowResults::checkForecast() 
{
	int state {};
	
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
	
	return state;
}

/***************************************************************************
**
** tFlowResults::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/

void tFlowResults::writeRestart(fstream & rStr) const
{
  BinaryWrite(rStr, limit);
  BinaryWrite(rStr, iimax);
  BinaryWrite(rStr, ribsOutput);
  BinaryWrite(rStr, writeFlag);
  BinaryWrite(rStr, count);
  for (int i = 0; i < limit; i++)
    BinaryWrite(rStr, fState[i]);

  for (int i = 0; i < limit; i++) {
    BinaryWrite(rStr, prr[i]);
    BinaryWrite(rStr, crr[i]);
    BinaryWrite(rStr, phydro[i]);
    BinaryWrite(rStr, mhydro[i]);
    BinaryWrite(rStr, HsrfRout[i]);
    BinaryWrite(rStr, SbsrfRout[i]);
    BinaryWrite(rStr, PsrfRout[i]);
    BinaryWrite(rStr, SatsrfRout[i]);
    BinaryWrite(rStr, max[i]);
    BinaryWrite(rStr, min[i]);
    BinaryWrite(rStr, msm[i]);
    BinaryWrite(rStr, msmRt[i]);
    BinaryWrite(rStr, msmU[i]);
    BinaryWrite(rStr, mgw[i]);
    BinaryWrite(rStr, met[i]);
    BinaryWrite(rStr, sat[i]);
    BinaryWrite(rStr, frac[i]);
    BinaryWrite(rStr, swe[i]);
    BinaryWrite(rStr, melt[i]);
    BinaryWrite(rStr, snsub[i]); // CJC2020
    BinaryWrite(rStr, snevap[i]); // CJC2020
    BinaryWrite(rStr, stC[i]);
    BinaryWrite(rStr, DUint[i]);
    BinaryWrite(rStr, slhf[i]);
    BinaryWrite(rStr, sshf[i]);
    BinaryWrite(rStr, sghf[i]);
    BinaryWrite(rStr, sphf[i]);
    BinaryWrite(rStr, srli[i]);
    BinaryWrite(rStr, srlo[i]);
    BinaryWrite(rStr, srsi[i]);
    BinaryWrite(rStr, intsn[i]);
    BinaryWrite(rStr, intsub[i]);
    BinaryWrite(rStr, intunl[i]);
    BinaryWrite(rStr, sca[i]);
    BinaryWrite(rStr, Perc[i]); //ASM Percolation option
    BinaryWrite(rStr, qunsat[i]); // CJC2025
  }
}

/***************************************************************************
**
** tFlowResults::readRestart() Function
**
***************************************************************************/

void tFlowResults::readRestart(fstream & rStr)
{
  BinaryRead(rStr, limit);
  BinaryRead(rStr, iimax);
  BinaryRead(rStr, ribsOutput);
  BinaryRead(rStr, writeFlag);
  BinaryRead(rStr, count);
  for (int i = 0; i < limit; i++)
    BinaryRead(rStr, fState[i]);

  for (int i = 0; i < limit; i++) {
    BinaryRead(rStr, prr[i]);
    BinaryRead(rStr, crr[i]);
    BinaryRead(rStr, phydro[i]);
    BinaryRead(rStr, mhydro[i]);
    BinaryRead(rStr, HsrfRout[i]);
    BinaryRead(rStr, SbsrfRout[i]);
    BinaryRead(rStr, PsrfRout[i]);
    BinaryRead(rStr, SatsrfRout[i]);
    BinaryRead(rStr, max[i]);
    BinaryRead(rStr, min[i]);
    BinaryRead(rStr, msm[i]);
    BinaryRead(rStr, msmRt[i]);
    BinaryRead(rStr, msmU[i]);
    BinaryRead(rStr, mgw[i]);
    BinaryRead(rStr, met[i]);
    BinaryRead(rStr, sat[i]);
    BinaryRead(rStr, frac[i]);
    BinaryRead(rStr, swe[i]);
    BinaryRead(rStr, melt[i]);
    BinaryRead(rStr, snsub[i]); // CJC2020
    BinaryRead(rStr, snevap[i]); // CJC2020
    BinaryRead(rStr, stC[i]);
    BinaryRead(rStr, DUint[i]);
    BinaryRead(rStr, slhf[i]);
    BinaryRead(rStr, sshf[i]);
    BinaryRead(rStr, sghf[i]);
    BinaryRead(rStr, sphf[i]);
    BinaryRead(rStr, srli[i]);
    BinaryRead(rStr, srlo[i]);
    BinaryRead(rStr, srsi[i]);
    BinaryRead(rStr, intsn[i]);
    BinaryRead(rStr, intsub[i]);
    BinaryRead(rStr, intunl[i]);
    BinaryRead(rStr, sca[i]);
    BinaryRead(rStr, Perc[i]); //ASM Percolationoption
    BinaryRead(rStr, qunsat[i]); // CJC2025
  }
}

//=========================================================================
//
//
//                         End of tFlowResults.h
//
//
//=========================================================================
