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
**  tStorm.cpp: Functions for class tStorm.
**
**  A tStorm object generates random storms assuming an exponential
**  distribution of rainfall intensity, storm duration, and time to the
**  next storm. 
**
**************************************************************************/

#include "src/tStorm/tStorm.h"
#include "src/Headers/globalFns.h"
#include "src/Headers/globalIO.h"
#include <cstring>

/**************************************************************************
**
**  Default tStorm():  Constructor for storms. The default constructor
**             assigns a value of unity to storm depth, duration,
**             and interstorm duration. (Note: this does not allow sinusoidal 
**             variation in means).
**
**************************************************************************/
tStorm::tStorm()
:
optStoch(0),
stdurMean(1.0),
istdurMean(1.0),
pMean(1.0),
p(1.0),
stdur(1.0),
istdur(1.0),
seed(0),      // default seed
endtm(1.0E9)
{
	
}

/**************************************************************************
**
**  tStorm() Constructor
**
**  Alternative constructor that reads parameters directly from tInputFile
**  Reads option for variable storms,  mean values for rainfall intensity, 
**  duration, and interstorm period, and a random seed.
**
**  Also reads an option for long-term sinusoidal variations in the mean
**  values, and if the option is selected, reads the relevant parameters.
**  Variables p0, stdur0, and istdur0 are the mean values of the means;
**  pdev, stdurdev, and istdurdev are the range of variation.
**
**  Modifications:
**   - 3/00 initialization now includes creation of ".storm" file for storm
**     history (GT)
**
**************************************************************************/
tStorm::tStorm( SimulationControl *simCtrPtr, tInputFile &infile )
{  
	simCtrl = simCtrPtr;
	
	//Default values
	optStoch = 0;
	stdurMean = 1.0;
	istdurMean = 1.0;
	pMean = 1.0;
	seed = 0;      
	endtm = 1.0E6;
	
	SetStormVariables( infile );
}

tStorm::~tStorm()
{  
	Cout<<"tStorm Object has been destroyed..."<<endl<<flush;
}

/**************************************************************************
**
**  tStorm::SetStormVariables()
**
**************************************************************************/
void tStorm::SetStormVariables(tInputFile &infile)
{ 
	// Read + set parameters for storm intensity, duration, and spacing
	optStoch = infile.ReadItem( optStoch, "STOCHASTICMODE" );
	
	if (optStoch) {
		pMean = infile.ReadItem( pMean, "PMEAN" );
		stdurMean = infile.ReadItem( stdurMean, "STDUR" );
		istdurMean = infile.ReadItem( istdurMean, "ISTDUR" );
		endtm = infile.ReadItem( endtm, "RUNTIME" );
		
		// Set the current values
		p = pMean;
		stdur = stdurMean;
		istdur = istdurMean;
		
		// Handle option for sinuidoil variation in means
		if ( optStoch == 3 ||  optStoch == 4 || optStoch == 5) {
			p0 = pMean;
			stdur0 = stdurMean;
			istdur0 = istdurMean;
			
			twoPiLam = (2.0*PI)/(infile.ReadItem( twoPiLam, "PERIOD" ));
			pdev = infile.ReadItem( pdev, "MAXPMEAN" ) - pMean;
			if ( pdev < 0.0 ) 
				cout << "\nWarning: MAXPMEAN < PMEAN !"<<endl;
			
			else if ( pdev > pMean ) 
				cout << "\nWarning: MINPMEAN < 0 !"<<endl;
			
			stdurdev = infile.ReadItem( stdurdev, "MAXSTDURMN" ) - stdurMean;
			if ( stdurdev < 0.0 ) 
				cout << "\nWarning: MAXSTDURMN < STDURMN !"<<endl;
			
			else if ( stdurdev > stdurMean ) 
				cout << "\nWarning: MINSTDURMN < 0 !"<<endl;
			
			istdurdev = infile.ReadItem( istdurdev, "MAXISTDURMN" ) - istdurMean;
			if ( istdurdev < 0.0 ) 
				cout << "\nWarning: MAXISTDURMN < ISTDURMN !"<<endl;
			else if ( istdurdev > istdurMean ) 
				cout << "\nWarning: MINISTDURMN < 0 !"<<endl;
		}
		
		// Option for Weather Generator
		else if (optStoch == 6) {}
		
		// Create a file for writing t
		if ( optStoch != 1 ) {
			char fname[200];
			infile.ReadItem( fname, "OUTHYDROFILENAME" );
			strcat( fname, ".storm" );
			stormfile.open( fname );
			if ( !stormfile.good() )
				cout<<"\nWarning: unable to create storm data file '" 
					<<fname<<"'\n";
			
			// Read and initialize seed for random number generation
			seed = infile.ReadItem( seed, "SEED" );
		}
		
		// ############################################################
		// ### ONLY when REAL rain drives all other stochastic vars ###
		// ############################################################
		
		RealRn = 0;
		
		if (RealRn) {
			
			int i;
			int cntRn = 0;
			int inrain = 0;
			double tempo;
			ifstream InpRn;
			int endTime = infile.ReadItem( endTime, "RUNTIME" );
			
			rainIN = new double [endTime];
			assert(rainIN != 0);
			InpRn.open("Input/albuq_hrs");
			if (!InpRn) {
				cout <<"File 'Input/albuq_hrs' not found!\nExiting..."<<endl;
				exit(2);
			}
			// Reading meteorological data
			for (i=0; i < endTime; i++) {
				//InpRn>>tempo>>rainIN[i]>>tempo>>tempo>>tempo>>tempo;
				InpRn>>rainIN[i];
				if (rainIN[i] > 0.0)  
					cntRn++;
			}
			stbeg = new int [cntRn];
			assert(stbeg != 0);
			stend = new int [cntRn];
			assert(stend != 0);
			for (i=0; i < cntRn; i++)
				stbeg[i] = stend[i] = 0;
			
			rid = 0;
			for (i=0; i < endTime; i++) {
				// Input file must start from rainfall! I.e., at i=0
				if (rainIN[i] > 0.0) {
					// Consider every rain hour separately (which allows for)
					// temporal variability WITHIN the storm in the input data 
					stbeg[rid] = i;
					stend[rid] = i+1;
					rid++;
				}
			}
			rid = 0;  // index of the very first rain
		}
	}
	// No Stochastic Rainfall used
	else {;}
	return;
}


/**************************************************************************
**
**  tStorm::GenerateStorm()
**
**  Generates a new storm by drawing new values of p, stdur, and istdur from
**  an exponential distribution and updating the random seed.
**    If the minp parameter is greater than zero, the function will keep
**  picking storms until it finds one with p>minp. The total elapsed time,
**  including the rejected storms and their associated interstorm periods,
**  is stored istdur.
**
**  Inputs:      minp -- minimum value of rainfall rate p (default 0)
**               mind -- minimum storm depth to produce runoff (default 0)
**               tm -- current time in simulation
**  Members updated:  p, stdur, istdur take on new random values (if optVar)
**                    pMean, stdurMean, istdurMean adjusted (if optSinVar)
**  Assumptions:  pMean > 0
**
**************************************************************************/
void tStorm::GenerateStorm( double tm, double minp, double mind )
{
	float gamvar;
	
	if ( optStoch ) {
		// If option for sinusoidal variation is on, adjust the means of pdfs
		// Also set values for current storm to the means, in case option for
		// random storms is off.
		if ( optStoch == 3 ||  optStoch == 4 || optStoch == 5) {
			double sinfn = sin( tm*twoPiLam );
			pMean = p0 + pdev*sinfn;
			stdurMean = stdur0 + stdurdev*sinfn;
			istdurMean = istdur0 + istdurdev*sinfn;
		}
		if ( optStoch != 1 ) {
			// If option for random storms is on, pick a storm at random. 
			stdur = 0.0;
			istdur = 0.0;
			do {
				
				// --- Simulate interstorm duration ---
				istdur += istdurMean*ExpDev( &seed ) + stdur;
				//istdur += random_expon( 1/istdurMean ) + stdur;
				
				// --- Simulate storm duration ---
				stdur = stdurMean*ExpDev( &seed );
				//stdur = expon( 1/stdurMean );
				
				do {
					//p = pMean*ExpDev( &seed );
					//p = random_expon( 1/pMean );
					
					// --- Simulate storm depth as dependent variable of 'stdur' ---
					gamvar = gengam(1.0/(stdurMean*pMean), stdur/stdurMean);
					// --- Get storm rate ---
					p = (double)gamvar/stdur;
					
				} while (p<=0.0);
				
			} while( (p<=minp || (p*stdur)<=mind) && (tm+istdur+stdur<endtm) );
			
			// ############################################################
			// ### ONLY when REAL rain drives all other stochastic vars ###
			// ############################################################
			if (RealRn) {
				p = rainIN[ stbeg[rid] ]; // The very first value of rainfall pulse
				stdur = stend[rid] - stbeg[rid];
				istdur = stbeg[rid+1] - stend[rid];
				if (istdur < 0)  // In the end of stbeg array
					istdur = 1000.0;
				rid++;
			}
			
			stormfile<<currSeasID<<"  \t"<<p<<"  \t"<<stdur<<"  \t"<< istdur<<endl;
		}
	}
	//No Stochastic Rainfall used
	else {;}
	
	return;
}

/**************************************************************************
**
**  tStorm::ExpDev:  Finds a random number with an exponential distribution
**                   (adapted from Numerical Recipes).
**
**************************************************************************/
double tStorm::ExpDev( long *idum ) 
{
	double dum;
	do
		dum = ran3( idum );
	while ( dum==0.0 );
	return -log(dum);
}

/**************************************************************************
**
**  tStorm::IsDifferentSeason: Checks if the seson has changed
**
**************************************************************************/
int tStorm::DefineSeason(int mon) 
{ 
	int seas = 0;
	for (int i=0; i < NumSeas; i++) {
		// Special case: e.g., November through January
		if (SeasMo[i][0] > SeasMo[i][1]) {
			if (mon >= SeasMo[i][0]) {
				if (mon <= 12) 
					seas = i+1;
			}
			else {
				if (mon <= SeasMo[i][1]) 
					seas = i+1;
			}
		}
		else {
			if (mon >= SeasMo[i][0] && mon <= SeasMo[i][1]) 
				seas = i+1;
		}
	}
	return seas;
}

/**************************************************************************
**
**  tStorm "set" routines: set various variables
**
**************************************************************************/

void tStorm::setRainrate(double value) { p = value; }

void tStorm::setSeasonID(int value) { currSeasID = value; }

void tStorm::updateSeasonVars() 
{ 
	pMean      = MRate[currSeasID-1];
	stdurMean  = MStmDr[currSeasID-1];
	istdurMean = MSplDr[currSeasID-1];
	//cout<<"\t\t%%% SEASON CHANGED: "<<currSeasID<<" %%%"<<endl;
	//cout<<"\tMRate: "<<pMean<<"  \tMStDur: "<<stdurMean<<"\tMSpDur: "<<istdurMean<<endl;
	return;
}

void tStorm::setSeasonMonth(int sID, int val1, int val2)
{ 
	assert(sID > 0 && sID < 13);
	SeasMo[sID-1][0] = val1;
	SeasMo[sID-1][1] = val2;
	return;
}

void tStorm::setSeasonPMean(int sID, double val) 
{ 
	assert(sID > 0 && sID < 13);
	MRate[sID-1] = val;
	return;
}

void tStorm::setSeasonStmDurMean(int sID, double val) 
{ 
	assert(sID > 0 && sID < 13);
	MStmDr[sID-1] = val;
	return;
}

void tStorm::setSeasonSplDurMean(int sID, double val) 
{ 
	assert(sID > 0 && sID < 13);
	MSplDr[sID-1] = val;
	return;
}

void tStorm::allocSeasonMemory(int val) 
{ 
	assert(val > 0 && val < 13);
	
	// Assign number of seasons variable
	NumSeas = val;
	
	SeasMo = new int* [val];
	for (int i=0; i < val; i++) {
		SeasMo[i] = new int[2];
		assert(SeasMo[i] != 0);
	}
	MStmDr = new double[val];
	assert(MStmDr != 0);
	
	MSplDr = new double[val];
	assert(MSplDr != 0);
	
	MRate = new double[val];
	assert(MRate != 0);
	
	currSeasID = 0;
	return;
}

/***************************************************************************
**
** tStorm::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/
void tStorm::writeRestart(fstream & rStr) const
{
  BinaryWrite(rStr, RealRn);
  BinaryWrite(rStr, rid);
  BinaryWrite(rStr, optStoch);
  BinaryWrite(rStr, stdurMean);
  BinaryWrite(rStr, istdurMean);
  BinaryWrite(rStr, pMean);
  BinaryWrite(rStr, p);
  BinaryWrite(rStr, stdur);
  BinaryWrite(rStr, istdur);
  BinaryWrite(rStr, p0);
  BinaryWrite(rStr, stdur0);
  BinaryWrite(rStr, istdur0);
  BinaryWrite(rStr, pdev);
  BinaryWrite(rStr, stdurdev);
  BinaryWrite(rStr, istdurdev);
  BinaryWrite(rStr, twoPiLam);
  BinaryWrite(rStr, seed);
  BinaryWrite(rStr, endtm);

  BinaryWrite(rStr, currSeasID);
  BinaryWrite(rStr, NumSeas);
}

/***************************************************************************
**
** tStorm::readRestart() Function
**
***************************************************************************/
void tStorm::readRestart(fstream & rStr)
{
  BinaryRead(rStr, RealRn);
  BinaryRead(rStr, rid);
  BinaryRead(rStr, optStoch);
  BinaryRead(rStr, stdurMean);
  BinaryRead(rStr, istdurMean);
  BinaryRead(rStr, pMean);
  BinaryRead(rStr, p);
  BinaryRead(rStr, stdur);
  BinaryRead(rStr, istdur);
  BinaryRead(rStr, p0);
  BinaryRead(rStr, stdur0);
  BinaryRead(rStr, istdur0);
  BinaryRead(rStr, pdev);
  BinaryRead(rStr, stdurdev);
  BinaryRead(rStr, istdurdev);
  BinaryRead(rStr, twoPiLam);
  BinaryRead(rStr, seed);
  BinaryRead(rStr, endtm);

  BinaryRead(rStr, currSeasID);
  BinaryRead(rStr, NumSeas);
}

/**************************************************************************
**
**  tStorm "get" routines: return various variables
**
**************************************************************************/

double tStorm::getStormDuration() const { return stdur; }
double tStorm::interstormDur()    const { return istdur; }
double tStorm::getRainrate()      const { return p; }
double tStorm::getMeanStormDur()  const {return stdurMean;}
double tStorm::getMeanInterstormDur() const {return istdurMean;}
double tStorm::getMeanPrecip()    const {return pMean;}
int tStorm::getoptStorm() const {return optStoch;}
int tStorm::getSeasonID() { return currSeasID; }

//=========================================================================
//
//
//			End of tStorm.cpp                   
//
//
//=========================================================================
