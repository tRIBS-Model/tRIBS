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
**  tRunTimer.cpp: Functions for class tRunTimer
**
***************************************************************************/

#include "src/tSimulator/tRunTimer.h"
#include "src/Headers/globalFns.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//                  Section 1: tRunTimer Constructors/Destructors
//
//
//=========================================================================

/*****************************************************************************
** 
**  tRunTimer::tRunTimer()
**
**  A default constructor
**  A constructor that sets the run duration and output interval 
**  (and if desired sets the option to print time steps to stdout; 
**  A constructor that reads run duration and output interval from a 
**  tInputFile object 
**
*****************************************************************************/

// Default Constructor
tRunTimer::tRunTimer()
{	
	currentTime = 0;
	endTime = 1;
	outputInterval = 1;
	SPOutputInterval = 0;
	nextOutputTime = 0;
	nextSPOutputTime = 0;
	tstep = 0;
}

tRunTimer::tRunTimer( double duration, int opint )
{
	currentTime = 0;
	endTime = duration;
	outputInterval = opint;
	SPOutputInterval = 0;
	nextOutputTime = 0;
	nextSPOutputTime = 0;
	tstep = 0;
}

tRunTimer::tRunTimer( tInputFile &infile )
{
	int  cum = 0;
	int tmp[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

	for (int k=0; k < 13; k++) {  // Initialization of arrays
		cum += tmp[k];
		days[k] = tmp[k]; 
		cumdays[k] = cum;
	}
	InitializeTimer( infile );
	
}

void tRunTimer::InitializeTimer( tInputFile &infile )
{
	char namos[50];
	
	// Start Time
	if (infile.IsItemIn( "STARTDATE" ))
		infile.ReadItem( namos, "STARTDATE" );
	else if (infile.IsItemIn( "STARTINGDATETIME" ))
		infile.ReadItem( namos, "STARTINGDATETIME" );
	else {
		cerr<<"\nThe start date is not specified"
		<<"\nCorrect the .in file and re-run the model"<<endl;
		exit(2);
	}

	minuteS = 0; // GMnSKY2008MLE

	sscanf(namos, "%02d/%02d/%04d/%02d/%02d",
		   &monthS,&dayS,&yearS,&hourS,&minuteS);
	
	if (fabs((float)minuteS) > 60)
		minuteS = 0;
	
	minute = minuteRn = minuteS;
	hour   = hourRn   = hourS;
	day    = dayRn    = dayS;
	month  = monthRn  = monthS;
	year   = yearRn   = yearS;
	
	// if (hour > 24 || day > days[month] || day<0 || month>12 || month<0) {
	if (hour > 24 || day > days[month] || day<1 || month>12 || month<1) { // SKYnGM2008LU: check also later if hour > 23, and maybe add hour <0
		cout << "\nError: Check hours & dates of input data\n";
		cout <<"Exiting Program...\n"<<endl;
		exit(2);
	}
	
	Cout<<"\nStart TIME: \tYEAR"<<"\tMONTH"<<"\tDAY"
		<<"\tHOUR"<<"\tMINUTE"<<endl;
	Cout<<"\t\t"<<yearS<<"\t"<<monthS<<"\t"<<dayS<<"\t"<<hourS<<"\t"<<minuteS<<endl;
	
	// End Time
	endTime = infile.ReadItem( endTime, "RUNTIME" ); 
	Cout<<"\nTotal TIME:\t "<<endTime<<" hours";
	
	currentTime = 0; 
	for (int i=0; i<(int)endTime; i++)
		Advance( 1.0 );
	
	Cout<<"\n\nEnd TIME: \tYEAR"<<"\tMONTH"<<"\tDAY"
	    <<"\tHOUR"<<"\tMINUTE"<<endl;
	Cout<<"\t\t"<<year<<"\t"<<month<<"\t"<<day<<"\t"<<hour<<"\t"<<minute<<endl;
	
	minute= minuteS; 
	hour  = hourS;
	day   = dayS;    
	month = monthS;
	year  = yearS;
	
	// Output time steps 
	
	outputInterval = infile.ReadItem( outputInterval, "OPINTRVL" );
	Cout<<"\nOutput INTERVAL: \t\t\t"<<outputInterval<<" hours"<<endl;
	
	SPOutputInterval = infile.ReadItem( SPOutputInterval, "SPOPINTRVL" );
	Cout<<"\nSpatial Output INTERVAL: \t\t"<<SPOutputInterval<<" hours"<<endl;
	
	// Rainfall time steps and increments
	
	dtRain = infile.ReadItem(dtRain, "RAININTRVL");
	Cout<<"Rainfall Input TIME STEP: \t\t"<<dtRain<<" hours"<<endl;
	
	if ( (dtRain > 1) && (ceil(dtRain) != floor(dtRain)) ) {
		cout<<"\nError: Rainfall time increment larger than 1 hour "
		<<"must be multiple of 1 hour, e.g. 2, 3, 5, 10, etc."<<endl;
		cout<<"Exiting program..."<<endl;
		exit(2);}
	
	else if ( (dtRain < 1) && (fmod(1.0, dtRain) != 0) ) {
		cout<<"\nError: Rainfall time increment less than 1 hour "
		//<<"must correspond to one of its multiples, e.g. 0.25, 0.5, etc."<<endl;
		<<"must correspond to one of its divisors, e.g. 0.25, 0.5, etc."<<endl; // SKYnGM2008LU
		cout<<"Exiting program..."<<endl;
		exit(2);}
	
	metstep = infile.ReadItem( metstep, "METSTEP");
	Cout <<"Meteorological Input TIME STEP: \t"<<metstep/60.<<" hours"<<endl;
	metstep /= 60.;
	
	// Computation Time Steps
	
	tstep = infile.ReadItem( tstep, "TIMESTEP");
	Cout<<"Unsaturated zone TIME STEP: \t\t"<<tstep<<" minutes"<<endl;
	tstep /= 60.;
	
	gwatstep = infile.ReadItem( gwatstep, "GWSTEP"); 
	Cout<<"Saturated Zone TIME STEP: \t\t"<<gwatstep<<" minutes"<<endl;
	gwatstep /= 60.;
	
	// Select minimum of rainfall, met time step for ET/I
	// if (dtRain > metstep)
	if (dtRain >= metstep) // SKY2008Snow 
		etistep = metstep;
	else
		etistep = dtRain;

	// Rainfall Forecast Times
	optForecast = infile.ReadItem(optForecast, "FORECASTMODE");
	if (optForecast != 0) {
		fTime = infile.ReadItem(fTime, "FORECASTTIME");
		fLength = infile.ReadItem(fLength, "FORECASTLENGTH");
		fLead = infile.ReadItem(fLead, "FORECASTLEADTIME");
	}
	
	// Adjustments to Time with Rainfall Input
	
	double help;
	help = 0.0;  
	
	startHour = help;
	currentTime = help;
	RainTime = currentTime;
	// Advance Rain Time
	addRainTime(); 
	
	endTime += help;
	nextOutputTime = help + outputInterval;
	nextSPOutputTime = help + SPOutputInterval;
	
	//Initialize met time as tstep to get first hour
	MetTime = EtITime = tstep;
	
	//Initialize storm time for stochastic rainfall
	StormTime_1 = StormTime = 0;
	
}

//=========================================================================
//
//
//                  Section 2: tRunTimer Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  tRunTimer::addRainTime()
**  
**  Increments RainTime by period between successive input rainfall files. 
**
*****************************************************************************/
void tRunTimer::addRainTime()
{
	double minutet;  //Local minute as double
	
	RainTime += dtRain; 
	correctCalendarTime(RainTime, dtRain, &minuteRn, 
						&hourRn, &dayRn, &monthRn, &yearRn);
	return;
}

/*****************************************************************************
**  
**  tRunTimer::addMetTime() and getMetTime()
**  
**  Increments MetTime and EtITime using flag
**
*****************************************************************************/
void tRunTimer::addMetTime(int flag)
{
	if (flag == 1)
		MetTime += metstep; 
	if (flag == 2)
		EtITime += etistep;
	return;
}

double tRunTimer::getMetTime(int flag) 
{
	if (flag == 1)
		return MetTime;
	if (flag == 2)
		return EtITime;
	return 0;
}

/*****************************************************************************
**  
**  tRunTimer::isGaugeTime(double dt)
**  
**  Checks to see if it is time for rain gauge data input. 
**
*****************************************************************************/
int tRunTimer::isGaugeTime(double dt)
{
	if((currentTime-tstep)/dt == floor((currentTime-tstep)/dt))
		return 1;
	else
		return 0;
}

/*****************************************************************************
**  
**  tRunTimer::Start(double start, double end )
**   
**  Sets the current time and run duration.
**
*****************************************************************************/
void tRunTimer::Start( double start, double end ) 
{
	assert( end > start );
	currentTime = start;
	endTime = end;
}

/*****************************************************************************
**  
**  tRunTimer::Advance( double dt )
**   
**  Increments the current time by dt. Returns 1 if there is still time
**  remaining, 0 if the time is up.
**
*****************************************************************************/
int tRunTimer::Advance( double dt ) 
{     
	currentTime += dt;
	correctCalendarTime(currentTime, dt, &minute, &hour, &day, &month, &year);
	return( currentTime < endTime );
}

/*****************************************************************************
**  
**  tRunTimer::correctCalendarTime()
**
**  If change of hour occurs, this function updates all other calendar vars
**  NOTE: The function is supposed to be called AFTER generic time 'cTime'
**  has been incremented with a generic time step 'dt'   
**
*****************************************************************************/
void tRunTimer::correctCalendarTime(double cTime, double dt, int *mi, 
									int *hr, int *dy, int *mo, int *yr)
{ 
	int cntHrs = 0;
	double minutet;  //Local minute as double
	
	minutet = cTime - floor(cTime);
	
	if (minutet == 0 && dt < 1) { //Correction for dt < 1
		(*mi) = 0;
		cntHrs = 1;
	}
	else if (minutet == 0 && dt >= 1) {
		(*mi) = 0;
		cntHrs = (int)(ceil(dt));
	}
	else if (minutet != 0 && dt < 1) {
		(*mi) = (int)(minutet*60);         // Convert value to minutes
		cntHrs = 0;
	}
	else if (minutet != 0 && dt >= 1) {
		(*mi) = (int)((dt - floor(dt))*60);
		cntHrs = (int)(floor(dt));
	}
	else {  // 1 HOUR increment
		(*mi) = 0;
		cntHrs = 1;
	}
	
	// Adjust calendar time for required # of hours
	for (int i = 0; i < cntHrs; i++) {
		
		// Increment hours
		(*hr)++;
		
		if ((*hr) == 24) {
			(*dy)++; 
			(*hr) = 0;
		}
		
		if ((*mo) == 2 && (*dy) > days[(*mo)]) { // February
			if ((*hr) == 0) {
				if ( ((*yr)%4) != 0 || ((*dy) > 29)) {
					(*mo)++; // After/Instead February 29th
					(*dy) = 1;
				}
			}
		}
		else if ((*dy) > days[(*mo)]) {
			(*mo)++;
			(*dy) = 1;
			if ((*mo) > 12) {
				(*yr)++;
				(*mo) = 1;
			}
		}
	}
	return;
}

/*****************************************************************************
**  
**  tRunTimer::CheckOutputTime()
**   
**  Checks to see whether it's time to write output yet.
**
*****************************************************************************/
int tRunTimer::CheckOutputTime() 
{	
	if ( currentTime >= nextOutputTime ) {
		nextOutputTime += outputInterval;
		return 1;
	}
	else return 0;
}

/*****************************************************************************
**  
**  tRunTimer::CheckSpatialOutputTime()
**   
**  Checks to see whether it's time to write output yet.
**
*****************************************************************************/
int tRunTimer::CheckSpatialOutputTime() 
{	
	if ( currentTime >= nextSPOutputTime ) {
		nextSPOutputTime += SPOutputInterval;
		return 1;
	}
	else return 0;
}

/*****************************************************************************
**  
**  tRunTimer::IsFinished()
**   
**  Checks to see whether the run is finished.
**
*****************************************************************************/
int tRunTimer::IsFinished() 
{	
	return( currentTime >= endTime );
}

//=========================================================================
//
//
//                  Section 3: tRunTimer Get/Set Functions
//
//
//=========================================================================
double tRunTimer::RemainingTime(double t) { 
	return (endTime - t); }

double tRunTimer::getOutputIntervalSec() {
	return outputInterval*3600.0; }

double tRunTimer::getSpatialOutputIntervalSec() {
	return SPOutputInterval*3600.0; }

int tRunTimer::getElapsedSteps(double timeelapsed) {  
	return (int)( timeelapsed/tstep ); }  

int tRunTimer::getElapsedETISteps(double timeelapsed) {  
	return (int)( timeelapsed/etistep ); }  

int tRunTimer::getElapsedMETSteps(double timeelapsed) {  
	return (int)( timeelapsed/metstep ); }  

int tRunTimer::getResStep(double time) {
	return ( (int)floor((time+currentTime)/outputInterval) ); }  

double tRunTimer::get_abs_hour(double t) {
	return (currentTime+t); }

int tRunTimer::getStepForSpecifiedDT(double time, double dtrf) {
	return ( (int)ceil(time/dtrf) ); } 

//=========================================================================
//
//
//                  Section 4: tRunTimer TimeStep Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  tRunTimer::res_time_mid(step, ihour, imin)
**   
**  Compute the time (hours, minutes) corresponding to the center of a 
**  time step. 
**
**  Algorithm:
**    obtain absolute hour of center of time step
**    get hours and minutes
**
*****************************************************************************/
void tRunTimer::res_time_mid(int step, int *ihour, int *imin) 
{ 
	double t_hour=((double)step+.5)*outputInterval + startHour;
	*ihour = (int)t_hour;
	*imin  = (int)((t_hour-(double)(*ihour))*60.+.5);
	return;
}

double tRunTimer::res_hour_mid(int it) 
{
	return startHour + ((double)it+.5)*outputInterval;
}

/*****************************************************************************
**  
**  tRunTimer::res_time_begin(step,ihour,imin)
**             res_hour_begin(int it)
**
**  Compute the time (hours, minutes) corresponding to the beginning of a 
**  time step. 
**
**  Algorithm:
**    obtain absolute hour of center of time step
**    get hours and minutes
**
*****************************************************************************/
void tRunTimer::res_time_begin(int step, int *ihour, int *imin) 
{
	double t_hour=(double)step*outputInterval + startHour;
	*ihour=(int)t_hour;
	*imin=(int)((t_hour-(double)(*ihour))*60.+.5);
	return;
}

double tRunTimer::res_hour_begin(int it) 
{
	return startHour + (double)it*outputInterval;
}

/*****************************************************************************
**  
**  tRunTimer::res_time_end(step,ihour,imin)
**   
**  Compute the time (hours, minutes) corresponding to the end of a 
**  time step. 
**
**  Algorithm:
**    obtain absolute hour of center of time step
**    get hours and minutes
**
*****************************************************************************/
void tRunTimer::res_time_end(int step, int *ihour, int *imin) 
{
	double t_hour=((double)step+1.0)*outputInterval + startHour;
	*ihour=(int)t_hour;
	*imin=(int)((t_hour-(double)(*ihour))*60.+.5);
	return;
}

double tRunTimer::res_hour_end(int it) 
{
	return startHour + (double)(it+1)*outputInterval;
}

/***************************************************************************
**
** tRunTimer::Get Rainfall Forecasting Functions
**
**
***************************************************************************/

double tRunTimer::getfTime() { return fTime; }

double tRunTimer::getfLength() { return fLength; }

double tRunTimer::getfLead() { return fLead; }

int tRunTimer::getoptForecast() { return optForecast; }

/***************************************************************************
**
** tRunTimer::Stochastic Storm Functions
**
**
***************************************************************************/
void tRunTimer::UpdateStorm( double update ) 
{
	StormTime_1 = StormTime;
	StormTime = StormTime + update;
	return;
}

double tRunTimer::getStormTime() 
{
	return StormTime;
}

double tRunTimer::getPrevStormTime() 
{
	return StormTime_1;
}

/***************************************************************************
**
** tWaterBalance::writeRestart() Function
** 
** Called from tSimulator during simulation loop
** 
***************************************************************************/
void tRunTimer::writeRestart(fstream & rStr) const
{
  for (int i = 0; i < 13; i++) {
    BinaryWrite(rStr, days[i]);
    BinaryWrite(rStr, cumdays[i]);
  } 
  BinaryWrite(rStr, minute);
  BinaryWrite(rStr, hour);
  BinaryWrite(rStr, day);
  BinaryWrite(rStr, month);
  BinaryWrite(rStr, year);
  BinaryWrite(rStr, minuteRn);
  BinaryWrite(rStr, hourRn);
  BinaryWrite(rStr, dayRn);
  BinaryWrite(rStr, monthRn);
  BinaryWrite(rStr, yearRn);
  BinaryWrite(rStr, minuteS);
  BinaryWrite(rStr, hourS);
  BinaryWrite(rStr, minS);
  BinaryWrite(rStr, dayS);
  BinaryWrite(rStr, monthS);
  BinaryWrite(rStr, yearS);
  BinaryWrite(rStr, optForecast);
  BinaryWrite(rStr, dtRain);
  BinaryWrite(rStr, startHour);
  BinaryWrite(rStr, currentTime);
  BinaryWrite(rStr, nextOutputTime);
  BinaryWrite(rStr, nextSPOutputTime);
  BinaryWrite(rStr, tstep);
  BinaryWrite(rStr, gwatstep);
  BinaryWrite(rStr, metstep);
  BinaryWrite(rStr, outputInterval);
  BinaryWrite(rStr, SPOutputInterval);
  BinaryWrite(rStr, RainTime);
  BinaryWrite(rStr, fTime);
  BinaryWrite(rStr, fLength);
  BinaryWrite(rStr, fLead);
  BinaryWrite(rStr, etistep);
  BinaryWrite(rStr, MetTime);
  BinaryWrite(rStr, EtITime);
  BinaryWrite(rStr, StormTime);
  BinaryWrite(rStr, StormTime_1);
}

/***************************************************************************
**
** tRunTimer::readRestart() Function
**
***************************************************************************/
void tRunTimer::readRestart(fstream & rStr)
{
  for (int i = 0; i < 13; i++) {
    BinaryRead(rStr, days[i]);
    BinaryRead(rStr, cumdays[i]);
  }
  BinaryRead(rStr, minute);
  BinaryRead(rStr, hour);
  BinaryRead(rStr, day);
  BinaryRead(rStr, month);
  BinaryRead(rStr, year);
  BinaryRead(rStr, minuteRn);
  BinaryRead(rStr, hourRn);
  BinaryRead(rStr, dayRn);
  BinaryRead(rStr, monthRn);
  BinaryRead(rStr, yearRn);
  BinaryRead(rStr, minuteS);
  BinaryRead(rStr, hourS);
  BinaryRead(rStr, minS);
  BinaryRead(rStr, dayS);
  BinaryRead(rStr, monthS);
  BinaryRead(rStr, yearS);
  BinaryRead(rStr, optForecast);
  BinaryRead(rStr, dtRain);
  BinaryRead(rStr, startHour);
  BinaryRead(rStr, currentTime);
  BinaryRead(rStr, nextOutputTime);
  BinaryRead(rStr, nextSPOutputTime);
  BinaryRead(rStr, tstep);
  BinaryRead(rStr, gwatstep);
  BinaryRead(rStr, metstep);
  BinaryRead(rStr, outputInterval);
  BinaryRead(rStr, SPOutputInterval);
  BinaryRead(rStr, RainTime);
  BinaryRead(rStr, fTime);
  BinaryRead(rStr, fLength);
  BinaryRead(rStr, fLead);
  BinaryRead(rStr, etistep);
  BinaryRead(rStr, MetTime);
  BinaryRead(rStr, EtITime);
  BinaryRead(rStr, StormTime);
  BinaryRead(rStr, StormTime_1);
}

//=========================================================================
//
//
//                          End of tRunTimer.cpp
//
//
//=========================================================================
