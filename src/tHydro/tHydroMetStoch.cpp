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
**  tHydroMetStoch.cpp: Functions for class tHydroMetStoch 
**
***************************************************************************/

#include "src/tHydro/tHydroMetStoch.h"

//=========================================================================
//
//
//           Section 1: tHydroMetStoch Constructors and Destructors
//
//
//=========================================================================
tHydroMetStoch::tHydroMetStoch()
{
	gmt = 0;
	latitude = 0.0;
	longitude = 0.0;
	basinLong = 0.0;
	basinLat = 0.0;
}

tHydroMetStoch::tHydroMetStoch(tMesh<tCNode> *gptr, tRunTimer *t, tInputFile &infile, 
							   tEvapoTrans *et, tRainfall *storm)
{
	gridPtr = gptr;
	timer   = t;
	etPtr   = et;
	rainPtr = storm;
	
	double tmpa[3][3] = { {0.567, 0.086, -0.002}, {0.253, 0.504, -0.05}, {-0.006, -0.039, 0.244} };
	double tmpb[3][3] = { {0.782, 0, 0}, {0.328, 0.637, 0}, {0.238, -0.341, 0.873} };
	A = new double* [3];
	assert(A != 0);
	B = new double* [3];
	assert(B != 0);
	
	for (int i=0; i<3; i++) {
		A[i] = new double[3];
		assert(A[i] != 0);
		B[i] = new double[3];
		assert(B[i] != 0);
		for (int j=0; j<3; j++) {
			A[i][j] = tmpa[i][j];
			B[i][j] = tmpb[i][j];
		}
	}
	
	// Initialize hydromet variables
	SetStochasticHydroMet( infile );
	
}

tHydroMetStoch::~tHydroMetStoch()
{
	for (int i=0; i<3; i++) {
		delete [] A[i];       
		delete [] B[i];
	}
	delete [] A;
	delete [] B;
	hout.close();
}


//=========================================================================
//
//
//                  Section 2: tHydroMetStoch Functions
//
//
//=========================================================================

/***************************************************************************
**
**  SetStochasticHydroMet
**
**  Sets basic variables of the class and provides initial estimates 
**
***************************************************************************/
void tHydroMetStoch::SetStochasticHydroMet(tInputFile &infile)
{
	int i;
	double randv[3];
	char WeatherFile[kMaxNameSize];
	
	// Initialize time
	minute = timer->minute; 
	hour   = timer->hour;
	day    = timer->day;    
	month  = timer->month;
	year   = timer->year;
	
	// Check if the start hour is == 0
	if (hour) {
		cerr<<"\n\tCan not run stochastic simulator for start hour > 0!"<<endl;
		cerr<<"\tAssign hour in STARTDATE to 0!"<<endl;
		cerr<<"\n\tExiting..."<<endl<<flush;
		exit(2);
	}
	
	Idir=Idir_vis=Idir_nir=Idif=Idif_vis=Idif_nir=0.0;
	
	// Read parameters from input file
	infile.ReadItem(WeatherFile, "WEATHERTABLENAME"); // input table
	
	// Read weather parameters from a file
	ReadWeatherParameters( WeatherFile );
	
	// Estimate expected value of annual precipitation
	// "JDay1" and "Jday2" specify the period for which the
	// storm parameters are representative, e.g. growing season
	double JDay1, JDay2;
	JDay1 = 0.;
	JDay2 = 365;
	if (rainPtr->getMeanInterstormDur()+rainPtr->getMeanStormDur() > 0) {
		AnnRain = (JDay2-JDay1)*24./(rainPtr->getMeanInterstormDur()+rainPtr->getMeanStormDur());
		AnnRain = AnnRain*rainPtr->getMeanStormDur()*rainPtr->getMeanPrecip();
	}
	else 
		AnnRain = 0.0;
	
	RainOpt = infile.ReadItem(RainOpt, "STOCHASTICMODE");
	if (RainOpt == 6) {
		rainPtr->setSeasonID( rainPtr->DefineSeason(timer->month) );
		rainPtr->updateSeasonVars();
	}
	
	// Simulate atmospheric pressure for the first DAY
	Patm_mean = 1000;  // default value
	Patm_std = 6.0;    // default value
	atmPress = EstimateAR1Var(Patm_mean, Patm_std, 0.85, Patm_mean);
	atmPress = 1013;   // [mb] 
	
	// Simulate wind speed for the first hour
	SetWindVars();
	windSpeed = WindSp_mean; // random_Weibull(WindSp_mean, WindSh_fact);
	
	// Simulate cloud cover for the first hour
	SetSkyVars();
	skyCover_1 = 10;
	skyCover = SimulateSkyCover(1);
	
	// Initialize temperature model
	TempOpt = 1;
	
	// Curtis and Eagleson [1982] air temperature model
	if (TempOpt == 1) {
		SetTemperatureVars();
	}

	// SKY2008Snow from AJR2007 starts here
	// =============================================================
	// Initialize air and dew point temperature estimates
	// Generate mean monthly temperature from N(Tmon_i, sgma^2) 
	else if (TempOpt == 2) {
		T_ro1 = 0.7;
		SetTemperatParameters();

		for (i=0; i<2; i++) {
			// Initialization estimates 
			if ( !i ) {
				// Initial estimates for some variables (td, etc. for day 1)
				td = tdnext = tmone;
				deltt = sgma; // HERE *2

				// Obtain intitial estimate of the residuals vector
				GetVectRandomNorm(randv, 3, 0.0, 1.0);
				MatrixVectorProduct(B, 3, 3, randv, hi);

				// Approximate expected values for expected tmax and tmin
				ExpectedTmaxTmin();

				// Use multivariate process model to get expected tmax and tmin
				// NOTE: tmaxe & tmine can be conditioned based on rainfall events
				SimulateTmaxTmin(tmaxe, sgma, tmine, sgma);
 
				// Define time shift for temperature peak in sin function
				SetPeakTemperatShift();
			}
			// Estimate diurnal air and dew point temperature cycles for day 1
			else
				SimulateAirDPTemperature();
		}
	}
	// SKY2008Snow from AJR2007 ends here

	// Open output file for generated variables 
	if (etPtr->simCtrl->Verbose_label == 'Y') {
#ifdef ALPHA_64
		hout.open("_hydrometout.dat");
#elif defined LINUX_32
		hout.open("_hydrometout.dat");
#elif defined WIN
		hout.open("_hydrometout.dat");
#else 
		hout.open("_hydrometout.dat");
#endif
	}
	
	if (!hout.good()) {
		cerr<<"\nFile not created!\nExiting Program..."<<endl;
		exit(2);
	}
	hout.setf(ios::fixed, ios::floatfield);
	
	return;
}

/***************************************************************************
**
**  SimulateHydrometVars
**
**  Provides estimates of basic hydrometeorologic variables for the current
**  hour. Either uses values defined previously or simulates diurnal cycles.
**
***************************************************************************/
void tHydroMetStoch::SimulateHydrometVars()
{  
	// ===================================
	// Update variables if it is mid-night
	// ===================================
	if (!timer->hour && hour != timer->hour) {
		ComputeDailyEpCld(0.0, 0.0);
		
		// Simulate Atmospheric pressure
		//atmPress = EstimateAR1Var(Patm_mean, Patm_std, 0.85, atmPress);
		atmPress = 1013.0; // [mb]
		
		// For temperature model make required f-n calls
		if (TempOpt == 1) {
			; // Empty for now
		}
		// SKY2008Snow from AJR2007
		else if (TempOpt == 2)
			SimulateAirDPTemperature();

	}
	
	// ===================================
	// Check if the month has changed...
	// ===================================
	if (month != timer->month) {
		
		if (CldParOpt == 2)
			SetCloudinessVars();
		
		if (TempOpt == 1) 
			SetTemperatureVars();

		// SKY2008Snow from AJR2007
		if (TempOpt == 2)
			SetTemperatParameters();

		if (RainOpt == 6) {
			rainPtr->setSeasonID( rainPtr->DefineSeason(timer->month) );
			rainPtr->updateSeasonVars();
		}
	}
	
	// ==============================================
	// Generate wind speed as an INDEPENDENT variable
	// Direction is NOT currently simulated
	// ==============================================
	// windSpeed = random_Weibull(WindSp_mean, WindSh_fact);
	windSpeed = SimulateWindCurtis();
	
	// ##### Special simulation option: no randomness in climate 
	if (etPtr->simCtrl->smooth_weather == 'Y')
		windSpeed = 3.0;
	
	// =============================================================
	// Simulate sky cover (or as fraction of 2D rainfall over basin) 
	// =============================================================
	skyCover_1 = skyCover;     // (t-1) cloudiness
	if (timer->getCurrentTime() <= 
		(timer->getStormTime()-rainPtr->interstormDur()) && rainPtr->getRainrate()) {
		skyCover = 10.0;
	}
	// ===========================================
	// Use a non-stationary model of Curtis [1982]
	// ===========================================
	else
		skyCover = SimulateSkyCover(1);
	
	// ##### Special simulation option: no randomness in climate 
	if (etPtr->simCtrl->smooth_weather == 'Y')
		skyCover = 0.0;
	
	// =====================================================
	// Simulate pressure and temperatures
	// Check if it is called for the first time
	// =====================================================
	if (!timer->hour && hour != timer->hour) {
		atmPress = EstimateAR1Var(Patm_mean, Patm_std, 0.85, atmPress);
		if (TempOpt == 1) {
			; // Empty for now
		}

		// SKY2008Snow from AJR2007
		else if (TempOpt == 2)
			SimulateAirDPTemperature();

	}
	
	// =====================================================
	// Simulate temperature / obtain from daily array
	// =====================================================
	if (TempOpt == 1) {
		airTemp_1 = airTemp; // (t-1) temperature
		airTemp = SimulateAirDPTemperatureCurtis();
		dewTemp = SimulateDewPointTemperature(timer->hour - etPtr->getDeltaT(), 
											  etPtr->getSunRiseHour()); 
	}
	
	// SKY2008Snow from AJR2007
	else if (TempOpt == 2) {
		airTemp = DiurnT[timer->hour];
		dewTemp = DiurnDP[timer->hour];
	}

	// =====================================================
	// Condition dew point temperature on rainfall occurence
	// CANCELLED: leads to very high std of Tdew
	// =====================================================
	if (timer->getCurrentTime() <= 
		(timer->getStormTime()-rainPtr->interstormDur()) &&
		rainPtr->getRainrate()) {
		//dewTemp = airTemp;    // ######
		//dewTemp_1 = dewTemp;
	}
	
	// Store current values of time variables in the object
	year   = timer->year;
	month  = timer->month;
	day    = timer->day;
	hour   = timer->hour;
	minute = timer->minute;
	
	return;
}

/***************************************************************************
**
**  OutputHydrometVars()
**
**  Outputs simulated hydrometeorological variables to a file
**  "_hydrometout.dat"
**
***************************************************************************/
void tHydroMetStoch::OutputHydrometVars()
{
	double tmp;
	char tempstr[200];
	
	//sprintf(tempstr,"%02d%02d%04d %02d", month, day, year, hour);
	
	// Output hydrometeorologic file if needed
	// Various format options are considered: DO NOT eliminate commented out!
	
	/*
		// For post-processing routines that need radiation fluxes
	 if (SinH < 0) {
		 hout<<tempstr<<" "<<etPtr->getinShortWave()<<"\t"<<etPtr->getinLongWave()<<"\t"
         <<airTemp<<"\t"<<dewTemp<<"\t"
         <<windSpeed<<"\t"<<rainPtr->getRainrate()<<"\t"<<skyCover<<"\t0.0\t0.0"<<endl;
	 }
	 else {
		 hout<<tempstr<<" "<<etPtr->getinShortWave()<<"\t"<<etPtr->getinLongWave()<<"\t"
         <<airTemp<<"\t"<<dewTemp<<"\t"
         <<windSpeed<<"\t"<<rainPtr->getRainrate()<<"\t"<<skyCover<<"\t"
		 <<Idir*SinH<<"\t"<<Idif<<endl;
	 }
	 */
	
	/*
		if (tNextRn < 1.0 && tNextRn > 0.0 && !rainPtr->getRainrate())
	 tmp = 1.0;
	 else 
	 tmp = rainPtr->getRainrate();
	 
	 hout<<tempstr<<" "<<SunH<<"\t"<<tmp<<"\t"<<skyCover<<"\t"
	 <<Idir<<"\t"<<Idir_vis<<"\t"<<Idir_nir<<"\t"
	 <<Idif<<"\t"<<Idif_vis<<"\t"<<Idif_nir<<"\t"
	 <<airTemp<<"\t"<<dewTemp<<"\t"<<windSpeed<<endl;
	 */
	
	// 'ParameterEstimation.cpp' input format
	/*
		sprintf(tempstr,"%02d%02d%04d%02d", month, day, year, hour);
	 hout<<tempstr<<"\t"
	 <<airTemp<<"\t"<<dewTemp<<"\t"<<windSpeed<<"\t"
	 <<rainPtr->getRainrate()<<"\t"<<skyCover<<endl;
	 */
	
	// The order of the first 10 variables corresponds to 
	// that of the output of 'ParameterEstimation.cpp' in 
	// order to match their location with the measured data
	/*
		if (SinH < 0)
     hout<<tempstr<<"\t"<<etPtr->julianDay()
	 <<"\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t"
	 <<skyCover<<"\t0.0\t0.0\t0.0\t0.0\t"
	 <<SinH<<"\t"<<SunH<<endl;  
	 else 
     hout<<tempstr<<"\t"<<etPtr->julianDay()<<"\t"
	 <<Io*SinH<<"\t"<<Io<<"\t"<<Idir<<"\t"
	 <<(Idir*SinH+Idif)<<"\t"<<Idir*SinH<<"\t"<<Idif<<"\t"
	 <<skyCover<<"\t"
	 <<Idir_vis*SinH<<"\t"<<Idir_nir*SinH<<"\t"
	 <<Idif_vis     <<"\t"<<Idif_nir<<"\t"
	 <<SinH<<"\t"<<SunH<<endl;  
	 */
	return;
}

/***************************************************************************
**
**  SimulateAirDPTemperatureCurtis
**
**  Provides estimates of basic hydrometeorologic variables for the current
**  hour. Either uses values defined previously or simulates diurnal cycles.
**
***************************************************************************/
double tHydroMetStoch::SimulateAirDPTemperatureCurtis()
{
	double Kbar, airT, tmp;
	double TL, TP, T12, T23, Gt_tp;
	
	// Compute daily model coefficients
	if (!timer->hour) 
		ComputeDailyTempCoeffs();
	
	// Convert standard time to local time [hour] 
	TL = timer->hour - etPtr->getDeltaT();
	
	// Beginning of simulation, t_prime [hour]
	TP = -etPtr->getDeltaT() - 1.0;
	T12 = 12. - etPtr->getDeltaT();
	
	// Attenuation factor due to cloud cover [-]
	// NOTE: 'skyCover' must be known
	Kbar = 1 - 0.65*skyCover*skyCover/100.;
	
	// Compute longwave radiation with vars of the previous time step
	qt_1 = inLongWave(getAirTemp_1(), getSkyCover_1()/10.);
	
	//qt_1 = etPtr->getinLongWave();
	
	// ## Factor to be consistent with the off-line routine ##
	qt_1 /= 1000.;
	
	tmp = (1-exp(-b1))*exp(b1*TL)/b1;
	
	I1 = K1*(exp(b1*TL)-exp(b1*TP));
	I4 = b4*tmp*qt_1 + I4_1;
	// The following components are currently DISABLED
	// They account for the effects of ground temperature,
	// wind speed and wind direction on air temperature
	I5 = 0.; // I5 = b5*tmp*grndTemp  + I5_1;
	I6 = 0.; // I6 = b6*tmp*windSpeed + I6_1;
	I7 = 0.; // I7 = b7*tmp*windDir   + I7_1;
	
	// Subtotal
	Gt_tp = I1 + I4 + I5 + I6 + I7;
	
	// Compute differently for different ranges
	// 1.) Before the sunrise -- no change in 'Gt_tp'
	if (TL > etPtr->getSunRiseHour()) {
		
		// 2.) First hour (or a fraction) after sunrise
		if (TL <= etPtr->getSunRiseHour() + 1.0)
			Function0(2, TL, etPtr->getSunRiseHour(), Kbar);
		
		// 3.) After sunrise and before noon
		else if (TL <= T12)
			Function0(3, TL, TL-1, Kbar);
		
		// 4.) First hour after local noon
		else if (TL <= T12 + 1.0)
			Function0(4, TL, TL-1, Kbar);
		
		// 5.) After local noon, before sunset
		else if (TL <= etPtr->getSunSetHour())
			Function0(5, TL, TL-1, Kbar);
		
		// 6.) First hour after local sunset
		else if (TL <= etPtr->getSunSetHour()+1)
			Function0(6, etPtr->getSunSetHour(), TL-1, Kbar);
		
		// 7.) After  sunset
		else if (TL > etPtr->getSunSetHour()+1 && TL < 23)
			;
		Gt_tp += (I2 + I3);
  }
	
	// --- Deterministic temperature component ---
	airT = K0*exp(-b1*(TL-TP)) + Gt_tp*exp(-b1*TL);
	
	// Store the temperature (its _deterministic_ component)
	if (timer->hour == 23)
		temp11pm = airT;
	
	// Add random component 
	Tdev_1 = Tdev;
	Tdev = EstimateAR1Var(MeanTDev, StdTDev, AutoCorTDev, Tdev_1);
	
	// ##### Special simulation option: no randomness in climate 
	if (etPtr->simCtrl->smooth_weather == 'Y')
		Tdev = 0.0;
	
	airT += Tdev;
	
	// Define the 'max' and 'min' values
	if (airT > Tmax)
		Tmax = airT;
	if (airT < Tmin)
		Tmin = airT; 
	
	// Store the 't-1' values
	I2_1 = I2;
	I3_1 = I3;
	I4_1 = I4;
	I5_1 = I5;
	I6_1 = I6;
	I7_1 = I7;
	
	// =========================================================
	/*  
		cout<<"\nHOUR = "<<timer->hour<<";\tAirT = "<<airT-Tdev
		<<";\tK0 term = "<<K0*exp(-b1*(TL-TP))
		<<";\tGt_tp term = "<<Gt_tp*exp(-b1*TL)<<endl<<flush;
	cout<<"\tClouds = "<<skyCover<<";\tShrtWv = "<<etPtr->getinShortWave()
		<<";\tLongWv = "<<etPtr->getinLongWave()<<endl<<flush;
	cout<<"K0 = "<<K0<<"\tK1 = "<<K1<<"\tK2 = "<<K2<<"\tK3 = "<<K3
		<<"\tK4 = "<<K4<<"\tK5 = "<<K5<<"\tK6 = "<<K6<<endl;
	cout<<"I1 = "<<I1<<";\tI2 = "<<I2<<";\tI3 = "<<I3
		<<";\tI4 = "<<I4<<endl<<flush;
	cout<<"Tdev = "<<Tdev<<"; Tdev_1 = "<<Tdev_1<<endl;
	*/
	// =========================================================
	
	
	// Store the tprime temperature (deterministic)
	if (timer->hour == 23) {
		I2=I3=I4=I2_1=I3_1=I4_1=I5_1=I6_1=I7_1=0.0;
	}
	
	return airT;
}

/***************************************************************************
**
**  inLongWave()
**
**  Computes incoming longwave radiation based on air temperature and clouds
**
***************************************************************************/
double tHydroMetStoch::inLongWave(double temp, double cld)
{
	double sigm, kCloud, airTempK, Ea, Rlin;
	sigm = 5.67e-8;
	kCloud = 1 + 0.17*pow(cld,2);
	airTempK = temp + 273.15;
	Rlin = kCloud*sigm*pow(airTempK,4);
	return Rlin;
}

/***************************************************************************
**
**  ComputeDailyTempCoeffs
**
**  Computes the required constants for a new day 
**  
***************************************************************************/
void tHydroMetStoch::ComputeDailyTempCoeffs()
{
	double b2p2, del, phi;
	
	// Define constants
	del = etPtr->getDeltaAngle();
	phi = etPtr->getPhiAngle();
	
	// Define coefficients
	b2p2 = b1*b1 + PI12*PI12;
	
	// Temperature at 23 pm of the previous day
	K0 = temp11pm; 
	
	// Shortwave radiation variables
	K1 = b0/b1;
	K2 = b2/b1*sin(del)*sin(phi);
	K3 = b1*b2*cos(del)*cos(phi)/b2p2;
	K4 = PI12*b2*cos(del)*cos(phi)/b2p2;
	K5 = PI12*PI12*b3*cos(del)*cos(phi)/b2p2;
	K6 = PI12*b1*b3*cos(del)*cos(phi)/b2p2;
	
	return;
}

/***************************************************************************
**
**  Function0, Function1, Function2, Function3
**
**  Provides estimates of variables that are used in the model of Curtis &
**  Eagleson [1982]. Their evaluation depend on time of the day.
**
***************************************************************************/
void tHydroMetStoch::Function0(int opt, double t1, double t2, double Kbar)
{
	// 'I2' does not change when (TL > Sunset)
	// 'I3' does not change when (TL > 12)
	// Thus previously computed values are used (they are stored 
	// as the data members of the 'tHydroMetStoch' class)
	
	// Calculation of I2
	if (opt >= 2 && opt <= 6) {
		I2  = Function1(t1, t2);
		I2 -= Function2(t1);
		I2 += Function2(t2);
		I2 *= Kbar;
		if (opt != 2)
			I2 += I2_1;
	}
	// Calculation of I3
	if (opt == 2 || opt == 3 || opt == 4) {
		if (opt == 4)
			I3  = Function3(12.0);
		else 
			I3  = Function3(t1);
		I3 -= Function3(t2);
		I3 *= Kbar; 
		if (opt != 2)
			I3 += I3_1;
	}
	return;
}

double tHydroMetStoch::Function1(double LocT, double SunR)
{
	return(K2*(exp(b1*LocT)-exp(b1*SunR)));
}

double tHydroMetStoch::Function2(double time)
{
	return(exp(b1*time)*(K3*cos(PI12*time) + K4*sin(PI12*time)));
}

double tHydroMetStoch::Function3(double time)
{
	return(exp(b1*time)*(K6*sin(PI12*time) - K5*cos(PI12*time)));
}

/***************************************************************************
**
**  SimulateSkyCover
**
**  Simulates cloudiness using a non-stationary model of Curtis and Eagleson
**  [1982]. While the mean cloudiness changes based on an exponential model
**  given information about the current and the following storm, the std
**  remains the same. The noise is introduced using the lag-1 Markov model.
**  Required parameters: mean cloudiness, std, lag-1 correlation coefficient,
**  rate of cloudiness growth/decay before/after a storm.
**
***************************************************************************/
double tHydroMetStoch::SimulateSkyCover(int opt)
{ 
	float betaVar;
	double mt, clouds, Pt;
	double denom, l1, l2, L1, L2;
	if (opt == 1) {
		tNextRn = timer->getStormTime()   - timer->getCurrentTime();
		tAftrRn = timer->getCurrentTime() - timer->getPrevStormTime() - 
			rainPtr->getStormDuration();
		// Get transition value first
		Pt = GetCloudTransitValue();
		if (Pt > 0.99)
			Pt = 1.0;
		
		// For 'fairweather' only
		if (fabs(Pt) > 0.99) {
			
			// Get initial limits
			l1 = M0 + (1-M0)*(1-Pt);
			l2 = 1 - l1;
			
			// Compute limits for the random term
			denom = sqrt(1-m_ro1*m_ro1)*StdSky*Pt;
			L1 = (-l1 - m_ro1*Pt*mt_1)/denom;
			L2 = ( l2 - m_ro1*Pt*mt_1)/denom;
			
			// Set appropriate parameters of the Beta d-n
			SetCondBetaPars( skyCover_1 );
			
			// Simulate Beta variable E [0, 1]
			betaVar = genbet(betaA, betaB);
			
			// Transform to the corresponding limits, i.e. E [L1, L2]
			betaVar = L1 + (L2-L1)*betaVar;
		}
		else 
			betaVar = random_normal(0.0,1.0);
		
		//mt = EstimateAR1Var(0.0, StdSky, m_ro1, mt_1);
		mt = m_ro1*mt_1 + sqrt(1-m_ro1*m_ro1)*StdSky*betaVar;
		
		mt_1 = mt;
		
		clouds = M0 + (1-M0)*(1-Pt) + mt*Pt;
		if (clouds > 1)
			clouds = 1.0;
		else if (clouds < 0) 
			clouds = 0.0;
		
		
		// Compute the cloud correction factors
		if (clouds > 0.0) {
			// Assume that the amount of low and moderately high
			// cloudiness is about (1 - Cloud transition value)
			// Cloud shape factor 'f0' (is around the same)
			Nlm = 1 - Pt;
			f0  = uniform(0.3, 1.0);
			
			// During the fairweather period, do not use 'Pt'
			if (fabs(Pt) > 0.99)
				Nlm = uniform(0.0, 1.0);
			
			if (timer->getCurrentTime() <= 
				(timer->getStormTime()-rainPtr->interstormDur()) &&
				rainPtr->getRainrate())
				Nlm = f0 = 1.0;
			
			// The rest of the cloudiness can be composed of cirrus
			Nci = uniform(0.0, (1.0-Nlm));
			
			// Scale back to the actual cloudiness
			Nlm *= clouds;
			Nci *= clouds;
		}
		
		// To scale to the values that are used 
		clouds *= 10.; 
	}
	return clouds;
}

/***************************************************************************
**
**  GetCloudTransitValue
**
**  Defines the P(t) value in the transition period
**
***************************************************************************/
double tHydroMetStoch::GetCloudTransitValue()
{ 
	return ((1-exp(-GammaSky*tAftrRn))*(1-exp(-GammaSky*tNextRn)));
}

/***************************************************************************
**
**  SetCondBetaPars
**
**  Sets parameters of the Beta distribution depending on cloudiness value
**
***************************************************************************/
void tHydroMetStoch::SetCondBetaPars(double cloud)
{
	int ind, mon;
	if (cloud-floor(cloud) > 0.5)
		ind = (int)(ceil(cloud));
	else
		ind = (int)(floor(cloud));
	
	mon = timer->month;
	
	// Select appropriate Beta parameters
	// cerr<<"cloud = "<<cloud<<";\t ind = "<<ind<<endl;
	betaA = EtA[mon-1][ind];
	betaB = EtB[mon-1][ind];
	return;
}

/***************************************************************************
**
**  SimulateWindCurtis
**
**  Estimates wind speed as an indpendent variables. Based on the model of
**  Curtis and Ealgeson [1982] that uses the lag-1 Markov model to transform
**  the simulated normal variables to a skewed distribution. Four parameters
**  are required: wind speed mean, standard deviation, lag-1 correlation
**  coefficicent, and skewness coefficient
**  NOTE: Curtis and Eagleson [1982] considered dependence of the wind speed
**        mean and std on time of the day -- NOT considered here (easy to 
**        implement though).
**
***************************************************************************/
double tHydroMetStoch::SimulateWindCurtis()
{ 
	int skfact = 0;
	double arv;
	double EpsT, windSp;  
	windSpeed_1 = windSpeed;
	
	arv = random_normal(0.0,1.0);
	
	// Curtis [1982]: To counteract the problem of sudden shifts in a generated
	// time series whose variate is skewed and has a high (e.g. > 0.8) lag-1
	// autocorrelation coefficient, restrict the usage of the tail of the N(0,1)
	// that causes the problem. By restricting excursions into the offending tail
	// to absolute values below 2.8, only 0.26 percent of the distribution is 
	// restricted.
	// 1.) If the skew is negative, restrict the negative tail of N(0,1)
	// 2.) If the skew is positive, restrict the positive tail of N(0,1)
	
	while (!skfact) {
		if ((w_skew < 0 && arv <= -2.8) || (w_skew > 0 && arv >= 2.8)) {
			arv = random_normal(0.0,1.0);
			skfact = 0;
		}
		else
			skfact = 1;
	}
	
	// Use Curtis and Eagleson [1982] model for skewed wind speed distribution
	EpsT  = (1+w_skew_gam*arv/6.-pow(w_skew_gam,2.)/36.);
	EpsT  = pow(EpsT,3.)-1.0;
	EpsT *= 2/w_skew_gam;
	
	// Markov lag-1 model accounting for the skewness of wind speed
	windSp =  WindSp_mean + w_ro1*(windSpeed_1-WindSp_mean);
	windSp += (EpsT*WindSp_std*(sqrt(1-w_ro1*w_ro1)));
	
	if (windSp <= 0.0)
		windSp = 0.1;
	return windSp; 
}

void tHydroMetStoch::SetWindVars()
{ 
	w_skew_gam = w_skew*(1-pow(w_ro1,3.))/pow((1-pow(w_ro1,2.)),1.5);
	return;
}

void tHydroMetStoch::SetSkyVars()
{ 
	// Rainfall is generated at the very beginning of simulation
	AvCld = AvCld_1 = 10.0; 
	mt_1    = 1.0;
	tNextRn = timer->getStormTime()   - timer->getPrevStormTime();
	tAftrRn = timer->getCurrentTime() - timer->getPrevStormTime();
	
	if (CldParOpt == 2)
		SetCloudinessVars();
	return;
}

void tHydroMetStoch::SetCloudinessVars()
{ 
	int mon;
	// If monthly parameter values are inputted
	// read them from corresponding arrays
	if (CldParOpt == 2) {
		mon      = timer->month;
		M0       =  Cldi[mon-1][0];
		StdSky   =  Cldi[mon-1][1];
		m_ro1    =  Cldi[mon-1][2];
		GammaSky =  Cldi[mon-1][3];
	}
	return;
}

void tHydroMetStoch::SetTemperatureVars()
{ 
	int mon = timer->month;
	
	// If there are monthly values of regression parameters
	// Read those that correspond to the current month
	if (TempOpt == 1 && TempParOpt == 2) {
		
		// Parameters of temperature deviations
		MeanTDev     = Tdevi[mon-1][0];
		StdTDev      = Tdevi[mon-1][1];
		AutoCorTDev  = Tdevi[mon-1][2];
		
		// Regression parameters
		b0 = Bi[mon-1][0];
		b1 = Bi[mon-1][1];
		b2 = Bi[mon-1][2];
		b3 = Bi[mon-1][3];
		b4 = Bi[mon-1][4];
		b5 = Bi[mon-1][5];
		b6 = Bi[mon-1][6];
		b7 = Bi[mon-1][7];                               
	}
	
	// Initialization of simulation variables 
	if (!timer->getCurrentTime()) {
		airTemp = airTemp_1 = dewTemp = dewTemp_1 = temp11pm;
		I1_1=I2_1=I3_1=I4_1=I5_1=I6_1=I7_1= 0.0;  
		Tdev_1 = Tdev = 0.;
		AvEnergy = AvEnergy_1 = 0.0;
		AvEp = AvEp_1 = 0.0;
		AvCld = AvCld_1 = 0.0;
		EF = 0.;
		Tmax = Tmax_1 = airTemp;
		Tmin = Tmin_1 = airTemp;
	}
	
	// Dew point model initialization
	AnnRain = AnnR[mon-1];
	
	return;
}

/***************************************************************************
**
**  SimulateDewPointTemperature
**
**  Computes dew point temperature based on the time of sunrise 
**  Arguments: LOCAL time and LOCAL Sun rise time
**
***************************************************************************/
double tHydroMetStoch::SimulateDewPointTemperature(double tlocal, 
												   double timeSunR)
{
	double dewT, Trange;
	// If the time is after midnight/ before the sun rise/ any other time
	if ((int)(floor(timeSunR) - floor(tlocal)) != 1) {
		dewT = dewTemp_1;
		if (dewT > airTemp)
			dewT = airTemp;
	}
	else {
		// At the hour preceding the hour of sun rise
		// --> approximate the temperature amplitude
		// with the previous day values
		if (timer->getCurrentTime() > 23.0)
			Trange = Tmax_1-Tmin_1;
		else
			Trange = Tmax-Tmin;
		
		// Call the function that adjusts the dew point 
		// temperature from the daily minimum value 
		dewT = ComputeDailyTdewKimball(airTemp, Trange);
		if (dewT > airTemp)
			dewT = airTemp;
	}
	dewTemp_1 = dewT;
	
	// ##### Special simulation option: no randomness in climate 
	if (etPtr->simCtrl->smooth_weather == 'Y')
		dewT = 12.787;
	
	// Now, convert 'dewT' to other units: 'mb' and relative humidity
	double latHeat;
	double rv = 461.5;
	double eo = 6.112;
	double to = 273.15;
	
	// The Latent Heat of Vaporization as a function of temperature [J/kg]
	latHeat = (597.3 - airTemp*0.57)*4186.8;
	
	// Vapor pressure from dew point temperature [mb]
	vaporPress = eo*exp( (latHeat/rv)*(1/to - 1/(dewT+to)) );
	
	// Compute relative humidity [-] 
	rHumidity = 100*vaporPress/(6.112*exp((17.67*airTemp)/(airTemp+243.5)));
	
	return dewT;
}

/***************************************************************************
**
**  ComputeDailyTdewKimball
**
**  Computes daily dew point temperature based on an empirical model 
**  of Kimball et al. [1997]. Values of EF, T max and Tmin are apprroximated
**  with the previous day values and scaled by a factor accounting for 
**  anticipated daily cloudiness.
**
***************************************************************************/
double tHydroMetStoch::ComputeDailyTdewKimball(double Tmind, double Trange)
{
	double TdewK, aa;
	double ExpClouds, delT1;
	
	// Compute daily EF using Ep from the previous day
	// and day length of the current day 
	ComputeEF(etPtr->getDayLength(), AvEp_1, AnnRain);
	
	// Adjust EF value by using a factor that
	// relates previous and current day cloudness
	if (timer->hour < etPtr->getSunRiseHour())
		delT1 = etPtr->getSunRiseHour() - timer->hour + 1.0;
	else if (timer->hour > etPtr->getSunSetHour())
		delT1 = etPtr->getSunRiseHour() + 24.0 - timer->hour;
	
	ExpClouds = ExpectedCloudiness(1, delT1, etPtr->getDayLength());
	cloudfactor = (1.0-0.65*pow(ExpClouds/10.,2.0))/
		(1.0-0.65*pow(AvCld_1/10.,2.0));
	//cerr<<"\t=> Current time = "<<timer->getCurrentTime()<<endl;
	//cerr<<"\t=> Prev Clouds = "<<AvCld_1<<"\t ExpClouds = "<<ExpClouds
	//    <<"\t Cloud factor = "<<cloudfactor<<endl;
	
	EF *= cloudfactor;
	Trange *= cloudfactor;
	
	// Use an empirical formula of Kimball et al. [1997]
	Tmind += 273.15;
	TdewK = (-0.127+1.121*(1.003-1.444*EF+12.312*EF*EF - 
						   32.766*pow(EF,3.0))+6.0E-4*Trange);
	
	//aa = TdewK - 273.15/(Tmind-273.15)*(1 - TdewK);
	//cerr<<"\t====> CoeffCorr = "<<TdewK<<"\tCorr ToC = "<<aa
	//    <<"\tEF = "<<EF<<"\tTrange = "<<Trange<<endl;
	
	// Leads to unrealistic values otherwise
	if (TdewK < 0.8)
		TdewK = 0.8;
	TdewK *= Tmind;
	
	if (TdewK < 200.0) {
		cout<<"WRONG:"<<endl;
		cout<<"\tPrev Clouds = "<<AvCld_1<<"\t Cloud factor = "<<cloudfactor<<endl;
		cout<<"\tExpClouds = "<<ExpClouds<<"\t Tmind = "<<Tmind<<endl;
		cout<<"\tTrange = "<<Trange<<"\tdelT1 = "<<delT1<<endl;
		cout<<"\tTdewK = "<<TdewK<<"\tEF = "<<EF<<endl;
	}
	
	// Maximum allowed change in dew temperature
	double dTcr = 10.0;
	aa = TdewK - 273.15 - dewTemp_1;
	if (fabs(aa) > dTcr && timer->getCurrentTime() > 23.0)
		(aa > 0.0) ? TdewK -= (aa-dTcr) : TdewK -= (aa+dTcr);
	
	return(TdewK-273.15);
}

/***************************************************************************
**
**  ComputeEF
**
**  Estimates evaporative ratio for the day as a function of Ep rate (drate)
**  day lenght (dlength), and mean annual rainfall
**
***************************************************************************/
void tHydroMetStoch::ComputeEF(double dlength, double drate, double arain)
{
	if (arain)
		EF = drate*dlength/arain;
	else 
		EF = 0.0;
	return;
}

/***************************************************************************
**
**  ComputeDailyAvailEnergy & ComputeDailyEpCld
**
**  Compute daily average values of available energy or potential
**  evaporation. Call is made from tEvapotrans::callEvapoPotential()
**  Available energy      - [W/m^2] 
**  Potential evaporation - mm/hour
**  Cloudiness            - ranges [0,10]
**
***************************************************************************/
void tHydroMetStoch::ComputeDailyAvailEnergy(double rad)
{
	if (timer->hour) {
		AvEnergy *= (timer->hour);
		AvEnergy += rad;
		AvEnergy /= ((timer->hour)+1);
		//cerr<<"AvEnergy = "<<AvEnergy<<"; Rad = "<<rad<<endl;
	}
	else {
		AvEnergy_1 = AvEnergy;
		AvEnergy   = rad;
		//cerr<<"\tAvEnergy = "<<AvEnergy_1<<endl;
	}
	return;
}

void tHydroMetStoch::ComputeDailyEpCld(double rate, double cld)
{
	if (timer->hour) {
		AvEp  *= HrCount;
		AvCld *= HrCount;
		AvEp  += rate;
		AvCld += cld;
		HrCount++;
		AvEp  /= HrCount;
		AvCld /= HrCount;
		//cerr<<"\t\tAvEp = "<<AvEp<<"; AvCld = "<<AvCld
		//	   <<"; Rate = "<<rate<<"; HrCount = "<<HrCount<<endl;
	}
	else {
		// cerr<<"\n\tHHHHH>>Tmax = "<<Tmax<<";\tTmin = "<<Tmin<<endl;
		// cerr<<"\tHHHHH>>AvEp   = "<<AvEp<<";\tAvCld   = "<<AvCld<<endl;
		// if (AvEp_1 > 0 && (Tmax_1-Tmin_1) > 0) 
		// hout<<cloudfactor<<" "
		//     <<AvEp/AvEp_1<<" "
		//     <<(Tmax-Tmin)/(Tmax_1-Tmin_1)<<endl;
		// cerr<<"\tFactor predict = "<<cloudfactor
		//     <<"\tFactor obs Ep = "<<AvEp/AvEp_1
		//     <<"\tFactor obs Tair = "<<(Tmax-Tmin)/(Tmax_1-Tmin_1)<<endl;
		if (HrCount) {
			AvEp_1  = AvEp;
			AvCld_1 = AvCld;
			Tmax_1  = Tmax;
			Tmin_1  = Tmin;
		}
		AvEp    = rate;
		AvCld   = cld;
		Tmax    = -9999;
		Tmin    =  9999;
		HrCount = 0;
		//cerr<<"\tHHHHH>>AvEp_1 = "<<AvEp_1<<";\tAvCld_1 = "<<AvCld_1<<endl;
	}
	return;
}

/***************************************************************************
**
**  ExpectedCloudiness
**
**  Estimates expected cloudness during the time of the day with positive
**  shortwave radiation flux. 'delT1' provides the time difference between
**  the time at which the function call is made and the time of Sun Rise.
**  'delT2' is the day length (shortwave > 0).
**  NOTE: Only mean values are used (non-stationary model is used).
**        Noise is not added.
**
***************************************************************************/
double tHydroMetStoch::ExpectedCloudiness(int opt, double delT1, double delT2)
{
	int cnt = 0;
	int flag = 0;
	int countR = 0;
	double clouds, cloudSUM, Pt;
	double tmpTC, tmpN, tmpA, mt, mtlag;
	
	if (opt == 1) {
		
		// To adjust the evaluation period to the period of
		// positive shortwave radiation 'delT1'
		tmpTC = timer->getCurrentTime();
		tmpTC += delT1;
		mtlag = mt_1;
		mt = 0.0;
		cloudSUM = 0;
		for (int i=0; i<(int)(ceil(delT2)); i++) {
			// If the current hour during previous storm...
			if (tmpTC <= (timer->getStormTime()-rainPtr->interstormDur()) && 
				rainPtr->getRainrate()) {
				clouds = 1.0;
			}
			else {
				// If coming from the PREVIOUS storm/interstorm period...
				if (!flag) {
					// If 'tmpTC' is between the coming storm time [getStormTime()]
					// and end of the previous storm [getPrevStormTime() + getStormDuration()]
					if (tmpTC < timer->getStormTime()) {
						tmpN = timer->getStormTime() - tmpTC;
						tmpA = tmpTC - timer->getPrevStormTime() - 
							rainPtr->getStormDuration();
						mt = EstimateAR1Var(0.0, StdSky, m_ro1, mtlag);
						mtlag= mt;
						// Non-stationary mean cloud cover model
						Pt = (1-exp(-GammaSky*tmpA))*(1-exp(-GammaSky*tmpN));
					}
					// If 'tmpTC' is larger than the time tag of the coming storm
					// --> we don't know the duration of that storm --> approximate
					// with the mean value and then exponentially decay cloudiness (flag=1)
					else if (tmpTC >= timer->getStormTime()) {
						if (countR < rainPtr->getMeanStormDur()-1) {
							clouds = 1.0;
							countR++;
						}
						else {
							clouds = 1.0;
							flag = 1;
						}          
					}
				}
				// If coming from the NEXT storm to the NEXT interstorm period...
				else {
					tmpA = (double)flag;
					Pt = 1 - (1-exp(-GammaSky*tmpA));
					flag++;
				}
				
				// Evaluate cloud cover
				clouds = M0 + (1-M0)*(1-Pt) + mt*Pt;
				
				if (clouds > 1)
					clouds = 1.0;
				else if (clouds < 0) 
					clouds = 0.0;
			}
			cloudSUM += clouds;
			tmpTC++;
		}
		cloudSUM /= ceil(delT2);
		clouds = cloudSUM*10.; // To scale to the values that are used 
	}
	return clouds;
}

/***************************************************************************
**
**  ReadWeatherParameters
**
**  Sets basic variables of the class 
**
***************************************************************************/
void tHydroMetStoch::ReadWeatherParameters(char *WeatherFile)
{
	int i,j;
	int topt;
	int numS;
	double t1, t2;
	char lineIn[300];
	ifstream Inp0(WeatherFile);
	if (!Inp0) {
		cout <<"File "<<WeatherFile<<" not found!!!"<<endl;
		cout <<"Exiting Program...\n\n"<<endl;
		exit(2);
	}
	// Read mean basin latitude and longitude
	Inp0.getline(lineIn, 300);
	if (etPtr->simCtrl->Verbose_label == 'Y') {
		cout<<"\n   Weather generator input parameters:"<< endl;
		cout<<"\n\t"<< lineIn << endl;
	}
	Inp0>>latitude>>longitude;
	Inp0.get();
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<<latitude<<" "<<longitude<<endl;
	
	// Read difference with GMT
	Inp0.getline(lineIn, 300);
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<< lineIn << endl;
	Inp0>>gmt;
	Inp0.get();
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<<gmt<<endl; 
	
	// Read wind speed par-s
	Inp0.getline(lineIn, 200);
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<< lineIn << endl;
	Inp0>>WindSp_mean>>WindSp_std>>w_ro1>>w_skew>>WindSh_fact;
	Inp0.get();
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<<WindSp_mean<<" "<<WindSp_std<<" "
			<<w_ro1<<" "<<w_skew<<" "<<WindSh_fact<<endl;
	
	// Read cloudiness model par-s
	Inp0.getline(lineIn, 200);
	if (etPtr->simCtrl->Verbose_label == 'Y') {
		cout<<"\t"<< lineIn << endl;
		cout<<"\t"<<"Option for input of cloudiness model parameters"<<endl;
	}
	Inp0>>CldParOpt;
	Inp0.get();
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<<CldParOpt<<endl;
	Inp0.getline(lineIn, 200);
	
	if (CldParOpt == 1) {
		Inp0>>M0>>StdSky>>m_ro1>>GammaSky;
		Inp0.get();
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<<M0<<" "<<StdSky<<" "<<m_ro1<<" "<<GammaSky<<endl;
	}
	else if (CldParOpt == 2) {
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t<<<Cloudiness parameters are inputted for each month >>>\n"<<endl;
		Cldi = new double* [12];
		for (i=0; i < 12; i++) {
			Cldi[i] = new double[4];
			assert(Cldi[i] != 0);
			for (j=0; j < 4; j++)
				Inp0>>Cldi[i][j];
			if (etPtr->simCtrl->Verbose_label == 'Y')
				PrintOutArray(Cldi[i], 4);
		}
		Inp0.get();
	}
	
	// Read in parameters of Beta distribution for each month
	Inp0.getline(lineIn, 200);
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<< lineIn << endl;
	Inp0.getline(lineIn, 200);
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<< lineIn << endl;
	for (i=0; i < 12; i++) {
		for (j=0; j < 11; j++) {
			Inp0>>EtA[i][j];
			if (EtA[i][j] <= 0) {
				cerr<<"Wrong input of beta-distribution parameters!"<<endl;
				cerr<<"Parameters must be > 0!\n Exiting..."<<endl;
				exit(2);
			}
		}
		if (etPtr->simCtrl->Verbose_label == 'Y')
			PrintOutArray(EtA[i], 11);
	}
	Inp0.get();
	Inp0.getline(lineIn, 200);
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<< lineIn << endl;
	for (i=0; i < 12; i++) {
		for (j=0; j < 11; j++) {
			Inp0>>EtB[i][j];
			if (EtB[i][j] <= 0) {
				cerr<<"Wrong input of beta-distribution parameters!"<<endl;
				cerr<<"Parameters must be > 0!\n Exiting..."<<endl;
				exit(2);
			}
		}
		if (etPtr->simCtrl->Verbose_label == 'Y')
			PrintOutArray(EtB[i], 11);
	}
	Inp0.get();
	
	
	// Read initial temperature
	Inp0.getline(lineIn, 200);
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<< lineIn << endl;
	Inp0>>temp11pm;
	Inp0.get();
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<<temp11pm<<endl;
	
	
	// Read temperature model par-s
	Inp0.getline(lineIn, 200);
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<<"Option for input of temperature model parameters"<<endl;
	Inp0>>TempParOpt;
	Inp0.get();
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<<TempParOpt<<endl;
	
	if (TempParOpt == 1) {
		Inp0.getline(lineIn, 200);
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<< lineIn << endl;
		Inp0>>MeanTDev>>StdTDev>>AutoCorTDev;
		Inp0.get();
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<<MeanTDev<<" "<<StdTDev<<" "<<AutoCorTDev<<endl;
		
		// Read temperature model regression
		Inp0.getline(lineIn, 200);
		Inp0.getline(lineIn, 200);
		Inp0>>b0>>b1>>b2>>b3>>b4>>b5>>b6>>b7;
		Inp0.get();
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<<b0<<" "<<b1<<" "<<b2<<" "<<b3<<" "
				<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<endl;
	}
	else if (TempParOpt == 2) {
		Inp0.getline(lineIn, 200);
		if (etPtr->simCtrl->Verbose_label == 'Y') {
			cout<<"\t"<< lineIn << endl;
			cout<<"\t<<<Temperature parameters are inputted for each month>>>\n"<<endl;
		}
		Tdevi = new double* [12];
		for (i=0; i < 12; i++) {
			Tdevi[i] = new double[3];
			assert(Tdevi[i] != 0);
			for (j=0; j < 3; j++)
				Inp0>>Tdevi[i][j];
			if (etPtr->simCtrl->Verbose_label == 'Y')
				PrintOutArray(Tdevi[i], 3);
		}
		Inp0.get();  
		
		Inp0.getline(lineIn, 200);
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<< lineIn << endl;
		Inp0.getline(lineIn, 200);
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<< lineIn << endl;
		
		Bi = new double* [12];
		for (i=0; i < 12; i++) {
			Bi[i] = new double[8];
			assert(Bi[i] != 0);
			for (j=0; j < 8; j++)
				Inp0>>Bi[i][j];
			if (etPtr->simCtrl->Verbose_label == 'Y')
				PrintOutArray(Bi[i], 8);
		} 
		Inp0.get();
	}
	Inp0.get();
	
	// Read in parameters for the dew temperature model
	Inp0.getline(lineIn, 200);
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<< lineIn << endl;
	for (i=0; i < 12; i++) {
		Inp0>>AnnR[i];
		if (AnnR[i] <= 0) {
			cerr<<"Wrong input for the dew point temperature model!"<<endl;
			cerr<<"Coefficients must be > 0!\n Exiting..."<<endl;
			exit(2);
		}
	}
	if (etPtr->simCtrl->Verbose_label == 'Y')
		PrintOutArray(AnnR, 12);
	Inp0.get();
	
	
	// Read in parameters for rainfall seasonality
	Inp0.getline(lineIn, 200);
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"\t"<< lineIn << endl;
	Inp0>>numS;
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"# of seasons: "<<numS<<endl; 
	
	// Allocate memory in tStorm
	rainPtr->allocSeasonMemory(numS);
	
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"Season months:";
	for (i=0; i < numS; i++) { 
		Inp0>>t1>>t2;
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<<(int)t1<<"-"<<(int)t2;
		rainPtr->setSeasonMonth(i+1, (int)t1, (int)t2);
	}
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<endl;
	
	
	Inp0.get();
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"[T_int]: ";
	for (i=0; i < numS; i++) { 
		Inp0>>t1;
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<<t1;
		rainPtr->setSeasonSplDurMean(i+1, t1);
	}
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<endl;
	
	
	Inp0.get();
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"[T_stm]: ";
	for (i=0; i < numS; i++) { 
		Inp0>>t1;
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<<t1;
		rainPtr->setSeasonStmDurMean(i+1, t1);
	}
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<endl;
	
	
	Inp0.get();
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<"[Rrate]: ";
	for (i=0; i < numS; i++) { 
		Inp0>>t1;
		if (etPtr->simCtrl->Verbose_label == 'Y')
			cout<<"\t"<<t1;
		rainPtr->setSeasonPMean(i+1, t1);
	}
	if (etPtr->simCtrl->Verbose_label == 'Y')
		cout<<endl;
	
	
	Inp0.get();
	Inp0.close();
	return;
}

/***************************************************************************
**
** Set() and Get() functions to identify location
**
***************************************************************************/
void tHydroMetStoch::setLat(double lat, int option)
{
	if(option == 1)
		latitude = lat;
	else
		basinLat = lat;
}

void tHydroMetStoch::setSunH(double val){
	SunH = val;
}

void tHydroMetStoch::setSinH(double val){
	SinH = val;
}

void tHydroMetStoch::setIo(double val){
	Io = val;
}

void tHydroMetStoch::setIdir(double val){
	Idir = val;
}

void tHydroMetStoch::setIdir_vis(double val){
	Idir_vis = val;
}

void tHydroMetStoch::setIdir_nir(double val){
	Idir_nir = val;
}

void tHydroMetStoch::setIdif(double val){
	Idif = val;
}

void tHydroMetStoch::setIdif_vis(double val){
	Idif_vis = val;
}

void tHydroMetStoch::setIdif_nir(double val){
	Idif_nir = val;
}

double tHydroMetStoch::getLat(int option){
	if(option==1)
		return latitude;
	else 
		return basinLat;
}

void tHydroMetStoch::setLong(double longit, int option){
	if(option == 1)
		longitude = longit;
	else
		basinLong = longit;
}

double tHydroMetStoch::getLong(int option){
	if(option == 1)
		return longitude;
	else
		return basinLong;
}

void tHydroMetStoch::setGmt(int gmtime){
	gmt = gmtime;
}

int tHydroMetStoch::getGmt(){
	return gmt;
}

/***************************************************************************
**
** Set() and Get() Functions for Hydrometeorologic Data
**
***************************************************************************/

int tHydroMetStoch::getYear(){
	return year;
}

int tHydroMetStoch::getMonth(){
	return month;
}

int tHydroMetStoch::getDay(){
	return day;
}

int tHydroMetStoch::getHour(){
	return hour;
}

double tHydroMetStoch::getAirTemp(){
	return airTemp;
}

double tHydroMetStoch::getAirTemp_1(){
	return airTemp_1;
}

double tHydroMetStoch::getDewTemp(){
	return dewTemp;
}

double tHydroMetStoch::getAtmPress(){
	return atmPress;
}

double tHydroMetStoch::getWindSpeed(){
	return windSpeed;
}

double tHydroMetStoch::getSkyCover(){
	return skyCover;
}

double tHydroMetStoch::getSkyCover_1(){
	return skyCover_1;
}

double tHydroMetStoch::getRHumidity(){
	return rHumidity;
}

double tHydroMetStoch::getSunH(){
	return SunH;
}

double tHydroMetStoch::getSinH(){
	return SinH;
}

double tHydroMetStoch::getIo(){
	return Io;
}

double tHydroMetStoch::getIdir(){
	return Idir;
}

double tHydroMetStoch::getIdir_vis(){
	return Idir_vis;
}

double tHydroMetStoch::getIdir_nir(){
	return Idir_nir;
}

double tHydroMetStoch::getIdif(){
	return Idif;
}

double tHydroMetStoch::getIdif_vis(){
	return Idif_vis;
}

double tHydroMetStoch::getIdif_nir(){
	return Idif_nir;
}

double tHydroMetStoch::getVaporPress(){
	return vaporPress;
}

double tHydroMetStoch::getSurfTemp(){
	return 9999.99;
}

double tHydroMetStoch::getNetRad(){
	return 9999.99;
}

double tHydroMetStoch::getPanEvap(){
	return 9999.99;
}

double tHydroMetStoch::getNlm() {return Nlm;}
double tHydroMetStoch::getNci() {return Nci;}
double tHydroMetStoch::getF0()  {return f0;}

// SKY2008Snow from AJR2007 starts here
//##########################################################################
// 
//  The following is the code that estimates air/ dew point temperatures
//  using a mixture of Richardson [1981] model and some modifications
//  to allow subdaily variations, i.e. using a sin() function to fit the
//  diurnal variations between simulated Tmax and Tmin. Time of Tmax and 
//  Tmin is simulated as random variables that defines the phase shift of
//  the sin() function. Day-to-day variations of the sin() are fitted
//  using a linear displacement vector that has increase in values closer
//  to the beginning of the following day (i.e midnight). 
//  NOTE: The model is NOT coupled to rainfall/shortwave/longwave models
//  and was used only for test purposes.
//  I guess I should state that this is an obsolete but still working
//  version usefull when other parts of the tRIBS model are developed. E.g.
//  vegetation, infiltration, routing, etc. NOT to be used for climatology.
// 
//##########################################################################

/***************************************************************************
**
**  SimulateAirDPTemperature
**
**  Simulates diurnal cycles of air (DiurnT) and dew point (DiurnDP) tempe-
**  rature based on autoregressive processes (assumed that certain values
**  for times (t) and (t-1) are known by the time of the function call) 
**
***************************************************************************/
void tHydroMetStoch::SimulateAirDPTemperature()
{
	double gradt0, gradt;
	double minDaily;
	double tempo[24];
	double *ct;

	// Store the values so that the NEXT day variables
	// become CURRENT day variables
	StoreCurrValues();

	// Use AR(1) model for mean daily temperature
	tdnext = EstimateAR1Var(tmone, sgma, T_ro1, td); // HERE *2

	// Estimate vars for a day following the current day 
	// (this will allow to obtain temperature diurnal cycle)
	// Obtain estimate of residuals vector
	EstimateResidualsVector();

	// Approximate expected values for expected tmax and tmin
	ExpectedTmaxTmin();

	// Use multivariate process model to get expected tmax and tmin
	// NOTE: tmaxe & tmine can be conditioned based on rainfall events
	SimulateTmaxTmin(tmaxe, sgma, tmine, sgma);

	// Define the temp. amplitude for the simulation day
	SetTemperatAmplitude(tmax, tmin);

	// Define time shift for temperature peak in sin function
	SetPeakTemperatShift();

	// Define displacement correction 'gradt0' & 'gradt'
	SetDisplaceCorrects(&gradt0, &gradt);

	// Define temperature displacement vector
	ct = new double[24-(int)tShift];
	SetDisplaceVector(ct, int(24-tShift), gradt0, gradt);   //cast changed, ERV
	AddDisplaceVector(ct, int(24-tShift), tempo);          //cast changed, ERV

	// Define the temp. amplitude for the _next_ day
	SetTemperatAmplitude(tmaxnext, tminnext);

	// Compute toC diurnal cycle accounting for peak time shift
	DiurnalCycleTemperature(tempo, &minDaily);

	// Compute Dew Point diurnal cycle
	DiurnalCycleDewPointTemperature(&minDaily);

	delete [] ct;

	return;
}

/***************************************************************************
**
**  EstimateResidualsVector
**
**  The function estimates a vector of residuals for multivariate auto-
**  regressive process
**
***************************************************************************/

void tHydroMetStoch::EstimateResidualsVector()
{
	double tmp1[3];
	double tmp2[3];
	double randv[3];
	GetVectRandomNorm(randv, 3, 0.0, 1.0);
	MatrixVectorProduct(A, 3, 3, hi_1, tmp1);
	MatrixVectorProduct(B, 3, 3, randv, tmp2);
	VectorSummation(tmp1, tmp2, 3, hi);
	return;
}

/***************************************************************************
**
**  ExpectedTMaxTMin
**
**  Estimates expected values of tmax and tmin. Can be conditioned by
**  precipitiation events occurrence
**
***************************************************************************/

void tHydroMetStoch::ExpectedTmaxTmin()
{ 
	tmaxe = tdnext + sgma;
	tmine = tdnext - sgma;
	return;
}

/***************************************************************************
**
**  SimulateTmaxTmin
**
**  Using estimates of std and mean for both variables, obtains an AR(1)
**  multivariate estimate of tmax and tmin (residuals vector hi must be 
**  computed by the time of the function call)
**
***************************************************************************/

void tHydroMetStoch::SimulateTmaxTmin(double maxe, double stdmax,
				      double mine, double stdmin)
{ 
	tmaxnext = hi[0]*stdmax + maxe;
	tminnext = hi[1]*stdmin + mine;
	return;
}

/***************************************************************************
**
**  SetDisplaceCorrects
**
**  Defines the absolute temperature displacements: beginning & end of day
**
***************************************************************************/

void tHydroMetStoch::SetDisplaceCorrects(double *gr0, double *gr1)
{
	// Define temperature gradient and vertical
	// displacement corrections to allow smooth transition
	if (abs(tShift - tShift_1) > 1E-3) 
		*gr0 = deltt/2*(sin(2*M_PI*(1-(tShift_1+1)/24)) - 
			sin(2*M_PI*(1-(tShift+1)/24)));
	else 
		*gr0 = 0;

	*gr1 = tdnext - td;
	return;
}

/***************************************************************************
**
**  SetDisplaceVector
**
**  Defines a vector of temperature displacements
**
***************************************************************************/

void tHydroMetStoch::SetDisplaceVector(double *ct, int n,
				       double gr0, double gr1)
{
	assert(n>0);
	for (int i = 0; i<n; i++)
  	ct[i] = (double)(i)*(gr1-gr0)/(n-1) + gr0;
	return;
}

/***************************************************************************
**
**  AddDisplaceVector
**
**  Adds displacement vector 'ct' to vector 'tempo' that makes first 
**  approximation of the diurnal cycle using sine function
**
***************************************************************************/

void tHydroMetStoch::AddDisplaceVector(double *ct, int n, double *tempo)
{
	assert(n>0);
	for (int i = 0; i<24; i++) {
		tempo[i] = (td + deltt/2.*sin(2*M_PI*(1-(i+1)/24.)));
		if (i >= tShift) 
			tempo[i] += ct[i-(int)tShift];
	}
	return;
}

/***************************************************************************
**
**  DiurnalCycleTemperature
**
**  Computes diurnal cycle of air temperature based on previously simulated
**  values of tmax and tmin as well peak time shift
**
***************************************************************************/

void tHydroMetStoch::DiurnalCycleTemperature(double *tempo, double *minDaily)
{
	*minDaily = 99999.;
	int k = 1;
	for (int i = 0; i<24; i++) {
		if (i <= 24-tShift-1) {
			// Use cycle fraction computed for that day 
			DiurnT[i] = tempo[i+(int)tShift];
		}
		else {
			// Use temperatures which you would compute for the next day
			// NOTE: deltt must be computed using tmaxnext and tminnext!
			DiurnT[i] = tdnext + deltt/2.*sin(2*M_PI*(1-k/24.));
			k = k+1;
         
			if (DiurnT[i] < (*minDaily)) 
				(*minDaily) = DiurnT[i];
     		}
  	}
  	return;
}

/***************************************************************************
**
**  DiurnalCycleDewPointTemperature
**
**  Computes diurnal cycle of dew point temperature based on previously 
**  simulated diurnal cycle of air temperature
**
***************************************************************************/

void tHydroMetStoch::DiurnalCycleDewPointTemperature(double *minDaily)
{
	for (int i = 0; i<24; i++) {
		DiurnDP[i] = (*minDaily);
		// Add some noise here..
		// Condition on rain: if Rain --> DiurnDP(l) = DiurnT(l);
		if (DiurnDP[i] > DiurnT[i])
			DiurnDP[i] = DiurnT[i];
	}
	return;
}

/***************************************************************************
**
**  SetPeakTemperatShift
**
**  Randomly chooses the phase shift for sin peak between 2 and 4 hours 
**  Assumed uniform distribution
**
***************************************************************************/

void tHydroMetStoch::SetPeakTemperatShift()
{
	// Define time shift for temperature peak for the 
	// simulation day: 2 to 4 hours
	tShift = uniform(2, 4);
	if (tShift - floor(tShift) > 0.5) 
		tShift = ceil(tShift);
	else 
		tShift = floor(tShift);
	return;
}

/***************************************************************************
**
**  SetTemperatAmplitude
**
**  Defines the temperature amplitude
**
***************************************************************************/

void tHydroMetStoch::SetTemperatAmplitude(double max, double min)
{
	deltt = max - min;
	if (deltt < 1) 
		deltt = deltt_1;
	return;
}

/***************************************************************************
**
**  SetTemperatParametersSetTemperatParameters
**
**  Defines monthly mean and std for air temperature
**
***************************************************************************/

void tHydroMetStoch::SetTemperatParameters()
{
	sgma = Stdmon[timer->month-1]; //HERE /2
	tmone = random_normal(0.0,1.0)*Stdmon[timer->month-1]+Tmon[timer->month-1];
	//cout<<"\tSimulated monthly toC and std values:"<<endl;
	//cout<<"\tTMON = "<<tmone<<"\tSGMA = "<<sgma<<endl;
	return;
}

/***************************************************************************
**
**  StoreCurrValues
**
**  Stores values for time t
**
***************************************************************************/

void tHydroMetStoch::StoreCurrValues()
{ 
	// Store as the current values
	hi_1[0] = hi[0];
	hi_1[1] = hi[1];

	td = tdnext;      // Will be the current day values next
	tmax = tmaxnext;
	tmin = tminnext;
	deltt_1 = deltt;
	tShift_1 = tShift; 
	return;
}

/***************************************************************************
**
**  GetVectRandomNorm
**
**  Simulates normally distributed values with par-rs N(m, std^2) and
**  fills the vector a of length n
**
***************************************************************************/

void tHydroMetStoch::GetVectRandomNorm(double *a, int n, double m, double std)
{ 
	for (int i=0; i<n; i++)
		a[i] = random_normal(m,std);
	return;
}
// SKY2008Snow from AJR2007 ends here

void tHydroMetStoch::PrintOutArray(double *a, int n)
{
	for (int i = 0; i<n; i++)  
		cout<<a[i]<<" ";
	cout<<endl;
}

/***************************************************************************
**
** tHydroMetStoch::writeRestart() Function
** 
** Called from tSimulator during simulation loop **
***************************************************************************/
                                                                                 void tHydroMetStoch::writeRestart(fstream & rStr) const                                     {
  BinaryWrite(rStr, gmt);
  BinaryWrite(rStr, latitude);
  BinaryWrite(rStr, longitude);
  BinaryWrite(rStr, basinLat);
  BinaryWrite(rStr, basinLong);
  for (int i = 0; i < 12; i++) {     BinaryWrite(rStr, Tmon[i]);
    BinaryWrite(rStr, Stdmon[i]);   }
  BinaryWrite(rStr, T_ro1);
  BinaryWrite(rStr, tmone);
  BinaryWrite(rStr, sgma);
  BinaryWrite(rStr, tmaxe);
  BinaryWrite(rStr, tmine);
  BinaryWrite(rStr, td);
  BinaryWrite(rStr, tmax);
  BinaryWrite(rStr, tmin);
  BinaryWrite(rStr, tdnext);
  BinaryWrite(rStr, tmaxnext);
  BinaryWrite(rStr, tminnext);
  BinaryWrite(rStr, deltt);
  BinaryWrite(rStr, deltt_1);
  BinaryWrite(rStr, tShift);
  BinaryWrite(rStr, tShift_1);
  for (int i = 0; i < 3; i++) {
    BinaryWrite(rStr, hi[i]);
    BinaryWrite(rStr, hi_1[i]);
  }
  for (int i = 0; i < 24; i++) {
     BinaryWrite(rStr, DiurnT[i]);
     BinaryWrite(rStr, DiurnDP[i]);
  }
  BinaryWrite(rStr, SunH);
  BinaryWrite(rStr, SinH);;
  BinaryWrite(rStr, Io);
  BinaryWrite(rStr, Idir);
  BinaryWrite(rStr, Idir_vis);
  BinaryWrite(rStr, Idir_nir);
  BinaryWrite(rStr, Idif);
  BinaryWrite(rStr, Idif_vis);
  BinaryWrite(rStr, Idif_nir);
  BinaryWrite(rStr, atmPress);
  BinaryWrite(rStr, Patm_mean);
  BinaryWrite(rStr, Patm_std);
  BinaryWrite(rStr, skyCover);
  BinaryWrite(rStr, skyCover_1);
  BinaryWrite(rStr, M0);
  BinaryWrite(rStr, StdSky);
  BinaryWrite(rStr, GammaSky);
  BinaryWrite(rStr, m_ro1);
  BinaryWrite(rStr, mt_1);
  BinaryWrite(rStr, tAftrRn);
  BinaryWrite(rStr, tNextRn);
  BinaryWrite(rStr, Nlm);
  BinaryWrite(rStr, Nci);
  BinaryWrite(rStr, f0);
  BinaryWrite(rStr, windSpeed);
  BinaryWrite(rStr, windSpeed_1);
  BinaryWrite(rStr, WindSp_mean);
  BinaryWrite(rStr, WindSp_std);
  BinaryWrite(rStr, WindSh_fact);
  BinaryWrite(rStr, w_ro1);
  BinaryWrite(rStr, w_skew);
  BinaryWrite(rStr, w_skew_gam);
  BinaryWrite(rStr, airTemp);
  BinaryWrite(rStr, airTemp_1);
  BinaryWrite(rStr, dewTemp);
  BinaryWrite(rStr, dewTemp_1);
  BinaryWrite(rStr, TempOpt);
  BinaryWrite(rStr, CldParOpt);
  BinaryWrite(rStr, TempParOpt);
  BinaryWrite(rStr, RainOpt);
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 4; j++)
      BinaryWrite(rStr, Cldi[i][j]);
    for (int j = 0; j < 3; j++)
      BinaryWrite(rStr, Tdevi[i][j]);
    for (int j = 0; j < 8; j++)
      BinaryWrite(rStr, Bi[i][j]);
    for (int j = 0; j < 11; j++) {
      BinaryWrite(rStr, EtA[i][j]);
      BinaryWrite(rStr, EtB[i][j]);
    }
    BinaryWrite(rStr, AnnR[i]);
  }
  BinaryWrite(rStr, betaA);
  BinaryWrite(rStr, betaB);
  BinaryWrite(rStr, K0);
  BinaryWrite(rStr, K1);
  BinaryWrite(rStr, K2);
  BinaryWrite(rStr, K3);
  BinaryWrite(rStr, K4);
  BinaryWrite(rStr, K5);
  BinaryWrite(rStr, K6);
  BinaryWrite(rStr, I1_1);
  BinaryWrite(rStr, I2_1);
  BinaryWrite(rStr, I3_1);
  BinaryWrite(rStr, I4_1);
  BinaryWrite(rStr, I5_1);
  BinaryWrite(rStr, I6_1);
  BinaryWrite(rStr, I7_1);
  BinaryWrite(rStr, temp11pm);
  BinaryWrite(rStr, qt_1);
  BinaryWrite(rStr, MeanTDev);
  BinaryWrite(rStr, StdTDev);
  BinaryWrite(rStr, AutoCorTDev);
  BinaryWrite(rStr, Tdev);
  BinaryWrite(rStr, Tdev_1);
  BinaryWrite(rStr, AnnRain);
  BinaryWrite(rStr, EF);
  BinaryWrite(rStr, AvEnergy);
  BinaryWrite(rStr, AvEnergy_1);
  BinaryWrite(rStr, AvEp);
  BinaryWrite(rStr, AvEp_1);
  BinaryWrite(rStr, AvCld);
  BinaryWrite(rStr, AvCld_1);
  BinaryWrite(rStr, Tmax);
  BinaryWrite(rStr, Tmin);
  BinaryWrite(rStr, Tmax_1);
  BinaryWrite(rStr, Tmin_1);
  BinaryWrite(rStr, cloudfactor);
  BinaryWrite(rStr, vaporPress);
  BinaryWrite(rStr, rHumidity);
  BinaryWrite(rStr, year);
  BinaryWrite(rStr, month);
  BinaryWrite(rStr, day);
  BinaryWrite(rStr, hour);
  BinaryWrite(rStr, minute);
  BinaryWrite(rStr, HrCount);

  BinaryWrite(rStr, Xm1);
  BinaryWrite(rStr, Xm2);
  BinaryWrite(rStr, Xa1);
  BinaryWrite(rStr, Xa2);
  for (int i = 0; i < 32; i++) {
    BinaryWrite(rStr, Xcg1[i]);
    BinaryWrite(rStr, Xcg2[i]);
    BinaryWrite(rStr, Xig1[i]);
    BinaryWrite(rStr, Xig2[i]);
    BinaryWrite(rStr, Xlg1[i]);
    BinaryWrite(rStr, Xlg2[i]);
    BinaryWrite(rStr, Xqanti[i]);
  }
  BinaryWrite(rStr, Xa1w);
  BinaryWrite(rStr, Xa2w);
  BinaryWrite(rStr, Xa1vw);
  BinaryWrite(rStr, Xa2vw);
}

/***************************************************************************
**
** HydroMetStoch::readRestart() Function
**
***************************************************************************/

void tHydroMetStoch::readRestart(fstream & rStr)
{
  BinaryRead(rStr, gmt);
  BinaryRead(rStr, latitude);
  BinaryRead(rStr, longitude);
  BinaryRead(rStr, basinLat);
  BinaryRead(rStr, basinLong);
  for (int i = 0; i < 12; i++) {
    BinaryRead(rStr, Tmon[i]);
    BinaryRead(rStr, Stdmon[i]);
  }
  BinaryRead(rStr, T_ro1);
  BinaryRead(rStr, tmone);
  BinaryRead(rStr, sgma);
  BinaryRead(rStr, tmaxe);
  BinaryRead(rStr, tmine);
  BinaryRead(rStr, td);
  BinaryRead(rStr, tmax);
  BinaryRead(rStr, tmin);
  BinaryRead(rStr, tdnext);
  BinaryRead(rStr, tmaxnext);
  BinaryRead(rStr, tminnext);
  BinaryRead(rStr, deltt);
  BinaryRead(rStr, deltt_1);
  BinaryRead(rStr, tShift);
  BinaryRead(rStr, tShift_1);
  for (int i = 0; i < 3; i++) {
    BinaryRead(rStr, hi[i]);
    BinaryRead(rStr, hi_1[i]);
  }
  for (int i = 0; i < 24; i++) {
     BinaryRead(rStr, DiurnT[i]);
     BinaryRead(rStr, DiurnDP[i]);
  }
  BinaryRead(rStr, SunH);
  BinaryRead(rStr, SinH);;
  BinaryRead(rStr, Io);
  BinaryRead(rStr, Idir);
  BinaryRead(rStr, Idir_vis);
  BinaryRead(rStr, Idir_nir);
  BinaryRead(rStr, Idif);
  BinaryRead(rStr, Idif_vis);
  BinaryRead(rStr, Idif_nir);
  BinaryRead(rStr, atmPress);
  BinaryRead(rStr, Patm_mean);
  BinaryRead(rStr, Patm_std);
  BinaryRead(rStr, skyCover);
  BinaryRead(rStr, skyCover_1);
  BinaryRead(rStr, M0);
  BinaryRead(rStr, StdSky);
  BinaryRead(rStr, GammaSky);
  BinaryRead(rStr, m_ro1);
  BinaryRead(rStr, mt_1);
  BinaryRead(rStr, tAftrRn);
  BinaryRead(rStr, tNextRn);
  BinaryRead(rStr, Nlm);
  BinaryRead(rStr, Nci);
  BinaryRead(rStr, f0);
  BinaryRead(rStr, windSpeed);
  BinaryRead(rStr, windSpeed_1);
  BinaryRead(rStr, WindSp_mean);
  BinaryRead(rStr, WindSp_std);
  BinaryRead(rStr, WindSh_fact);
  BinaryRead(rStr, w_ro1);
  BinaryRead(rStr, w_skew);
  BinaryRead(rStr, w_skew_gam);
  BinaryRead(rStr, airTemp);
  BinaryRead(rStr, airTemp_1);
  BinaryRead(rStr, dewTemp);
  BinaryRead(rStr, dewTemp_1);
  BinaryRead(rStr, TempOpt);
  BinaryRead(rStr, CldParOpt);
  BinaryRead(rStr, TempParOpt);
  BinaryRead(rStr, RainOpt);
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 4; j++)
      BinaryRead(rStr, Cldi[i][j]);
    for (int j = 0; j < 3; j++)
      BinaryRead(rStr, Tdevi[i][j]);
    for (int j = 0; j < 8; j++)
      BinaryRead(rStr, Bi[i][j]);
    for (int j = 0; j < 11; j++) {
      BinaryRead(rStr, EtA[i][j]);
      BinaryRead(rStr, EtB[i][j]);
    }
    BinaryRead(rStr, AnnR[i]);
  }
  BinaryRead(rStr, betaA);
  BinaryRead(rStr, betaB);
  BinaryRead(rStr, K0);
  BinaryRead(rStr, K1);
  BinaryRead(rStr, K2);
  BinaryRead(rStr, K3);
  BinaryRead(rStr, K4);
  BinaryRead(rStr, K5);
  BinaryRead(rStr, K6);
  BinaryRead(rStr, I1_1);
  BinaryRead(rStr, I2_1);
  BinaryRead(rStr, I3_1);
  BinaryRead(rStr, I4_1);
  BinaryRead(rStr, I5_1);
  BinaryRead(rStr, I6_1);
  BinaryRead(rStr, I7_1);
  BinaryRead(rStr, temp11pm);
  BinaryRead(rStr, qt_1);
  BinaryRead(rStr, MeanTDev);
  BinaryRead(rStr, StdTDev);
  BinaryRead(rStr, AutoCorTDev);
  BinaryRead(rStr, Tdev);
  BinaryRead(rStr, Tdev_1);
  BinaryRead(rStr, AnnRain);
  BinaryRead(rStr, EF);
  BinaryRead(rStr, AvEnergy);
  BinaryRead(rStr, AvEnergy_1);
  BinaryRead(rStr, AvEp);
  BinaryRead(rStr, AvEp_1);
  BinaryRead(rStr, AvCld);
  BinaryRead(rStr, AvCld_1);
  BinaryRead(rStr, Tmax);
  BinaryRead(rStr, Tmin);
  BinaryRead(rStr, Tmax_1);
  BinaryRead(rStr, Tmin_1);
  BinaryRead(rStr, cloudfactor);
  BinaryRead(rStr, vaporPress);
  BinaryRead(rStr, rHumidity);
  BinaryRead(rStr, year);
  BinaryRead(rStr, month);
  BinaryRead(rStr, day);
  BinaryRead(rStr, hour);
  BinaryRead(rStr, minute);
  BinaryRead(rStr, HrCount);

  BinaryRead(rStr, Xm1);
  BinaryRead(rStr, Xm2);
  BinaryRead(rStr, Xa1);
  BinaryRead(rStr, Xa2);
  for (int i = 0; i < 32; i++) {
    BinaryRead(rStr, Xcg1[i]);
    BinaryRead(rStr, Xcg2[i]);
    BinaryRead(rStr, Xig1[i]);
    BinaryRead(rStr, Xig2[i]);
    BinaryRead(rStr, Xlg1[i]);
    BinaryRead(rStr, Xlg2[i]);
    BinaryRead(rStr, Xqanti[i]);
  }
  BinaryRead(rStr, Xa1w);
  BinaryRead(rStr, Xa2w);
  BinaryRead(rStr, Xa1vw);
  BinaryRead(rStr, Xa2vw);
}

//=========================================================================
//
//
//                         End of tHydroMetStoch.cpp
//
//
//=========================================================================
