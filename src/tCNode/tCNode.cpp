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
**  tCNode.cpp: Functions for derived class tCNode (see tCNode.h)
**
***************************************************************************/

#include "src/tCNode/tCNode.h"

//=========================================================================
//
//
//                  Section 1: tCNode Constructors/Destructors
//
//
//=========================================================================

// Default Constructor

tCNode::tCNode() :tNode() 
{
	alpha = 0.0; 
	NwtOld = NwtNew = 0.0;
	MuOld = MuNew = 0.0; 
	MiOld = MiNew = 0.0;
	NtOld = NtNew = 0.0;
	NfOld = NfNew = 0.0;
	RuOld = RuNew = 0.0;
	RiOld = RiNew = 0.0;
	Qin = Qout = Qpin = Qpout = 0.0;
	Rain = 0.0;
	srf_hr=srf=hsrf=esrf=psrf=satsrf=rsrf=sbsrf=cumsrf=0.0; //added cumsrf CJC2021
	flowedge = 0;
	gwchange = 0.0;
	traveltime = hillpath = streampath = 0.0;
	intstorm = 0.0; 
	Interception=CumIntercept = 0.0;
	NetPrecipitation = 0.0;
	// CanStorage = CanopyStorage = 0.0; // SKYnGM2008LU
	CanStorage = CanopyStorVol = 0.0; // SKYnGM2008LU 
	UnSaturatedStorage = SaturatedStorage = 0.0;
	Recharge = UnSatFlowIn = UnSatFlowOut = 0.0;
	PotEvaporation=ActEvaporation=EvapWetCanopy = 0.0;
	EvapDryCanopy=EvapSoil=EvapoTranspiration = 0.0;
	SoilMoisture = SoilMoistureSC = SoilMoistureUNSC = 0.0;
	RootMoisture = RootMoistureSC = 0.0;
	StormLength = 0;
	nVerts = -999;
	VertsX = 0;
	VertsY = 0;
	bndEdge1 = 0;
	bndEdge2 = 0;
	ContrArea = Curvature = BasinArea = 0.0;
	BedrockDepth = 0.0;
	Transmissivity = 0.0;
	AirTemp = DewTemp = SurfTemp = NetRad = SoilTemp = 0.0;
	VapPress = RelHumid = WindSpeed = SkyCover = AirPressure = 0.0;
	GridET = LongRadIn = LongRadOut = 0.0;
	ShortRadIn = ShortRadIn_dir = ShortRadIn_dif = 0.0;
	ShortAbsbVeg = ShortAbsbSoi = ShortRadSlope = 0.0;
	gFlux = hFlux = lFlux = Gnod = 0.0;
	QgwIn = QgwOut = 0.0;
	Hlevel = Qstrm = Width = Roughness = RunOn = 0.0;
	StreamPtr = 0;
	TimeInd = 0;
	Qeff = 0;
	FlowVelocity = 0.0;
	tracer = flood = soiID = LandUse = satOccur = 0;
   Reach = -1;
	hsrfOccur=psrfOccur=satsrfOccur=sbsrfOccur=0.0;
	percOccur=0.0; //ASM perc in 00i
	avPerc=0.0; //ASM perc in 00i
	RechDisch = Aspect = VegFraction = 0.0;
	VegFraction = 0.0; 
	AvSoilMoisture = AvEvapFract = AvET = 0.0;
	
	xC=yC=-1;

	// ASM 2/10/2017
	ChannelPerc = 0.0;
	Ft = 0.0;

	// SKYnGM2008LU
	//added for land use grid AJR 2007
	LandUseAlb = ThroughFall = VegHeight = StomRes = 0.0; 

	// SKYnGM2008LU
	IntercepCoeff = CanFieldCap = DrainCoeff = 0.0;  
	DrainExpPar = OptTransmCoeff = LeafAI = 0.0;	
	CanStorParam = 0.0;
  EvapThresh = TransThresh = 0.0; // CJC2025

	// SKYnGM2008LU
	LandUseAlbInPrevGrid = LandUseAlbInUntilGrid = ThroughFallInPrevGrid = ThroughFallInUntilGrid = 0.0;
	VegHeightInPrevGrid = VegHeightInUntilGrid = StomResInPrevGrid = StomResInUntilGrid = 0.0;
	VegFractionInPrevGrid = VegFractionInUntilGrid = CanStorParamInPrevGrid = CanStorParamInUntilGrid = 0.0;
	IntercepCoeffInPrevGrid = IntercepCoeffInUntilGrid = CanFieldCapInPrevGrid = CanFieldCapInUntilGrid = 0.0;
	DrainCoeffInPrevGrid = DrainCoeffInUntilGrid = DrainExpParInPrevGrid = DrainExpParInUntilGrid = 0.0;
	OptTransmCoeffInPrevGrid = OptTransmCoeffInUntilGrid = LeafAIInPrevGrid = LeafAIInUntilGrid = 0.0;
  EvapThreshInPrevGrid = EvapThreshInUntilGrid = TransThreshInPrevGrid = TransThreshInUntilGrid = 0.0; // CJC2025

	// SKYnGM2008LU
	AvCanStorParam = AvIntercepCoeff = AvThroughFall = AvCanFieldCap = 0.0;
	AvDrainCoeff = AvDrainExpPar = AvLandUseAlb = AvVegHeight = 0.0;
	AvOptTransmCoeff = AvStomRes = AvVegFraction = AvLeafAI = 0.0;
  AvEvapThresh = AvTransThresh = 0.0; // CJC2025
	
	// SKY2008Snow from AJR2007
	//snowpack -- RINEHART 2007 @ NMT
	liqWEq = iceWEq =  liqRoute = dU = 0.0;
	snTemperC = 0.0;
	crAge = densAge = ETage = 0.0;
	cumLHF = cumMelt = 0.0;
	snSub = snEvap = 0.0; // Initialize snowpack cumulative sublimation & evaporation CJC2020
	cumSnSub = cumSnEvap = cumTotEvap = cumBarEvap = 0.0; // Initialize snowpack cumulative sublimation & evaporation CJC2020
	snLHF = snSHF = snGHF = snPHF = snRLin = snRLout = snRSin = 0.0;
	Unode = Uerror = 0.0;
	cumUError = 0.0;
	cumHrsSun = cumHrsSnow = 0.0;
	cumLHF = cumSHF = cumPHF = cumRLin = cumRLout = cumRSin = cumGHF = cumMelt = 0.0;  
	persTime = persTimeTemp = peakSWE = peakSWEtemp = initPackTime = peakPackTime = 0.0;
	initPackTimeTemp = 0.0;  
	//snowintercept -- RINEHART 2007 @ NMT
	intSWEq = 0.0;
	intSnUnload = 0.0;
	intSub = intPrec = 0.0;
	cumIntSub = cumIntUnl = cumHrsSun = 0.0;  
	//shelter members -- RINEHART 2007 @ NMT
	horizonAngle0000 = 0.0;
	horizonAngle0225 = 0.0;
	horizonAngle0450 = 0.0;
	horizonAngle0675 = 0.0;
	horizonAngle0900 = 0.0;
	horizonAngle1125 = 0.0;
	horizonAngle1350 = 0.0;
	horizonAngle1575 = 0.0;
	horizonAngle1800 = 0.0;
	horizonAngle2025 = 0.0;
	horizonAngle2250 = 0.0;
	horizonAngle2475 = 0.0;
	horizonAngle2700 = 0.0;
	horizonAngle2925 = 0.0;
	horizonAngle3150 = 0.0;
	horizonAngle3375 = 0.0; 
	sfact = 0.0;
	lfact = 0.0;

    // updates to beta_func
    soil_cutoff = 0.0;
    root_cutoff = 0.0;

}

// Input File Constructor 

tCNode::tCNode(tInputFile &infile) :tNode() {
	alpha = 0.0; 
	NwtOld = NwtNew = 0.0;
	MuOld = MuNew = 0.0; 
	MiOld = MiNew = 0.0;
	NtOld = NtNew = 0.0;
	NfOld = NfNew = 0.0;
	RuOld = RuNew = 0.0;
	RiOld = RiNew = 0.0;
	Qin = Qout = Qpin = Qpout = 0.0;
	Rain = 0.0;
	srf_hr=srf=hsrf=esrf=psrf=satsrf=rsrf=sbsrf=cumsrf=0.0; //added cumsrf CJC2021
	flowedge = 0;
	gwchange = 0.0;
	traveltime = hillpath = streampath = 0.0;
	intstorm = 0.0; 
	Interception=CumIntercept = 0.0;
	NetPrecipitation = 0.0;
	// CanStorage = CanopyStorage = 0.0; // SKYnGM2008LU
	CanStorage = CanopyStorVol = 0.0; // SKYnGM2008LU
	UnSaturatedStorage = SaturatedStorage = 0.0;
	Recharge = UnSatFlowIn = UnSatFlowOut = 0.0;
	PotEvaporation=ActEvaporation=EvapWetCanopy = 0.0;
	EvapDryCanopy=EvapSoil=EvapoTranspiration = 0.0;
	SoilMoisture = SoilMoistureSC = SoilMoistureUNSC = 0.0;
	RootMoisture = RootMoistureSC = 0.0;
	StormLength = 0;
	nVerts = -999;
	VertsX = 0;
	VertsY = 0;
	bndEdge1 = 0;
	bndEdge2 = 0;
	ContrArea = Curvature = BasinArea = 0.0;
	BedrockDepth = 0.0;
	Transmissivity = 0.0;
	AirTemp = DewTemp = SurfTemp = NetRad = SoilTemp = 0.0;
	VapPress = RelHumid = WindSpeed = SkyCover = AirPressure = 0.0;
	GridET = ShortRadIn = LongRadIn = LongRadOut = 0.0;
	ShortRadIn = ShortRadIn_dir = ShortRadIn_dif = 0.0;
	ShortAbsbVeg = ShortAbsbSoi = ShortRadSlope = 0.0;
	gFlux = hFlux = lFlux = Gnod = 0.0;
	QgwIn = QgwOut = 0.0;
	Hlevel = Qstrm = Width = Roughness = RunOn = 0.0;
	StreamPtr = 0;
	TimeInd = 0;
	Qeff = 0;
	FlowVelocity = 0.0;
	tracer = flood = soiID = LandUse = satOccur = 0;
   Reach = -1;
	hsrfOccur=psrfOccur=satsrfOccur=sbsrfOccur=0.0;
	percOccur=0.0; //ASM perc in 00i
	avPerc=0.0; //ASM perc in 00i
	RechDisch = Aspect = VegFraction = 0.0;
	VegFraction = 0.0; 
	AvSoilMoisture = AvEvapFract = AvET = 0.0;
	
	xC=yC=-1;

	// ASM 2/10/2017
	ChannelPerc = 0.0;
	Ft = 0.0;

	// SKYnGM2008LU
	//added for land use grid AJR 2007
	LandUseAlb = ThroughFall = VegHeight = StomRes = 0.0; 

	// SKYnGM2008LU
	IntercepCoeff = CanFieldCap = DrainCoeff = 0.0; 
	DrainExpPar = OptTransmCoeff = LeafAI = 0.0;
        CanStorParam = 0.0;	
  EvapThresh = TransThresh = 0.0; // CJC2025

	// SKYnGM2008LU
	LandUseAlbInPrevGrid = LandUseAlbInUntilGrid = ThroughFallInPrevGrid = ThroughFallInUntilGrid = 0.0;
	VegHeightInPrevGrid = VegHeightInUntilGrid = StomResInPrevGrid = StomResInUntilGrid = 0.0;
	VegFractionInPrevGrid = VegFractionInUntilGrid = CanStorParamInPrevGrid = CanStorParamInUntilGrid = 0.0;
	IntercepCoeffInPrevGrid = IntercepCoeffInUntilGrid = CanFieldCapInPrevGrid = CanFieldCapInUntilGrid = 0.0;
	DrainCoeffInPrevGrid = DrainCoeffInUntilGrid = DrainExpParInPrevGrid = DrainExpParInUntilGrid = 0.0;
	OptTransmCoeffInPrevGrid = OptTransmCoeffInUntilGrid = LeafAIInPrevGrid = LeafAIInUntilGrid = 0.0;
  EvapThreshInPrevGrid = EvapThreshInUntilGrid = TransThreshInPrevGrid = TransThreshInUntilGrid = 0.0; // CJC2025

	// SKYnGM2008LU
	AvCanStorParam = AvIntercepCoeff = AvThroughFall = AvCanFieldCap = 0.0;
	AvDrainCoeff = AvDrainExpPar = AvLandUseAlb = AvVegHeight = 0.0;
	AvOptTransmCoeff = AvStomRes = AvVegFraction = AvLeafAI = 0.0;
  AvEvapThresh = AvTransThresh = 0.0; // CJC2025

	// SKY2008Snow from AJR2007
	//snowpack -- RINEHART 2007 @ NMT
	liqWEq = iceWEq =  liqRoute = dU = 0.0;
	snTemperC = 0.0;
	crAge = densAge = ETage = 0.0;
	cumLHF = cumMelt = 0.0;
	snSub = snEvap = 0.0; // Initialize snowpack cumulative sublimation & evaporation CJC2020
	cumSnSub = cumSnEvap = cumTotEvap = cumBarEvap = 0.0; // Initialize snowpack cumulative sublimation & evaporation CJC2020
	snLHF = snSHF = snGHF = snPHF = snRLin = snRLout = snRSin = 0.0;
	Unode = Uerror = 0.0;
	cumUError = 0.0;
	cumLHF = cumSHF = cumPHF = cumRLin = cumRLout = cumRSin = cumGHF = 0.0;
	cumMelt = 0.0;
	cumHrsSun = cumHrsSnow = 0.0;
	persTime = persTimeTemp = peakSWE = peakSWEtemp = initPackTime = peakPackTime = 0.0;
	initPackTimeTemp = 0.0;
	//snowintercept -- RINEHART 2007 @ NMT
	intSWEq = 0.0;
	intSnUnload = 0.0;
	intSub = intPrec = 0.0;
	cumIntSub = cumIntUnl = 0.0;
	//shelter -- RINEHART 2007 @ NMT
	horizonAngle0000 = 0.0;
	horizonAngle0225 = 0.0;
	horizonAngle0450 = 0.0;
	horizonAngle0675 = 0.0;
	horizonAngle0900 = 0.0;
	horizonAngle1125 = 0.0;
	horizonAngle1350 = 0.0;
	horizonAngle1575 = 0.0;
	horizonAngle1800 = 0.0;
	horizonAngle2025 = 0.0;
	horizonAngle2250 = 0.0;
	horizonAngle2475 = 0.0;
	horizonAngle2700 = 0.0;
	horizonAngle2925 = 0.0;
	horizonAngle3150 = 0.0;
	horizonAngle3375 = 0.0; 
	sfact = 0.0;
	lfact = 0.0;
    
    // Added by Giuseppe Mascaro in 2016 to allow ingestion of soil grids
    Ks = 0.0;
    ThetaS = 0.0;
    ThetaR = 0.0;
    PoreSize = 0.0;
    AirEBubPres = 0.0;
    DecayF = 0.0;
    SatAnRatio = 0.0;
    UnsatAnRatio = 0.0;
    Porosity = 0.0;
    VolHeatCond = 0.0;
    SoilHeatCap = 0.0;

}

tCNode::~tCNode() {//}
	
	// GMnSKY2008MLE to fix memory leaks
	deleteVertArrays(); 
	//deleteDataStack(); //WR not necessary as TimeInd and Qeff are now smart pointers
}

//=========================================================================
//
//
//                  Section 2: tCNode Get/Set/Add Functions
//
//
//=========================================================================

// Get Functions 

double tCNode::getNwtOld() { return NwtOld; }
double tCNode::getNwtNew() { return NwtNew; }
double tCNode::getMuOld()  { return MuOld; }
double tCNode::getMuNew()  { return MuNew; }
double tCNode::getMiOld()  { return MiOld; }
double tCNode::getMiNew()  { return MiNew; }
double tCNode::getNtOld()  { return NtOld; }
double tCNode::getNtNew()  { return NtNew; }
double tCNode::getNfOld()  { return NfOld; }
double tCNode::getNfNew()  { return NfNew; }
double tCNode::getRuOld()  { return RuOld; }
double tCNode::getRuNew()  { return RuNew; }
double tCNode::getRiOld()  { return RiOld; }
double tCNode::getRiNew()  { return RiNew; }
double tCNode::getQin()    { return Qin;  }
double tCNode::getQpin()   { return Qpin; }
double tCNode::getQout()   { return Qout; }
double tCNode::getQpout()  { return Qpout; }
double tCNode::getRain()   { return Rain;  }
double tCNode::getSrf_Hr() { return srf_hr;  }
double tCNode::getSrf()    { return srf;  }
double tCNode::getCumSrf() { return cumsrf;  } // added CJC2021
double tCNode::getHsrf()   { return hsrf; }
double tCNode::getPsrf()   { return psrf; }
double tCNode::getSatsrf() { return satsrf; }
double tCNode::getSbsrf()  { return sbsrf;  }
double tCNode::getRsrf()   { return rsrf; }
double tCNode::getEsrf()   { return esrf; }
double tCNode::getRunOn()  { return RunOn; }
double tCNode::getTTime()  { return traveltime; }
double tCNode::getHillPath()   { return hillpath; }
double tCNode::getStreamPath() { return streampath; }
double tCNode::getGwaterChng() { 
    // Sort groundwater flux contributions and sum
    // This allows groundwater to be calculated the same
    // for all nodes. This is important for the parallel
    // runs. SMM - 09232008
    if (gwc.size() > 0) {
        if (gwc.size() > 1) gwc.sort();
        list<double>::iterator id;
        for (id = gwc.begin(); id != gwc.end(); ++id) {
            gwchange += (*id);
        }
        gwc.clear();
    }
    return gwchange;
}
//SMM added 10132008
list<double>& tCNode::getGwaterChngList() { return gwc; }
double tCNode::getIntStormVar()   { return intstorm; }
double tCNode::getInterceptLoss() { return Interception; }
double tCNode::getCanStorage()    { return CanStorage; }
double tCNode::getPotEvap()       { return PotEvaporation; }
double tCNode::getActEvap()       { return ActEvaporation; }
double tCNode::getNetPrecipitation() { return NetPrecipitation; }
double tCNode::getStormLength()   { return StormLength; }
double tCNode::getCumIntercept()  { return CumIntercept; }
double tCNode::getEvapoTrans()    { return EvapoTranspiration; }
double tCNode::getEvapWetCanopy() { return EvapWetCanopy; }
double tCNode::getEvapDryCanopy() { return EvapDryCanopy; }
double tCNode::getEvapSoil()      { return EvapSoil; }
double tCNode::getSoilMoisture()     { return SoilMoisture; }
double tCNode::getSoilMoistureSC()   { return SoilMoistureSC; }
double tCNode::getSoilMoistureUNSC() { return SoilMoistureUNSC; }
int    tCNode::getReach()            { return Reach; }
double tCNode::getRootMoisture()     { return RootMoisture; }
double tCNode::getRootMoistureSC()   { return RootMoistureSC; }
double tCNode::getTransmiss()  { return Transmissivity; }
double tCNode::getAirPressure(){ return AirPressure; }
double tCNode::getDewTemp()    { return DewTemp;  }
double tCNode::getRelHumid()   { return RelHumid; }
double tCNode::getVapPressure(){ return VapPress; }
double tCNode::getSkyCover()   { return SkyCover;  }
double tCNode::getWindSpeed()  { return WindSpeed; }
double tCNode::getAirTemp()    { return AirTemp;  }
double tCNode::getSurfTemp()   { return SurfTemp; }
double tCNode::getSoilTemp()   { return SoilTemp; }
double tCNode::getNetRad()     { return NetRad;}
double tCNode::getGridET()     { return GridET;}

double tCNode::getShortRadIn()     { return ShortRadIn; }
double tCNode::getShortRadIn_dir() { return ShortRadIn_dir; }
double tCNode::getShortRadIn_dif() { return ShortRadIn_dif; }
double tCNode::getShortAbsbVeg()   { return ShortAbsbVeg; }
double tCNode::getShortAbsbSoi()   { return ShortAbsbSoi; }
double tCNode::getShortRadSlope()  { return ShortRadSlope; }
double tCNode::getLongRadIn()      { return LongRadIn;  }
double tCNode::getLongRadOut()     { return LongRadOut; }

double tCNode::getGFlux() { return gFlux; }
double tCNode::getHFlux() { return hFlux; }
double tCNode::getLFlux() { return lFlux; }
double tCNode::getGnod()  { return Gnod; }

double tCNode::getChannelPerc() { return ChannelPerc; } //ASM 2/10/2017
double tCNode::getFt() {return Ft; } //ASM
double tCNode::getPercOccur() { return percOccur; } //ASM
double tCNode::getavPerc() { return avPerc; } //ASM

double tCNode::getQgwIn() { return QgwIn; }
double tCNode::getQgwOut(){ return QgwOut; }

double tCNode::getBedrockDepth() { return BedrockDepth; }
double tCNode::getContrArea() { return ContrArea; }
double tCNode::getCurvature() { return Curvature; }
double tCNode::getBasinArea() { return BasinArea; }

double tCNode::getHlevel() { return Hlevel; }
double tCNode::getQstrm()  { return Qstrm; }
double tCNode::getChannelWidth() { return Width; }
double tCNode::getRoughness()    { return Roughness; }
double tCNode::getFlowVelocity() { return FlowVelocity; }

//double tCNode::getCanopyStorage()     { return CanopyStorage; }
// SKYnGM2008LU: getCanopyStorage() now replaced by getCanopyStorVol() below.
double tCNode::getCanopyStorVol()     { return CanopyStorVol; }

double tCNode::getUnSaturatedStorage(){ return UnSaturatedStorage; }
double tCNode::getSaturatedStorage()  { return SaturatedStorage; }
double tCNode::getRecharge()          { return Recharge; }
double tCNode::getUnSatFlowIn()       { return UnSatFlowIn; }
double tCNode::getUnSatFlowOut()      { return UnSatFlowOut; }
double tCNode::getAspect()            { return Aspect; } 
double tCNode::getVegFraction()       { return VegFraction; } 
double tCNode::getAvSoilMoisture()    { return AvSoilMoisture; } 
double tCNode::getAvEvapFract()       { return AvEvapFract; } 
double tCNode::getAvET()              { return AvET; } 

// SKYnGM2008LU
//added for land use grid AJR 2007
double tCNode::getLandUseAlb() {return LandUseAlb;}
double tCNode::getThroughFall() {return ThroughFall;}
double tCNode::getVegHeight() {return VegHeight;}
double tCNode::getStomRes() {return StomRes;}

// SKYnGM2008LU
double tCNode::getIntercepCoeff() {return IntercepCoeff;}
double tCNode::getCanFieldCap() {return CanFieldCap;}
double tCNode::getDrainCoeff() {return DrainCoeff;}
double tCNode::getDrainExpPar() {return DrainExpPar;}
double tCNode::getOptTransmCoeff() {return OptTransmCoeff;}
double tCNode::getLeafAI() {return LeafAI;}
double tCNode::getCanStorParam() {return CanStorParam;}
// CJC2025: New Parameters
double tCNode::getEvapThresh() {return EvapThresh;}
double tCNode::getTransThresh() {return TransThresh;}

// SKYnGM2008LU
double tCNode::getLandUseAlbInPrevGrid() {return LandUseAlbInPrevGrid;}
double tCNode::getLandUseAlbInUntilGrid() {return LandUseAlbInUntilGrid;}
double tCNode::getThroughFallInPrevGrid() {return ThroughFallInPrevGrid;}
double tCNode::getThroughFallInUntilGrid() {return ThroughFallInUntilGrid;}
double tCNode::getVegHeightInPrevGrid() {return VegHeightInPrevGrid;}
double tCNode::getVegHeightInUntilGrid() {return VegHeightInUntilGrid;}
double tCNode::getStomResInPrevGrid() {return StomResInPrevGrid;}
double tCNode::getStomResInUntilGrid() {return StomResInUntilGrid;}
double tCNode::getVegFractionInPrevGrid() {return VegFractionInPrevGrid;}
double tCNode::getVegFractionInUntilGrid() {return VegFractionInUntilGrid;}
double tCNode::getCanStorParamInPrevGrid() {return CanStorParamInPrevGrid;}
double tCNode::getCanStorParamInUntilGrid() {return CanStorParamInUntilGrid;}
double tCNode::getIntercepCoeffInPrevGrid() {return IntercepCoeffInPrevGrid;}
double tCNode::getIntercepCoeffInUntilGrid() {return IntercepCoeffInUntilGrid;}
double tCNode::getCanFieldCapInPrevGrid() {return CanFieldCapInPrevGrid;}
double tCNode::getCanFieldCapInUntilGrid() {return CanFieldCapInUntilGrid;}
double tCNode::getDrainCoeffInPrevGrid() {return DrainCoeffInPrevGrid;}
double tCNode::getDrainCoeffInUntilGrid() {return DrainCoeffInUntilGrid;}
double tCNode::getDrainExpParInPrevGrid() {return DrainExpParInPrevGrid;}
double tCNode::getDrainExpParInUntilGrid() {return DrainExpParInUntilGrid;}
double tCNode::getOptTransmCoeffInPrevGrid() {return OptTransmCoeffInPrevGrid;}
double tCNode::getOptTransmCoeffInUntilGrid() {return OptTransmCoeffInUntilGrid;}
double tCNode::getLeafAIInPrevGrid() {return LeafAIInPrevGrid;}
double tCNode::getLeafAIInUntilGrid() {return LeafAIInUntilGrid;}
// CJC2025: New Parameters
double tCNode::getEvapThreshInPrevGrid() {return EvapThreshInPrevGrid;} 
double tCNode::getEvapThreshInUntilGrid() {return EvapThreshInUntilGrid;} 
double tCNode::getTransThreshInPrevGrid() {return TransThreshInPrevGrid;} 
double tCNode::getTransThreshInUntilGrid() {return TransThreshInUntilGrid;} 

// SKYnGM2008LU
double tCNode::getAvCanStorParam() {return AvCanStorParam;}
double tCNode::getAvIntercepCoeff() {return AvIntercepCoeff;}
double tCNode::getAvThroughFall() {return AvThroughFall;}
double tCNode::getAvCanFieldCap() {return AvCanFieldCap;}
double tCNode::getAvDrainCoeff() {return AvDrainCoeff;}
double tCNode::getAvDrainExpPar() {return AvDrainExpPar;}
double tCNode::getAvLandUseAlb() {return AvLandUseAlb;}
double tCNode::getAvVegHeight() {return AvVegHeight;}
double tCNode::getAvOptTransmCoeff() {return AvOptTransmCoeff;}
double tCNode::getAvStomRes() {return AvStomRes;}
double tCNode::getAvVegFraction() {return AvVegFraction;}
double tCNode::getAvLeafAI() {return AvLeafAI;}
// CJC2025: New Parameters
double tCNode::getAvEvapThresh() {return AvEvapThresh;}
double tCNode::getAvTransThresh() {return AvTransThresh;}

int    tCNode::getFloodStatus()  { return flood; }
int    tCNode::getSoilID()       { return soiID; }
int    tCNode::getLandUse()      { return LandUse; }
tEdge * tCNode::getFlowEdg()     { return flowedge; }
tCNode * tCNode::getStreamNode() { return StreamPtr; }

tList< int >    * tCNode::getTimeIndList() { return TimeInd.get(); }
tList< double > * tCNode::getQeffList()    { return Qeff.get(); }

// SKY2008Snow from AJR2007
// snowpack -- RINEHART 2007 @ NMT
double tCNode::getLiqWE() {return liqWEq;}//state
double tCNode::getIceWE() {return iceWEq;}//state
double tCNode::getDU() {return dU;}//flux
double tCNode::getLiqRouted() {return liqRoute;}//flux
double tCNode::getSnTempC() {return snTemperC;}//state
double tCNode::getCrustAge() {return crAge;}//state
double tCNode::getDensityAge() {return densAge;}//state
double tCNode::getEvapoTransAge() {return ETage;}//flux
double tCNode::getSnLHF() {return snLHF;}//flux
double tCNode::getSnSHF() {return snSHF;}//flux
double tCNode::getSnSub() {return snSub;}//flux // CJC2020
double tCNode::getSnEvap() {return snEvap;}//flux // CJC2020
double tCNode::getSnGHF() {return snGHF;}//flux
double tCNode::getSnPHF() {return snPHF;}//flux
double tCNode::getSnRLin() {return snRLin;}//flux
double tCNode::getSnRLout() {return snRLout;}//flux
double tCNode::getSnRSin() {return snRSin;}//flux
double tCNode::getUnode() {return Unode;}//flux
double tCNode::getUerror() {return Uerror;}//flux
double tCNode::getCumLHF() {return cumLHF;}//flux
double tCNode::getCumSnSub() {return cumSnSub;}//flux // CJC2020
double tCNode::getCumSnEvap() {return cumSnEvap;}//flux // CJC2020
double tCNode::getCumTotEvap() {return cumTotEvap;}//flux // CJC2020
double tCNode::getCumBarEvap() {return cumBarEvap;}//flux // CJC2020 
double tCNode::getCumMelt() {return cumMelt;}//flux
double tCNode::getCumSHF() {return cumSHF;}//flux
double tCNode::getCumPHF() {return cumPHF;}//flux
double tCNode::getCumRLin() {return cumRLin;}//flux
double tCNode::getCumRLout() {return cumRLout;}//flux
double tCNode::getCumRSin() {return cumRSin;}//flux
double tCNode::getCumGHF() {return cumGHF;}//flux
double tCNode::getCumHrsSun() {return cumHrsSun;}//flux
double tCNode::getCumHrsSnow() {return cumHrsSnow;}//flux
double tCNode::getCumUerror() {return cumUError;}//flux
double tCNode::getPersTimeMax() {return persTime;}
double tCNode::getPersTime() {return persTimeTemp;}
double tCNode::getPeakSWE() {return peakSWE;}
double tCNode::getPeakSWETemp() {return peakSWEtemp;}
double tCNode::getInitPackTime() {return initPackTime;}
double tCNode::getInitPackTimeTemp() {return initPackTimeTemp;}
double tCNode::getPeakPackTime() {return peakPackTime;}
// snowintercept -- RINEHART 2007 @ NMT
double tCNode::getIntSWE() {return intSWEq;}//state
double tCNode::getIntPrec() {return intPrec;}//flux 
double tCNode::getIntSnUnload() {return intSnUnload;}//flux
double tCNode::getIntSub() {return intSub;}//flux
double tCNode::getCumIntSub() {return cumIntSub;}//flux
double tCNode::getCumIntUnl() {return cumIntUnl;}//flux
// shelter -- RINEHART 2007 @ NMT
double tCNode::getHorAngle0000() {return horizonAngle0000;}//initialized state
double tCNode::getHorAngle0225() {return horizonAngle0225;}//initialized state
double tCNode::getHorAngle0450() {return horizonAngle0450;}//initialized state
double tCNode::getHorAngle0675() {return horizonAngle0675;}//initialized state
double tCNode::getHorAngle0900() {return horizonAngle0900;}//initialized state
double tCNode::getHorAngle1125() {return horizonAngle1125;}//initialized state
double tCNode::getHorAngle1350() {return horizonAngle1350;}//initialized state
double tCNode::getHorAngle1575() {return horizonAngle1575;}//initialized state
double tCNode::getHorAngle1800() {return horizonAngle1800;}//initialized state
double tCNode::getHorAngle2025() {return horizonAngle2025;}//initialized state
double tCNode::getHorAngle2250() {return horizonAngle2250;}//initialized state
double tCNode::getHorAngle2475() {return horizonAngle2475;}//initialized state
double tCNode::getHorAngle2700() {return horizonAngle2700;}//initialized state
double tCNode::getHorAngle2925() {return horizonAngle2925;}//initialized state
double tCNode::getHorAngle3150() {return horizonAngle3150;}//initialized state
double tCNode::getHorAngle3375() {return horizonAngle3375;}//initialized state
double tCNode::getSheltFact() {return sfact;}//derived
double tCNode::getLandFact() {return lfact;}//derived

// Added by Giuseppe Mascaro in 2016 to allow ingestion of soil grids
double tCNode::getKs() {return Ks;}
double tCNode::getThetaS() {return ThetaS;}
double tCNode::getThetaR() {return ThetaR;}
double tCNode::getPoreSize() {return PoreSize;}
double tCNode::getAirEBubPres() {return AirEBubPres;}
double tCNode::getDecayF() {return DecayF;}
double tCNode::getSatAnRatio() {return SatAnRatio;}
double tCNode::getUnsatAnRatio() {return UnsatAnRatio;}
double tCNode::getPorosity() {return Porosity;}
double tCNode::getVolHeatCond() {return VolHeatCond;}
double tCNode::getSoilHeatCap() {return SoilHeatCap;}

//beta update
double tCNode::getSoilCutoff() {return soil_cutoff;}
double tCNode::getRootCutoff() {return root_cutoff;}

// Set Functions

void tCNode::setMuOld(double value)  { MuOld = value; }
void tCNode::setMuNew(double value)  { MuNew = value; }
void tCNode::setRiOld(double value)  { RiOld = value; }
void tCNode::setRiNew(double value)  { RiNew = value; }
void tCNode::setRuOld(double value)  { RuOld = value; }
void tCNode::setRuNew(double value)  { RuNew = value; }
void tCNode::setNfOld(double value)  { NfOld = value; }
void tCNode::setNfNew(double value)  { NfNew = value; }
void tCNode::setNtOld(double value)  { NtOld = value; }
void tCNode::setNtNew(double value)  { NtNew = value; }
void tCNode::setNwtOld(double value) { NwtOld = value; }
void tCNode::setNwtNew(double value) { NwtNew = value; }
void tCNode::setMiOld(double value)  { MiOld = value; }
void tCNode::setMiNew(double value)  { MiNew = value; }
void tCNode::setQpout(double value)  { Qpout = value; }
void tCNode::setQpin(double value)   { Qpin = value; }
void tCNode::setRain(double value)   { Rain = value; }
void tCNode::setSrf_Hr(double value) { srf_hr = value; }
void tCNode::setsrf(double value)    { srf = value; }
void tCNode::sethsrf(double value)   { hsrf = value; }
void tCNode::setesrf(double value)   { esrf = value; }
void tCNode::setpsrf(double value)   { psrf = value; }
void tCNode::setsatsrf(double value) { satsrf = value; }
void tCNode::setsbsrf(double value)  { sbsrf = value; }
void tCNode::setRunOn(double value)  { RunOn = value; }
void tCNode::setFlowEdg(tEdge * edgs) { flowedge = edgs; }
void tCNode::setGwaterChng(double value)  { 
   gwchange = value; 
   gwc.clear(); //SMM - 09232008, empty list of groundwater contributions
}
void tCNode::setFloodStatus( int status )  { flood = status; }
void tCNode::setIntStormVar( double value ) { intstorm = value; }
void tCNode::setTTime(double value)  { traveltime = value; } 
void tCNode::setSoilID(int ids)      { soiID = ids; }
void tCNode::setHillPath(double l)   { hillpath = l; }
void tCNode::setStreamPath(double l) { streampath = l; }
void tCNode::setLandUse(int landuse) { LandUse = landuse; }
void tCNode::setReach(int r) { Reach = r; }
void tCNode::setInterceptLoss(double interception) {  
	Interception = interception; }
void tCNode::setNetPrecipitation(double netprecip) {  
	NetPrecipitation = netprecip; }
void tCNode::setCanStorage(double canstore) { CanStorage = canstore; }
void tCNode::setPotEvap(double potEvap) { PotEvaporation = potEvap; }
void tCNode::setActEvap(double actEvap) { ActEvaporation = actEvap; }
void tCNode::setTransmiss(double Tr)    { Transmissivity = Tr; }
void tCNode::setStormLength(int option, double time){
	if (option == 0)
		StormLength = 0.0;
	else if (option == 1)
		StormLength+= time;
}
void tCNode::setCumIntercept(double cum) {CumIntercept = cum;}
void tCNode::setEvapoTrans(double evapoTrans) {
	EvapoTranspiration = evapoTrans;}
void tCNode::setEvapWetCanopy(double evapWetCanopy) {
	EvapWetCanopy = evapWetCanopy;}
void tCNode::setEvapDryCanopy(double evapDryCanopy) {
	EvapDryCanopy = evapDryCanopy;}
void tCNode::setEvapSoil(double evapSoil) {EvapSoil = evapSoil;}
void tCNode::setSoilMoisture(double soilMoisture) {
	SoilMoisture = soilMoisture;}
void tCNode::setSoilMoistureSC(double soilMoisture) {
	SoilMoistureSC = soilMoisture;}
void tCNode::setSoilMoistureUNSC(double soilMoisture) {
	SoilMoistureUNSC = soilMoisture;}
void tCNode::setRootMoisture(double rootMoisture) {
	RootMoisture = rootMoisture;}
void tCNode::setRootMoistureSC(double rootMoisture) {
	RootMoistureSC = rootMoisture;}
void tCNode::setAirPressure(double airpress){ AirPressure = airpress;}
void tCNode::setDewTemp(double dewtemp)     { DewTemp = dewtemp;}
void tCNode::setRelHumid(double relhumid)   { RelHumid = relhumid;}
void tCNode::setVapPressure(double vappres) { VapPress = vappres;}
void tCNode::setSkyCover(double skycover)   { SkyCover = skycover;}
void tCNode::setWindSpeed(double windspeed) { WindSpeed = windspeed;}
void tCNode::setAirTemp(double airtemp)     { AirTemp = airtemp;}
void tCNode::setSurfTemp(double surftemp)   { SurfTemp = surftemp; }
void tCNode::setSoilTemp(double soiltemp)   { SoilTemp = soiltemp; }
void tCNode::setNetRad(double netrad)       { NetRad = netrad; }
void tCNode::setGridET(double gridet)       { GridET = gridet; }

void tCNode::setShortRadIn(double val)     { ShortRadIn = val; }
void tCNode::setShortRadIn_dir(double val) { ShortRadIn_dir = val; }
void tCNode::setShortRadIn_dif(double val) { ShortRadIn_dif = val; }
void tCNode::setShortAbsbVeg(double val)   { ShortAbsbVeg = val; } 
void tCNode::setShortAbsbSoi(double val)   { ShortAbsbSoi = val; }
void tCNode::setShortRadSlope(double val)  { ShortRadSlope = val; }
void tCNode::setLongRadIn(double longin)   { LongRadIn = longin; }
void tCNode::setLongRadOut(double longout) { LongRadOut = longout; }

void tCNode::setGFlux(double gF){ gFlux = gF; }
void tCNode::setHFlux(double hF){ hFlux = hF; }
void tCNode::setLFlux(double lF){ lFlux = lF; }
void tCNode::setGnod(double Go){ Gnod = Go; }

void tCNode::setChannelPerc(double cP){ ChannelPerc = cP; } //ASM 2/10/2017
void tCNode::setFt(double Ft_prime){ Ft = Ft_prime; } //ASM

void tCNode::setBedrockDepth(double brock) { BedrockDepth = brock; }
void tCNode::setContrArea(double value)    { ContrArea = value; }
void tCNode::setCurvature(double value)    { Curvature = value; }
void tCNode::setBasinArea(double barea)    { BasinArea = barea; }

void tCNode::setHlevel(double value)       { Hlevel = value; }
void tCNode::setQstrm(double value)        { Qstrm = value; }
void tCNode::setChannelWidth(double value) { Width = value; }
void tCNode::setRoughness(double value)    { Roughness = value; }
void tCNode::setFlowVelocity(double value) { FlowVelocity = value; }

//void tCNode::setCanopyStorage(double value)      { CanopyStorage = value; }
// SKYnGM2008LU: setCanopyStorage(double) now replaced by setCanopyStorVol() below.
void tCNode::setCanopyStorVol(double value) { CanopyStorVol = value; }

void tCNode::setUnSaturatedStorage(double value) { UnSaturatedStorage = value; }
void tCNode::setSaturatedStorage(double value)   { SaturatedStorage = value; }
void tCNode::setRecharge(double value)           { Recharge = value; }
void tCNode::setUnSatFlowIn(double value)        { UnSatFlowIn = value; }
void tCNode::setUnSatFlowOut(double value)       { UnSatFlowOut = value; }

void tCNode::setAspect(double value)   { Aspect = value; }
void tCNode::setVegFraction(double value) { VegFraction = value; }
void tCNode::setAvSoilMoisture(double value) { AvSoilMoisture = value; }
void tCNode::setAvEvapFract(double value)    { AvEvapFract = value; }
void tCNode::setAvET(double value)           { AvET = value; }

// SKY2008Snow from AJR2007
//snowpack -- RINEHART 2007 @ NMT
void tCNode::setLiqWE(double lwe) {liqWEq = lwe;}//state
void tCNode::setIceWE(double iwe) {iceWEq = iwe;}//state
void tCNode::setDU(double value) {dU = value;}//flux
void tCNode::setLiqRouted(double lr) {liqRoute = lr;}//flux
void tCNode::setSnTempC(double temperature) {snTemperC = temperature;}//state
void tCNode::setCrustAge(double ca) {crAge = ca;}//state
void tCNode::setDensityAge(double da) {densAge = da;}//state
void tCNode::setEvapoTransAge(double ea) {ETage = ea;}//flux
void tCNode::setSnLHF(double value) {snLHF = value;}//flux
void tCNode::setSnSHF(double value) {snSHF = value;}//flux
void tCNode::setSnSub(double value) {snSub = value;}//flux // Snowpack sublimation CJC2020
void tCNode::setSnEvap(double value) {snEvap = value;}//flux // Snowpack evaporation CJC2020
void tCNode::setSnGHF(double value) {snGHF = value;}//flux
void tCNode::setSnPHF(double value) {snPHF = value;}//flux
void tCNode::setSnRLout(double value) {snRLout = value;}
void tCNode::setSnRLin(double value) {snRLin = value;}//flux
void tCNode::setSnRSin(double value) {snRSin = value;}//flux
void tCNode::setUnode(double value) {Unode = value;}//flux
void tCNode::setUerror(double value) {Uerror = value;}//flux
void tCNode::setPersTimeMax(double value) {persTime = value;}
void tCNode::setPersTime(double value) {persTimeTemp = value;}
void tCNode::setPeakSWE(double value) {peakSWE = value;}
void tCNode::setPeakSWETemp(double value) {peakSWEtemp = value;}
void tCNode::setInitPackTime(double value) {initPackTime = value;}
void tCNode::setInitPackTimeTemp(double value) {initPackTimeTemp = value;}
void tCNode::setPeakPackTime(double value) {peakPackTime = value;}
//snowintercept -- RINEHART 2007 @ NMT
void tCNode::setIntSWE(double swe) {intSWEq = swe;}//state
void tCNode::setIntPrec(double value) {intPrec = value;}//flux
void tCNode::setIntSnUnload(double we) {intSnUnload = we;}//flux
void tCNode::setIntSub(double is) {intSub = is;}//flux
//shelter -- RINEHART 2007 @ NMT
void tCNode::setHorAngle0000(double angle) {horizonAngle0000 = angle;}//initialized state
void tCNode::setHorAngle0225(double angle) {horizonAngle0225 = angle;}//initialized state
void tCNode::setHorAngle0450(double angle) {horizonAngle0450 = angle;}//initialized state
void tCNode::setHorAngle0675(double angle) {horizonAngle0675 = angle;}//initialized state
void tCNode::setHorAngle0900(double angle) {horizonAngle0900 = angle;}//initialized state
void tCNode::setHorAngle1125(double angle) {horizonAngle1125 = angle;}//initialized state
void tCNode::setHorAngle1350(double angle) {horizonAngle1350 = angle;}//initialized state
void tCNode::setHorAngle1575(double angle) {horizonAngle1575 = angle;}//initialized state
void tCNode::setHorAngle1800(double angle) {horizonAngle1800 = angle;}//initialized state
void tCNode::setHorAngle2025(double angle) {horizonAngle2025 = angle;}//initialized state
void tCNode::setHorAngle2250(double angle) {horizonAngle2250 = angle;}//initialized state
void tCNode::setHorAngle2475(double angle) {horizonAngle2475 = angle;}//initialized state
void tCNode::setHorAngle2700(double angle) {horizonAngle2700 = angle;}//initialized state
void tCNode::setHorAngle2925(double angle) {horizonAngle2925 = angle;}//initialized state
void tCNode::setHorAngle3150(double angle) {horizonAngle3150 = angle;}//initialized state
void tCNode::setHorAngle3375(double angle) {horizonAngle3375 = angle;}//initialized state
void tCNode::setSheltFact(double sf) {sfact = sf;}//derived
void tCNode::setLandFact(double lf) {lfact = lf;}//derived


void tCNode::setStreamNode(tCNode *cn) { StreamPtr = cn; } 

// SKYnGM2008LU
//added for landuse grid AJR 2007
void tCNode::setLandUseAlb(double value) { LandUseAlb = value; }
void tCNode::setThroughFall(double value) { ThroughFall = value; }
void tCNode::setVegHeight(double value) { VegHeight = value; }
void tCNode::setStomRes(double value) {StomRes = value; }

// SKYnGM2008LU
void tCNode::setIntercepCoeff(double value) { IntercepCoeff = value; }
void tCNode::setCanFieldCap(double value) { CanFieldCap = value; }
void tCNode::setDrainCoeff(double value) { DrainCoeff = value; }
void tCNode::setDrainExpPar(double value) {DrainExpPar = value; }
void tCNode::setOptTransmCoeff(double value) { OptTransmCoeff = value; }
void tCNode::setLeafAI(double value) { LeafAI = value; }
void tCNode::setCanStorParam(double value) { CanStorParam = value; }
// CJC2025: New Parameters
void tCNode::setEvapThresh(double value) { EvapThresh = value; }
void tCNode::setTransThresh(double value) { TransThresh = value; }

// SKYnGM2008LU
void tCNode::setLandUseAlbInPrevGrid(double value) { LandUseAlbInPrevGrid = value; }
void tCNode::setLandUseAlbInUntilGrid(double value) { LandUseAlbInUntilGrid = value; }
void tCNode::setThroughFallInPrevGrid(double value) { ThroughFallInPrevGrid = value; }
void tCNode::setThroughFallInUntilGrid(double value) { ThroughFallInUntilGrid = value; }
void tCNode::setVegHeightInPrevGrid(double value) { VegHeightInPrevGrid = value; }
void tCNode::setVegHeightInUntilGrid(double value) { VegHeightInUntilGrid = value; }
void tCNode::setStomResInPrevGrid(double value) { StomResInPrevGrid = value; }
void tCNode::setStomResInUntilGrid(double value) { StomResInUntilGrid = value; }
void tCNode::setVegFractionInPrevGrid(double value) { VegFractionInPrevGrid = value; }
void tCNode::setVegFractionInUntilGrid(double value) { VegFractionInUntilGrid = value; }
void tCNode::setCanStorParamInPrevGrid(double value) { CanStorParamInPrevGrid = value; }
void tCNode::setCanStorParamInUntilGrid(double value) { CanStorParamInUntilGrid = value; }
void tCNode::setIntercepCoeffInPrevGrid(double value) { IntercepCoeffInPrevGrid = value; }
void tCNode::setIntercepCoeffInUntilGrid(double value) { IntercepCoeffInUntilGrid = value; }
void tCNode::setCanFieldCapInPrevGrid(double value) { CanFieldCapInPrevGrid = value; }
void tCNode::setCanFieldCapInUntilGrid(double value) { CanFieldCapInUntilGrid = value; }
void tCNode::setDrainCoeffInPrevGrid(double value) { DrainCoeffInPrevGrid = value; }
void tCNode::setDrainCoeffInUntilGrid(double value) { DrainCoeffInUntilGrid = value; }
void tCNode::setDrainExpParInPrevGrid(double value) { DrainExpParInPrevGrid = value; }
void tCNode::setDrainExpParInUntilGrid(double value) { DrainExpParInUntilGrid = value; }
void tCNode::setOptTransmCoeffInPrevGrid(double value) { OptTransmCoeffInPrevGrid = value; }
void tCNode::setOptTransmCoeffInUntilGrid(double value) { OptTransmCoeffInUntilGrid = value; }
void tCNode::setLeafAIInPrevGrid(double value) { LeafAIInPrevGrid = value; }
void tCNode::setLeafAIInUntilGrid(double value) { LeafAIInUntilGrid = value; }
// CJC2025: New Parameters
void tCNode::setEvapThreshInPrevGrid(double value) { EvapThreshInPrevGrid = value; }
void tCNode::setEvapThreshInUntilGrid(double value) { EvapThreshInUntilGrid = value; }
void tCNode::setTransThreshInPrevGrid(double value) { TransThreshInPrevGrid = value; }
void tCNode::setTransThreshInUntilGrid(double value) { TransThreshInUntilGrid = value; }

// SKYnGM2008LU
void tCNode::setAvCanStorParam(double value) { AvCanStorParam = value; }
void tCNode::setAvIntercepCoeff(double value) { AvIntercepCoeff = value; }
void tCNode::setAvThroughFall(double value) { AvThroughFall = value; }
void tCNode::setAvCanFieldCap(double value) { AvCanFieldCap = value; }
void tCNode::setAvDrainCoeff(double value) { AvDrainCoeff = value; }
void tCNode::setAvDrainExpPar(double value) { AvDrainExpPar = value; }
void tCNode::setAvLandUseAlb(double value) { AvLandUseAlb = value; }
void tCNode::setAvVegHeight(double value) { AvVegHeight = value; }
void tCNode::setAvOptTransmCoeff(double value) { AvOptTransmCoeff = value; }
void tCNode::setAvStomRes(double value) { AvStomRes = value; }
void tCNode::setAvVegFraction(double value) { AvVegFraction = value; }
void tCNode::setAvLeafAI(double value) { AvLeafAI = value; }
// CJC2025: New Parameters
void tCNode::setAvEvapThresh(double value) { AvEvapThresh = value; }
void tCNode::setAvTransThresh(double value) { AvTransThresh = value; }

// Added by Giuseppe Mascaro in 2016 to allow ingestion of soil grids
void tCNode::setKs(double value) { Ks = value;}
void tCNode::setThetaS(double value) { ThetaS = value;}
void tCNode::setThetaR(double value) { ThetaR = value;}
void tCNode::setPoreSize(double value) { PoreSize = value;}
void tCNode::setAirEBubPres(double value) { AirEBubPres = value;}
void tCNode::setDecayF(double value) { DecayF = value;}
void tCNode::setSatAnRatio(double value) { SatAnRatio = value;}
void tCNode::setUnsatAnRatio(double value) { UnsatAnRatio = value;}
void tCNode::setPorosity(double value) { Porosity = value;}
void tCNode::setVolHeatCond(double value) { VolHeatCond = value;}
void tCNode::setSoilHeatCap(double value) { SoilHeatCap = value;}



//beta update
void tCNode::setSoilCutoff(double value) {soil_cutoff = value;}
void tCNode::setRootCutoff(double value) {root_cutoff = value;}



// Add Functions

void tCNode::addGwaterChng(double value) { 
    // Add groundwater contribution to list for later
    // sorting and summing. SMM - 09232008
    gwc.push_back(value);
}
void tCNode::addQpin(double value)  { Qpin += value; }
void tCNode::addTTime(double value) { traveltime += value; } 
void tCNode::addIntStormVar(double value) { intstorm += value; }
void tCNode::addContrArea(double value)   { ContrArea += value; }
void tCNode::addQgwIn(double value) { QgwIn += value; }
void tCNode::addQgwOut(double value){ QgwOut += value; }
void tCNode::addQstrm(double value) { Qstrm += value; }
void tCNode::addAvSoilMoisture(double value) { AvSoilMoisture += value; }
void tCNode::addSrf_Hr(double value) { srf_hr += value; }
void tCNode::addCumSrf(double value) { cumsrf += value; } // added CJC2021

// SKY2008Snow from AJR2007
//snow, snow interception and sheltering -- RINEHART 2007 @ NMT
void tCNode::addLatHF(double value) { cumLHF += value;}
void tCNode::addSHF(double value) {cumSHF += value;}
void tCNode::addSnSub(double value) { cumSnSub += value;} // Snowpack sublimation CJC2020
void tCNode::addSnEvap(double value) {cumSnEvap += value;} // Snowpack evaporation CJC2020
void tCNode::addTotEvap(double value) {cumTotEvap += value;} // total ET CJC2020
void tCNode::addBarEvap(double value) {cumBarEvap += value;} // bare soil evap CJC2020
void tCNode::addPHF(double value) {cumPHF += value;}
void tCNode::addRLin(double value) {cumRLin += value;}
void tCNode::addRLout(double value) {cumRLout += value;}
void tCNode::addRSin(double value) {cumRSin += value;}
void tCNode::addGHF(double value) {cumGHF += value;}
void tCNode::addMelt(double value) { cumMelt += value;}
void tCNode::addIntSub(double value) { cumIntSub += value;}
void tCNode::addIntUnl(double value) { cumIntUnl += value;}
void tCNode::addCumHrsSun(double value) { cumHrsSun += value;}
void tCNode::addCumUerror(double value) {cumUError += value;}
void tCNode::addCumHrsSnow(double value) {cumHrsSnow += value;}

//=========================================================================
//
//
//                  Section 3: tCNode Tracer Sorting
//
//
//=========================================================================

/**************************************************************************
**
**  tCNode:: Tracer-sorting routines:
**
**  These routines are utilities that are used in sorting the nodes
**  according to their position within the drainage network. The main
**  sorting algorithm is implemented in tFlowNet::SortNodesByNetOrder().
**  The sorting method works by introducing a "tracer" at each point,
**  then allowing the tracers to iteratively cascade downstream. At each
**  step any nodes not containing tracers are moved to the back of the
**  list. The result is a list sorted in upstream-to-downstream order.
**
**  These utilities do the following:
**    ActivateSortTracer -- injects a single tracer at a node
**    AddTracer -- adds a tracer to a node (ignored if node is a bdy)
**    MoveSortTracerDownstream -- removes a tracer and sends it to the
**                                downstream neighbor (unless the node is
**                                a sink; then the tracer just vanishes)
**    NoMoreTracers -- reports whether there are any tracers left here
**
**************************************************************************/

void tCNode::setTracer(int cnt)
{ 
	tracer = cnt; 
}

int tCNode::getTracer()
{ 
	return tracer; 
}

void tCNode::ActivateSortTracer()
{ 
	tracer = 1; 
}

void tCNode::DeactivateTracer()
{
	tracer = 0;
}

void tCNode::MoveSortTracerDownstream()
{
	tracer--;
	getDownstrmNbr()->AddTracer();
}

void tCNode::AddTracer()
{
	if((boundary==0)||(boundary==3)) tracer++;
}

int tCNode::NoMoreTracers()
{
	return( tracer==0 );
}

tCNode * tCNode::getDownstrmNbr()
{
	if( flowedge == 0 ) return 0;
	return (tCNode *)flowedge->getDestinationPtrNC();     
}

void tCNode::deleteVertArrays()
{
	delete [] VertsX;
	delete [] VertsY;
	nVerts = -999; 
}

void tCNode::allocVertArrays(int n)
{
	nVerts = n;
	VertsX = new double [ nVerts ];
	assert(VertsX != 0);
	VertsY = new double [ nVerts ];
	assert(VertsY != 0);
	return;
}

void tCNode::allocDataStack()
{
    //WR debug convert to smart pointers
	//TimeInd = new tList< int >;
	//Qeff = new tList< double >;
    TimeInd = std::make_shared<tList<int> >();
    Qeff = std::make_shared<tList<double> >();

	assert(TimeInd != 0);
	assert(Qeff != 0);
	return; 
}

void tCNode::deleteDataStack()
{
/*SMM - 08132008, this code is currently causing an error
	if (TimeInd != 0) {
		if (!(TimeInd->isEmpty())) 
			TimeInd->Flush();
		delete TimeInd;
	}
	if (Qeff != 0) {
		if (!(Qeff->isEmpty())) 
			Qeff->Flush();
		delete Qeff;
	}
	return;
*/
}

//=========================================================================
//
//
//                  Section 4: tCNode Associated Voronoi Geometry
//
//
//=========================================================================

/**************************************************************************
**
**  tCNode:: Get Centroid of associated voronoi polygon
**
**************************************************************************/
double tCNode::getCentroidX()
{
	if(xC == -1){
		double areaT=0.0;
	
		tArray<double> xy(2);
		int      nPoints;
		int cnt=0;
		tEdge *firstedg;
		tEdge *curedg;
		
		firstedg = getFlowEdg();
		if (!firstedg)
			firstedg = getEdg();
		curedg = firstedg->getCCWEdg();
		cnt++;
		while (curedg != firstedg) {
			curedg = curedg->getCCWEdg();
			cnt++;
		}

		nPoints = cnt;
		
		double *vXs = new double [cnt+1];
		double *vYs = new double [cnt+1];

		int iv = 0;
		firstedg = getFlowEdg(); 
	
		xy = firstedg->getRVtx();
		vXs[iv] = xy[0];
		vYs[iv] = xy[1];
		iv++;
		curedg = firstedg->getCCWEdg();
		while (curedg != firstedg) {
			xy = curedg->getRVtx();
			if (xy[0] == vXs[iv-1] && xy[1] == vYs[iv-1]) {
				iv--;         //If points coincide
				nPoints--; //just skip it...
			} 
			else { 
				vXs[iv] = xy[0];
				vYs[iv] = xy[1];
			}
			iv++;
			curedg = curedg->getCCWEdg();
		}
		vXs[iv] = vXs[0];
		vYs[iv] = vYs[0];

		cnt = polyCentroid(vXs, vYs, nPoints, 
							&xC, &yC, &areaT);
							
		delete [] vXs;
		delete [] vYs;
							
	}
	return xC;
}

double tCNode::getCentroidY()
{
	if(yC == -1){
		double areaT=0.0;
	
		tArray<double> xy(2);
		int      nPoints;
		int cnt=0;
		tEdge *firstedg;
		tEdge *curedg;
		
		firstedg = getFlowEdg();
		if (!firstedg)
			firstedg = getEdg();
		curedg = firstedg->getCCWEdg();
		cnt++;
		while (curedg != firstedg) {
			curedg = curedg->getCCWEdg();
			cnt++;
		}

		nPoints = cnt;
		
		double *vXs = new double [cnt+1];
		double *vYs = new double [cnt+1];

		int iv = 0;
		firstedg = getFlowEdg(); 
	
		xy = firstedg->getRVtx();
		vXs[iv] = xy[0];
		vYs[iv] = xy[1];
		iv++;
		curedg = firstedg->getCCWEdg();
		while (curedg != firstedg) {
			xy = curedg->getRVtx();
			if (xy[0] == vXs[iv-1] && xy[1] == vYs[iv-1]) {
				iv--;         //If points coincide
				nPoints--; //just skip it...
			} 
			else { 
				vXs[iv] = xy[0];
				vYs[iv] = xy[1];
			}
			iv++;
			curedg = curedg->getCCWEdg();
		}
		vXs[iv] = vXs[0];
		vYs[iv] = vYs[0];

		cnt = polyCentroid(vXs, vYs, nPoints, 
							&xC, &yC, &areaT);
							
		delete [] vXs;
		delete [] vYs;
							
	}
	return yC;
}

int tCNode::polyCentroid(double x[], double y[], int n,
						double *xCentroid, double *yCentroid, double *area) 
{
	int i, j;
	double ai, atmp = 0.0, xtmp = 0.0, ytmp = 0.0;
	if (n < 3) return 1;
	for (i = n-1, j = 0; j < n; i = j, j++)
    {
		ai = x[i] * y[j] - x[j] * y[i];
		atmp += ai;
		xtmp += (x[j] + x[i]) * ai;
		ytmp += (y[j] + y[i]) * ai;
    }
	*area = atmp / 2.0;
	if (atmp != 0.0)
    {
		*xCentroid =  xtmp / (3.0 * atmp);
		*yCentroid =  ytmp / (3.0 * atmp);
		return 0;
    }
	return 2;
}

/***************************************************************************
**
** tCNode::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/

void tCNode::writeRestart(fstream& rStr) const
{
  BinaryWrite(rStr, srf_hr);
  BinaryWrite(rStr, cumsrf); //added CJC2021
  BinaryWrite(rStr, RunOn);
  BinaryWrite(rStr, VapPress);
  BinaryWrite(rStr, ShortRadIn_dir);
  BinaryWrite(rStr, ShortRadIn_dif);
  BinaryWrite(rStr, ShortAbsbVeg);
  BinaryWrite(rStr, ShortAbsbSoi);
  BinaryWrite(rStr, Gnod);
  BinaryWrite(rStr, VegFraction);

  BinaryWrite(rStr, tracer);
  BinaryWrite(rStr, flood);
  BinaryWrite(rStr, soiID);
  BinaryWrite(rStr, LandUse);
  BinaryWrite(rStr, Reach);
  BinaryWrite(rStr, Qin);
  BinaryWrite(rStr, Qout);
  BinaryWrite(rStr, traveltime);
  BinaryWrite(rStr, hillpath);
  BinaryWrite(rStr, streampath);
  BinaryWrite(rStr, GridET);
  BinaryWrite(rStr, ContrArea);
  BinaryWrite(rStr, Curvature);
  BinaryWrite(rStr, BasinArea);
  BinaryWrite(rStr, BedrockDepth);
  BinaryWrite(rStr, hFlux);
  BinaryWrite(rStr, QgwIn);
  BinaryWrite(rStr, QgwOut);
  BinaryWrite(rStr, Aspect);
  BinaryWrite(rStr, Width);
  BinaryWrite(rStr, Roughness);

  BinaryWrite(rStr, satOccur);
  BinaryWrite(rStr, hsrfOccur);
  BinaryWrite(rStr, percOccur); //ASM
  BinaryWrite(rStr, avPerc);  //ASM
  BinaryWrite(rStr, psrfOccur);
  BinaryWrite(rStr, satsrfOccur);
  BinaryWrite(rStr, sbsrfOccur);
  BinaryWrite(rStr, RechDisch);
  BinaryWrite(rStr, NwtOld);
  BinaryWrite(rStr, MuOld);
  BinaryWrite(rStr, MiOld);
  BinaryWrite(rStr, NtOld);
  BinaryWrite(rStr, NfOld);
  BinaryWrite(rStr, RuOld);
  BinaryWrite(rStr, RiOld);
  BinaryWrite(rStr, Qpout);
  BinaryWrite(rStr, Rain);
  BinaryWrite(rStr, intstorm);
  BinaryWrite(rStr, Interception);
  BinaryWrite(rStr, NetPrecipitation);
  BinaryWrite(rStr, CanStorage);
  BinaryWrite(rStr, PotEvaporation);
  BinaryWrite(rStr, ActEvaporation);
  BinaryWrite(rStr, StormLength);
  BinaryWrite(rStr, CumIntercept);
  BinaryWrite(rStr, EvapWetCanopy);
  BinaryWrite(rStr, EvapDryCanopy);
  BinaryWrite(rStr, EvapSoil);
  BinaryWrite(rStr, EvapoTranspiration);
  BinaryWrite(rStr, SoilMoisture);
  BinaryWrite(rStr, SoilMoistureSC);
  BinaryWrite(rStr, SoilMoistureUNSC);
  BinaryWrite(rStr, RootMoisture);
  BinaryWrite(rStr, RootMoistureSC);
  BinaryWrite(rStr, Transmissivity);
  BinaryWrite(rStr, AirTemp);
  BinaryWrite(rStr, DewTemp);
  BinaryWrite(rStr, RelHumid);
  BinaryWrite(rStr, SurfTemp);
  BinaryWrite(rStr, SoilTemp);
  BinaryWrite(rStr, WindSpeed);
  BinaryWrite(rStr, SkyCover);
  BinaryWrite(rStr, NetRad);
  BinaryWrite(rStr, AirPressure);
  BinaryWrite(rStr, ShortRadIn);
  BinaryWrite(rStr, LongRadIn);
  BinaryWrite(rStr, LongRadOut);
  BinaryWrite(rStr, gFlux);
  BinaryWrite(rStr, lFlux);
  BinaryWrite(rStr, AvSoilMoisture);
  BinaryWrite(rStr, AvEvapFract);
  BinaryWrite(rStr, AvET);
  BinaryWrite(rStr, Hlevel);
  BinaryWrite(rStr, Qstrm);
  BinaryWrite(rStr, FlowVelocity);
  BinaryWrite(rStr, CanopyStorVol);
  BinaryWrite(rStr, UnSaturatedStorage);
  BinaryWrite(rStr, SaturatedStorage);
  BinaryWrite(rStr, Recharge);
  BinaryWrite(rStr, UnSatFlowIn);
  BinaryWrite(rStr, UnSatFlowOut);
  BinaryWrite(rStr, ChannelPerc); //ASM 2/10/2017
  BinaryWrite(rStr, Ft); //ASM

  BinaryWrite(rStr, liqWEq); // Snowpack
  BinaryWrite(rStr, iceWEq);
  BinaryWrite(rStr, dU);
  BinaryWrite(rStr, Unode);
  BinaryWrite(rStr, Uerror);
  BinaryWrite(rStr, cumUError);
  BinaryWrite(rStr, liqRoute);
  BinaryWrite(rStr, snTemperC);
  BinaryWrite(rStr, crAge);
  BinaryWrite(rStr, densAge);
  BinaryWrite(rStr, ETage);
  BinaryWrite(rStr, snLHF);
  BinaryWrite(rStr, snSHF);
  BinaryWrite(rStr, snGHF);
  BinaryWrite(rStr, snPHF);
  BinaryWrite(rStr, snRLin);
  BinaryWrite(rStr, snRLout);
  BinaryWrite(rStr, snRSin);
  BinaryWrite(rStr, persTime);
  BinaryWrite(rStr, persTimeTemp);
  BinaryWrite(rStr, peakSWE);
  BinaryWrite(rStr, peakSWEtemp);
  BinaryWrite(rStr, initPackTime);
  BinaryWrite(rStr, initPackTimeTemp);
  BinaryWrite(rStr, peakPackTime);

  BinaryWrite(rStr, intSWEq); // Snow intercept
  BinaryWrite(rStr, intSnUnload);
  BinaryWrite(rStr, intSub);
  BinaryWrite(rStr, intPrec);

  BinaryWrite(rStr, horizonAngle0000); // shelter
  BinaryWrite(rStr, horizonAngle0225);
  BinaryWrite(rStr, horizonAngle0450);
  BinaryWrite(rStr, horizonAngle0675);
  BinaryWrite(rStr, horizonAngle0900);
  BinaryWrite(rStr, horizonAngle1125);
  BinaryWrite(rStr, horizonAngle1350);
  BinaryWrite(rStr, horizonAngle1575);
  BinaryWrite(rStr, horizonAngle1800);
  BinaryWrite(rStr, horizonAngle2025);
  BinaryWrite(rStr, horizonAngle2250);
  BinaryWrite(rStr, horizonAngle2475);
  BinaryWrite(rStr, horizonAngle2700);
  BinaryWrite(rStr, horizonAngle2925);
  BinaryWrite(rStr, horizonAngle3150);
  BinaryWrite(rStr, horizonAngle3375);
  BinaryWrite(rStr, sfact);
  BinaryWrite(rStr, lfact);
  BinaryWrite(rStr, VegFraction);

  BinaryWrite(rStr, LandUseAlb); // Additions for landuse grids
  BinaryWrite(rStr, ThroughFall);
  BinaryWrite(rStr, VegHeight);
  BinaryWrite(rStr, StomRes);
  BinaryWrite(rStr, IntercepCoeff);
  BinaryWrite(rStr, CanFieldCap);
  BinaryWrite(rStr, DrainCoeff);
  BinaryWrite(rStr, DrainExpPar);
  BinaryWrite(rStr, OptTransmCoeff);
  BinaryWrite(rStr, LeafAI);
  BinaryWrite(rStr, CanStorParam);
  // CJC2025: New Parameters
  BinaryWrite(rStr, EvapThresh);
  BinaryWrite(rStr, TransThresh);

  BinaryWrite(rStr, LandUseAlbInPrevGrid);
  BinaryWrite(rStr, LandUseAlbInUntilGrid);
  BinaryWrite(rStr, ThroughFallInPrevGrid);
  BinaryWrite(rStr, ThroughFallInUntilGrid);
  BinaryWrite(rStr, VegHeightInPrevGrid);
  BinaryWrite(rStr, VegHeightInUntilGrid);
  BinaryWrite(rStr, StomResInPrevGrid);
  BinaryWrite(rStr, StomResInUntilGrid);
  BinaryWrite(rStr, VegFractionInPrevGrid);
  BinaryWrite(rStr, VegFractionInUntilGrid);
  BinaryWrite(rStr, CanStorParamInPrevGrid);
  BinaryWrite(rStr, CanStorParamInUntilGrid);
  BinaryWrite(rStr, IntercepCoeffInPrevGrid);
  BinaryWrite(rStr, IntercepCoeffInUntilGrid);
  BinaryWrite(rStr, CanFieldCapInPrevGrid);
  BinaryWrite(rStr, CanFieldCapInUntilGrid);
  BinaryWrite(rStr, DrainCoeffInPrevGrid);
  BinaryWrite(rStr, DrainCoeffInUntilGrid);
  BinaryWrite(rStr, DrainExpParInPrevGrid);
  BinaryWrite(rStr, DrainExpParInUntilGrid);
  BinaryWrite(rStr, OptTransmCoeffInPrevGrid);
  BinaryWrite(rStr, OptTransmCoeffInUntilGrid);
  BinaryWrite(rStr, LeafAIInPrevGrid);
  BinaryWrite(rStr, LeafAIInUntilGrid);
  // CJC2025: New Parameters
  BinaryWrite(rStr, EvapThreshInPrevGrid);
  BinaryWrite(rStr, EvapThreshInUntilGrid);
  BinaryWrite(rStr, TransThreshInPrevGrid);
  BinaryWrite(rStr, TransThreshInUntilGrid);

  BinaryWrite(rStr, AvCanStorParam);
  BinaryWrite(rStr, AvIntercepCoeff);
  BinaryWrite(rStr, AvThroughFall);
  BinaryWrite(rStr, AvCanFieldCap);
  BinaryWrite(rStr, AvDrainCoeff);
  BinaryWrite(rStr, AvDrainExpPar);
  BinaryWrite(rStr, AvLandUseAlb);
  BinaryWrite(rStr, AvVegHeight);
  BinaryWrite(rStr, AvOptTransmCoeff);
  BinaryWrite(rStr, AvStomRes);
  BinaryWrite(rStr, AvVegFraction);
  BinaryWrite(rStr, AvLeafAI);
  // CJC2025: New Parameters
  BinaryWrite(rStr, AvEvapThresh);
  BinaryWrite(rStr, AvTransThresh);

  BinaryWrite(rStr, cumHrsSun); // Snow
  BinaryWrite(rStr, cumLHF);
  BinaryWrite(rStr, cumSHF);
  BinaryWrite(rStr, cumSnSub); // Write snowpack sublimation to restart file CJC2020
  BinaryWrite(rStr, cumSnEvap); // Write snowpack evaporation to restart file CJC2020
  BinaryWrite(rStr, cumTotEvap); // Write snowpack evaporation to restart file CJC2020
  BinaryWrite(rStr, cumBarEvap); // Write snowpack evaporation to restart file CJC2020
  BinaryWrite(rStr, cumPHF);
  BinaryWrite(rStr, cumRLin);
  BinaryWrite(rStr, cumRLout);
  BinaryWrite(rStr, cumRSin);
  BinaryWrite(rStr, cumGHF);
  BinaryWrite(rStr, cumMelt);
  BinaryWrite(rStr, cumIntSub);
  BinaryWrite(rStr, cumIntUnl);
  BinaryWrite(rStr, cumHrsSnow);
  
  // Added by Giuseppe Mascaro in 2016 to allow ingestion of soil grids
  BinaryWrite(rStr, Ks); 
  BinaryWrite(rStr, ThetaS);
  BinaryWrite(rStr, ThetaR);
  BinaryWrite(rStr, PoreSize);
  BinaryWrite(rStr, AirEBubPres);
  BinaryWrite(rStr, DecayF);
  BinaryWrite(rStr, SatAnRatio);
  BinaryWrite(rStr, UnsatAnRatio);
  BinaryWrite(rStr, Porosity);
  BinaryWrite(rStr, VolHeatCond);
  BinaryWrite(rStr, SoilHeatCap);

  int size;
  if (TimeInd != 0) {
    size = TimeInd->getSize();
    BinaryWrite(rStr, size);
    tListIter< int > IndIter;
    IndIter.Reset(*TimeInd);
    IndIter.First();
    while ( !(IndIter.AtEnd())) {
      BinaryWrite(rStr, IndIter.DatRef());
      IndIter.Next();
    }
  } else {
    size = 0;
    BinaryWrite(rStr, size);
  }

  if (Qeff != 0) {
    size = Qeff->getSize();
    BinaryWrite(rStr, size);
    tListIter< double > QIter;
    QIter.Reset(*Qeff);
    QIter.First();
    while ( !(QIter.AtEnd())) {
      BinaryWrite(rStr, QIter.DatRef());
      QIter.Next();
    }
  } else {
    size = 0;
    BinaryWrite(rStr, size);
  }
}

/***************************************************************************
**
** tCNode::readRestart() Function
**
***************************************************************************/

void tCNode::readRestart(fstream& rStr)
{
  BinaryRead(rStr, srf_hr);
  BinaryRead(rStr, cumsrf); //added CJC2021
  BinaryRead(rStr, RunOn);
  BinaryRead(rStr, VapPress);
  BinaryRead(rStr, ShortRadIn_dir);
  BinaryRead(rStr, ShortRadIn_dif);
  BinaryRead(rStr, ShortAbsbVeg);
  BinaryRead(rStr, ShortAbsbSoi);
  BinaryRead(rStr, Gnod);
  BinaryRead(rStr, VegFraction);

  BinaryRead(rStr, tracer);
  BinaryRead(rStr, flood);
  BinaryRead(rStr, soiID);
  BinaryRead(rStr, LandUse);
  BinaryRead(rStr, Reach);
  BinaryRead(rStr, Qin);
  BinaryRead(rStr, Qout);
  BinaryRead(rStr, traveltime);
  BinaryRead(rStr, hillpath);
  BinaryRead(rStr, streampath);
  BinaryRead(rStr, GridET);
  BinaryRead(rStr, ContrArea);
  BinaryRead(rStr, Curvature);
  BinaryRead(rStr, BasinArea);
  BinaryRead(rStr, BedrockDepth);
  BinaryRead(rStr, hFlux);
  BinaryRead(rStr, QgwIn);
  BinaryRead(rStr, QgwOut);
  BinaryRead(rStr, Aspect);
  BinaryRead(rStr, Width);
  BinaryRead(rStr, Roughness);

  BinaryRead(rStr, satOccur);
  BinaryRead(rStr, hsrfOccur);
  BinaryRead(rStr, percOccur); //ASM
  BinaryRead(rStr, avPerc); //ASM
  BinaryRead(rStr, psrfOccur);
  BinaryRead(rStr, satsrfOccur);
  BinaryRead(rStr, sbsrfOccur);
  BinaryRead(rStr, RechDisch);
  BinaryRead(rStr, NwtOld);
  BinaryRead(rStr, MuOld);
  BinaryRead(rStr, MiOld);
  BinaryRead(rStr, NtOld);
  BinaryRead(rStr, NfOld);
  BinaryRead(rStr, RuOld);
  BinaryRead(rStr, RiOld);
  BinaryRead(rStr, Qpout);
  BinaryRead(rStr, Rain);
  BinaryRead(rStr, intstorm);
  BinaryRead(rStr, Interception);
  BinaryRead(rStr, NetPrecipitation);
  BinaryRead(rStr, CanStorage);
  BinaryRead(rStr, PotEvaporation);
  BinaryRead(rStr, ActEvaporation);
  BinaryRead(rStr, StormLength);
  BinaryRead(rStr, CumIntercept);
  BinaryRead(rStr, EvapWetCanopy);
  BinaryRead(rStr, EvapDryCanopy);
  BinaryRead(rStr, EvapSoil);
  BinaryRead(rStr, EvapoTranspiration);
  BinaryRead(rStr, SoilMoisture);
  BinaryRead(rStr, SoilMoistureSC);
  BinaryRead(rStr, SoilMoistureUNSC);
  BinaryRead(rStr, RootMoisture);
  BinaryRead(rStr, RootMoistureSC);
  BinaryRead(rStr, Transmissivity);
  BinaryRead(rStr, AirTemp);
  BinaryRead(rStr, DewTemp);
  BinaryRead(rStr, RelHumid);
  BinaryRead(rStr, SurfTemp);
  BinaryRead(rStr, SoilTemp);
  BinaryRead(rStr, WindSpeed);
  BinaryRead(rStr, SkyCover);
  BinaryRead(rStr, NetRad);
  BinaryRead(rStr, AirPressure);
  BinaryRead(rStr, ShortRadIn);
  BinaryRead(rStr, LongRadIn);
  BinaryRead(rStr, LongRadOut);
  BinaryRead(rStr, gFlux);
  BinaryRead(rStr, lFlux);
  BinaryRead(rStr, AvSoilMoisture);
  BinaryRead(rStr, AvEvapFract);
  BinaryRead(rStr, AvET);
  BinaryRead(rStr, Hlevel);
  BinaryRead(rStr, Qstrm);
  BinaryRead(rStr, FlowVelocity);
  BinaryRead(rStr, CanopyStorVol);
  BinaryRead(rStr, UnSaturatedStorage);
  BinaryRead(rStr, SaturatedStorage);
  BinaryRead(rStr, Recharge);
  BinaryRead(rStr, UnSatFlowIn);
  BinaryRead(rStr, UnSatFlowOut);
  BinaryRead(rStr, ChannelPerc); //ASM 2/10/2017
  BinaryRead(rStr, Ft); //ASM

  BinaryRead(rStr, liqWEq); // Snowpack
  BinaryRead(rStr, iceWEq);
  BinaryRead(rStr, dU);
  BinaryRead(rStr, Unode);
  BinaryRead(rStr, Uerror);
  BinaryRead(rStr, cumUError);
  BinaryRead(rStr, liqRoute);
  BinaryRead(rStr, snTemperC);
  BinaryRead(rStr, crAge);
  BinaryRead(rStr, densAge);
  BinaryRead(rStr, ETage);
  BinaryRead(rStr, snLHF);
  BinaryRead(rStr, snSHF);
  BinaryRead(rStr, snGHF);
  BinaryRead(rStr, snPHF);
  BinaryRead(rStr, snRLin);
  BinaryRead(rStr, snRLout);
  BinaryRead(rStr, snRSin);
  BinaryRead(rStr, persTime);
  BinaryRead(rStr, persTimeTemp);
  BinaryRead(rStr, peakSWE);
  BinaryRead(rStr, peakSWEtemp);
  BinaryRead(rStr, initPackTime);
  BinaryRead(rStr, initPackTimeTemp);
  BinaryRead(rStr, peakPackTime);

  BinaryRead(rStr, intSWEq); // Snow intercept
  BinaryRead(rStr, intSnUnload);
  BinaryRead(rStr, intSub);
  BinaryRead(rStr, intPrec);

  BinaryRead(rStr, horizonAngle0000); // shelter
  BinaryRead(rStr, horizonAngle0225);
  BinaryRead(rStr, horizonAngle0450);
  BinaryRead(rStr, horizonAngle0675);
  BinaryRead(rStr, horizonAngle0900);
  BinaryRead(rStr, horizonAngle1125);
  BinaryRead(rStr, horizonAngle1350);
  BinaryRead(rStr, horizonAngle1575);
  BinaryRead(rStr, horizonAngle1800);
  BinaryRead(rStr, horizonAngle2025);
  BinaryRead(rStr, horizonAngle2250);
  BinaryRead(rStr, horizonAngle2475);
  BinaryRead(rStr, horizonAngle2700);
  BinaryRead(rStr, horizonAngle2925);
  BinaryRead(rStr, horizonAngle3150);
  BinaryRead(rStr, horizonAngle3375);
  BinaryRead(rStr, sfact);
  BinaryRead(rStr, lfact);
  BinaryRead(rStr, VegFraction);

  BinaryRead(rStr, LandUseAlb); // Additions for landuse grids
  BinaryRead(rStr, ThroughFall);
  BinaryRead(rStr, VegHeight);
  BinaryRead(rStr, StomRes);
  BinaryRead(rStr, IntercepCoeff);
  BinaryRead(rStr, CanFieldCap);
  BinaryRead(rStr, DrainCoeff);
  BinaryRead(rStr, DrainExpPar);
  BinaryRead(rStr, OptTransmCoeff);
  BinaryRead(rStr, LeafAI);
  BinaryRead(rStr, CanStorParam);
  // CJC2025: New Parameters
  BinaryRead(rStr, EvapThresh);
  BinaryRead(rStr, TransThresh);

  BinaryRead(rStr, LandUseAlbInPrevGrid);
  BinaryRead(rStr, LandUseAlbInUntilGrid);
  BinaryRead(rStr, ThroughFallInPrevGrid);
  BinaryRead(rStr, ThroughFallInUntilGrid);
  BinaryRead(rStr, VegHeightInPrevGrid);
  BinaryRead(rStr, VegHeightInUntilGrid);
  BinaryRead(rStr, StomResInPrevGrid);
  BinaryRead(rStr, StomResInUntilGrid);
  BinaryRead(rStr, VegFractionInPrevGrid);
  BinaryRead(rStr, VegFractionInUntilGrid);
  BinaryRead(rStr, CanStorParamInPrevGrid);
  BinaryRead(rStr, CanStorParamInUntilGrid);
  BinaryRead(rStr, IntercepCoeffInPrevGrid);
  BinaryRead(rStr, IntercepCoeffInUntilGrid);
  BinaryRead(rStr, CanFieldCapInPrevGrid);
  BinaryRead(rStr, CanFieldCapInUntilGrid);
  BinaryRead(rStr, DrainCoeffInPrevGrid);
  BinaryRead(rStr, DrainCoeffInUntilGrid);
  BinaryRead(rStr, DrainExpParInPrevGrid);
  BinaryRead(rStr, DrainExpParInUntilGrid);
  BinaryRead(rStr, OptTransmCoeffInPrevGrid);
  BinaryRead(rStr, OptTransmCoeffInUntilGrid);
  BinaryRead(rStr, LeafAIInPrevGrid);
  BinaryRead(rStr, LeafAIInUntilGrid);
  // CJC2025: New Parameters
  BinaryRead(rStr, EvapThreshInPrevGrid);
  BinaryRead(rStr, EvapThreshInUntilGrid);
  BinaryRead(rStr, TransThreshInPrevGrid);
  BinaryRead(rStr, TransThreshInUntilGrid);

  BinaryRead(rStr, AvCanStorParam);
  BinaryRead(rStr, AvIntercepCoeff);
  BinaryRead(rStr, AvThroughFall);
  BinaryRead(rStr, AvCanFieldCap);
  BinaryRead(rStr, AvDrainCoeff);
  BinaryRead(rStr, AvDrainExpPar);
  BinaryRead(rStr, AvLandUseAlb);
  BinaryRead(rStr, AvVegHeight);
  BinaryRead(rStr, AvOptTransmCoeff);
  BinaryRead(rStr, AvStomRes);
  BinaryRead(rStr, AvVegFraction);
  BinaryRead(rStr, AvLeafAI);
  // CJC2025: New Parameters
  BinaryRead(rStr, AvEvapThresh);
  BinaryRead(rStr, AvTransThresh);

  BinaryRead(rStr, cumHrsSun); // Snow
  BinaryRead(rStr, cumLHF);
  BinaryRead(rStr, cumSHF);
  BinaryRead(rStr, cumSnSub); // Read snowpack sublimation from restart file CJC2020
  BinaryRead(rStr, cumSnEvap); // Read snowpack evaporation from restart file CJC2020
  BinaryRead(rStr, cumTotEvap); // Read total evaporation from restart file CJC2020
  BinaryRead(rStr, cumBarEvap); // Read soil evaporation from restart file CJC2020
  BinaryRead(rStr, cumPHF);
  BinaryRead(rStr, cumRLin);
  BinaryRead(rStr, cumRLout);
  BinaryRead(rStr, cumRSin);
  BinaryRead(rStr, cumGHF);
  BinaryRead(rStr, cumMelt);
  BinaryRead(rStr, cumIntSub);
  BinaryRead(rStr, cumIntUnl);
  BinaryRead(rStr, cumHrsSnow);
  
  // Added by Giuseppe Mascaro in 2016 to allow ingestion of soil grids
  BinaryRead(rStr, Ks); 
  BinaryRead(rStr, ThetaS);
  BinaryRead(rStr, ThetaR);
  BinaryRead(rStr, PoreSize);
  BinaryRead(rStr, AirEBubPres);
  BinaryRead(rStr, DecayF);
  BinaryRead(rStr, SatAnRatio);
  BinaryRead(rStr, UnsatAnRatio);
  BinaryRead(rStr, Porosity);
  BinaryRead(rStr, VolHeatCond);
  BinaryRead(rStr, SoilHeatCap);

  int size;
  int timeInd;
  BinaryRead(rStr, size);
  for (int i = 0; i < size; i++) {
    BinaryRead(rStr, timeInd);
    TimeInd->insertAtBack(timeInd);
  }

  BinaryRead(rStr, size);
  double qeff;
  for (int i = 0; i < size; i++) {
    BinaryRead(rStr, qeff);
    Qeff->insertAtBack(qeff);
  }
}

/***************************************************************************
**
** tCNode::printVariables() Function
**
***************************************************************************/

void tCNode::printVariables()
{
  cout << " satOccur " << satOccur << endl;
  cout << " tracer " << tracer << endl;
  cout << " flood " << flood << endl;
  cout << " soiID " << soiID << endl;
  cout << " LandUse " << LandUse << endl;
  cout << " Reach " << Reach << endl;
  cout << " hsrfOccur " << hsrfOccur << endl;
  cout << " percOccur " << percOccur <<endl; //ASM
  cout << " avPerc " << avPerc <<endl; //ASM
  cout << " psrfOccur " << psrfOccur << endl;
  cout << " satsrfOccur " << satsrfOccur << endl;
  cout << " sbsrfOccur " << sbsrfOccur << endl;
  cout << " RechDisch " << RechDisch << endl;
  cout << " NwtOld " << NwtOld << endl;
  cout << " MuOld " << MuOld << endl;
  cout << " MiOld " << MiOld << endl;
  cout << " NtOld " << NtOld << endl;
  cout << " NfOld " << NfOld << endl;
  cout << " RuOld " << RuOld << endl;
  cout << " RiOld " << RiOld << endl;
  cout << " Qin " << Qin << endl;
  cout << " Qout " << Qout << endl;
  cout << " Qpout " << Qpout << endl;
  cout << " Rain " << Rain << endl;
  cout << " intstorm " << intstorm << endl;
  cout << " traveltime " << traveltime << endl;
  cout << " hillpath " << hillpath << endl;
  cout << " streampath " << streampath << endl;
  cout << " Interception " << Interception << endl;
  cout << " NetPrecipitation " << NetPrecipitation << endl;
  cout << " CanStorage " << CanStorage << endl;
  cout << " PotEvaporation " << PotEvaporation << endl;
  cout << " ActEvaporation " << ActEvaporation << endl;
  cout << " StormLength " << StormLength << endl;
  cout << " CumIntercept " << CumIntercept << endl;
  cout << " EvapWetCanopy " << EvapWetCanopy << endl;
  cout << " EvapDryCanopy " << EvapDryCanopy << endl;
  cout << " EvapSoil " << EvapSoil << endl;
  cout << " EvapoTranspiration " << EvapoTranspiration << endl;
  cout << " SoilMoisture " << SoilMoisture << endl;
  cout << " SoilMoistureSC " << SoilMoistureSC << endl;
  cout << " SoilMoistureUNSC " << SoilMoistureUNSC << endl;
  cout << " RootMoisture " << RootMoisture << endl;
  cout << " RootMoistureSC " << RootMoistureSC << endl;
  cout << " Transmissivity " << Transmissivity << endl;
  cout << " AirTemp " << AirTemp << endl;
  cout << " DewTemp " << DewTemp << endl;
  cout << " RelHumid " << RelHumid << endl;
  cout << " SurfTemp " << SurfTemp << endl;
  cout << " SoilTemp " << SoilTemp << endl;
  cout << " WindSpeed " << WindSpeed << endl;
  cout << " SkyCover " << SkyCover << endl;
  cout << " NetRad " << NetRad << endl;
  cout << " AirPressure " << AirPressure << endl;
  cout << " GridET " << GridET << endl;
  cout << " ContrArea " << ContrArea << endl;
  cout << " Curvature " << Curvature << endl;
  cout << " BasinArea " << BasinArea << endl;
  cout << " BedrockDepth " << BedrockDepth << endl;
  cout << " ShortRadIn " << ShortRadIn << endl;
  cout << " LongRadIn " << LongRadIn << endl;
  cout << " LongRadOut " << LongRadOut << endl;
  cout << " gFlux " << gFlux << endl;
  cout << " hFlux " << hFlux << endl;
  cout << " lFlux " << lFlux << endl;
  cout << " QgwIn " << QgwIn << endl;
  cout << " QgwOut " << QgwOut << endl;
  cout << " Aspect " << Aspect << endl;
  cout << " AvSoilMoisture " << AvSoilMoisture << endl;
  cout << " AvEvapFract " << AvEvapFract << endl;
  cout << " AvET " << AvET << endl;
  cout << " Hlevel " << Hlevel << endl;
  cout << " Qstrm " << Qstrm << endl;
  cout << " Width " << Width << endl;
  cout << " Roughness " << Roughness << endl;
  cout << " FlowVelocity " << FlowVelocity << endl;
  cout << " CanopyStorVol " << CanopyStorVol << endl;
  cout << " UnSaturatedStorage " << UnSaturatedStorage << endl;
  cout << " SaturatedStorage " << SaturatedStorage << endl;
  cout << " Recharge " << Recharge << endl;
  cout << " UnSatFlowIn " << UnSatFlowIn << endl;
  cout << " UnSatFlowOut " << UnSatFlowOut << endl;
  cout << " ChannelPerc " << ChannelPerc << endl; //ASM 2/10/2017
  cout << " Ft " << Ft << endl; //ASM

  cout << " liqWEq " << liqWEq; // Snowpack
  cout << " iceWEq " << iceWEq;
  cout << " dU " << dU;
  cout << " Unode " << Unode;
  cout << " Uerror " << Uerror;
  cout << " cumUError " << cumUError;
  cout << " liqRoute " << liqRoute;
  cout << " snTemperC " << snTemperC;
  cout << " crAge " << crAge;
  cout << " densAge " << densAge;
  cout << " ETage " << ETage;
  cout << " snLHF " << snLHF;
  cout << " snSHF " << snSHF;
  cout << " snGHF " << snGHF;
  cout << " snPHF " << snPHF ;
  cout << " snRLin " << snRLin;
  cout << " snRLout " << snRLout ;
  cout << " snRSin " << snRSin;
  cout << " persTime " << persTime;
  cout << " persTimeTemp " << persTimeTemp;
  cout << " peakSWE " << peakSWE;
  cout << " peakSWEtemp " << peakSWEtemp;
  cout << " initPackTime " << initPackTime;
  cout << " initPackTimeTemp " << initPackTimeTemp ;
  cout << " peakPackTime " << peakPackTime;

  cout << " intSWEq " << intSWEq; // Snow intercept
  cout << " intSnUnload " << intSnUnload;
  cout << " intSub " << intSub;
  cout << " intPrec " << intPrec;

  cout << " horizonAngle0000 " << horizonAngle0000; // shelter
  cout << " horizonAngle0225 " << horizonAngle0225;
  cout << " horizonAngle0450 " << horizonAngle0450;
  cout << " horizonAngle0675 " << horizonAngle0675;
  cout << " horizonAngle0900 " << horizonAngle0900;
  cout << " horizonAngle1125 " << horizonAngle1125;
  cout << " horizonAngle1350 " << horizonAngle1350;
  cout << " horizonAngle1575 " << horizonAngle1575;
  cout << " horizonAngle1800 " << horizonAngle1800;
  cout << " horizonAngle2025 " << horizonAngle2025;
  cout << " horizonAngle2250 " << horizonAngle2250;
  cout << " horizonAngle2475 " << horizonAngle2475;
  cout << " horizonAngle2700 " << horizonAngle2700;
  cout << " horizonAngle2925 " << horizonAngle2925;
  cout << " horizonAngle3150 " << horizonAngle3150;
  cout << " horizonAngle3375 " << horizonAngle3375;
  cout << " sfact " << sfact;
  cout << " lfact " << lfact;
  cout << " VegFraction " << VegFraction;

  cout << " LandUseAlb " << LandUseAlb; // Additions for landuse grids
  cout << " ThroughFall " << ThroughFall;
  cout << " VegHeight " << VegHeight;
  cout << " StomRes " << StomRes;
  cout << " IntercepCoeff " << IntercepCoeff;
  cout << " CanFieldCap " << CanFieldCap;
  cout << " DrainCoeff " << DrainCoeff;
  cout << " DrainExpPar " << DrainExpPar;
  cout << " OptTransmCoeff " << OptTransmCoeff;
  cout << " LeafAI " << LeafAI;
  cout << " CanStorParam " << CanStorParam;
  // CJC2025: New Parameters
  cout << " EvapThresh " << EvapThresh;
  cout << " TransThresh " << TransThresh;

  cout << " LandUseAlbInPrevGrid " << LandUseAlbInPrevGrid;
  cout << " LandUseAlbInUntilGrid " << LandUseAlbInUntilGrid;
  cout << " ThroughFallInPrevGrid " << ThroughFallInPrevGrid;
  cout << " ThroughFallInUntilGrid " << ThroughFallInUntilGrid;
  cout << " VegHeightInPrevGrid " << VegHeightInPrevGrid;
  cout << " VegHeightInUntilGrid " << VegHeightInUntilGrid;
  cout << " StomResInPrevGrid " << StomResInPrevGrid;
  cout << " StomResInUntilGrid " << StomResInUntilGrid;
  cout << " VegFractionInPrevGrid " << VegFractionInPrevGrid;
  cout << " VegFractionInUntilGrid " << VegFractionInUntilGrid;
  cout << " CanStorParamInPrevGrid " << CanStorParamInPrevGrid;
  cout << " CanStorParamInUntilGrid " << CanStorParamInUntilGrid;
  cout << " IntercepCoeffInPrevGrid " << IntercepCoeffInPrevGrid;
  cout << " IntercepCoeffInUntilGrid " << IntercepCoeffInUntilGrid;
  cout << " CanFieldCapInPrevGrid " << CanFieldCapInPrevGrid;
  cout << " CanFieldCapInUntilGrid " << CanFieldCapInUntilGrid;
  cout << " DrainCoeffInPrevGrid " << DrainCoeffInPrevGrid;
  cout << " DrainCoeffInUntilGrid " << DrainCoeffInUntilGrid;
  cout << " DrainExpParInPrevGrid " << DrainExpParInPrevGrid;
  cout << " DrainExpParInUntilGrid " << DrainExpParInUntilGrid;
  cout << " OptTransmCoeffInPrevGrid " << OptTransmCoeffInPrevGrid;
  cout << " OptTransmCoeffInUntilGrid " << OptTransmCoeffInUntilGrid;
  cout << " LeafAIInPrevGrid " << LeafAIInPrevGrid;
  cout << " LeafAIInUntilGrid " << LeafAIInUntilGrid;
  // CJC2025: New Parameters
  cout << " EvapThreshInPrevGrid " << EvapThreshInPrevGrid;
  cout << " EvapThreshInUntilGrid " << EvapThreshInUntilGrid;
  cout << " TransThreshInPrevGrid " << TransThreshInPrevGrid;
  cout << " TransThreshInUntilGrid " << TransThreshInUntilGrid;

  cout << " AvCanStorParam " << AvCanStorParam;
  cout << " AvIntercepCoeff " << AvIntercepCoeff;
  cout << " AvThroughFall " << AvThroughFall;
  cout << " AvCanFieldCap " << AvCanFieldCap;
  cout << " AvDrainCoeff " << AvDrainCoeff;
  cout << " AvDrainExpPar " << AvDrainExpPar;
  cout << " AvLandUseAlb " << AvLandUseAlb;
  cout << " AvVegHeight " << AvVegHeight;
  cout << " AvOptTransmCoeff " << AvOptTransmCoeff;
  cout << " AvStomRes " << AvStomRes;
  cout << " AvVegFraction " << AvVegFraction;
  cout << " AvLeafAI " << AvLeafAI;
  // CJC2025: New Parameters
  cout << " AvEvapThresh " << AvEvapThresh;
  cout << " AvTransThresh " << AvTransThresh;

  cout << " cumHrsSun " << cumHrsSun; // Snow
  cout << " cumLHF " << cumLHF;
  cout << " cumSHF " << cumSHF;
  cout << " cumSnSub " << cumSnSub; // CJC2020
  cout << " cumSnEvap " << cumSnEvap; // CJC2020
  cout << " cumPHF " << cumPHF;
  cout << " cumRLin " << cumRLin;
  cout << " cumRLout " << cumRLout;
  cout << " cumRSin " << cumRSin;
  cout << " cumGHF " << cumGHF;
  cout << " cumMelt " << cumMelt;
  cout << " cumIntSub " << cumIntSub;
  cout << " cumIntUnl " << cumIntUnl;
  cout << " cumHrsSnow " << cumHrsSnow;
    
  // Added by Giuseppe Mascaro in 2016 to allow ingestion of soil grids
  cout << " Ks " << Ks; 
  cout << " ThetaS " << ThetaS;
  cout << " ThetaR " << ThetaR;
  cout << " PoreSize " << PoreSize;
  cout << " AirEBubPres " << AirEBubPres;
  cout << " DecayF " << DecayF;
  cout << " SatAnRatio " << SatAnRatio;
  cout << " UnsatAnRatio " << UnsatAnRatio;
  cout << " Porosity " << Porosity;
  cout << " VolHeatCond " << VolHeatCond;
  cout << " SoilHeatCap " << SoilHeatCap;
}


//=========================================================================
//
//
//                 	   End of tCNode.cpp
//
//
//=========================================================================
