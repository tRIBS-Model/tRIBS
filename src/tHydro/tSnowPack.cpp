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
**  tSnowPack.cpp:   Function file for tSnowPack class (see tSnowPack.h)
**
***************************************************************************/

#include "src/tHydro/tSnowPack.h"
#include "src/Headers/globalIO.h"

//===========================================================================
//
//		Section 1: Constructor and Initialization Routines
//
//===========================================================================

//---------------------------------------------------------------------------
//
//  			tSnowPack() Constructor and Destructor
//
//  Construct the tSnowPack object. This class inherits from tEvapoTrans but
//  does not require anymore inputed parameters than tEvapoTrans. Additional
//  parameters are initialized using SetSnowVariables and SetSnowPack variables.
//  The logic in having two set-functions is to allow the easy construction of
//  other snow-type classes that require similar information as all other snow
//  classes (e.g., thermal properties, energy balance information, SWE, ...).
//
//---------------------------------------------------------------------------

tSnowPack::tSnowPack() {

}

tSnowPack::tSnowPack(SimulationControl *simCtrPtr, tMesh<tCNode> *gridRef,
                     tInputFile &infile, tRunTimer *t, tResample *resamp,
                     tHydroModel *hydro, tRainfall *storm)
        : tEvapoTrans(simCtrPtr, gridRef, infile, t, resamp, hydro, storm) {

    gridPtr = gridRef;
    //timerET = t;
    timer = t; // SKY2008Snow
    simCtrl = simCtrPtr;

    //set variables
    SetSnowVariables(infile);
    SetSnowPackVariables(infile);
    SetSnowInterceptVariables();

}

tSnowPack::~tSnowPack() {
    Cout << "tSnowPack Object has been destroyed..." << endl;
}

//---------------------------------------------------------------------------
//
//	tSnowPack::SetSnowPackVariables()
//
//	Auxiliary function used in robust constructor to initialize snow 
//	variables
//
//---------------------------------------------------------------------------

void tSnowPack::SetSnowPackVariables(tInputFile &infile) {

    //parameters
    //minSnTemp = infile.ReadItem(minSnTemp, "MINSNTEMP");
    //snliqfrac = infile.ReadItem(snliqfrac, "SNLIQFRAC"); // Added by CJC 2020
    hillAlbedoOption = infile.ReadItem(hillAlbedoOption, "HILLALBOPT");
    ETAge = 0.0;
    compactParam = 0.3;
    //rhoSnFreshkg = 100.0;
    snOnOff = 0.0;
}

void tSnowPack::SetSnowInterceptVariables() {
    Qcs = Ce = I = psiS = 0.0;
    Imax = prec = LAI = 0.0;
    RH = D = rhoVap = Omega = 0.0;
    Sh = Nu = Re = 0.0;
    kc = 0.010; //-
    iceRad = 500e-6; //m
    Mwater = 18.01; //kg/kmol
    R = 8313.0; //J/kmol K
    RdryAir = 287.0; //J/kg K
    nu = 1.3e-5; //m^2/s
    KtAtm = 0.024; //J/msK WR updated value from Liston and Elder
    esatIce = 0.0;
    beta = 0.9;
    acoefficient = 0.0;
    Lm = 0.0;
}

void tSnowPack::SetSnowVariables(tInputFile &infile) {
    

    //time steps
    timeStepm = infile.ReadItem(timeStepm, "METSTEP");
    timeSteph = timeStepm / 60.0;
    timeSteps = 60.0 * timeStepm;
    minutelyTimeStep = 0.0;

    //state variables
    liqWE = iceWE = snWE = 0.0;
    liqWEm = iceWEm = snWEm = 0.0;
    Utot = Usn = Uwat = 0.0;
    liqWatCont = 0.0;
    liqTempC = iceTempC = snTempC = 0.0;
    liqTempK = iceTempK = snTempK = 0.0;
    crustAge = 0.0;
    albedo = 0.88;
    canWE = 0.0;
    Iold = Utotold = 0;

    //fluxes
    H = L = G = Prec = Rn = 0.0;
    snPrec = liqPrec = 0.0;
    snPrecm = liqPrecm = 0.0;
    snPrecmm = liqPrecmm = 0.0;
    vapPressSmb = vapPresskSPa = 0.0;
    snSub = snEvap = liqRoute = 0.0;
    RSin = 0;

    //density
    rholiqcgs = 1.0;
    rhoicecgs = 0.92;
    rhosncgs = 0.1;
    rholiqkg = 1000.0; //kg/m^3
    rhoicekg = 920.0; //kg/m^3
    rhoAir = 1.3;

    //thermal properties
    cpicekJ = 2.1;
    cpwaterkJ = 4.190;
    cpairkJ = 1.006;
    latFreezekJ = 334.0;
    latVapkJ = 2470.0;
    latSubkJ = latFreezekJ + latVapkJ;

    //output variables
    snDepth = snDepthm = 0.0;
    snOnOff = 0.0;
    peakSnWE = peakSnWEtemp = 0.0;
    persMax = persMaxtemp = 0.0;
    inittime = peaktime = 0.0;

    //conversions
    naughttokilo = 1e-3;
    kilotonaught = 1e3;
    cgsRHOtomks = 1e3;
    mksRHOtocgs = 1e-3;
    naughttocm = 100.0; // Used convert m to cm
    cmtonaught = 0.01; // Used convert cm to m
    ctom = 10.0;
    mtoc = 0.1;

    // Set hardcoded default snow parameters
    rhoSnFreshkg    = 60.0;   // kg/m3
    snliqfrac       = 0.04;    // 4% irreducible water content
    kSatRef         = 0.005;   // m/s
    minSnTemp       = -30.0;   // C
    canopyWindAtten = 0.5;
    snRoughness     = 0.001;   // m
    snInitialAlbedo = 0.85;
    snLambdaDry     = 0.96;
    snLambdaWet     = 0.82;
    snMinAlbedo     = 0.4;
    albResetThresh  = 5.0;     // 5mm to reset albedo
    optPartition    = 0;       // Default to Linear Transition
    maxTempWB       = 5.0;     // Default for Wet Bulb
    minTempRain     = 0.0;     // Default for Linear
    maxTempSnow     = 2.0;     // Default for Linear

    // Read Input Parameters from snow parameter file
    // Read the path to the snow parameter file from the main input file
    char snowParamFile[kName]; 

    // If user defined snow parameter file in input then use that
    if (infile.IsItemIn("SNOWFILENAME")) {

        infile.ReadItem(snowParamFile, "SNOWFILENAME"); // the existense of this file has already been check in tPreProcess.cpp
        // Create a new tInputFile object specifically for the snow parameters
        tInputFile snowIn(snowParamFile); 

        // Read the parameters using keywords defined in file.
        // Density
        rhoSnFreshkg = snowIn.ReadItem(rhoSnFreshkg, "FRESH_SNOW_DENSITY");

        // General
        snliqfrac = snowIn.ReadItem(snliqfrac,  "IRREDUCIBLE_SAT"); 
        kSatRef   = snowIn.ReadItem(kSatRef,    "K_SAT_REF");
        minSnTemp = snowIn.ReadItem(minSnTemp,  "MIN_SNOW_TEMP");
        canopyWindAtten = snowIn.ReadItem(canopyWindAtten, "CANOPY_WIND_ATTENUATION");
        snRoughness = snowIn.ReadItem(snRoughness, "ROUGHNESS_LENGTH");

        // Albedo
        snInitialAlbedo = snowIn.ReadItem(snInitialAlbedo, "ALBEDO_FRESH");
        snLambdaDry     = snowIn.ReadItem(snLambdaDry,     "ALBEDO_DECAY_DRY");
        snLambdaWet     = snowIn.ReadItem(snLambdaWet,     "ALBEDO_DECAY_WET");
        snMinAlbedo     = snowIn.ReadItem(snMinAlbedo,     "ALBEDO_MIN");
        albResetThresh  = snowIn.ReadItem(albResetThresh,  "ALBEDO_RESET_THRESHOLD");

        // Precipitation Partitioning
        optPartition = snowIn.ReadItem(optPartition, "OPTPRECPARTITION");
        if (optPartition == 0) { // Wet-bulb method
            // Check input actually exists
            if (!snowIn.IsItemIn("MAX_WETBULB_TEMP")) {
                cout << "\nERROR: You selected Precipitation Partition Option 0 (Wet Bulb)," << endl;
                cout << "       but the parameter 'MAX_WETBULB_TEMP' is missing from the file." << endl;
                exit(1);
            }

            maxTempWB = snowIn.ReadItem(maxTempWB, "MAX_WETBULB_TEMP");
        } else if (optPartition == 1) { // Linear transition method
            // Check input actually exists
            if (!snowIn.IsItemIn("MIN_TEMP_RAIN")) {
                cout << "\nERROR: You selected Precipitation Partition Option 1 (Linear Transition)," << endl;
                cout << "       but the parameter 'MIN_TEMP_RAIN' is missing from the file." << endl;
                exit(1);
            } else if (!snowIn.IsItemIn("MAX_TEMP_SNOW")) {
                cout << "\nERROR: You selected Precipitation Partition Option 1 (Linear Transition)," << endl;
                cout << "       but the parameter 'MAX_TEMP_SNOW' is missing from the file." << endl;
                exit(1);
            }

            minTempRain = snowIn.ReadItem(minTempRain, "MIN_TEMP_RAIN");
            maxTempSnow = snowIn.ReadItem(maxTempSnow, "MAX_TEMP_SNOW");
        } else { // Invalid option
            cout << "\nPrecipitation Partition Option " << optPartition;
            cout << " not valid." << endl;
            cout << "\tPlease use :" << endl;
            cout << "\t\t(0) for Wet Bulb Threshold" << endl;
            cout << "\t\t(1) for Linear Transition" << endl;
            cout << "\nExiting Program..." << endl << endl;
            exit(1);
        }
    }

    // Initialize the density state variable to the parameter value
    rhosnkg = rhoSnFreshkg; 
}

//---------------------------------------------------------------------------
//
//			tSnowPack::callSnowPack()
//
//  Calling function for snow pack dymamics. This function has the dual role
//  of representing both snow physics and evapotranspiration in the warm 
//  landscape. This is because we have to check to see if there is snow either
//  on the ground or coming down from the canopy or atmosphere before we can
//  tell if if we need to deal with snow or not.
//
//  The structure of this function is similar to that of 
//  tEvapoTrans::callEvapoPotential and tEvapoTrans::callEvapoTrans. It begins
//  by setting the sun variables and preparing the stochastic weather 
//  simulator if necessary. Then, the mesh list is initialized and the loop
//  through the node list begins.
//
//  At this point, the slope, aspect and elevation of every node is computed,
//  the appropriate sheltering algorithm is implemented, the appropriate met
//  and rainfall data is called, the land-surface is reinitialized, and the 
//  current state of the land-surface is found (snow or no-snow). If lapse
//  rates are implemented, then the air temperature and/or the rainfall is 
//  adjusted for elevation.
//
//  At this point, the actual physics begin. If the interception scheme is
//  on and there is canopy, then tSnowIntercept::callSnowIntercept() is
//  called and the snow-canopy interaction is accounted for.
//
//  Then existence of snow pack, snow precipitation or snow throughfall is 
//  checked for. If it does not exist, then tEvapoTrans::callEvapoPotential 
//  and tEvapoTrans::callEvapoTrans are recreated. RMK: THIS ASSUMES THAT 
//  ETISTEP AND METSTEP ARE THE SAME, WHICH HAS BEEN DONE IMPLICITLY THROUGH
//  MOST OF THE EXISTING MODEL.
//
//  If it does exist, or if there is a positive mass flux, we enter the actual
//  snow-physics portion of the code. Here, we begin by dealing with the 
//  atmospheric snow mass balance (adding/subtracting precip, throughfall,
//  turbulent latent heat flux). We have to make the distinction between a 
//  developed and undeveloped snow pack in order to accurately represent the
//  heat flux. If there is no snow on the ground, then there is, in a sense,
//  no precipitation heat flux. The pack's state is initialized to that of the
//  air temperature and precipitation phase masses. Otherwise, we deal with 
//  both the mass and energy precipitation fluxes, and initialize the energy
//  of the pack using the updated mass and the previous temperature.
//
//  RMK: IF SOMEONE REWORKS THIS PORTION OF THE MODEL, RESEQUENCING THE MB
//  AND EB WOULD BE A GOOD IDEA, IN ORDER TO ELIMINATE SOME OF THE ERROR THAT
//  WILL COME FROM THE CURRENT SEQUENCING. UNFORTUNATELY, I COULDN'T FIGURE
//  OUT A WAY TO HAVE THE EB OCCUR DURING THE ACCUMULATION PERIOD WITHOUT THIS
//  SEQUENCING.
//
//  We are now discussing the energy balance portion of the code. If there is
//  no snow left after the initial mass balance, then we do not actually 
//  compute the energy balance. This is to ensure stability of temperatures.
//
//  If there is snow left, then we initialize the internal energy of the pack,
//  compute the change in energy, add the change in energy, and compute the
//  new state of the model. If the internal energy is greater than 0J/m^2, then
//  there is liquid water in the pack. Otherwise, we compute the temperature
//  of the snow pack (which will be less than 0).
//
//  RMK: WE ASSUME A SINGLE EQUIVALENT TEMPERATURE IS REPRESENTATIVE, BUT WE
//  KNOW THAT IT IS NOT.
//
//  If there the amount of liquid mass is greater than 40% of the solid mass,
//  then the excess water is routed out of the pack and is treated as 
//  precipitation in tHydroModel::UnSaturatedZone(). That's it. It seems so
//  simple now that its done.
//
//  09 July 2007 -- AJR @ New Mexico Tech
//
//
//				    Pertinent references: Wigmosta et al (1994)
//							  Tarboton and Luce (1996)
//							  Tuteja et al (1996)
//							  Marks et al (1999)
//							  Anderson (1976)
//							  Jordan (1991)
//				    
//---------------------------------------------------------------------------
// SKY2008Snow, AJR2008
void tSnowPack::callSnowPack(tIntercept *Intercept, int flag) {

    tCNode *cNode;
    tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());

    int count = 0;
    int cnt = 0;
    double EP = 0.0; //double tmp = 0.0;
    double SkyC = 0.0; //double tmpC = 0.0;
    double vegHeight = 0;

    // SKY2008Snow, AJR2008
    //  metHour = metStep;
    //  etHour = etStep;

    if (simCtrl->Verbose_label == 'Y') {
        cout << "\nSnowPack Routine Call ..." << endl;
    }

    // Set time, sun, and meteorological variables -- AJR2008, SKY2008Snow
    SetEnvironment();

    // Resample Meteorological Grids, if option
    if (metdataOption == 2) {
        resampleGrids(timer);
    }


    if (luOption == 1) { // resampling Land Use grids done here, i.e., dynamic case
        if (AtFirstTimeStepLUFlag) {
            initialLUGridAssignment();
            AtFirstTimeStepLUFlag = 0;
        } else {
            LUGridAssignment();
        }
    }

    // BEGIN LOOP THROUGH NODES
    cNode = nodeIter.FirstP();
    while (nodeIter.IsActive()) {
        double precip = 0.0;

        // CJC2025 this block was overwriting land use parameters when it is not needed
        /* landPtr->setLandPtr(cNode->getLandUse());
        cNode->setCanStorParam(landPtr->getLandProp(1));
        cNode->setIntercepCoeff(landPtr->getLandProp(2));
        cNode->setThroughFall(landPtr->getLandProp(3));
        cNode->setCanFieldCap(landPtr->getLandProp(4));
        cNode->setDrainCoeff(landPtr->getLandProp(5));
        cNode->setDrainExpPar(landPtr->getLandProp(6));
        cNode->setLandUseAlb(landPtr->getLandProp(7));
        cNode->setVegHeight(landPtr->getLandProp(8));
        cNode->setOptTransmCoeff(landPtr->getLandProp(9));
        cNode->setStomRes(landPtr->getLandProp(10));
        cNode->setVegFraction(landPtr->getLandProp(11));
        cNode->setLeafAI(landPtr->getLandProp(12)); */

        // Elapsed MET steps from the beginning, used for averaging dynamic LU grid values below over time for integ. output
        auto te = (double) timer->getElapsedMETSteps(timer->getCurrentTime());
        integratedLUVars(cNode, te);

        ID = cNode->getID();

        //Get Rainfall
        rain = cNode->getRain(); // get new rainfall

        //Set Elevation, Slope and Aspect
        slope = fabs(atan(cNode->getFlowEdg()->getSlope()));
        aspect = cNode->getAspect();
        elevation = cNode->getZ();
        snOnOff = 0.0;

        checkShelter(cNode);

        // Load the full set of default parameters from the table.
        // This provides the fallback values for any parameters not supplied as grids.
        setCoeffs(cNode);

        // If using gridded data (dynamic or static), overwrite the defaults.
        if (luOption == 1 || luOption == 2) {

            // For DYNAMIC grids (Option 1), update/interpolate node values first.
            if (luOption == 1) {
                if (luInterpOption == 1) {
                    interpolateLUGrids(cNode);
                } else {
                    constantLUGrids(cNode);
                }
            }
            // For STATIC grids (Option 2), values are already on the node from initialization, so no action is needed here.
            // For both gridded options, load the parameters from the node into the local 'coeff' variables.
            newLUGridData(cNode);
        }

        //updates meteorological variables if not in stochastic mode
        if (metdataOption == 1) {
            thisStation = assignedStation[count];
            newHydroMetData(hourlyTimeStep); //read in met data from station file -- inherited function

            if (fabs(skyCover-9999.99)<1.0E-3){ // assumed if first value is 9999.99 the rest are
                skycover_flag =1;
            }

        } else if (metdataOption == 2) {
            newHydroMetGridData(cNode); // set up and get appropriate data -- inherited function

            if (fabs(skyCover-9999.99)<1.0E-3){ // work around since nodata from grids only set to 9999.99 once in tvariannt
                skycover_flag =1;
            }
        }

        // Set the observed values to the node:
        // they will be required by other function calls
        vPress = vaporPress();


        // Check/modify cloud cover values
        if (skycover_flag == 1) {
            skyCover = compSkyCover();
        }

        cNode->setAirTemp(airTemp); // celsius
        cNode->setDewTemp(dewTemp);
        cNode->setRelHumid(rHumidity);
        cNode->setVapPressure(vPress);
        cNode->setWindSpeed(windSpeed);
        cNode->setAirPressure(atmPress);
        cNode->setSkyCover(skyCover);
        cNode->setShortRadIn(inShortR);

        //Set Soil/Surface Temperature
        if (hourlyTimeStep == 0) {
            cNode->setSoilTemp(Tlo - 273.15);
            cNode->setSurfTemp(Tso - 273.15);
        }

        if (Ioption == 0) {
            cNode->setNetPrecipitation(rain);
        }

        //Call Beta functions
        betaFunc(cNode); // inherited
        betaFuncT(cNode); // inherited

        // Get Soil/Surface Temperature --WR debug 01032024 same setup as callEvapPotential.
         Tso = cNode->getSurfTemp() + 273.15;
         Tlo = cNode->getSoilTemp() + 273.15;

        //get the necessary information from tCNode for snow model
        getFrNodeSnP(cNode);

        // ensure routed liquid is reset
        cNode->setLiqRouted(0.0);

        snUnload = 0.0;
        canWE = cNode->getIntSWE();

        //No Snow on ground and canopy and not snowing
        if ((snWE <= 1e-4) && (rain * snowFracCalc() <= 5e-2) && rholiqkg * cmtonaught * (cNode->getIntSWE()) < 1e-3) {

            // Following block of code mirrors callEvapoPotential and callEvapoTrans in tEvapoTrans as no snow occurs at any
            // level of the system: i.e. snowpack, canopy, or snowing. —refactored by WR 6/21/23

            //Calculate the Potential and Actual Evaporation
            if (evapotransOption == 1) {
                EvapPenmanMonteith(cNode); // SKY2008Snow
            } else if (evapotransOption == 2) {
                EvapPan();
            } else {
                cout << "\nEvapotranspiration Option " << evapotransOption;
                cout << " not valid." << endl;
                cout << "\tPlease use :" << endl;
                cout << "\t\t(1) for Penman-Monteith Method" << endl;
                cout << "\t\t(2) for Pan Evaporation Measurements" << endl;
                cout << "Exiting Program...\n\n" << endl;
                exit(1);
            }
            // Set
            setToNode(cNode);
            ComputeETComponents(Intercept, cNode, count, flag);

            snTempC = 0.0; // reinitialize snTemp
            snWE = 0.0; // reinitialize snWE
            snSub = 0.0; // No sublimation occurs CJC2020
            snEvap = 0.0; // No evaporation occurs CJC2020
            dUint = RLin = RLout = RSin = H = L = G = Prec = 0.0; //reinitialize energy terms
            ETAge = ETAge + timeStepm;
            liqRoute = 0.0;
            iceWE = 0.0;
            liqWE = 0.0;
            Utot = 0.0;
            Usn = 0.0;
            Uwat = 0.0;
            snTempC = 0.0;
            crustAge = 0.0;

            cNode->setIntSub(0.0);
            setToNodeSnP(cNode);
        }//end no-snow

        else //condtions include some combination of snowpack and snow in canopy, and snowing, rain on snow, or no precip
        {

            // Implement interception schemes for snow, refactored WR 6/21/23
            if (Ioption && Intercept->IsThereCanopy(cNode)) { // && coeffV>0
                callSnowIntercept(cNode, Intercept, count);
                snUnload = cNode->getIntSnUnload(); //calculated in callSnowIntercept() units in cm
                snCanWE = cNode->getIntSWE();//units in cm
            } else {
                cNode->setNetPrecipitation(rain);
            }

            precip = cNode->getNetPrecipitation(); //units in mm
            precip += snUnload * ctom; // units in mm

            //change mass (volume) quantities to correct units (kJ, m, C, s)
            iceWE = iceWE * cmtonaught; // mm to m
            liqWE = liqWE * cmtonaught; // mm to m
            snWE = iceWE + liqWE; // mm to m

            // Here precip is being set by net precipitation which is set from callSnowIntercept
            // and represents throughfall + unloaded snow  from snow interception (so scaled by coeffV) and precipitation
            // that falls on the no vegetated fraction of the cell (1-coeffV). refactored WR 6/21/23

            // BEGIN Calculate current snow depth for use in the turbulent heat flux calculations and output.
            // We don't yet know the new state of the new mass balance variables so we will project what they will be
            // in order to calculate the new snow depth.

            // Calculate Incoming Mass (temporary variables)
            // precip is in mm. Convert to meters (0.001).
            double newSnowM = (precip * snowFracCalc()) * 0.001; 
            double newRainM = (precip * (1.0 - snowFracCalc())) * 0.001;

            // Calculate projected Mass States
            // Approximation of what the pack will look like after this timestep.
            // We use these only to calculate depth/roughness right now.
            double projIceWE = iceWE + newSnowM;       // Total Ice Matrix
            double projTotWE = snWE + newSnowM + newRainM; // Total Weight (for compaction)

            // Weighted Average Mixing
            // We have to do this step because of the single layer snow model.
            // Only mix if we actually have new solid snow.
            if (newSnowM > 1e-6) {
                double freshRho = freshDensityCalc(airTemp); // e.g. 80 kg/m3
                
                double oldIceKg = iceWE * 1000.0;
                double newIceKg = newSnowM * 1000.0;
                
                // Conservation of Volume
                double vol_old = oldIceKg / rhosnkg;
                double vol_new = newIceKg / freshRho;
                
                // New Mixed Density
                rhosnkg = (oldIceKg + newIceKg) / (vol_old + vol_new);
            }

            // Apply Compaction 
            // We pass projTotWE because the weight of the liqwe contributes to squashing the snow.
            if (projTotWE > 1e-6) {
                rhosnkg = densityCompactionCalc(rhosnkg, projTotWE * 1000.0, snTempC, timeSteps);
            }
            
            // Safety Clamps on Density
            if (rhosnkg < 50.0) rhosnkg = 50.0;
            if (rhosnkg > 600.0) rhosnkg = 600.0; // Max density for dry snow

            // Calculate Final Depth
            // We divide by rhosnkg so liquid water doesn't artificially inflate depth.
            snDepthm = (projIceWE * 1000.0) / rhosnkg;

            // END Calculate current snow depth

            //account for veg height
            if (coeffH == 0) {
                vegHeight = 0.1;
            } else {
                vegHeight = coeffH;
            }

            if (airTemp > 0.0) {
                snTempC = 0.0;
            } else {
                snTempC = airTemp;
            }

            if (snWE < 1e-5) {

                //no precipitation heat flux, as it is totally accounted for in the snow pack energy
                //   initialization
                phfOnOff = 0.0;

                //reinitialize crust age
                crustAge = 0.0;


                //snowMB
                // evaporate liquid from ripe pack snWE +=
                // note evaporation/sublimation flux can be negative or positive, with latter representing condensation/deposition
                if (snTempC == 0.0) {
                    updateRipeSnowPack(precip);
                }
                    //sublimate solid from frozen pack
                else {
                    updateSolidSnowPack(precip);
                }

                snWE = iceWE + liqWE; // unit in meters here
                snSub *= naughttocm; // m to cm
                snEvap *= naughttocm; // m to cm

                //set other fluxes
                L = H = G = Prec = Utotold = 0.0;

                // We assume the density of the ice matrix (rhosnkg) did not change during melt, 
                // only the amount of ice changed. 
                if (rhosnkg > 0) {
                        snDepthm = (iceWE * 1000.0) / rhosnkg;
                } else {
                        snDepthm = 0.0;
                }

                // Calculate the Max Liquid Capacity (in meters)
                // We treat snliqfrac as the Volumetric Irreducible Saturation.
                // Capacity = Depth * Porosity * snliqfrac
                double maxPoreStorageM = snDepthm * porosityCalc(rhosnkg) * snliqfrac; 

                //put in routing bucket
                if (liqWE > maxPoreStorageM) { 
                    //there is enough water left over
                    if (liqWE != snWE) {
                        // Calculate drainage using Darcy's Law approximation
                        liqRoute = drainageCalc(liqWE, snDepthm, rhosnkg, snliqfrac, timeSteps);
                        liqWE = liqWE - liqRoute;
                        snWE = liqWE + iceWE;
                    }
                        //there is no more pack
                    else {
                        liqRoute = snWE;
                        liqWE = 0.0;
                        iceWE = 0.0;
                        snWE = 0.0;
                    }
                }
                if (liqWE < 0) {    // caused by snVap, XYT2023
                    liqWE = 0;
                    snWE = iceWE;
                }

                //initialize and record energy balance
                Utot = dUint = iceWE * rhoicekg * cpicekJ * snTempC // Changed to use rhoicekg CJC 2020
                               + liqWE * rholiqkg * latFreezekJ;

            } else {

                //account for precipitation heat flux
                phfOnOff = 1.0;

                //reset crust age if snowing out
                if (precip * snowFracCalc() > albResetThresh) {
                    crustAge = 0.0;
                }

                //snowMB
                //ripe pack -- evaporate water
                if (snTempC == 0.0) {
                    updateRipeSnowPack(precip);
                }
                    //sublimate solid from frozen pack
                else {
                    updateSolidSnowPack(precip);
                }

                snWE = iceWE + liqWE;// unit in meters here
                snSub *= naughttocm;// m to cm
                snEvap *= naughttocm;// m to cm

                //snowEB
                ETAge = 0.0;
                L = latentHFCalc(resFactCalc());
                Prec = precipitationHFCalc();

                //if there is no snow left at this point, then bail out of
                //energy balance.
                if ((snWE <= 5e-6) || (snTempC < -800.0)) {
                    liqRoute = 0.0;
                    snWE = 0.0;
                    iceWE = 0.0;
                    liqWE = 0.0;
                    Utot = 0.0;
                    Usn = 0.0;
                    Uwat = 0.0;
                    snTempC = 0.0;
                    crustAge = 0.0;
                } else {

                    //find initial state of energy
                    Utot = Utotold = iceWE * rhoicekg * cpicekJ * snTempC + liqWE * rholiqkg *
                                                                            latFreezekJ; // I am pretty sure this should be rhoicekg CJC 2020
                    //adjust albedo for age
                    albedo = agingAlbedo();

                    //calculate dU
                    snowEB(ID, cNode); // AJR2008, SKY2008Snow

                    //check for balance
                    Uerr = (Utot - Utotold) - dUint;

                    if (Utot < 0.0) { //frozen pack -- change temperature
                        Usn = Utot;// all energy in solid phase
                        Uwat = 0.0;// no energy in liquid phase
                        liqWatCont = 0.0;// no liquid content
                        liqWE = 0.0;// no liq WE
                        liqTempC = 0.0; // reset liq temp to default

                        iceWE = snWE;

                        //calculate sn temperature, modified THM 2012
                        if (iceWE < 0.1 && iceWE > 0) {
                            iceTempC = Usn / (cpicekJ * rhoicekg * 0.1);
                        } else {
                            iceTempC = Usn / (cpicekJ * rhoicekg * snWE);
                        }

                        //adjust to minimum snow temperature
                        //	RMK: THIS IS A KLUGE THAT IS NECESSARY B/C OF THE
                        //	ONE-LAYER ASSUMPTION. IT IS ALSO REQUIRED B/C OF THE
                        //	SIMPLISTIC WAY WE MODEL GROUND HEAT FLUX. SOMEONE
                        //	NEEDS TO INCORPORATE MULTIPLE LAYERS.
                        if (iceTempC <= minSnTemp) {
                            iceTempC = minSnTemp;
                        }

                        //set pack temperature to ice temperature
                        snTempC = iceTempC;

                    }//end -- frozen pack

                    else {//melt
                        Uwat = Utot;
                        Usn = 0.0;

                        liqWE = Uwat / (latFreezekJ * rholiqkg);

                        //make sure that there is enough SWE in the pack for the melt
                        if (liqWE >= snWE) {
                            liqWE = snWE; // this is here because the liqWE += term above
                        }                  //  can result in liqWE > snWE
                        //assign water equivalents
                        iceWE = snWE - liqWE; 

                        // We assume the density of the ice matrix (rhosnkg) did not change during melt, 
                        // only the amount of ice changed. 
                        if (rhosnkg > 0) {
                                snDepthm = (iceWE * 1000.0) / rhosnkg;
                        } else {
                                snDepthm = 0.0;
                        }

                        // Calculate the Max Liquid Capacity (in meters)
                        // We treat 'snliqfrac' as the Volumetric Irreducible Saturation.
                        // Capacity = Depth * Porosity * snliqfrac
                        double maxPoreStorageM = snDepthm * porosityCalc(rhosnkg) * snliqfrac; 

                        //put in routing bucket
                        if (liqWE > maxPoreStorageM) { 
                            //there is enough water left over
                            if (liqWE != snWE) {
                                // Calculate drainage using Darcy's Law approximation
                                liqRoute = drainageCalc(liqWE, snDepthm, rhosnkg, snliqfrac, timeSteps);
                                liqWE = liqWE - liqRoute;
                                snWE = liqWE + iceWE;
                            }
                                //there is no more pack
                            else {
                                liqRoute = snWE;
                                liqWE = 0.0;
                                iceWE = 0.0;
                                snWE = 0.0;
                            }
                        }
                        //set temperatures to 0 Celsius
                        snTempC = 0.0;
                        iceTempC = 0.0;
                        liqTempC = 0.0;
                    }//end -- melt
                }//end -- snow left after initial mass decrement
            }//end -- existing pack at beginning of time step

            //make sure that we still have snow
            if (snWE <= 5e-6) {
                liqRoute += snWE;
                snWE = 0.0;
                liqWE = 0.0;
                iceWE = 0.0;
                crustAge = 0.0;
                Utot = 0.0;
                Usn = 0.0;
                Uwat = 0.0;
            } else {
                crustAge += timeSteph;
                snOnOff = 1.0;
            }

            //mass balance leaves >= 0 snow, then prepare for output in cm
            snWE = naughttocm * snWE; // m to cm
            iceWE = naughttocm * iceWE; // m to cm
            liqWE = naughttocm * liqWE; // m to cm
            liqRoute = naughttocm * liqRoute; // m to cm


            // Set ET variables equal to zero due to snowpack
            // ET variables are set to zero for canopy when snow in canopy, see tSnowIntercept
            cNode->setEvapSoil(0.0);
            cNode->setActEvap(0.0);
            // This maybe redundant.This setEvapoTrans is called in ComputeETComponets, which
            // is only called when no snow is in the system or in CallSnowIntercept, if there is
            // no snow in the canopy.
            cNode->setEvapoTrans(cNode->getEvapWetCanopy() + cNode->getEvapDryCanopy());



            setToNodeSnP(cNode);
            //setToNode(cNode); // WR 01032024this also being set in callSnowIntercept, may be source of variation in AtmPress?

        }//end yes-snow

        //get the next node information for while loop
        cNode = nodeIter.NextP();
        count++;

    }//end while-nodes
    timeCount++;
    oldTimeStep = hourlyTimeStep;
    hourlyTimeStep++;
}

//---------------------------------------------------------------------------
//
//		tSnowIntercept::callSnowIntercept()
//
//    Calls the physical algorithms from tSnow::callSnowPack(). Some of
//    tIntercept::callIntercept() is implemented for the case when there is
//    no snow.
//
//---------------------------------------------------------------------------

void tSnowPack::callSnowIntercept(tCNode *node, tIntercept *interceptModel, int count) {
    double CanStorage;
    double subFrac, unlFrac, precip, Isnow, throughfall, scover;// SKY2008Snow, AJR2008
    int flag;
    flag = 1;
    CanStorage = node->getCanStorage();

    //set meteorolgical conditions
    rHumidity = node->getRelHumid();
    vPress = node->getVapPressure();
    dewTemp = node->getDewTemp();
    windSpeed = node->getWindSpeed();
    atmPress = node->getAirPressure();
    inShortR = node->getShortRadIn();
    airTemp = node->getAirTemp();
    airTempK = CtoK(airTemp);

    precip = node->getRain(); // precip should not be scaled by vegfrac CJC 2025
    // VF scaling here does not match tIntercept logic
    LAI = coeffLAI;

    //reinitialize snow interception model
    Iold = rholiqkg * cmtonaught * (node->getIntSWE());
    Qcs = 0.0;
    Lm = 0.0;

    if ((precip * snowFracCalc() < 1e-4) && (Iold < 1e-3)) {
        //The below code block account for the case where there is no snow in canopy and it's not snowing
        // but could be raining (i.e. rain on snow). Basically this emulates callEvapPotential and callEvapoTrans,
        // This is necessary to simulate potential evaporation and subsequently evaporation from the canopy. Terms that are related
        // to soil are reset to zero on the node.

        I = Iold = Lm = Qcs = 0.0;

        //Calculate the Potential and Actual Evaporation
        if (evapotransOption == 1) {
            EvapPenmanMonteith(node); // call to get EvapPot, but energy balance for soil es
        } else if (evapotransOption == 2) {
            EvapPan();
        } else {
            cout << "\nEvapotranspiration Option " << evapotransOption;
            cout << " not valid." << endl;
            cout << "\tPlease use :" << endl;
            cout << "\t\t(1) for Penman-Monteith Method" << endl;
            cout << "\t\t(2) for Pan Evaporation Measurements" << endl;
            cout << "Exiting Program...\n\n" << endl;
            exit(1);
        }


        // Set ground element-scale fluxes to zero and surface and soil to snow temp
        node->setNetRad(0.0);
        node->setGFlux(0.0);
        node->setHFlux(0.0);
        node->setLFlux(0.0);
        node->setLongRadOut(0.0);
        node->setSurfTemp(node->getSnTempC()); //assume surface temp == snow temp, note if snWE < 1e-4 snow intercept should not be called
        node->setSoilTemp(node->getSnTempC()); //this should be updated with n-layer snow model

        // follows structure of
        setToNode(node);

        // Set actual evaporation to 0 since snow pack exists and soil evaporation is function of actual evap
        node->setActEvap(0.0);

        ComputeETComponents(interceptModel, node, count, 1);

        //Set canopy snow components to 0
        node->setIntSWE(0);
        node->setIntSnUnload(0);
        node->setIntSub(0);
        node->setIntPrec(0);
        node->setEvapSoil(0.0);


    }//end -- no snow

        //snowing with or without snow in canopy
    else {


        //albedo = 0.8; WR debug this is set elsewhere and should be double checked

        //Note no actual unit conversion is necessary from mm to kg/m^2\
        //Here mm values (i.e. canopy storage, precip) are assumed to be converted to kg/m^2

        //Add CanStorage to I_old and reset to node canopy storage to zero
        if (CanStorage > 1e-5) {
            Iold += CanStorage; //CanStorage has been scaled by coeffV, through scaling of precip
            node->setCanStorage(0.0);
            flag = 0; // WR refactor, initial setup did not compute Qcs and Lm on first event of snow
        }

        //maximum mass of snow stored in canopy (kg/m^2)
        Imax = 4.4 * LAI;

        //compute new intercepted snow (kg/m^2)
        Isnow = 0.7 * (Imax - Iold) * (1 - exp(-precip / Imax));
        I = Iold + Isnow;

        //precip minus intercepted snow (i.e. throughfall)
        throughfall = precip - Isnow;

        //if there was old snow, sublimate and unload
        if (Iold > 0.0 && flag) {
            computeSub();
            computeUnload();
        } else {
            Qcs = 0.0;//sublimation term
            Lm = 0.0;//unloading term
        }

        I += Qcs -
             Lm; //I == interception (kg), Qcs == sublimation (kg) (sign computed), Lm == unloading (computed positive) (kg)

        // SKY2008Snow based on AJR2008's recommendation starts here (water balance now preserved)
        if (I < 0.0) {

            if (Qcs < 0.0) {
                subFrac = fabs(Qcs) / (fabs(Qcs) + Lm);
                unlFrac = Lm / (fabs(Qcs) + Lm);
            } else {
                subFrac = 0.0;
                unlFrac = 1.0;
            }
            Qcs -= I * subFrac;
            Lm += I * unlFrac;
            I = 0.0;
        }


        //adjust amount of snow in canopy
        Iold = I; //WR debug moved to below catch for I<0

        // SKY2008Snow based on AJR2008's recommendation ends here

        // Set adjusted fluxes and states to node
        // Flux variables represent the flux averaged over the entire voronoi cell (scaled by veg fraction)
        // State variables represent teh state of the canopy itself (un-scaled by veg fraction)
        node->setIntSWE(naughttocm * (1 / rholiqkg) * I); //length units in cm
        node->setIntSnUnload(naughttocm * (1 / rholiqkg) * Lm * coeffV);
        node->setIntSub(naughttocm * (1 / rholiqkg) * Qcs * coeffV);
        node->addIntSub(naughttocm * (1 / rholiqkg) * Qcs * coeffV);
        node->addIntUnl(naughttocm * (1 / rholiqkg) * Lm * coeffV);
        node->setIntPrec(Isnow * (1 / rholiqkg) * naughttocm);
        // Rain rate for the _ENTIRE_ cell:
        double total_net_precip = (throughfall * coeffV) + (node->getRain() * (1 - coeffV));
        node->setNetPrecipitation(total_net_precip);
        // note mm and kg/m^2 requires no conversion

        //set wet and dry evap to 0 when snow in canopy
        node->setEvapWetCanopy(0.0);
        node->setEvapDryCanopy(0.0);
        node->setEvapSoil(0.0);
        node->setEvapoTrans(0.0);
        node->setPotEvap(0.0);
    }//end -- snow exists

    return;
}

/****************************************************************************
**
**		      tSnowPack -- Interact w/ tCNode
**
** Functions to update snowpack in callsnowpack
**
****************************************************************************/

void tSnowPack::updateRipeSnowPack(double precip) {
    //liq WE update
    snEvap = (1 - coeffV) * latentHFCalc(resFactCalc()) * timeSteps /
             (rholiqkg * latVapkJ); // units should be in meters
    liqWE += cmtonaught * ((mtoc * precip) * (1 - snowFracCalc())) * timeSteps /
             3600; // Removed snUnload term CJC2020

    if (liqWE + snEvap <= 0) {
        snEvap = liqWE;
        liqWE = 0;
    } else {
        liqWE += snEvap;
    }

    //solid WE update
    iceWE += cmtonaught * (mtoc * (precip * snowFracCalc())) * timeSteps / 3600;
    snSub = 0.0; // No sublimation occurs CJC2020
}

void tSnowPack::updateSolidSnowPack(double precip) {
    //liq WE update
    liqWE += cmtonaught * (mtoc * (precip * (1 - snowFracCalc()))) * timeSteps /
             3600; // Removed snUnload term CJC2020
    snEvap = 0.0; // No evaporation occurs CJC2020


    //ice WE update
    snSub = (1 - coeffV) * latentHFCalc(resFactCalc()) * timeSteps /
            (rholiqkg * latSubkJ);// units should be in meters
    iceWE += cmtonaught * (mtoc * (precip * snowFracCalc())) * timeSteps / 3600;

    if (iceWE + snSub <= 0) {
        snSub = iceWE;
        iceWE = 0;
    } else {
        iceWE += snSub;
    }
}
/****************************************************************************
**
**		      tSnowPack -- Interact w/ tCNode
**
**	This section of functions do the bulk on interacting w/ tCNode for
**	the tSnowPack class.
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//			      tSnowPack::getFrNodeSnP()
//    
//	Get last time step's state variable information from tCNode object.
//
//---------------------------------------------------------------------------

void tSnowPack::getFrNodeSnP(tCNode *node) {

    liqWE = node->getLiqWE(); //cm
    iceWE = node->getIceWE(); //cm
    liqRoute = 0.0; //cm
    snWE = liqWE + iceWE; //cm

    // Deal with snow density
    rhosnkg = node->getRhoSn();
    // If there is no snow, density should be initialized to fresh snow
    if (snWE <= 1e-5) {
        rhosnkg = rhoSnFreshkg; // Uses the parameter you read from the file
    }

    //deal with the snow temperatures
    if (snWE > 1e-4) {

        snTempC = node->getSnTempC(); //Celsius
        iceTempC = snTempC; //Celsius
        liqTempC = 0.0; //default -- if this is 0.0, then the other should have been
        //		 set to 0.0 during the last time step.

    } else if (airTemp > 0.0) {
        //in case it is snowing out and the air temperature is greater than 0.0 (this is
        //actually a possibility). This also assumes that there is no pack.
        snTempC = 0.0;
        iceTempC = 0.0;
        liqTempC = 0.0;
    } else if (airTemp <= 0.0) {
        //if the air temperature is less than zero and there is no pack, assume
        //that the temperature of the phases are all the air temp.
        snTempC = airTemp;
        iceTempC = airTemp;
        liqTempC = airTemp;
    }

    crustAge = node->getCrustAge();
    ETAge = node->getEvapoTransAge();

    persMax = node->getPersTimeMax();
    persMaxtemp = node->getPersTime();

    peakSnWE = node->getPeakSWE();
    peakSnWEtemp = node->getPeakSWETemp();

    inittime = node->getInitPackTime();
    inittimeTemp = node->getInitPackTimeTemp();
    peaktime = node->getPeakPackTime();

}

//---------------------------------------------------------------------------
//
//			  tSnowPack::setToNodeSnP()
//
//	Sets the state variables, fluxes and outputs to the node.
//
//---------------------------------------------------------------------------

void tSnowPack::setToNodeSnP(tCNode *node) {

    //state variables
    node->setLiqWE(liqWE);
    node->setIceWE(iceWE);
    node->setSnTempC(snTempC);
    node->setCrustAge(crustAge);
    node->setEvapoTransAge(ETAge);
    node->setSnSub(snSub); // scaled by non-vegetated area
    node->setSnEvap(snEvap); // scaled by non-vegetated area

    //mass flux
    node->setLiqRouted(liqRoute);

    //energy fluxes, changes
    node->setDU(dUint);
    node->setUnode(Utot);

    //outputs
    node->setUerror(Uerr);

    //times and peaks
    if ((iceWE + liqWE) <= 1e-5) {
        persMaxtemp = 0.0;
        inittimeTemp = 0.0;
    }

    if (snWE > 1e-5) {
        persMaxtemp += 1.0;

        if (persMaxtemp > persMax) {
            persMax = persMaxtemp;

            if (persMaxtemp - 1 < 1)
                inittimeTemp = hourlyTimeStep;

        }
    }

    if (snWE >= peakSnWE) {
        peakSnWE = snWE;
        peaktime = hourlyTimeStep;
        inittime = inittimeTemp;
    }

    // Save the updated density
    node->setRhoSn(rhosnkg);

    node->setPersTimeMax(persMax);
    node->setPersTime(persMaxtemp);
    node->setPeakSWE(peakSnWE);
    node->setPeakPackTime(double(peaktime));
    node->setInitPackTime(double(inittime));
    node->setInitPackTimeTemp(double(inittimeTemp));
    // Convert meters to centimeters for the output/node storage
    node->setSnDepth(snDepthm * 100.0); 

    //cumulative outputs
    node->addLatHF(L);
    node->addSnSub(snSub); // cumulative snow sublimation CJC2020
    node->addSnEvap(snEvap); // cumulative snow evaporation CJC2020
    node->addMelt(liqRoute);
    node->addSHF(H * timeSteps);
    node->addPHF(Prec * timeSteps);
    node->addGHF(G * timeSteps);
    node->addRLin(RLin * timeSteps);
    node->addRLout(RLout * timeSteps);
    node->addRSin(RSin * timeSteps);
    node->addCumUerror(Uerr * timeSteps);
    node->addCumHrsSnow(snOnOff);

    //reset fluxes to zero
    L = H = Prec = G = RLin = RLout = RSin = dUint = 0.0;
    liqRoute = 0.0;
}


/****************************************************************************
**
**		  tSnowPack -- Internal Physical Routines
**
**	Deals with purely internal changes of the pack (not at surface)
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//			  tSnowPack::freshDensityCalc()
//
//	  Calculates the density of newly fallen snow based on air temperature.
//	  This uses the empirical relationship derived for continental snowpacks.
//
//							 (Hedstrom-Pomeroy 1998)
//
//---------------------------------------------------------------------------

double tSnowPack::freshDensityCalc(double airTempC) {

    // Returns density in kg/m^3
    double rho_fresh = rhoSnFreshkg + 51.25 * exp(airTempC / 2.59);
    
    // Safety clamps from https://doi.org/10.1175/1520-0477(2000)081%3C1577:DOFFSI%3E2.3.CO;2
    if (rho_fresh < 50.0) rho_fresh = 50.0;
    if (rho_fresh > 200.0) rho_fresh = 200.0;
    
    return rho_fresh;
}

//---------------------------------------------------------------------------
//
//			  tSnowPack::densityCompactionCalc()
//
//	  Calculates the change in bulk snow density due to two processes:
//	  1. Destructive Metamorphism (crystal settling/rounding)
//	  2. Overburden Pressure (compaction due to weight)
//
//						(Jordan 1991 / Anderson 1976)
//
//---------------------------------------------------------------------------

double tSnowPack::densityCompactionCalc(double rho, double SWE_kgm2, double snTempC, double dt_sec) {

    // Constants from Jordan (1991) / Anderson (1976)
    // currentDensity in kg/m^3
    // currentSnWE_kgm2 in kg/m^2
    // snTempC in Celsius
    // dt_sec is timestep in seconds
    
    double density = rho;
    double T_kelvin = snTempC + 273.15;
    double T_diff = 273.15 - T_kelvin;
    double gravity = 9.81;
    
    if (snTempC < -20.0) return rho; // Too cold for significant compaction

    // Destructive Metamorphism (Settling)
    // Constants converted from hourly (0.01) to secondly (2.77e-6)
    double c1 = 2.777e-6; // s^-1 
    double c3 = 0.04;     // K^-1
    double c4 = 0.046;    // m^3/kg
    
    // Fractional rate (s^-1)
    double fractional_meta = c1 * exp(-c3 * T_diff) * exp(-c4 * (density - 100.0));
    
    // Overburden Pressure (Compaction)
    // Pressure at midpoint = 0.5 * SWE
    double P = 0.5 * SWE_kgm2 * gravity;
    
    // Viscosity (Pa s) - Jordan (1991)
    double eta0 = 3.6e6; 
    double a_eta = 0.08; 
    double b_eta = 0.021; 
    double eta = eta0 * exp(a_eta * T_diff + b_eta * density);
    
    double fractional_overburden = P / eta; // s^-1
    
    // Integration
    double total_rate = fractional_meta + fractional_overburden;
    
    // Use exponential growth: rho(t) = rho0 * exp(rate * dt)
    double new_density = density * exp(total_rate * dt_sec);
    
    // Safety Clamps
    if (new_density > 600.0) new_density = 600.0;
    
    return new_density;
}

//---------------------------------------------------------------------------
//
//			  tSnowPack::drainageCalc()
//
//	  Calculates the vertical liquid water flux (drainage) out of the 
//	  snowpack using a Darcy's Law approximation for gravity drainage.
//	  
//	  The function treats 'irredSat' as the Volumetric Irreducible 
//	  Water Saturation (fraction of pore space holding immobile water).
//	  Flux only occurs when liquid water exceeds this holding capacity.
//
//						(Colbeck 1972 / Jordan 1991)
//
//---------------------------------------------------------------------------

double tSnowPack::drainageCalc(double liqWE_m, double snDepth_m, double rho, double irredSat, double dt_sec) {
    
    // Sanity Check
    if (liqWE_m <= 0.0 || snDepth_m <= 0.0) return 0.0;

    // Calculate Porosity
    double porosity = porosityCalc(rho);

    // Calculate Effective Saturation
    double pore_volume = snDepth_m * porosity;
    double S_absolute = liqWE_m / pore_volume;
    double S_eff = 0.0;

    // Only drain water above the irreducible holding capacity
    if (S_absolute > irredSat) {
        S_eff = (S_absolute - irredSat) / (1.0 - irredSat);
    } else {
        return 0.0; // Water is held by capillary forces
    }

    // Calculate Flux (Darcy's Law approximation)
    // Flux = K_sat * (S_eff)^3, 3 is a theoretical value from Colbeck
    double drainage_flux_mps = kSatRef * S_eff * S_eff * S_eff;

    // Since timestep is hardcoded to 1 hour, Flux (m/hr) = Volume (m)
    double route_vol = drainage_flux_mps * dt_sec;

    // Mass Balance Limit
    // Cannot drain more water than is "mobile" (Liquid - Holding Capacity)
    double mobile_water = liqWE_m - (pore_volume * irredSat);
    
    // Clamp to 0 just in case floating point math goes slightly negative
    if (mobile_water < 0.0) mobile_water = 0.0;
    if (route_vol > mobile_water) {
        route_vol = mobile_water;
    }
    
    return route_vol;
}

//---------------------------------------------------------------------------
//
//        tSnowPack::porosityCalc()
//
//	  Calculate bulk porosity based on current snow density. Assumes 
//	  a pure ice grain density of 917 kg/m^3.
//	  
//	  Includes a safety clamp (min porosity 0.05) to prevent division 
//	  by zero errors in drainage calculations for highly compacted ice.
//
//---------------------------------------------------------------------------
double tSnowPack::porosityCalc(double rho) {

    // Calculate theoretical porosity
    double por = 1.0 - (rho / rhoicekg);

    // Safety Clamp
    // Prevent divide by zero in dense ice
    if (por < 0.05) { 
        por = 0.05; 
    }

    return por;
}

/****************************************************************************
**
**		    tSnowPack -- EB Functions (basic calcs)
**
**	  This set of functions forms the basis of the energy balance 
**	  calculations later. This includes all of the heat fluxes (H, L, 
**	  P, G, Rn) and the ancillary information needed to calculate these.
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//			  tSnowPack::latentHFCalc()
//
//	  Calculates the turbulent latent heat flux given a given aerodynamic
//	  resistivity.  There are two separate conditions. The first is 
//	  sublimation from dry snow and the second is evaporation of the 
//	  liquid phase. Other than that, it is straight forward (REF).
//
//	  RMK: DOES NOT TAKE INTO ACCOUNT ATMOSPHERIC INSTABILITY (DOES NOT
//	  INCORPORATE RICHARD'S NUMBER).
//
//							  (Dingman 2002)
//
//---------------------------------------------------------------------------

double tSnowPack::latentHFCalc(double Kaero) {

    double lhf;
    double snTemp_corr = snTempC;

    // If snowpack has liquid then make sure temp isn't below 0C
    if (liqWE > 1e-5) {
        snTemp_corr = 0.0;
    }

    // If the snow surface is at or above the melting point, the process is evaporation.
    if (snTemp_corr == 0.0) {
        // Use latent heat of vaporization.
        lhf = (latVapkJ * 0.622 * rhoAir * Kaero * (vPress - 6.111) / atmPress); //evaporation by THM 2012
    // If the snow surface is frozen, the process is sublimation.
    } else {
        // Use latent heat of sublimation.
        // Saturation Vapor Pressure over ICE
        lhf = (latSubkJ * 0.622 * rhoAir * Kaero * (vPress - 6.112 * exp((21.87 * snTemp_corr) / (snTemp_corr + 265.5))) /
               atmPress); //sublimation by THM 2012
    }
    return lhf;
}

//----------------------------------------------------------------------------
//
//			    tSnowPack::sensibleHFCalc()
//
//	  Calculates the turbulent sensible heat from from the surface of the
//	  snow pack given a aerodynamic resistivity.
//
//	  RMK: DOES NOT ACCOUNT FOR ATMOSPHERIC INSTABILITY (DOES NOT
//	  INCORPORATE RICHARDS NUMBER).
//
//							    (Dingman 2002)
//
//----------------------------------------------------------------------------

double tSnowPack::sensibleHFCalc(double Kaero) {

    double shf;
    double snTemp_corr = snTempC;

    // If snowpack has liquid then make sure temp isn't below 0C
    if (liqWE > 1e-5) {
        snTemp_corr = 0.0;
    }

    shf = (rhoAir * cpairkJ * Kaero * ((airTemp + 273.15) - (snTemp_corr + 273.15)));
    return shf;
}

//-----------------------------------------------------------------------------
//
//			    tSnowPack::snowFracCalc()
//
//	  Partitions precipitation between solid and liquid phase based on 
//	  air temperature . This is linear interpolation b/t two
//	  temperatures. Ice can exist above freezing and liquid can exist below
//	  freezing, but there are certain minimun and maximum temperatures at 
//	  which the precipitation becomes all ice and all liquid, respectively.
//
//							(Wigmosta et al 1994)
//
//-----------------------------------------------------------------------------

double tSnowPack::snowFracCalc() {

    double snowfrac;
    double Tw, RH, Ta, f1;
    Ta = airTemp;
    RH = rHumidity;
    if (optPartition == 0) { // Wet-bulb partitioning

        // Calculate wet-Bulb Temperature according to Stull (2011) https://doi.org/10.1175/JAMC-D-11-0143.1
        Tw = Ta*atan(0.151977*pow(RH + 8.313659,0.5)) + atan(Ta + RH) - atan(RH - 1.676331) +
            0.00391838*pow(RH,1.5)*atan(0.023101*RH) - 4.686035; // in degC

        // Calculate  snowfall fraction according to Wang et al. (2019) https://doi.org/10.1029/2019GL085722
        f1 = 1 + 6.99*pow(10,-5)*exp(2*(Tw + 3.97));

        // Wet-bulb temperature > 5C corresponds to Ta = 10C and RH = 50%, assuming no snowfall
        if ( Tw > maxTempWB ) {
            snowfrac = 0;
        } else {
            snowfrac = 1/f1;
        }
    } else { // Linear ramp partitioning

        if (Ta <= minTempRain) {
            snowfrac = 1; // all ice
        } else if (Ta >= maxTempSnow) {
            snowfrac = 0; // all liquid
        } else if ((Ta >= minTempRain) && (Ta <= maxTempSnow)) {
            snowfrac = (maxTempSnow - Ta) / (maxTempSnow - minTempRain); //mixture
        }
    }
    return snowfrac;
}

//------------------------------------------------------------------------------
//
//			tSnowPack::precipitationHFCalc()
//			
//	  Finds the amount of energy advected into the pack from precipitation.
//	  It has two cases:
//
//		Ta > 0	  ==>	Tsol = 0,   Tliq = Ta
//		Ta <= 0	  ==>	Tsol = Ta,  Tliq = 0
//
//	  In other words, if the air temperature is greater than zero, the ice
//	  just froze and the liquid is at ambient temperature, but if the air
//	  temperature is less than 0, then the ice is at the ambient temperature
//	  and the liquid is freezing. This difference is displayed in the phase
//	  of the precipitation has the heat capacity and temperature incorportated
//	  into the calculation.
//
//------------------------------------------------------------------------------

double tSnowPack::precipitationHFCalc() {
    double phf = 0;
    double frac;

    //  frac = snowFracCalc();
    snPrec = (snowFracCalc() * (rain + ctom * snUnload)) * mtoc; //convert from mm to cm
    liqPrec = ((1 - snowFracCalc()) * (rain + ctom * snUnload)) * mtoc; //convert from mm to cm

    if (airTemp > 0) {
        phf = (cmtonaught * snPrec * 0 * rholiqkg * cpicekJ +
               cmtonaught * liqPrec * (latFreezekJ + airTemp * rholiqkg * cpwaterkJ)) / 3600;

    } else if (airTemp <= 0) {
        phf = (cmtonaught * snPrec * airTemp * rholiqkg * cpicekJ +
               cmtonaught * liqPrec * latFreezekJ * rholiqkg) / 3600;

    }

    return phf;
}

//-----------------------------------------------------------------------------------
//
//				tSnowPack::agingAlbedo()
//
//	Returns the effective albedo of the surface for a given age of the snow
//	surface. This mainly is in response to crystal structure changes but also
//	heuristically deals with increased amounts of incorporated dust. It has 
//	two sets of curves. The first is during the accumulation period, assumed
//	to have predominantly dry snow; and the second is for during the melt period,
//	assumed to have wet snow.
//
//							      Wigmosta et al. (1994)
//							      CRREL Review of Albedo
//								  Modification.
//
//-----------------------------------------------------------------------------------

double tSnowPack::agingAlbedo() {
    double alb;

    if (liqWE < 1.0E-5 && snTempC < -0.1) {
        // Dry snow aging
        alb = snInitialAlbedo * pow(snLambdaDry, pow(crustAge / 24.0, 0.58));
    } else {
        // Wet snow aging
        alb = snInitialAlbedo * pow(snLambdaWet, pow(crustAge / 24.0, 0.46));
    }

    // Enforce the minimum albedo floor
    if (alb < snMinAlbedo) {
        alb = snMinAlbedo;
    }

    return alb;
}

//-----------------------------------------------------------------------------------
//
//				tSnowPack::resFactCalc()
//
//	  Returns the effective aerodynamic resistance for use in the calculations
//	  of turbulent heat flux. This accounts for the change in effective 
//	  vegetation roughness height b/c of snow depth. The rest of the function is
//	  taken directly from tEvapoTrans.
//
//-----------------------------------------------------------------------------------

double tSnowPack::resFactCalc() {
    const double vonKarm = 0.41;
    const double g = 9.81;
    double windSpeedC, windSpeedS; // windspeed 2 meters above canopy & snow/surface JB2025 @ASU
    double rf, ra;
    double vegHeight, vegFrac, vegBare;
    double zm, zom, zov, d, rav, ras;
	double Ri_cr = 0.1;  // Critical Richardson number, 0.2 is the most aggresive

    windSpeedC = (windSpeed == 0.0 || fabs(windSpeed - 9999.99) < 1e-3) ? 0.01 : windSpeed;

    // Vegetation height adjusted for snow
    vegHeight = (coeffH <= 0) ? 0.1 : coeffH;
    vegHeight = (vegHeight > snDepthm) ? vegHeight - snDepthm : 0.1;

    vegBare = 0.4;
    vegFrac = coeffV;

    // Compute below canopy windspeed at snow surface following equation Moreno et al. (2016) CJC 2020
    if (snDepthm < coeffH) {
        windSpeedS = windSpeedC * exp(-canopyWindAtten * coeffLAI * (1.0 - (snDepthm / coeffH)));
    } else {
        windSpeedS = windSpeedC;  // No canopy attenuation
    }
    
    // Aerodynamic resistance for vegetation
    zm = 2.0 + vegHeight;
    zom = 0.123 * vegHeight;
    zov = 0.0123 * vegHeight;
    d = 0.67 * vegHeight;

    rav = log((zm - d) / zom) * log((zm - d) / zov) / (windSpeedC * pow(vonKarm, 2)); //Uses canopy level wind speed

    // Aerodynamic resistance for snow surface
    zom = snRoughness;
    zov = snRoughness / 10.0; // 10% of the momentum roughness like with the original code
    vegBare = snRoughness / 0.123; // Back calculate effective obstacle height using same relation as original code
    zm = 2.0 + vegBare;
    //zom = 0.123 * vegBare;
    //zov = 0.0123 * vegBare;
    d = 0.67 * vegBare;

    ras = log((zm - d) / zom) * log((zm - d) / zov) / (windSpeedS * pow(vonKarm, 2)); // Uses snow/surface level wind

    // Weighted resistance
    ra = (1 - vegFrac) * ras + vegFrac * rav;

    // If snowpack has liquid then make sure snow temp isn't below 0C
    double Ts; // Snow temperature in K
    if (liqWE > 1e-5) {
        Ts = 273.15;
    } else {
        Ts = snTempC + 273.15; 
    }

    double Ta = airTemp + 273.15;    // Air temp in K
    double T_avg = 0.5 * (Ta + Ts);  // Mean temp
    // Effective vertical distance between reference air temperature (2 m AGL) and snow surface.
    // If snow depth exceeds 2 m, zm_eff becomes negative or zero, which is non-physical.
    // Clamp to minimum of 0.1 m to preserve numerical stability in RiB calculation.
    double zm_eff = std::max(0.1, 2.0 - snDepthm);
    double z0 = 0.123 * vegHeight;    // Roughness length, same as zom above

    // Considering atmospheric stability, Andreadis et al. (2009), JB2025 @ASU
    // https://doi.org/10.1029/2008WR007042

    // Bulk Richardson Number (Eq. 17)
    // Use windSpeedS (attenuated wind at snow surface) in RiB calc.
    // Though Ta is from 2m AGL, we prioritize the snow–air interface.
    // This matches the aerodynamic resistance formulation and maintains internal consistency.
    double RiB = (g * zm_eff * (Ta - Ts)) / (T_avg * windSpeedC* windSpeedC);

    // Compute upper limit Ri_u (Eq. 24)
    double Ri_u = 1.0 / (log(zm_eff / z0) + 5.0);

    // Stability correction factor
    double C_stab;
    if (RiB < 0.0) {
        // Unstable conditions (Eq. 19)
        C_stab = 1.0 / pow(1.0 - 16.0 * RiB, 0.5);
    } else if (RiB <= Ri_u) {
        // Stable but not extreme (Eq. 25)
        C_stab = 1.0 / pow(1.0 - (RiB / Ri_cr), 2.0);
    } else {
        // Very stable, capped (Eq. 26)
        C_stab = 1.0 / pow(1.0 - (Ri_u / Ri_cr), 2.0);
    }

	// Clamp the stability correction factor to a reasonable range
	// As RiB approaches Ri_cr, C_stab approaches infinity, so we limit it
	if (C_stab > 15.0) { C_stab = 15.0; } // kaero can at most be 15 times the original value

    // Apply correction to aerodynamic resistance
    ra *= C_stab;

    // Convert resistance to conductance
    rf = 1.0 / ra; // Otherwise known as kaero

    return rf;
}


/****************************************************************************
**
**		      tSnowIntercept -- Physical Routines
**
**	Functions that compute changes internal to the canopy for
**	tSnowIntercept::callSnowIntercept. A loading function should probably
**	be implemented in order to fully modulate the algorithm.
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//		tSnowIntercept::computeSub()
//
//    Uses timplementation to compute the snow sublimated
//    from the canopy after the snow has existed in the canopy for more than a
//    singe time step. All notation is taken from Liston and Sturm (2006).
//
//----------------------------------------------------------------------------

void tSnowPack::computeSub() {

    //compute incoming shortwave radiation
    // inShortR = inShortWaveCan();// Currently commented out since inShortWaveCan is non-functional.
    // InshortR was being set to zero here, because Isw = Iv * (1.0 - albedo), and Iv was zero. Also not clear if InShortR
    // was properly being reset as it's replaced by Isw which is consists of multiple components that need to be further evaluated.
    // Lastly it appears that albedo is already called below, making this redundant  Isw = Iv * (1.0 - albedo). This function
    // needs to be refactored or maybe deprecated.

    //compute effective incident shortwave radiation on snow crystal
    Sp = PI * pow(iceRad, 2.0) * (1 - albedo) * inShortR;//check units--check (W)//WR debug change 0.8 to albedo

    //Find coefficient for changing windspeed
    acoefficient = beta * coeffLAI;

    //find windspeed
    if (windSpeed == 0.0) {
        windSpeedC = 0.1; // WR 01032024 switched to windSpeedC since that is what is set to node.
    }
    windSpeedC = windSpeed * exp(-acoefficient * 0.4);// WR 01032024 switched to windSpeedC since that is what is set to node.

    //Calculate Reynolds number
    Re = 2 * iceRad * windSpeedC / nu;

    //Calculate Sherwood number
    Sh = 1.79 + 0.606 * pow(Re, 0.5);

    //Calculate Nusselt number
    Nu = Sh;

    //calculate saturated vapor pressure of at ice interface
    esatIce = 611.15 * exp(22.452 * (airTempK - 273.16) / (airTempK - 0.61)); //check units--check

    //calculate density of vapor
    rhoVap = 0.622 * esatIce / (RdryAir * airTempK);

    //compute vapor diffusivity
    D = 2.06e-5 * pow(airTempK / 273.0, 1.75);

    //Place holder in algorithm
    Omega = (1 / (KtAtm * airTempK * Nu)) * (1000 * latSubkJ * Mwater / (R * airTempK) - 1);//check units--check

    //find change of mass of ice crystal with respect to time
    dmdt = (2.0 * PI * iceRad * (rHumidityC / 100 - 1) - Sp * Omega) / //WR 01032024 switched to rHumidtyC since that is what is set to node.
           (1000 * latSubkJ * Omega + (1 / (D * rhoVap * Sh)));//1000 conversion from KJ to J

    //relative sublimation from ice sphere
    psiS = dmdt / ((4.0 / 3.0) * PI * rhoicekg * iceRad * iceRad * iceRad);

    //canopy exposure coefficient
    Ce = kc * pow(I / Imax, -0.4);

    //compute total sublimated snow during timestep
    Qcs = Ce * I * psiS * timeSteps;
}


//---------------------------------------------------------------------------
//
//			tSnowIntercept::computeUnload()
//
//	Compute the amount of unloading during a timestep according to
//	Liston and Elder (2006). This is basically a degree day approach.
//
//----------------------------------------------------------------------------

void tSnowPack::computeUnload() {

    //find if over critical temperature
    if (airTempK >= 273.16) {
        Lm = 5.8e-5 * (airTempK - 273.16) * timeSteps;//unload
    } else {
        Lm = 0.0;//do not unload
    }

    //RMK: Lm IS AN INTERNAL VARIABLE AND DOES NOT NEED TO BE RETURNED TO THE
    //	 CALLING FUNCTION.


}


double tSnowPack::inShortWaveSn(tCNode *cNode) {
    double Is, N, Iv, Isw, Ir;
    double v, cosi, scover;
    double RadGlobClr;

    //Remaining variables in DirectDiffuse from v3 -- AJR2008, SKY2008Snow
    //  So, entire function changed to match inShortWave in tEvapoTrans
/*  double h0, m, pp0, Dh0ref, h0ref, drm, Tlinke;
  double TnTLK, Fdh0, A1p, A1, A2, A3;
  double pi = 4*atan(1.0);*/

    Ic = Is = Id = Ir = Ids = Ics = Isw = Iv = 0.0;

    // Elevation, Slope and Aspect have been set before

    SunHour = 0.0; //Rinehart 2007 -- initialize whether we see sun or not to NO

    if (alphaD > 0.0) {

        elevation = cNode->getZ(); //SMM 10142008
        DirectDiffuse(elevation);  // SKY2008Snow, AJR2007

        // because inShortR is incoming solar radiation--sky cover is already factored in--WR

        //		// Cloud cover information
        //		if (fabs(skyCover-9999.99) < 1.0E-3) {
        //			skyCover = compSkyCover();//ADDED BY RINEHART 2007 @ NMT
        //			scover = skyCover;
        //		}
        //		else
        //			scover = skyCover;
        //        skyCoverC = scover;


        N = 0; //10.0;

        // If observations (for a horizontal surface) exist -
        // use them, at least in an approximate manner
        if (tsOption > 1) {
            RadGlobClr = (inShortR / (1.0 - 0.65 * pow(N, 2.0)));
            Ic = Ic / (Ic * sinAlpha + Id) * RadGlobClr;
            Id = RadGlobClr - Ic * sinAlpha;
        }

        // 1) Slope aspect
        //Account for the aspect and slope of the element
        //Estimate 'cosi' and compare it with the Sun position
        //  'cosi' = cos(i), where 'i' is the angle between
        //  the sun beam and the normal to the slope surface
        //
        //Rinehart 2007 @ New Mexico Tech
        //
        //	We have incorporated sheltering options. Option 3 is
        //	no topographic shading. Option 0, the default, is
        //	local topographic shading. Option 1 is incorporation
        //	of horizon angles in calc of SV and LV. Option 2 is
        //	the total integration of local and remote sheltering.
        //
        //	Here, if any sheltering is turned on, then we calculate
        //	the local controls of slope and aspect.
        //
        //	RMK: SLOPE AND ASPECT ARE CALCULATED FROM THE FLOW EDGE
        //	AND ARE IN RADIANS.

        if (shelterOption < 4) {//CHANGED IN 2008

            cosi = 1.0 * (cos(slope) * sinAlpha + sin(slope) * cos(asin(sinAlpha)) * cos(sunaz - aspect));

            if (cosi >= 0.0) {
                Ics = Ic * cosi;
                SunHour = 1.0; //YES SUN
            } else {
                Ics = 0.0;
                SunHour = 0.0; //NO SUN
            }
        } else {
            Ics = 1.0 * Ic;
            SunHour = 1.0;
        }

        if ((shelterOption == 2) || (shelterOption == 1)) {//CHANGED IN 2008

            Ics *= aboveHorizon(ID); //check to see if we can see the sun (aboveHorizon() in tEvapoTrans)
            SunHour *= aboveHorizon(ID);
        }

        //2) Horizon factor for diffuse radiation?
        //Rinehart 2007 @ New Mexico Tech
        //
        //	See comment above about sheltering options.

        if ((shelterOption > 0) && (shelterOption < 3)) {
            v = shelterFactorGlobal; //incorporate remote sheltering
        } else if (shelterOption == 0 || shelterOption == 3) { //CHANGED 2008
            v = 0.5 * (1 + cos(slope)); //local sheltering
        } else {
            v = 1.0; // no sheltering
        }

        Ids = Id * v;

        // 3) Account for cloud cover
        Is = (1.0 - 0.65 * pow(N, 2)) * (Ics + Ids);

        //Reflected from surrounded sites radiation
        //
        //Modified by Rinehart 2007 @ New Mexico Tech
        //
        if (hillAlbedoOption == 0) {
            hillalbedo = albedo;
        } else if (hillAlbedoOption == 1) {
            hillalbedo = coeffAl;
        } else if (hillAlbedoOption == 2) {

            if (snCanWE == 0)
                hillalbedo = coeffV * coeffAl + (1 - coeffV) * albedo;
            else
                hillalbedo = albedo;
        }

        if (shelterOption == 0) {
            //local
            Ir = hillalbedo * Is * (1 - cos(slope)) * 0.5;
            Is += Ir;
        } else if ((shelterOption > 1) && (shelterOption < 4)) { //CHANGED IN 2008
            //remote
            Ir = hillalbedo * Is * (0.5 * (1 + cos(slope)) - shelterFactorGlobal); //CHANGED IN 2008
            landRefGlobal = 0.5 * (1 + cos(slope)) - shelterFactorGlobal;
            Is += Ir;

        } else { //CHANGED IN 2008
            Ir = 0.0;
            Is = Is;
        }

        // Account for vegetation
        if ((evapotransOption == 1) || (snowOption)) {
            //Iv = Is * coeffKt * coeffV + Is * (1.0 - coeffV);
            Iv = Is * exp((coeffKt - 1) * coeffLAI) * coeffV +
                 Is * (1.0 - coeffV); // Changed to use Beer-Lambert following Moreno et al. (2016) CJC 2020
        } else
            Iv = Is;

        // Account for albedo
        //Modified by Rinehart 2007 @ New Mexico Tech
        //	This is actually the main difference b/t the tEvapoTrans::inShortWave
        //	and this function. We no longer see the land surface and have
        //	calculated albedo as a function of surface age earlier in the
        //	algorithm.

        Isw = Iv * (1.0 - albedo);

    } //end -- alphaD > 0
    else {
        Ic = Is = N = Iv = Isw = Id = Ids = Ics = Ir = 0.0;
    } // end -- alphaD <= 0

    // Set shortwave variables to the node (partition is approximate)
    if (tsOption > 1) {
        cNode->setShortRadIn(inShortR); //or set(Is), they must be equal
    } else {
        cNode->setShortRadIn(Isw);
    }
    cNode->setShortRadIn_dir(Ics * (1.0 - 0.65 * pow(N, 2.0))); // TODO: should these be set as values above the canopy?
    cNode->setShortRadIn_dif((Ids + Ir) * (1.0 - 0.65 * pow(N, 2.0)));//AJR2008, SKY2008Snow

    return Isw;
}

double tSnowPack::inShortWaveCan() {
    double Is, N, Iv, Isw, Ir;
    double v, cosi, scover;
    double RadGlobClr;

    // WR refactor 8-31-2023 this is a almost the same as inShortWave, but returns Isw before
    // accounting for the effects of optical transmission through the canopy. There is
    // certainly a cleaner way to do this, but for now this will have to do. Note this was originally in tSnowIntercept
    // which was removed because it was mostly redundant. See original tSnowIntercept code/documentation for information
    // on the original intent of this function.

    Ic = Is = Id = Ir = Ids = Ics = Isw = Iv = 0.0;

    // Elevation, Slope and Aspect have been set before

    SunHour = 0.0; //Rinehart 2007 -- initialize whether we see sun or not to NO

    if (alphaD > 0.0) {

        DirectDiffuse(elevation);  // SKY2008Snow, AJR2007

        // because inShortR is incoming solar radiation--sky cover is already factored in--WR

        //		// Cloud cover information
        //		if (fabs(skyCover-9999.99) < 1.0E-3) {
        //			skyCover = compSkyCover();//ADDED BY RINEHART 2007 @ NMT
        //			scover = skyCover;
        //		}
        //		else
        //			scover = skyCover;
        //        skyCoverC = scover;


        N = 0; //10.0;

        // If observations (for a horizontal surface) exist -
        // use them, at least in an approximate manner
        if (tsOption > 1) {
            RadGlobClr = (inShortR / (1.0 - 0.65 * pow(N, 2.0)));
            Ic = Ic / (Ic * sinAlpha + Id) * RadGlobClr;
            Id = RadGlobClr - Ic * sinAlpha;
        }

        // 1) Slope aspect
        //Account for the aspect and slope of the element
        //Estimate 'cosi' and compare it with the Sun position
        //  'cosi' = cos(i), where 'i' is the angle between
        //  the sun beam and the normal to the slope surface
        //
        //Rinehart 2007 @ New Mexico Tech
        //
        //	We have incorporated sheltering options. Option 3 is
        //	no topographic shading. Option 0, the default, is
        //	local topographic shading. Option 1 is incorporation
        //	of horizon angles in calc of SV and LV. Option 2 is
        //	the total integration of local and remote sheltering.
        //
        //	Here, if any sheltering is turned on, then we calculate
        //	the local controls of slope and aspect.
        //
        //	RMK: SLOPE AND ASPECT ARE CALCULATED FROM THE FLOW EDGE
        //	AND ARE IN RADIANS.

        if (shelterOption < 4) {//CHANGED IN 2008

            cosi = 1.0 * (cos(slope) * sinAlpha + sin(slope) * cos(asin(sinAlpha)) * cos(sunaz - aspect));

            if (cosi >= 0.0) {
                Ics = Ic * cosi;
                SunHour = 1.0; //YES SUN
            } else {
                Ics = 0.0;
                SunHour = 0.0; //NO SUN
            }
        } else {
            Ics = 1.0 * Ic;
            SunHour = 1.0;
        }

        if ((shelterOption == 2) || (shelterOption == 1)) {//CHANGED IN 2008

            Ics *= aboveHorizon(ID); //check to see if we can see the sun (aboveHorizon() in tEvapoTrans)
            SunHour *= aboveHorizon(ID);
        }

        //2) Horizon factor for diffuse radiation?
        //Rinehart 2007 @ New Mexico Tech
        //
        //	See comment above about sheltering options.

        if ((shelterOption > 0) && (shelterOption < 3)) {
            v = shelterFactorGlobal; //incorporate remote sheltering
        } else if (shelterOption == 0 || shelterOption == 3) { //CHANGED 2008
            v = 0.5 * (1 + cos(slope)); //local sheltering
        } else {
            v = 1.0; // no sheltering
        }

        Ids = Id * v;

        // 3) Account for cloud cover
        Is = (1.0 - 0.65 * pow(N, 2)) * (Ics + Ids);

        //Reflected from surrounded sites radiation
        //
        //Modified by Rinehart 2007 @ New Mexico Tech
        //
        if (hillAlbedoOption == 0) {
            hillalbedo = albedo;
        } else if (hillAlbedoOption == 1) {
            hillalbedo = coeffAl;
        } else if (hillAlbedoOption == 2) {

            if (snCanWE == 0)
                hillalbedo = coeffV * coeffAl + (1 - coeffV) * albedo;
            else
                hillalbedo = albedo;
        }

        if (shelterOption == 0) {
            //local
            Ir = hillalbedo * Is * (1 - cos(slope)) * 0.5;
            Is += Ir;
        } else if ((shelterOption > 1) && (shelterOption < 4)) { //CHANGED IN 2008
            //remote
            Ir = hillalbedo * Is * (0.5 * (1 + cos(slope)) - shelterFactorGlobal); //CHANGED IN 2008
            landRefGlobal = 0.5 * (1 + cos(slope)) - shelterFactorGlobal;
            Is += Ir;

        } else { //CHANGED IN 2008
            Ir = 0.0;
            Is = Is;
        }
        Iv = Is; // added 6/3/2024 otherwise function returns zero, but still not fixed since it called in computeSub which
        // then also factors for albedo. So maybe line below needs to be commented out?

        // Account for albedo
        Isw = Iv * (1.0 - albedo);

    } //end -- alphaD > 0
    else {
        Ic = Is = N = Iv = Isw = Id = Ids = Ics = Ir = 0.0;
    } // end -- alphaD <= 0

    return Isw;
}

//----------------------------------------------------------------------------
//
//			      tSnowPack::emmisSn()
//
//	Calculates the emmissivity of snow. Similar to the thermal properties
//	above, this is here for incorporation into a multilayer model later.
//
//----------------------------------------------------------------------------

double tSnowPack::emmisSn() {

    double emiss = 0.9;
    return emiss;
}

/************************************************************************************
**
**			tSnowPack -- Energy Balance Functions 
**
**	This set of functions forms the computational heart of the code. It calculates
**	the change in energy in the snow pack over a single time step.
**
**	This set of functions will need to be expanded for a multilayer model.
**	
**	Units are all kW/m^2
**
************************************************************************************/

//-----------------------------------------------------------------------------------
//
//				  tSnowPack::snowEB()
//
//	Calculates the change in energy over a single time step for the single-layer
//	model given all of the fluxes.
//
//-----------------------------------------------------------------------------------

void tSnowPack::snowEB(int nodeID, tCNode *node) {

    double sigma(5.67e-8);
    double v1;

    if (shelterOption > 0 && shelterOption < 3)
        v1 = shelterFactorGlobal;
    else
        v1 = 1;

    //convert temperature
    snTempK = CtoK(snTempC);

    //set up resistance
    resFact = resFactCalc();

    //turbulent heat fluxes
    H = sensibleHFCalc(resFact);
    L = latentHFCalc(resFact);

    //precipitation heat flux
    Prec = phfOnOff * precipitationHFCalc();

    //atmospheric heat flux
    RSin = naughttokilo * inShortWaveSn(node); // AJR2008, SKY2008Snow
    RLin = naughttokilo * inLongWave(node); // AJR2008, SKY2008Snow
    RLout = -naughttokilo * v1 * emmisSn() * sigma * pow(snTempK, 4.0);


    //set up for output
    inShortR = kilotonaught * RSin;
    inLongR = kilotonaught * RLin;
    outLongR = kilotonaught * RLout;

    // SKY2008Snow, AJR 2008
    //	Set the non-snow fluxes to zero
    node->setHFlux(0.0);
    node->setLFlux(0.0);
    node->setLongRadIn(0.0);
    node->setLongRadOut(0.0);
    node->setShortAbsbVeg(0.0);
    node->setShortAbsbSoi(0.0);

    // Set the snow fluxes
    node->setSnLHF(L * kilotonaught);
    node->setSnSHF(H * kilotonaught);
    node->setSnPHF(Prec * kilotonaught);
    node->setSnGHF(G * kilotonaught);
    node->setSnRLin(RLin * kilotonaught);
    node->setSnRLout(RLout * kilotonaught);
    node->setSnRSin(RSin * kilotonaught);


    Rn = RSin + RLin + RLout;
    G = 0;

    //calculate total dU over given timestep
    dUint = (H + Rn + L + Prec + G) * timeSteps; // Changed *3600 to *timeSteps CJC2020

    //find new energy state of snow
    Utot += dUint;
}

/*****************************************************************************
**
**			  tSnowPack I/O Functions
**
**	This set of functions deal with the I/O of tSnowPack that is not already
**	derived from tEvapoTrans.
**
*****************************************************************************/

int tSnowPack::getSnowOpt() {
    return snowOption; //found in tEvapoTrans construction
}


/*****************************************************************************
**
**			  tSnowPack Conversion Functions
**
**	This set of functions deal with the I/O of tSnowPack that is not already
**	derived from tEvapoTrans.
**
*****************************************************************************/

//---------------------------------------------------------------------------
//
//				tSnowPack::CtoK()
//	
//	Convert Celsius to Kelvin
//	
//---------------------------------------------------------------------------

double tSnowPack::CtoK(double temperature) {

    return (temperature + 273.15);

}

//---------------------------------------------------------------------------
//
//				tSnowPack::KtoC()
//
//	Convert Kelvin to Celsius
//---------------------------------------------------------------------------

double tSnowPack::KtoC(double temperature) {

    return (temperature - 273.15);

}

/***************************************************************************
**
** tSnowPack::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/
void tSnowPack::writeRestart(fstream &rStr) const {
    BinaryWrite(rStr, hillAlbedoOption);
    BinaryWrite(rStr, rhosnkg);
    BinaryWrite(rStr, rainTemp);
    BinaryWrite(rStr, ETAge);

    BinaryWrite(rStr, timeSteph);
    BinaryWrite(rStr, timeSteps);
    BinaryWrite(rStr, timeStepm);
    BinaryWrite(rStr, minutelyTimeStep);

    BinaryWrite(rStr, liqWE);
    BinaryWrite(rStr, iceWE);
    BinaryWrite(rStr, snWE);
    BinaryWrite(rStr, canWE);
    BinaryWrite(rStr, liqRoute);
    BinaryWrite(rStr, liqWEm);
    BinaryWrite(rStr, iceWEm);
    BinaryWrite(rStr, snWEm);
    BinaryWrite(rStr, Utot);
    BinaryWrite(rStr, Usn);
    BinaryWrite(rStr, Uwat);
    BinaryWrite(rStr, Utotold);
    BinaryWrite(rStr, liqWatCont);
    BinaryWrite(rStr, liqTempC);
    BinaryWrite(rStr, iceTempC);
    BinaryWrite(rStr, snTempC);
    BinaryWrite(rStr, liqTempK);
    BinaryWrite(rStr, iceTempK);
    BinaryWrite(rStr, snTempK);
    BinaryWrite(rStr, crustAge);

    BinaryWrite(rStr, H);
    BinaryWrite(rStr, L);
    BinaryWrite(rStr, G);
    BinaryWrite(rStr, Prec);
    BinaryWrite(rStr, Rn);
    BinaryWrite(rStr, dUint);
    BinaryWrite(rStr, RLin);
    BinaryWrite(rStr, RLout);
    BinaryWrite(rStr, RSin);
    BinaryWrite(rStr, Uerr);
    BinaryWrite(rStr, snPrec);
    BinaryWrite(rStr, liqPrec);
    BinaryWrite(rStr, snPrecm);
    BinaryWrite(rStr, liqPrecm);
    BinaryWrite(rStr, snPrecmm);
    BinaryWrite(rStr, liqPrecmm);
    BinaryWrite(rStr, snUnload);
    BinaryWrite(rStr, snCanWE);
    BinaryWrite(rStr, vapPressSmb);
    BinaryWrite(rStr, vapPresskSPa);

    BinaryWrite(rStr, rholiqcgs);
    BinaryWrite(rStr, rhoicecgs);
    BinaryWrite(rStr, rhosncgs);
    BinaryWrite(rStr, rholiqkg);
    BinaryWrite(rStr, rhoicekg);
    BinaryWrite(rStr, rhosnkg);
    BinaryWrite(rStr, rhoAir);
    BinaryWrite(rStr, phfOnOff);

    BinaryWrite(rStr, cpsnowkJ);
    BinaryWrite(rStr, cpicekJ);
    BinaryWrite(rStr, cpwaterkJ);
    BinaryWrite(rStr, cpairkJ);
    BinaryWrite(rStr, latFreezekJ);
    BinaryWrite(rStr, latVapkJ);
    BinaryWrite(rStr, latSubkJ);

    BinaryWrite(rStr, resFact);
    BinaryWrite(rStr, albedo);
    BinaryWrite(rStr, hillalbedo);
    BinaryWrite(rStr, compactParam);
    BinaryWrite(rStr, rhoSnFreshkg);
    BinaryWrite(rStr, minSnTemp);

    BinaryWrite(rStr, snDepth);
    BinaryWrite(rStr, snDepthm);
    BinaryWrite(rStr, snOnOff);
    BinaryWrite(rStr, peakSnWE);
    BinaryWrite(rStr, peakSnWEtemp);
    BinaryWrite(rStr, persMax);
    BinaryWrite(rStr, persMaxtemp);
    BinaryWrite(rStr, inittime);
    BinaryWrite(rStr, inittimeTemp);
    BinaryWrite(rStr, peaktime);

    BinaryWrite(rStr, naughttokilo);
    BinaryWrite(rStr, kilotonaught);
    BinaryWrite(rStr, cgsRHOtomks);
    BinaryWrite(rStr, mksRHOtocgs);
    BinaryWrite(rStr, naughttocm);
    BinaryWrite(rStr, cmtonaught);
    BinaryWrite(rStr, ctom);
    BinaryWrite(rStr, mtoc);

    BinaryWrite(rStr, Qcs);
    BinaryWrite(rStr, Ce);
    BinaryWrite(rStr, I);
    BinaryWrite(rStr, Iold);
    BinaryWrite(rStr, psiS);
    BinaryWrite(rStr, Imax);
    BinaryWrite(rStr, prec);
    BinaryWrite(rStr, LAI);
    BinaryWrite(rStr, kc);
    BinaryWrite(rStr, iceRad);
    BinaryWrite(rStr, dmdt);
    BinaryWrite(rStr, Omega);
    BinaryWrite(rStr, Sp);
    BinaryWrite(rStr, RH);
    BinaryWrite(rStr, D);
    BinaryWrite(rStr, rhoVap);
    BinaryWrite(rStr, Sh);
    BinaryWrite(rStr, Nu);
    BinaryWrite(rStr, Re);
    BinaryWrite(rStr, KtAtm);
    BinaryWrite(rStr, Ta);
    BinaryWrite(rStr, Mwater);
    BinaryWrite(rStr, R);
    BinaryWrite(rStr, RdryAir);
    BinaryWrite(rStr, esatIce);
    BinaryWrite(rStr, nu);
    BinaryWrite(rStr, beta);
    BinaryWrite(rStr, acoefficient);
    BinaryWrite(rStr, Lm);
    BinaryWrite(rStr, airTempK);
    BinaryWrite(rStr, effPrecip);

    tEvapoTrans::writeRestart(rStr);
}

/***************************************************************************
**
** tSnowPack::readRestart() Function
**
***************************************************************************/
void tSnowPack::readRestart(fstream &rStr) {
    BinaryRead(rStr, hillAlbedoOption);
    BinaryRead(rStr, rhosnkg);
    BinaryRead(rStr, rainTemp);
    BinaryRead(rStr, ETAge);

    BinaryRead(rStr, timeSteph);
    BinaryRead(rStr, timeSteps);
    BinaryRead(rStr, timeStepm);
    BinaryRead(rStr, minutelyTimeStep);

    BinaryRead(rStr, liqWE);
    BinaryRead(rStr, iceWE);
    BinaryRead(rStr, snWE);
    BinaryRead(rStr, canWE);
    BinaryRead(rStr, liqRoute);
    BinaryRead(rStr, liqWEm);
    BinaryRead(rStr, iceWEm);
    BinaryRead(rStr, snWEm);
    BinaryRead(rStr, Utot);
    BinaryRead(rStr, Usn);
    BinaryRead(rStr, Uwat);
    BinaryRead(rStr, Utotold);
    BinaryRead(rStr, liqWatCont);
    BinaryRead(rStr, liqTempC);
    BinaryRead(rStr, iceTempC);
    BinaryRead(rStr, snTempC);
    BinaryRead(rStr, liqTempK);
    BinaryRead(rStr, iceTempK);
    BinaryRead(rStr, snTempK);
    BinaryRead(rStr, crustAge);

    BinaryRead(rStr, H);
    BinaryRead(rStr, L);
    BinaryRead(rStr, G);
    BinaryRead(rStr, Prec);
    BinaryRead(rStr, Rn);
    BinaryRead(rStr, dUint);
    BinaryRead(rStr, RLin);
    BinaryRead(rStr, RLout);
    BinaryRead(rStr, RSin);
    BinaryRead(rStr, Uerr);
    BinaryRead(rStr, snPrec);
    BinaryRead(rStr, liqPrec);
    BinaryRead(rStr, snPrecm);
    BinaryRead(rStr, liqPrecm);
    BinaryRead(rStr, snPrecmm);
    BinaryRead(rStr, liqPrecmm);
    BinaryRead(rStr, snUnload);
    BinaryRead(rStr, snCanWE);
    BinaryRead(rStr, vapPressSmb);
    BinaryRead(rStr, vapPresskSPa);

    BinaryRead(rStr, rholiqcgs);
    BinaryRead(rStr, rhoicecgs);
    BinaryRead(rStr, rhosncgs);
    BinaryRead(rStr, rholiqkg);
    BinaryRead(rStr, rhoicekg);
    BinaryRead(rStr, rhosnkg);
    BinaryRead(rStr, rhoAir);
    BinaryRead(rStr, phfOnOff);

    BinaryRead(rStr, cpsnowkJ);
    BinaryRead(rStr, cpicekJ);
    BinaryRead(rStr, cpwaterkJ);
    BinaryRead(rStr, cpairkJ);
    BinaryRead(rStr, latFreezekJ);
    BinaryRead(rStr, latVapkJ);
    BinaryRead(rStr, latSubkJ);

    BinaryRead(rStr, resFact);
    BinaryRead(rStr, albedo);
    BinaryRead(rStr, hillalbedo);
    BinaryRead(rStr, compactParam);
    BinaryRead(rStr, rhoSnFreshkg);
    BinaryRead(rStr, minSnTemp);

    BinaryRead(rStr, snDepth);
    BinaryRead(rStr, snDepthm);
    BinaryRead(rStr, snOnOff);
    BinaryRead(rStr, peakSnWE);
    BinaryRead(rStr, peakSnWEtemp);
    BinaryRead(rStr, persMax);
    BinaryRead(rStr, persMaxtemp);
    BinaryRead(rStr, inittime);
    BinaryRead(rStr, inittimeTemp);
    BinaryRead(rStr, peaktime);

    BinaryRead(rStr, naughttokilo);
    BinaryRead(rStr, kilotonaught);
    BinaryRead(rStr, cgsRHOtomks);
    BinaryRead(rStr, mksRHOtocgs);
    BinaryRead(rStr, naughttocm);
    BinaryRead(rStr, cmtonaught);
    BinaryRead(rStr, ctom);
    BinaryRead(rStr, mtoc);

    BinaryRead(rStr, Qcs);
    BinaryRead(rStr, Ce);
    BinaryRead(rStr, I);
    BinaryRead(rStr, Iold);
    BinaryRead(rStr, psiS);
    BinaryRead(rStr, Imax);
    BinaryRead(rStr, prec);
    BinaryRead(rStr, LAI);
    BinaryRead(rStr, kc);
    BinaryRead(rStr, iceRad);
    BinaryRead(rStr, dmdt);
    BinaryRead(rStr, Omega);
    BinaryRead(rStr, Sp);
    BinaryRead(rStr, RH);
    BinaryRead(rStr, D);
    BinaryRead(rStr, rhoVap);
    BinaryRead(rStr, Sh);
    BinaryRead(rStr, Nu);
    BinaryRead(rStr, Re);
    BinaryRead(rStr, KtAtm);
    BinaryRead(rStr, Ta);
    BinaryRead(rStr, Mwater);
    BinaryRead(rStr, R);
    BinaryRead(rStr, RdryAir);
    BinaryRead(rStr, esatIce);
    BinaryRead(rStr, nu);
    BinaryRead(rStr, beta);
    BinaryRead(rStr, acoefficient);
    BinaryRead(rStr, Lm);
    BinaryRead(rStr, airTempK);
    BinaryRead(rStr, effPrecip);

    tEvapoTrans::readRestart(rStr);
}

void tSnowPack::checkShelter(tCNode *cNode) {
    if ((shelterOption > 0) && (shelterOption < 4)) {//CHANGED 2008
        int tempIndex;
        for (tempIndex = 0; tempIndex < 16; tempIndex++) {
            switch (tempIndex) {
                case 0:
                    ha2700 = cNode->getHorAngle2700();
                    break;
                case 1:
                    ha2925 = cNode->getHorAngle2925();
                    break;
                case 2:
                    ha3150 = cNode->getHorAngle3150();
                    break;
                case 3:
                    ha3375 = cNode->getHorAngle3375();
                    break;
                case 4:
                    ha0000 = cNode->getHorAngle0000();
                    break;
                case 5:
                    ha0225 = cNode->getHorAngle0225();
                    break;
                case 6:
                    ha0450 = cNode->getHorAngle0450();
                    break;
                case 7:
                    ha0675 = cNode->getHorAngle0675();
                    break;
                case 8:
                    ha0900 = cNode->getHorAngle0900();
                    break;
                case 9:
                    ha1125 = cNode->getHorAngle1125();
                    break;
                case 10:
                    ha1350 = cNode->getHorAngle1350();
                    break;
                case 11:
                    ha1575 = cNode->getHorAngle1575();
                    break;
                case 12:
                    ha1800 = cNode->getHorAngle1800();
                    break;
                case 13:
                    ha2025 = cNode->getHorAngle2025();
                    break;
                case 14:
                    ha2250 = cNode->getHorAngle2250();
                    break;
                case 15:
                    ha2475 = cNode->getHorAngle2475();
                    break;
                default:
                    cout << "\nCheck tempInd -- did not exist or assign" << endl;
            }//end-switch
        }    //end-for

        shelterFactorGlobal = cNode->getSheltFact(); //computed in tShelter

    } else if (shelterOption == 0) { // local sheltering for factor only
        shelterFactorGlobal = 0.5 * (1 + cos(slope)); //computed here for output purposes
        cNode->setSheltFact(shelterFactorGlobal); // SKYnGM2008LU
    } else {
        shelterFactorGlobal = 1; //no sheltering
    }
    return;
}


/******************************************************************************
**
**		    tSnowPack -- END OF TSNOWPACK.CPP
**
******************************************************************************/
