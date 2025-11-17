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
**  tVariant.cpp: Functions for class tVariant (see tVariant.h) 
**
***************************************************************************/

#include "src/tRasTin/tVariant.h"

//=========================================================================
//
//
//                  Section 1: tVariant Constructors and Destructors
//
//
//=========================================================================

tVariant::tVariant()
{
	gridPtr = 0;
}

tVariant::tVariant(tMesh<tCNode> *gridRef, tResample *resamp) 
{
	gridPtr = gridRef;
	respPtr = resamp; 
}

//=========================================================================
//
//
//                  Section 2: tVariant Functions
//
//
//=========================================================================

/***************************************************************************
**
** setFileNames() Function
**
***************************************************************************/
void tVariant::setFileNames(char *file, char *ext)
{
	for(int ct=0;ct<kName;ct++) {
		inputName[ct] = file[ct];
	}
    for (int dt = 0; dt < kMaxExt; dt++)
        extension[dt] = ext[dt]; //TODO: Sanitizer warning Heap-buffer-overflow on address -WR, 2nd thought, maybe because file doesnt exist
    return;
}

/***************************************************************************
**
** newVariable() Function
**
** Char* argument used to identify the variant parameter of interest
** Now considers the land use parameters in addition to the meteorological 
** parameters in tEvapoTrans (SKYnGM2008LU). Could be expanded to include 
** other variables.
**
***************************************************************************/
void tVariant::newVariable(char *param)
{
	int id;
	tCNode * cn;
	tMeshListIter<tCNode> nodeIter( gridPtr->getNodeList() );
	
	id = 0;
	cn = nodeIter.FirstP();
	while( nodeIter.IsActive() ) {
		if (strcmp(param,"PA")==0)
			cn->setAirPressure( 0.0 );
		else if (strcmp(param,"TD")==0)
			cn->setDewTemp( 0.0 );
		else if (strcmp(param,"XC")==0)
			cn->setSkyCover( 0.0 );
		else if (strcmp(param,"US")==0)
			cn->setWindSpeed( 0.0 );
		else if (strcmp(param,"TA")==0)
			cn->setAirTemp( 0.0 );
		else if (strcmp(param,"TS")==0)
			cn->setSurfTemp( 0.0 );
		else if (strcmp(param,"NR")==0)
			cn->setNetRad( 0.0 );
		else if (strcmp(param,"ET")==0)
			cn->setGridET( 0.0 );
		else if (strcmp(param,"RH")==0)
			cn->setRelHumid( 0.0 );
		else if (strcmp(param,"VP")==0)
			cn->setVapPressure( 0.0 );
		else if (strcmp(param,"IS")==0) //E.R.V 3/6/2012
			cn->setShortRadIn( 0.0 );

		cn = nodeIter.NextP();
		id++;
	}
	return;
}

/***************************************************************************
**
** updateVariable() Function
**
** Char* argument used to identify the variant parameter of interest
** Resamples the current grid and assigns to tCNode
**
***************************************************************************/
void tVariant::updateVariable(char *param)
{
	int id;
	tCNode * cn;
	tMeshListIter<tCNode> nodeIter( gridPtr->getNodeList() );
	double *resample;
	
	// Resamples input ASCII grid
	resample = respPtr->doIt(fileIn, 1); 
	
	id = 0;
	cn = nodeIter.FirstP();
	while( nodeIter.IsActive() ) {
		if (strcmp(param,"PA")==0)
			cn->setAirPressure(resample[id]);
		else if (strcmp(param,"TD")==0)
			cn->setDewTemp( resample[id] );
		else if (strcmp(param,"XC")==0)
			cn->setSkyCover( resample[id] );
		else if (strcmp(param,"US")==0)
			cn->setWindSpeed(resample[id] );
		else if (strcmp(param,"TA")==0)
			cn->setAirTemp(resample[id] );
		else if (strcmp(param,"TS")==0)
			cn->setSurfTemp( resample[id] );
		else if (strcmp(param,"NR")==0)
			cn->setNetRad( resample[id] );
		else if (strcmp(param,"RH")==0)
			cn->setRelHumid( resample[id] );
		else if (strcmp(param,"VP")==0)
			cn->setVapPressure( resample[id] );
		else if (strcmp(param,"IS")==0)  //E.R.V 3/6/2012
			cn->setShortRadIn( resample[id] );
		else if (strcmp(param,"ET")==0)
			cn->setGridET( resample[id] );

		cn = nodeIter.NextP();
		id++; 
	}
	return;
}

// SKYnGM2008LU
/***************************************************************************
**
** updateLUVarOfPrevGrid() Function
**
** First char* argument used to identify the variant parameter of interest
** Second char* argument used is the input grid file that is resampled
** Resamples the grid previous to the model timestep and assigns to tCNode
** Usually useful only at time step 1, since for any subsequent time steps,
** when a new 'previous' grid time instant is crossed, its values are 
** obtained simply by setting the new 'previous' grid values to the old 
** 'until' grid values in updateLUVarOfBothGrids
** Also useful when final 'previous' grid time instant is crossed
***************************************************************************/
void tVariant::updateLUVarOfPrevGrid(const char *param, char *GridFileName)
{
	int id;
	tCNode * cn;
	tMeshListIter<tCNode> nodeIter( gridPtr->getNodeList() );
	double *resample;
	
	// Resamples input ASCII grid
	resample = respPtr->doIt(GridFileName, 1); 
	
	id = 0;
	cn = nodeIter.FirstP();
	while( nodeIter.IsActive() ) {
		if (strcmp(param,"AL") == 0) { 
			cn->setLandUseAlbInPrevGrid( resample[id] );
			cn->setLandUseAlb( resample[id] ); 
		}
		else if (strcmp(param,"TF") == 0) {
			cn->setThroughFallInPrevGrid( resample[id] );
			cn->setThroughFall( resample[id] );
		}
		else if (strcmp(param,"VH") == 0) {
			cn->setVegHeightInPrevGrid( resample[id] );
			cn->setVegHeight( resample[id] );
		}
		else if (strcmp(param,"SR") == 0) {
			cn->setStomResInPrevGrid( resample[id] );
			cn->setStomRes( resample[id] );
		}
		else if (strcmp(param,"VF") == 0) {
			cn->setVegFractionInPrevGrid( resample[id] );
			cn->setVegFraction( resample[id] );
		}
		else if (strcmp(param,"CS") == 0) {
			cn->setCanStorParamInPrevGrid( resample[id] );
			cn->setCanStorParam( resample[id] );
		}
		else if (strcmp(param,"IC") == 0) {
			cn->setIntercepCoeffInPrevGrid( resample[id] );
			cn->setIntercepCoeff( resample[id] );
		}
		else if (strcmp(param,"CC") == 0) {
			cn->setCanFieldCapInPrevGrid( resample[id] );
			cn->setCanFieldCap( resample[id] );
		}
		else if (strcmp(param,"DC") == 0) {
			cn->setDrainCoeffInPrevGrid( resample[id] );
			cn->setDrainCoeff( resample[id] );
		}
		else if (strcmp(param,"DE") == 0) {
			cn->setDrainExpParInPrevGrid( resample[id] );
			cn->setDrainExpPar( resample[id] );
		}
		else if (strcmp(param,"OT") == 0) {
			cn->setOptTransmCoeffInPrevGrid( resample[id] );
			cn->setOptTransmCoeff( resample[id] );
		}
		else if (strcmp(param,"LA") == 0) {
			cn->setLeafAIInPrevGrid( resample[id] );
			cn->setLeafAI( resample[id] );
		}
		// CJC2025: New parameters
		else if (strcmp(param,"SE") == 0) {
			cn->setEvapThreshInPrevGrid( resample[id] );
			cn->setEvapThresh( resample[id] );
		}
		else if (strcmp(param,"ST") == 0) {
			cn->setTransThreshInPrevGrid( resample[id] );
			cn->setTransThresh( resample[id] );
		}

		cn = nodeIter.NextP();
		id++; 
	}
	return;
}

// SKYnGM2008LU
/***************************************************************************
**
** updateLUVarOfBothGrids() Function
**
** First char* argument used to identify the variant parameter of interest
** Second char* argument used is the input grid file that is resampled
** Resamples the grid after and including the model timestep and assigns to tCNode
** when a new 'previous' grid time instant is crossed.
** Before doing so, it sets the new 'previous' grid values to old 'until' grid values
***************************************************************************/
void tVariant::updateLUVarOfBothGrids(const char *param, char *GridFileName)
{
	int id;
	tCNode * cn;
	tMeshListIter<tCNode> nodeIter( gridPtr->getNodeList() );
	double *resample;
       
	// Resamples input ASCII grid
	resample = respPtr->doIt(GridFileName, 1); 

	id = 0;
	cn = nodeIter.FirstP();
	while( nodeIter.IsActive() ) {
		if (strcmp(param,"AL") == 0) { 
			cn->setLandUseAlbInPrevGrid( cn->getLandUseAlbInUntilGrid() );
			cn->setLandUseAlb(cn->getLandUseAlbInUntilGrid() ); 
			cn->setLandUseAlbInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"TF") == 0) {
			cn->setThroughFallInPrevGrid( cn->getThroughFallInUntilGrid() );
			cn->setThroughFall( cn->getThroughFallInUntilGrid() );
			cn->setThroughFallInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"VH") == 0) {
			cn->setVegHeightInPrevGrid( cn->getVegHeightInUntilGrid() );
			cn->setVegHeight( cn->getVegHeightInUntilGrid() );
			cn->setVegHeightInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"SR") == 0) {
			cn->setStomResInPrevGrid( cn->getStomResInUntilGrid() );
			cn->setStomRes( cn->getStomResInUntilGrid() );
			cn->setStomResInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"VF") == 0) {
			cn->setVegFractionInPrevGrid( cn->getVegFractionInUntilGrid() );
			cn->setVegFraction( cn->getVegFractionInUntilGrid() );
			cn->setVegFractionInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"CS") == 0) {
			cn->setCanStorParamInPrevGrid( cn->getCanStorParamInUntilGrid() );
			cn->setCanStorParam( cn->getCanStorParamInUntilGrid() );
			cn->setCanStorParamInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"IC") == 0) {
			cn->setIntercepCoeffInPrevGrid( cn->getIntercepCoeffInUntilGrid() );
			cn->setIntercepCoeff( cn->getIntercepCoeffInUntilGrid() );
			cn->setIntercepCoeffInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"CC") == 0) {
			cn->setCanFieldCapInPrevGrid( cn->getCanFieldCapInUntilGrid() );
			cn->setCanFieldCap( cn->getCanFieldCapInUntilGrid() );
			cn->setCanFieldCapInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"DC") == 0) {
			cn->setDrainCoeffInPrevGrid( cn->getDrainCoeffInUntilGrid() );
			cn->setDrainCoeff( cn->getDrainCoeffInUntilGrid() );
			cn->setDrainCoeffInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"DE") == 0) {
			cn->setDrainExpParInPrevGrid( cn->getDrainExpParInUntilGrid() );
			cn->setDrainExpPar( cn->getDrainExpParInUntilGrid() );
			cn->setDrainExpParInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"OT") == 0) {
			cn->setOptTransmCoeffInPrevGrid( cn->getOptTransmCoeffInUntilGrid() );
			cn->setOptTransmCoeff( cn->getOptTransmCoeffInUntilGrid() );
			cn->setOptTransmCoeffInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"LA") == 0) {
			cn->setLeafAIInPrevGrid( cn->getLeafAIInUntilGrid() );
			cn->setLeafAI( cn->getLeafAIInUntilGrid() );
			cn->setLeafAIInUntilGrid( resample[id] );
		}
		// CJC2025: New parameters
		else if (strcmp(param,"SE") == 0) {
			cn->setEvapThreshInPrevGrid( cn->getEvapThreshInUntilGrid() );
			cn->setEvapThresh( cn->getEvapThreshInUntilGrid() );
			cn->setEvapThreshInUntilGrid( resample[id] );
		}
		else if (strcmp(param,"ST") == 0) {
			cn->setTransThreshInPrevGrid( cn->getTransThreshInUntilGrid() );
			cn->setTransThresh( cn->getTransThreshInUntilGrid() );
			cn->setTransThreshInUntilGrid( resample[id] );
		}


		cn = nodeIter.NextP();
		id++; 
	}
	return;
}

/***************************************************************************
**
** composeFileName() Function
**
***************************************************************************/
int tVariant::composeFileName(tRunTimer *t)
{
	//if (infile)
	if (infile.is_open()) // SKYnGM2008LU
		infile.close();
	
	if (t->getMetStep() < t->getRainDT()){   //If 'minute' is NOT equal to '0'
		snprintf(fileIn, sizeof(fileIn), "%s%02d%02d%04d%02d%02d.%s", inputName,
				t->month, t->day, t->year, t->hour, t->minute, extension);
	}
	else{                    //If 'minute' IS equal to '0'
		snprintf(fileIn, sizeof(fileIn), "%s%02d%02d%04d%02d.%s", inputName,
				t->month, t->day, t->year, t->hour, extension);}
	
	infile.open(fileIn);
	//if (!infile)
	if (!infile.is_open()) // SKYnGM2008LU
		return 0;
	else 
		return 1;
}

/***************************************************************************
**
** noData() Function
**
** Assigns a no_data value 9999.99 to the tCNode member variable 
** with no ascii grid input.
**
***************************************************************************/
void tVariant::noData(char *param)
{
	int id;
	tCNode * cn;
	tMeshListIter<tCNode> nodeIter( gridPtr->getNodeList() );
	
	id = 0;
	cn = nodeIter.FirstP();
	while( nodeIter.IsActive() ) {
		if (strcmp(param,"PA")==0)
			cn->setAirPressure( 9999.99 );
		else if (strcmp(param,"TD")==0)
			cn->setDewTemp( 9999.99 );
		else if (strcmp(param,"XC")==0)  //E.R.V. 3/26/2012
			cn->setSkyCover( 9999.99 );
		else if (strcmp(param,"US")==0)
			cn->setWindSpeed( 9999.99 );
		else if (strcmp(param,"TA")==0)
			cn->setAirTemp( 9999.99 );
		else if (strcmp(param,"TS")==0)
			cn->setSurfTemp( 9999.99 );
		else if (strcmp(param,"NR")==0)
			cn->setNetRad( 9999.99 );
		else if (strcmp(param,"RH")==0)
			cn->setRelHumid( 9999.99 );
		else if (strcmp(param,"ET")==0)
			cn->setGridET( 9999.99 );
		else if (strcmp(param,"VP")==0)
			cn->setVapPressure( 9999.99 );
		else if (strcmp(param,"IS")==0)  //E.R.V 3/6/2012
			cn->setShortRadIn( 9999.99 );

		cn = nodeIter.NextP();
		id++;
	}
	return;
}

// Get Functions
char * tVariant::getInputName() { return inputName; }
char * tVariant::getExtension() { return extension; }

//=========================================================================
//
//
//                        End of tVariant.cpp
//
//
//=========================================================================
