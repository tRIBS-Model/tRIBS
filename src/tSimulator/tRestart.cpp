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
**  tRestart.cpp: Functions for class tRestart (see tRestart.h)
**
***************************************************************************/

#include <cassert>

#include "src/tSimulator/tRestart.h"

/*************************************************************************
**
** Constructor
**
*************************************************************************/

template< class tSubNode >
tRestart<tSubNode>::tRestart(
	tRunTimer* t,
	tMesh<tSubNode>* m,
	tKinemat* f,
	tWaterBalance* b,
	tHydroModel* h,
	tRainfall* r,
	tEvapoTrans* e,
	tIntercept* i,
   tSnowPack* s
   )
{
  timer = t;
  mesh = m;
  flow = f;
  balance = b;
  hydro = h;
  rainfall = r;
  evap = e;
  intercept = i;
  snowpack = s;
}

/*************************************************************************
**
** Write restart information for all controlled objects
**
*************************************************************************/

template< class tSubNode >
void tRestart<tSubNode>::writeRestart(fstream & rStr)
{
  timer->writeRestart(rStr);
  flow->writeRestart(rStr);
  balance->writeRestart(rStr);
  hydro->writeRestart(rStr);
  rainfall->writeRestart(rStr);
  intercept->readRestart(rStr);
  mesh->writeRestart(rStr);
	// Giuseppe DEBUG Restart 2012 - START 
	// I have introduced an IF that checks whether
	// if the snow module is on. If not, the relative variables 
	// are not saved in the binary Restart files.
	if (snowpack->getSnowOpt() != 0){
		snowpack->writeRestart(rStr);
	}
    else{
        evap->writeRestart(rStr);
    }// Giuseppe DEBUG Restart 2012 - END
	

}

/*************************************************************************
**
** Read restart information for all controlled objects
**
*************************************************************************/

template< class tSubNode >
void tRestart<tSubNode>::readRestart(fstream & rStr)
{
	timer->readRestart(rStr);
	flow->readRestart(rStr);
	balance->readRestart(rStr);
	hydro->readRestart(rStr);
	rainfall->readRestart(rStr);
	intercept->readRestart(rStr);
	mesh->readRestart(rStr);
	// Giuseppe DEBUG Restart 2012 - START 
	// I have introduced an IF that checks whether
	// if the snow module is on. This is to be consistent
	// with the writeRestart function.
	if (snowpack->getSnowOpt() != 0){	
		snowpack->readRestart(rStr);
	}
    else{
        evap->readRestart(rStr);
    }// Giuseppe DEBUG Restart 2012 - END// Giuseppe DEBUG Restart 2012 - END
}

//=========================================================================
//
//
//                        End of tRestart.cpp
//
//
//=========================================================================
