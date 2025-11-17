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
**  tWaterBalance.cpp:   Function file for tWaterBalance Class 
**                      (see tWaterBalance.h)
**
***************************************************************************/

#include "src/tHydro/tWaterBalance.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//                  Section 1: tWaterBalance Constructors/Destructors
//
//
//=========================================================================
tWaterBalance::tWaterBalance()
{
	gridPtr = 0;
	simCtrl = 0;
}

tWaterBalance::tWaterBalance(SimulationControl *simCtrPtr,tMesh<tCNode> *gridRef,
							 tInputFile &infile )
{ 
	gridPtr = gridRef;
	simCtrl = simCtrPtr;
	SetWaterBalance(infile);
}

tWaterBalance::~tWaterBalance()
{
	DeleteWaterBalance();
	Cout<<"tWaterBalance Object has been destroyed..."<<endl;
}

void tWaterBalance::SetWaterBalance(tInputFile &infile)
{
	initializeVariables(); 
	
	metStep = infile.ReadItem(metStep, "METSTEP");
	unsStep = infile.ReadItem(unsStep, "TIMESTEP");
	satStep = infile.ReadItem(satStep, "GWSTEP");
	
	return;
}

void tWaterBalance::DeleteWaterBalance()
{ 
	delete [] BasinStorages; 
	return;
}

//=========================================================================
//
//
//                  Section 2: tWaterBalance Functions
//
//
//=========================================================================


/***************************************************************************
**
** tWaterBalance::initializeVariables() Function
**
**
***************************************************************************/
void tWaterBalance::initializeVariables()
{
	BasinStorages = new double[6];
	for (int ct=0;ct<6;ct++)
		BasinStorages[ct] = 0.0;
	return;
}


//=========================================================================
//
//
//                  Section 2: Node Balance Functions
//
//
//=========================================================================


/***************************************************************************
**
** tWaterBalance::CanopyBalance() Function
**
** This function calculates the canopy storage for each computational node
** based on the amount of intercepted water and the evaporation loss from
** this storage.
**
**   DelI = Int - E
**      DelI = Change in Storage (mm/hr)
**      Int = Canopy Interception Rate (mm/hr)
**      E = Wet Canopy Evaporation (mm/hr)
** 
**   S = S+ DelI*dt*A/1000
**      S = Canopy Storage (m3)
**      dt = Time interval of meteorological scheme (hr)
**      A = Voronoi area (m2)
**
***************************************************************************/
void tWaterBalance::CanopyBalance()
{
	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	cNode = nodeIter.FirstP();
	
	double DelI, Int, E, dt, A, CS;
	
	dt = metStep/60.0;
	
	while(nodeIter.IsActive()){
		E = cNode->getEvapWetCanopy();
		Int = cNode->getInterceptLoss();
		A = cNode->getVArea();
		
		DelI = Int - E;
		
		// SKYnGM2008LU
		//cNode->setCanopyStorage(CS);

		CS = cNode->getCanopyStorVol() + DelI*dt*A/1000.0; // WR 12192023: Should this be scaled by vegetated fraction of cell?
		cNode->setCanopyStorVol(CS);

		cNode = nodeIter.NextP();
	}
	return;
}

/***************************************************************************
**
** tWaterBalance::UnSaturatedBalance() Function
**
** This function calculates the moisture storage in the unsaturated zone for 
** each computational node based on the amount of runoff, the evaporation loss,
** the infiltration and the lateral unsaturated flows.
**
**   DelU = Inf + QUin - QUout - Rech - ET - Run
**      DelU = Change in Storage (mm/hr)
**      Inf = Infiltration Rate (mm/hr) 
**      QUin = Lateral Flow In (mm/hr)
**      QUout = Lateral Flow Out (mm/hr)
**      ET = Evapotranspiration (mm/hr)
**      Rech = Recharge rate (mm/hr)
**      Run = Runoff from Unsaturated Zone (mm/hr)
** 
**   S = S+ DelU*dt*A/1000
**      S = Unsaturated Storage (m3)
**      dt = Time interval of unsaturated scheme (hr)
**      A = Voronoi area (m2)
**
***************************************************************************/
void tWaterBalance::UnSaturatedBalance()
{
	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	cNode = nodeIter.FirstP();
	
	double DelU, Inf, ET, Rech, QUin, QUout, Run, A, dt, USS, Melt, dMu;
	
	dt = unsStep/60.0;
	
	while(nodeIter.IsActive()){
		QUin = cNode->getUnSatFlowIn();
		QUout = cNode->getUnSatFlowOut();
		Rech = cNode->getRecharge();
		Run = cNode->getSrf()/unsStep;
		A = cNode->getVArea();
		Inf = cNode->getNetPrecipitation();
        Melt = cNode->getLiqRouted()*10.0;//WR 12192023: to mm, but implicitly mm/hr as thats total amount melted in 1 hr

        if(cNode->getLiqWE() + cNode->getIceWE() > 1e-4){
            Inf = Melt; //WR 12192023: snow on the ground Inf set to Melt
        }
        else{
            Inf = Inf+Melt; //WR 12192023: otherwise combined
        }


		ET = cNode->getEvapSoil()+cNode->getEvapDryCanopy();
		
		DelU = Inf + QUin - QUout - Rech - ET - Run;
		USS =  cNode->getUnSaturatedStorage() + DelU*dt*A/1000.0;
		cNode->setUnSaturatedStorage(USS);

        //WR 12192023: put check to error out if that chaning in total moisture above the water table (Mu) varies from the DelU by specified amount
        dMu = cNode->getMuNew()-cNode->getMuOld();
//        if (fabs(dMu-DelU*dt) > 10){
//            cerr<<"Change in total moisture above the water table, exceeds combined lateral and vertical fluxes by 1% of 1 mm."<<endl;
//        }



		cNode = nodeIter.NextP();
	}
	return;
}

/***************************************************************************
**
** tWaterBalance::SaturatedBalance() Function
**
** This function calculates the moisture storage in the saturated zone for 
** each computational node based on the amount of exfiltration, the recharge
** from the unsaturated zone and the lateral saturated flows.
**
**   DelG = Rech + QSTin - QSTout - Ex  [m3/hr]
**      DelG = Change in Storage (m3/hr)
**      Exf = Exfiltration rate (mm/hr) 
**      QSTin = Lateral Flow In (m3/hr)
**      QSTout = Lateral Flow Out (m3/hr)
**      Rech = Recharge rate (mm/hr)
** 
**   S = S+ DelG*dt*A/1000
**      S = Saturated Storage (m3)
**      dt = Time interval of saturated scheme (hr)
**      A = Voronoi area (m2)
** 
** Note: Exfiltration given by subsurface runoff
**
***************************************************************************/
void tWaterBalance::SaturatedBalance()
{
	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	cNode = nodeIter.FirstP();
	
	double DelG, Rech, QSTin, QSTout, Exf, A, dt, SS;
	
	dt = satStep/60.0;
	
	while(nodeIter.IsActive()){
		QSTin = cNode->getQgwIn() * 1.0E-9; //convert from mm3/hr to m3/hr
		QSTout = cNode->getQgwOut() * 1.0E-9;
		Rech = cNode->getRecharge();
		Exf = cNode->getSrf()/satStep;
		A =  cNode->getVArea();
		
		// Convert Rech, Exf from mm/hr to m3/hr
		DelG = QSTin - QSTout + (Rech - Exf)*A/1000.0;  
		SS =  cNode->getSaturatedStorage() + DelG*dt;
		cNode->setSaturatedStorage(SS);
		
		cNode = nodeIter.NextP();
	}
	return;
}

//=========================================================================
//
//
//                  Section 4: Watershed Balance Functions
//
//
//=========================================================================

/***************************************************************************
**
** tWaterBalance::BasinStorage Function
**
** This function computes the total basin storage for each of the land
** surface stores. 
**
** BasinRainfall (m3/hr) = Rainfall(mm/hr) * Area (m2)/1000.0
** BasinEvaporation (m3/hr) = (CanopyE + SurfaceE) (mm/hr) * Area (m2)/1000.0
** BasinUnsaturated (m3)
** BasinSaturated (m3)
**
***************************************************************************/
void tWaterBalance::BasinStorage( double time )
{
	tCNode * cNode;
	tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());
	cNode = nodeIter.FirstP();
	
	double BasinCanopy, BasinUnSaturated, BasinSaturated,BasinUnSaturated_old, BasinSaturated_old;
	double BasinRainfall, BasinEvaporation, BasinRunoff;
    double Balance;
	
	BasinCanopy = BasinUnSaturated = BasinRunoff = 0.0;
	BasinSaturated = BasinRainfall = BasinEvaporation = 0.0;
    BasinSaturated_old = BasinSaturated_old = 0.0;

    Balance = 0;
	
	while(nodeIter.IsActive()){
        //WR 12192023: set stores to snapshots in time, not cumulative values, (i.e. removed +=)
		BasinCanopy = cNode->getCanStorage()*(cNode->getVArea()/1000.0);//WR 12192023:  Needs to be scaled by veg fract?
		BasinRainfall += cNode->getRain()*(cNode->getVArea()/1000.0) * unsStep/60.0;
		BasinEvaporation+= cNode->getEvapoTrans()*(cNode->getVArea()/1000.0)* unsStep/60.0;//WR 12192023: All ready scaled with vegetation fraciton
		BasinRunoff += cNode->getSrf()*(cNode->getVArea()/1000.0);
        BasinUnSaturated_old = cNode->getMuOld()*(cNode->getVArea()/1000.0);
        BasinSaturated_old = cNode->getThetaS()*(cNode->getBedrockDepth() - cNode->getNwtOld())*(cNode->getVArea()/1000.0);
        BasinUnSaturated = cNode->getMuNew()*(cNode->getVArea()/1000.0);
		BasinSaturated = cNode->getThetaS()*(cNode->getBedrockDepth() - cNode->getNwtNew())*(cNode->getVArea()/1000.0); //WR 12192023:  Needs to be scaled with porosity correct?
		cNode = nodeIter.NextP();
	}

//    while(nodeIter.IsActive()){ //  in mm
//        //WR 12192023: set stores to snapshots in time, not cumulative values, (i.e. removed +=)
//        BasinCanopy = cNode->getCanStorage();//WR 12192023:  Needs to be scaled by veg fract?
//        BasinRainfall += cNode->getRain()* unsStep/60.0;
//        BasinEvaporation+= cNode->getEvapoTrans()* unsStep/60.0;//WR 12192023: All ready scaled with vegetation fraciton
//        BasinRunoff += cNode->getSrf();
//        BasinUnSaturated_old = cNode->getMuOld();
//        BasinSaturated_old = cNode->getThetaS()*(cNode->getBedrockDepth() - cNode->getNwtOld());
//        BasinUnSaturated = cNode->getMuNew();
//        BasinSaturated = cNode->getThetaS()*(cNode->getBedrockDepth() - cNode->getNwtNew()); //WR 12192023:  Needs to be scaled with porosity correct?
//        cNode = nodeIter.NextP();
//    }

    Balance = BasinRainfall-BasinEvaporation-BasinRunoff-(BasinSaturated-BasinSaturated_old+BasinUnSaturated-BasinUnSaturated_old+BasinCanopy-BasinStorages[2]);

    //WR 12192023: Set statement to check water balance error is above some acceptable threshold, in this case set to 1 m
//    if (fabs(Balance) > 0.1){
//        cerr << "Water balance error is above 1m^3: " << std::setprecision(5) << Balance << endl;
//    }

    BasinStorages[0] += BasinRainfall;
	BasinStorages[1] += BasinEvaporation;
	BasinStorages[2] = BasinCanopy;
	BasinStorages[3] = BasinUnSaturated;
	BasinStorages[4] = BasinSaturated;
	BasinStorages[5] += BasinRunoff;
	
	if(time == unsStep/60.0){
		cout<<"\n\tInitial Estimates";
		this->Print(BasinStorages);
	}
	return;
}

//=========================================================================
//
//
//                  Section 5: Output Functions
//
//
//=========================================================================
void tWaterBalance::Print(double * array)
{
	cout<<"\n-----------------------------------";
	cout<<"\n tRIBS Water Balance Calculations";
	cout<<"\n\n Basin Rainfall = " <<array[0]<<" (m3)";
	cout<<"\n Basin Evaporation = "<<array[1]<<" (m3)";
	cout<<"\n Basin Runoff = "<<array[5]<<" (m3)";
	cout<<"\n Basin Canopy Storage = "<<array[2]<<" (m3)";
	cout<<"\n Basin Unsaturated Zone Storage = "<<array[3]<<" (m3)";
	cout<<"\n Basin Saturated Zone Storage = "<<array[4]<<" (m3)";
	cout<<"\n-----------------------------------\n\n\n";
	return;
}

/***************************************************************************
**
** tWaterBalance::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/
void tWaterBalance::writeRestart(fstream & rStr) const
{
  BinaryWrite(rStr, metStep);
  BinaryWrite(rStr, unsStep);
  BinaryWrite(rStr, satStep);
  BinaryWrite(rStr, finalTime);
  for (int i = 0; i < 6; i++)
    BinaryWrite(rStr, BasinStorages[i]);
}

/***************************************************************************
**
** tWaterBalance::readRestart() Function
**
***************************************************************************/
void tWaterBalance::readRestart(fstream & rStr)
{
  BinaryRead(rStr, metStep);
  BinaryRead(rStr, unsStep);
  BinaryRead(rStr, satStep);
  BinaryRead(rStr, finalTime);
  for (int i = 0; i < 6; i++)
    BinaryRead(rStr, BasinStorages[i]);
}

//=========================================================================
//
//
//                       End of tWaterBalance.cpp
//
//
//=========================================================================
