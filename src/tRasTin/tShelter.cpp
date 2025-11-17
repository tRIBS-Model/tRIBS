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
**  tShelter.cpp:   Functions for tShelter classes (see tShelter.h)
**
**
***************************************************************************/

#include "src/tRasTin/tShelter.h"
#include "src/Headers/globalIO.h"

//==========================================================================
//
//		tShelter Constuctor/Destructor
//
//==========================================================================

tShelter::tShelter(SimulationControl *simCtrPtr, tMesh<tCNode> * gridPtr,
		    tInputFile &infile) : tResample(simCtrPtr , gridPtr) 

{ 

  int i(0), j(0), k(0), l(0);//looping variables
  int angleHolder(0);
  int delRow(0), delCol(0);
  int delXInt(0), delYInt(0);//variables for navigating through grid
  double delXReal(0.), delYReal(0.), remX(0.), remY(0.);//structural navigation 
  double tanAngle(0.), cosAngle(0.);			//variables
  double angle(0);//angle of transform to loop through
  double absTan, absCos;
  
  tCNode *cn;
  tMeshListIter< tCNode > niter ( mew->getNodeList() );
  tPtrList< tCNode >     NodesLst;
  tPtrListIter< tCNode > NodesIter( NodesLst );
  NodesLst.Flush();

  //get input information
  radSheltOpt = infile.ReadItem(radSheltOpt,"OPTRADSHELT");

  //get the file path to the DEM
  infile.ReadItem(GridInPath,"DEMFILE");

  angleDiv = 22.5;

//  cout << "Shelt option: " << radSheltOpt << endl;

  
  //derive HA maps
  if ( (radSheltOpt > 0) && (radSheltOpt < 4) ) {//CHANGED IN 2008
  
  cout << "read DEM" << endl;
  readInputGrid(GridInPath);


  // cout << "initialize grid" << endl;
  tempGrid = new double* [NR]; 
  assert(tempGrid != 0);
  
  for (i=0; i < NR; i++)  {
    tempGrid[i] = new double [MR];
    assert(gridIn[i] != 0);
    for (j=0; j < MR; j++)    {
      tempGrid[i][j] = gridIn[i][j];
    }
  }

  //i,j are for the address of the point of interest
  //k is for looping through the angles
  //
  //cout << "navigation variable loop" << endl;
  for (int tempInd = 0; tempInd < 16; tempInd = tempInd+1) {

    //cout << "navigation variable loop"<< tempInd << endl;

    angle = 22.5*double(tempInd);

    //find navigation variables and initialize properly
    //
    //	  This "pixelates" the navigation of the grid so that
    //	  we can work on angles.
    //
    //	  If the code is to use a dynamic discretization of the
    //	  azimuthal angles, then this will need to be automated.
    //	  But that is what held me up for so long last time.
    
    switch (tempInd) {
      case 0:
	//0.0 degrees
	delXReal = 1.0;
	delYReal = 0.0;
	break;
      case 1:
	//22.5 degrees
	delXReal = 1.0;
	delYReal = 0.414213562;
	break;
      case 2:
	//45.0 degrees
	delXReal = 1.0;
	delYReal = 1.0;
	break;
      case 3:
	//67.5 degrees
	delXReal = 0.414213562;
	delYReal = 1.0;
      case 4:
	//90.0 degrees
	delXReal = 0.0;
	delYReal = 1.0;
	break;
      case 5:
	//112.5 degrees
	delXReal = -0.414213562;
	delYReal = 1.0;
	break;
      case 6:
	//135.0 degrees
	delXReal = -1.0;
	delYReal = 1.0;
	break;
      case 7:
	//157.5 degrees
	delXReal = -1.0;
	delYReal = 0.414213562;
	break;
      case 8:
	//180.0 degrees
	delXReal = -1.0;
	delYReal = 0.0;
	break;
      case 9:
	//202.5 degrees
	delXReal = -1.0;
	delYReal = -0.414213562;
	break;
      case 10:
	//225.0 degrees
	delXReal = -1.0;
	delYReal = -1.0;
	break;
      case 11:
	//247.5 degrees
	delXReal = -0.414213562;
	delYReal = -1.0;
	break;
      case 12:
	//270.0 degrees
	delXReal = 0.0;
	delYReal = -1.0;
	break;
      case 13:
	//292.5 degrees
	delXReal = 0.414213562;
	delYReal = -1.0;
	break;
      case 14:
	//315.0 degrees
	delXReal = 1.0;
	delYReal = -1.0;
	break;
      case 15:
	//337.5 degrees
	delXReal = 1.0;
	delYReal = -0.414213562;
	break;
      default:
	cout << "\nSomething is screwed up in the for loop of tShelter alg" << endl;
    }
    
    //find trig values
    tanAngle = tan(3.1416*angle/180);
    cosAngle = cos(3.1416*angle/180);
    absTan = fabs(3.1416*tanAngle/180);
    absCos = fabs(3.1416*cosAngle/180);

    //initialize index stepping variables and remainder variables
    delXInt = 0;
    delYInt = 0;
    remX = 0;
    remY = 0;

//    cout << "go to unique (or starting) point in the grid" << endl;
    for(i = 0; i < NR; i++) {

      for(j = 0; j < MR; j++) {
	

	//ALL INITIALIZATIONS ARE FOR WHILE-LOOP
	
	//initialize elevation
	initZ = tempGrid[i][j];

	//intitialize HA to lowest value
  	gridIn[i][j] = 0.0;

	//initialize location of where HA occurs
	k = i;
	l = j;
	
	//initialize state variables
	maxTan = 0.0;
	maxRow = i;
	maxCol = j;
	maxZ = initZ;
	
	//initialize remainders and integer steps
	remX = delXReal;
	remY = delYReal;
	
	//set step in x-dir
	if (fabs(remX) >= 1.0) {
	  if (remX >= 0.0) {
	    delXInt = int(floor(remX));
	    remX = remX - 1.0;
	  }
	  else if (remX < 0.0) {
	    delXInt = int(ceil(remX));
	    remX = remX + 1.0;
	  }
	}

	//set integer step in y-dir
	if (fabs(remY) >= 1.0) {
	  if (remY >= 0.0) {
	    delYInt = int(floor(remY));//go down
	    remY = remY - 1.0;//find difference and redo remainder
	  }
	  else if (remY < 0.0) {
	    delYInt = int(ceil(remY));//go up
	    remY = remY + 1.0;//find difference
	  }
	
	}
	
	//initialize to next location away from i,j
	k += delYInt;
        l += delXInt;

	while (( (k >= 0)&&(k < NR)) && ( (l >= 0)&&(l < MR) )) {
	  
	  //find difference locations for length calculations
	  delRow = k - i;
	  delCol = l - j;

	  tempZ = tempGrid[k][l];
	  tempTan = (tempZ - initZ)/(
		  pow(pow(dR,2.0)*
		  (pow(delRow,2.0) + pow(delCol,2.0)), 0.5));

	  //compare to previous maxTan
	  if (tempTan > maxTan) {
	    maxTan = tempTan;
	    maxRow = k;
	    maxCol = l;
	  }

	  //get the next delXInt, delYInt, remX, remY
	  remX += delXReal;
	  remY += delYReal;

	  //find integer x-step if needed
	  if (fabs(remX) >= 1.0) {
	    if (remX >= 0.0) {
	      delXInt = int(floor(remX));
	      remX = remX - 1.0;
	    }
	    else {
	      delXInt = int(ceil(remX));
	      remX = remX + 1.0;
	    }
	  }

	  //find integer y-step if needed
	  if (fabs(remY) >= 1.0) {
	    if (remY >= 0.0) {
	      delYInt = int(floor(remY));
	      remY = remY - 1.0;
	    }
	    else {
	      delYInt = int(ceil(remY));
	      remY = remY + 1.0;
	    }
	  }	    

	  //change location
	  k += delYInt;
	  l += delXInt;
	  
	}//end of line of interest (while-loop)	

	//set to HA grid	
	if (maxTan > 1e-5) {
		gridIn[i][j] = fabs(atan(maxTan));
      	}
	else {
		gridIn[i][j] = 0.0;
	}
      }//end of j-loop
    }//end of i-loop


    //cout << "resample temporary grid to node (we need to set up the node function)" << endl;
    //	  This is mostly taken from tResample
    //
    //	  Resample to polygons
    int counter = 0;
    for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {  
      eta -> initializeVCell(simCtrl, this, vXs[counter], vYs[counter], nPoints[counter]);

      varFromGrid[counter] = dummy;
      varFromGrid[counter] = eta->convertToVoronoiFormat(1);
      eta -> DestrtvCell();
      
      counter++;
    }
    
    //cout << "set HA to tCNode objects" << endl;
    counter = 0;
    for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
      
      switch ( tempInd ) {
	case 0:
	  cn->setHorAngle2700(varFromGrid[counter]);
	  break;
	case 1:
	  cn->setHorAngle2925(varFromGrid[counter]);
	  break;
	case 2:
	  cn->setHorAngle3150(varFromGrid[counter]);
	  break;
	case 3:
	  cn->setHorAngle3375(varFromGrid[counter]);
	  break;
	case 4:
	  cn->setHorAngle0000(varFromGrid[counter]);
	  break;
	case 5:
	  cn->setHorAngle0225(varFromGrid[counter]);
	  break;
	case 6:
	  cn->setHorAngle0450(varFromGrid[counter]);
	  break;
	case 7:
	  cn->setHorAngle0675(varFromGrid[counter]);
	  break;
	case 8:
	  cn->setHorAngle0900(varFromGrid[counter]);
	  break;
	case 9:
	  cn->setHorAngle1125(varFromGrid[counter]);
	  break;
	case 10:
	  cn->setHorAngle1350(varFromGrid[counter]);
	  break;
	case 11:
	  cn->setHorAngle1575(varFromGrid[counter]);
	  break;
	case 12:
	  cn->setHorAngle1800(varFromGrid[counter]);
	  break;
	case 13:
	  cn->setHorAngle2025(varFromGrid[counter]);
	  break;
	case 14:
	  cn->setHorAngle2250(varFromGrid[counter]);
	  break;
	case 15:
	  cn->setHorAngle2475(varFromGrid[counter]);
	  break;
	default:
	  cout << "\nCheck tempInd -- did not exist or assign" << endl;
      }//end-switch
	
      counter++;
    }

  //cout << "loop through next angle" << endl;  
  }//loop through next angle

  
  //set up shelter and lv factors
  for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {
    
    //compute slope and aspect for polygon
    elevation = cn->getZ();
    slope = fabs(atan(cn->getFlowEdg()->getSlope()));
    aspect = cn->getAspect();

    //compute sv
    sv = 0.0;
    //integrate
    for (int tempInd(0); tempInd < 16; tempInd++) {
      
      switch ( tempInd ) {
	 case 0:
	    horAngle = cn->getHorAngle0000();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 -horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 1:
	    horAngle = cn->getHorAngle0225();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle- sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 2:
	    horAngle = cn->getHorAngle0450();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 3:
	    horAngle = cn->getHorAngle0675();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 -horAngle))) *
			(3.1416/8);
	    break;
	  case 4:
	    horAngle = cn->getHorAngle0900();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 5:
	    horAngle = cn->getHorAngle1125();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 -horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 6:
	    horAngle = cn->getHorAngle1350();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 -horAngle))) *
			(3.1416/8);
	    break;
	  case 7:
	    horAngle = cn->getHorAngle1575();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 8:
	    horAngle = cn->getHorAngle1800();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 9:
	    horAngle = cn->getHorAngle2025();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 10:
	    horAngle = cn->getHorAngle2250();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 -horAngle))) *
			(3.1416/8);
	    break;
	  case 11:
	    horAngle = cn->getHorAngle2475();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 12:
	    horAngle = cn->getHorAngle2700();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 13:
	    horAngle = cn->getHorAngle2925();
	    sv += (cos(atan(slope))*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(atan(slope))*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 14:
	    horAngle = cn->getHorAngle3150();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  case 15:
	    horAngle = cn->getHorAngle3375();
	    sv += (cos(slope)*sin(3.1416/2 - horAngle)*sin(3.1416/2 - horAngle) +
			sin(slope)*cos(tempInd*22.5*3.1416/180 - aspect)) *
			0.5*(3.1416/2 - horAngle - sin(2*(3.1416/2 - horAngle))) *
			(3.1416/8);
	    break;
	  default:
	    cout << "\ncheck for-loop in sheltFactorFunc()" << endl;
      }//switch
    }//for-loop
    sv /= (2*3.1416);
    cn->setSheltFact(sv);

  }//nodes-loop

  }//shelter-on loop

  else {

    //set everything to 0 and calculate other places.
    sv = 0.0;
    lv = 0.0;
    horAngle = 0.0;
    for (cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP()) {

	cn->setSheltFact(0.0);
        cn->setHorAngle0000(0.0);
        cn->setHorAngle0225(0.0);
        cn->setHorAngle0450(0.0);
        cn->setHorAngle0675(0.0);
        cn->setHorAngle0900(0.0);
        cn->setHorAngle1125(0.0);
        cn->setHorAngle1350(0.0);
        cn->setHorAngle1575(0.0);
        cn->setHorAngle1800(0.0);
	cn->setHorAngle2025(0.0);
	cn->setHorAngle2250(0.0);
	cn->setHorAngle2475(0.0);
	cn->setHorAngle2700(0.0);
	cn->setHorAngle2925(0.0);
	cn->setHorAngle3150(0.0);
	cn->setHorAngle3375(0.0);

    }
  }
 
  //cout << "delete statements at the end of the shelter constructor" << endl;
  if (radSheltOpt > 0 && radSheltOpt < 4) {//CHANGED IN 2008
    for (int i = 0; i<NR; i++)
	    delete [] tempGrid[i];
    delete [] tempGrid;

    for (int i=0; i<NR; i++) {
       delete [] gridIn[i];
    }
    delete [] gridIn;
//    delete [] GridInPath;
	delete [] coorXG;
	delete [] coorYG;


  }

}//end of constructor

tShelter::~tShelter() {
  Cout << "tShelter Object Destroyed..." << endl;
}

