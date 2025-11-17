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
**  tKinemat.cpp: Functions for class tKinemat (see tKinemat.h)
**           A Finite-Element Kinematic Wave Routing Algorithm
**
***************************************************************************/

#include "src/tFlowNet/tKinemat.h"
#include "src/Headers/globalIO.h"

#ifdef PARALLEL_TRIBS
#include "src/tGraph/tGraph.h"
#include "src/tParallel/tParallel.h"
#endif

//=========================================================================
//
//
//                  Section 1: tKinemat Constructors/Destructors
//
//
//=========================================================================

/****************************************************************************
**
**  tKinemat::tKinemat()
**
**  Constructor for tRIBS model use
**
*****************************************************************************/
tKinemat::tKinemat(SimulationControl *sPtr, tMesh<tCNode> *gridRef, tInputFile &infile,
                   tRunTimer *timptr)
        : ais(NULL), bis(NULL), his(NULL), reis(NULL), siis(NULL), rifis(NULL), // what is the : doing here?
          sumis(NULL), C(NULL), Y1(NULL), Y2(NULL), Y3(NULL), clis(NULL),
          tFlowNet(sPtr, gridRef, infile, timptr) {

    simCtrl = sPtr;

    char fullName1[kMaxNameSize + 20];
    char fullName2[kMaxNameSize + 20];

    n = m = m1 = 0;
    cHead = cOutlet = nullptr;
    TimeSteps = 0;    // Time steps elapsed
    qit = Qin = H0 = Qout = dt = 0;

    ChannelConduc = TransientConduc = reis1 = Pchannel = Preach = 0.0; //ASM 2/8/2017
    CountLimit = Count = 0; //ASM
    PsiB = PoreInd = 0.0;//ASM
    NodeLoss = clis = NULL; // ASM 2/9/2017 clis stands for channel loss and is supposed to mimic the ais, bis etc.

    dtReff = 0.5;   // Hour, for lateral influx time increment

    kincoef = infile.ReadItem(kincoef, "KINEMVELCOEF");
    Roughness = infile.ReadItem(Roughness, "CHANNELROUGHNESS");

    percolationOption = infile.ReadItem(percolationOption, "OPTPERCOLATION"); // ASM 2/9/2017
    if(percolationOption==1) { //TODO WR 10062023 should not need if statement but currently CHANPOREINDEX fails assert in read function
        ChannelConduc = infile.ReadItem(ChannelConduc, "CHANNELCONDUCTIVITY"); // ASM
        TransientConduc = infile.ReadItem(TransientConduc, "TRANSIENTCONDUCTIVITY"); //ASM
        TransientTime = infile.ReadItem(TransientTime, "TRANSIENTTIME"); //ASM
        channelPorosity = infile.ReadItem(channelPorosity, "CHANNELPOROSITY"); // ASM
        ChanWidth = infile.ReadItem(ChanWidth, "CHANNELWIDTH"); //ASM temporary fix
        PoreInd = infile.ReadItem(PoreInd, "CHANPOREINDEX");//ASM
        PsiB = infile.ReadItem(PsiB, "CHANPSIB");//ASM
        //PsiB = PsiB/1000.; //ASM convert to m
    }
    IntStormMax = infile.ReadItem(IntStormMax, "INTSTORMMAX"); //ASM

    Cout << "\nChannel Characteristics:" << endl << endl;
    Cout << "Kinematic velocity coefficient: " << kincoef << endl;
    Cout << "Roughness coefficient: \t\t" << Roughness << endl;

    Width = infile.ReadItem(Width, "CHANNELWIDTHCOEFF");
    if (Width <= 0.0) {
        Width = infile.ReadItem(Width, "CHANNELWIDTH");
        Cout << "Channel width: \t\t\t" << Width << " m" << endl;
    } else {
        Cout << "Channel width: \t\t\tVariable" << endl;
        Width = -999.;
        AssignChannelWidths(infile);
    }

    // Allocate memory for stream reach outlet levels
    OutletHlev = new double[NodesLstO.getSize()]; // # of stream outlets
    assert(OutletHlev != 0);

    for (int i = 0; i < NodesLstO.getSize(); i++)
        OutletHlev[i] = 0.0; // Initialization of outlet levels

    // Open file for model run control
    infile.ReadItem(fullName1, "OUTHYDROFILENAME");
    strcat(fullName1, ".cntrl");

#ifdef PARALLEL_TRIBS
                                                                                                                            // Add processor extension if running in parallel
   char procex[10];
   snprintf( procex,sizeof(procex),".%-d", tParallel::getMyProc());//WR--09192023: 'sprintf' is deprecated: This function is provided for compatibility reasons only.
   strcat(fullName1, procex);
#endif

    ControlOut.open(fullName1);

    if (!ControlOut.good()) {
        cout << "\nWarning: Simulation control file not created... "
             << "\nExiting program..." << endl << flush;
        exit(2);
    }
    ControlOut.setf(ios::fixed, ios::floatfield);

    // Open file for model streamflow at the OutletNode
#ifndef PARALLEL_TRIBS
    // When running sequentially, open this file now
    infile.ReadItem(fullName2, "OUTHYDROFILENAME");
    strcat(fullName2, "_Outlet.qout");
    theOFStream.open(fullName2);

    if (!theOFStream.good()) {
        cout << "\nWarning: Output file not created.... "
             << "\nExiting Program..." << endl << flush;
        exit(2);
    }
    if (simCtrl->Header_label=='Y')
        theOFStream << "1-Time,hr\t " << "2-Qstrm,m3/s\t" << "3-Hlev,m" << "\n";
#endif

    // Allocate memory for stacks in stream nodes
    tCNode *cn;
    tMeshListIter<tCNode> nodIter(gridPtr->getNodeList());
    for (cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP()) {
        if (cn->getBoundaryFlag() == kStream) {
            cn->allocDataStack();
        }
    }
    OutletNode->allocDataStack();

    /******* Edits by JECR 2015 Start ******/  
    optres = 0;
    tempVariable = 0.0;
    resTimeStep = 0.0;
    resRunTime = 0.0;
    optres = ResReadItem.IterReadItem(infile, tempVariable, "OPTRESERVOIR");
    resTimeStep = ResReadItem.IterReadItem(infile, tempVariable, "TIMESTEP");
    resRunTime = ResReadItem.IterReadItem(infile, tempVariable, "RUNTIME");
    if (optres == 1) {
        initialize_values(infile, resTimeStep);
    } else {
        /*// TODO: debug
     * The way the reservoir setup is currently LevelPool is created in the tKinemat.h so that functions that can use LevelPool can later be used.
     * In the process the call creates tReservoir, which then calls pointers to two instances of tResData with reservoirType not having an address (least as far as I cant tell).
     * Later when LevelPool is destructed it also trys to destruct pointers from reservoirType--which as far as I can tell don't actually point to anything. In this case I set both pointers to null.
     * I think this works, but may not as I imagine, hence the todo above -WR
     */
        LevelPool.reservoirTypes = nullptr;
        LevelPool.reservoirNodes = nullptr;
    }
    /******** Edits by JECR 2015 End *******/

}
//
// /****************************************************************************
//**
//**  tKinemat::tKinemat()
//**
//**  Constructor for testing purposes
//**
//*****************************************************************************/
//tKinemat::tKinemat() : ais(NULL), bis(NULL), his(NULL), reis(NULL),
//siis(NULL), rifis(NULL), sumis(NULL), C(NULL),
//Y1(NULL), Y2(NULL), Y3(NULL)
//{
//	ControlOut.open("h-cntr.stream");
//
//	if ( !ControlOut.good() ) {
//		cout<<"\nWarning: Simulation control file not created... "
//		<<"\nExiting program..."<<endl<<flush;
//		exit(2);
//	}
//
//	GeomtFile.open("artif_chann.dat");
//
//	if (!GeomtFile) {
//		cout<<"\nError: File artif_chann.dat not found!\nExiting Program..."<<endl;
//		exit(2);
//	}
//
//	theOFStream.open("_Outlet.qout");
//	if ( !theOFStream.good() ) {
//		cout<<"\nWarning: Output file not created... "
//		<<"\nExiting program..."<<endl<<flush;
//		exit(2);
//	}
//	if (simCtrl->Header_label == 'Y')
//		theOFStream<<"1-Time,hr\t "<<"2-Qstrm,m3/s\t"<<"3-Hlev,m"<<"\n";
//}
//
// /****************************************************************************
//**
//**  tKinemat::tKinemat()
//**
//**  Constructor for testing purposes
//**
//*****************************************************************************/
//tKinemat::tKinemat(char *argv[]) : ais(NULL), bis(NULL), his(NULL),
//reis(NULL), siis(NULL), rifis(NULL),
//sumis(NULL), C(NULL), Y1(NULL),
//Y2(NULL), Y3(NULL)
//{
//	GeomtFile.open( argv[1] );
//	if (!GeomtFile) {
//		cout<<"\nError: File "<<argv[1]<<" not found!\nExiting Program..."<<endl;
//		exit(2);
//	}
//
//	theOFStream.open(argv[2]);
//	if ( !theOFStream.good() ) {
//		cout<<"\nWarning: Output file not created... "
//		<<"\nExiting program..."<<endl<<flush;
//		exit(2);
//	}
//	if (simCtrl->Header_label == 'Y')
//		theOFStream<<"1-Time,hr\t "<<"2-Qstrm,m3/s\t"<<"3-Hlev,m"<<endl<<flush;
//
//	ControlOut.open("h_cntr.stream");
//
//	if ( !ControlOut.good() ) {
//		cout<<"\nWarning: Simulation control file not created... "
//		<<"\nExiting program..."<<endl<<flush;
//		exit(2);
//	}
//
//	n = m = m1 = 0;
//}


/****************************************************************************
**
**  tKinemat::~tKinemat()
**
**  Destructor
**
*****************************************************************************/
tKinemat::~tKinemat() {
    FreeMemory();
    delete[] OutletHlev;
#ifdef PARALLEL_TRIBS
                                                                                                                            // The theOFStream only exists on the
  // processor containing the last reach
  if (tGraph::hasLastReach())
#endif
    theOFStream.close();
    ControlOut.close();
    GeomtFile.close();

    // Deallocate memory for stacks in stream nodes //WR debug
    tCNode *cn;
    tMeshListIter<tCNode> nodIter(gridPtr->getNodeList());
    for (cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP()) {
        if (cn->getBoundaryFlag() == kStream) {
            cn->deleteDataStack();
        }
    }
    OutletNode->deleteDataStack();

    Cout << "tKinemat Object has been destroyed..." << endl << flush;
}

/*****************************************************************************
**
**  tKinemat::UpdateForNewRun
**
**  Used to update data members of tKinemat and tCNode when a new simulation
**  run is to be carried out (option -ON of the command line)
**
*****************************************************************************/
void tKinemat::UpdateForNewRun(tInputFile &infile, int keep) {
    char fullName1[kMaxNameSize + 20];
    char fullName2[kMaxNameSize + 20];

    tCNode *cn;
    tList<int> *TimeInd;  // Ptr to Time steps stack
    tList<double> *Qeff;     // Ptr to Qeff stack
    tMeshListIter<tCNode> nodIter(gridPtr->getNodeList());

    n = m = m1 = 0;
    cHead = cOutlet = NULL;
    TimeSteps = 0;    // Time steps elapsed
    qit = Qin = H0 = Qout = 0.0;
    ChannelConduc = TransientConduc = ReachLoss = 0.0; // ASM 2/8/2017

    kincoef = infile.ReadItem(kincoef, "KINEMVELCOEF");
    Roughness = infile.ReadItem(Roughness, "CHANNELROUGHNESS");

    Cout << "\nChannel Characteristics:" << endl << endl;
    Cout << "Kinematic velocity coefficient: " << kincoef << endl;
    Cout << "Roughness coefficient: \t\t" << Roughness << endl;

    // Depending on whether a uniform or variable
    // channel width is used -- different options
    Width = infile.ReadItem(Width, "CHANNELWIDTHCOEFF");
    if (Width <= 0.0) {
        Width = infile.ReadItem(Width, "CHANNELWIDTH");
        Cout << "Channel width: \t\t\t" << Width << " m" << endl;
    } else {
        Cout << "Channel width: \t\t\t is variable" << endl;
        Width = -999.;
        AssignChannelWidths(infile);
    }

    // Close the file and then re-open it
    ControlOut.close();

    // Open file for model run control
    infile.ReadItem(fullName1, "OUTHYDROFILENAME");
    strcat(fullName1, ".control");
    ControlOut.open(fullName1);

    if (!ControlOut.good()) {
        cout << "\nWarning: Simulation control file not created... "
             << "\nExiting program..." << endl << flush;
        exit(2);
    }
    ControlOut.setf(ios::fixed, ios::floatfield);

    // Close the file and then re-open it
    theOFStream.close();

    // Open file for model streamflow at the OutletNode
    infile.ReadItem(fullName2, "OUTHYDROFILENAME");
    strcat(fullName2, "_Outlet.qout");
    theOFStream.open(fullName2);

    if (!theOFStream.good()) {
        cout << "\nWarning: Output file not created... "
             << "\nExiting program..." << endl << flush;
        exit(2);
    }
    if (simCtrl->Header_label=='Y')
        theOFStream << "1-Time,hr\t " << "2-Qstrm,m3/s\t" << "3-Hlev,m" << "\n";


    // If a decision made to keep the state vars --
    // dont change anything, keep vars from previous run
    if (!keep) {
        for (int i = 0; i < NodesLstO.getSize(); i++)
            OutletHlev[i] = 0.0; // Initialization of outlet levels

        for (cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP()) {

            // Initialize the stream variables to '0'
            if (cn->getBoundaryFlag() == kStream) {
                cn->setHlevel(0.0);
                cn->setQstrm(0.0);
                cn->percOccur = 0.0; //ASM set initial percOccur to 0

                // Flush memory for stacks in stream nodes
                TimeInd = cn->getTimeIndList();
                Qeff = cn->getQeffList();
                TimeInd->Flush();
                Qeff->Flush();
            }
        }
        // Do the same for the basin outlet
        OutletNode->setHlevel(0.0);
        OutletNode->setQstrm(0.0);
        TimeInd = OutletNode->getTimeIndList();
        Qeff = OutletNode->getQeffList();
        TimeInd->Flush();
        Qeff->Flush();
    }
    return;
}

//=========================================================================
//
//
//                  Section 2: tKinemat Functions
//
//
//=========================================================================


/*************************************************************************
**
**  AssignChannelWidths
**
**  The function loops through the active node list to find stream nodes.
**  The procedure assigns the width of the channel using either
**  geomorphological functions or measured cross sections. In the
**  latter case, the procedure interpolates between the measured values
**  and uses geomorphological relationships for channels where there are
**  no measurements.
**
*************************************************************************/
void tKinemat::AssignChannelWidths(tInputFile &infile) {
    int i, cnt, flag, numXsections, wOption;
    char baseName[kMaxNameSize];
    double value, tempo, AccLgth, dW;
    double WCoeff, WExpnt;

    // Read in the parameter values for geomorph relations (power law)
    WCoeff = infile.ReadItem(WCoeff, "CHANNELWIDTHCOEFF");
    WExpnt = infile.ReadItem(WExpnt, "CHANNELWIDTHEXPNT");

    int *nodeList;
    double *widthList;
    tCNode *cn, *cnn, *cmove, *cprev, *cUpstream;
    tEdge *firstedg, *curedg;
    tMeshListIter<tCNode> nodIter(gridPtr->getNodeList());

    // Option 8 version looped on active nodes only because tGraph update
    // was called after this code.  With meshbuilder, reach nodes not in
    // this partition are moved to the back before this, so we have to
    // loop on all nodes MESHB
    int option = infile.ReadItem(option, "OPTMESHINPUT");

    if (option == 9) {
        for (cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP()) {
            if (cn->getBoundaryFlag() == kStream) {
                value = cn->getContrArea() * 1.0E-6;  // to give units of km^2
                value = WCoeff * pow(value, WExpnt);
                cn->setChannelWidth(value);
                cn->ActivateSortTracer(); // tracer is assigned to '1'
            }
        }
    } else {
        for (cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP()) {
            if (cn->getBoundaryFlag() == kStream) {
                value = cn->getContrArea() * 1.0E-6;  // to give units of km^2
                value = WCoeff * pow(value, WExpnt);
                cn->setChannelWidth(value);
                cn->ActivateSortTracer(); // tracer is assigned to '1'
            }
        }
    }
    value = OutletNode->getContrArea() * 1.0E-6;
    value = WCoeff * pow(value, WExpnt);
    OutletNode->setChannelWidth(value);

    // Now, try to access a file that contains cross section measurements.
    // If it does not exist use geomorphological relationships only
    infile.ReadItem(baseName, "CHANNELWIDTHFILE");   // output basename
    if (strlen(baseName) == 0) {
        Cout << "\nChannel widths determined from geomorphic relations..." << endl;
        return;
    }

    ifstream InWidthFile(baseName);
    if (!InWidthFile) {
        Cout << "\nChannel widths determined from geomorphic relations..." << endl;
        return;
    }

    // Read in width interpolation option
    // '0' - interpolation between measured and computed
    //       width values is enabled
    // '1' - interpolation between measured and computed
    //       width values is disabled (interpolation between
    //       measured width values only)

    wOption = infile.ReadItem(wOption, "WIDTHINTERPOLATION");

    InWidthFile >> numXsections;

    nodeList = new int[numXsections];
    widthList = new double[numXsections];

    // Complicated if X-s and Y-s inputted, easier to deal with node ID-s
    for (i = 0; i < numXsections; i++) {
        InWidthFile >> nodeList[i] >> widthList[i];
    }

    // Put all the corresponding nodes in the stack
    tPtrList<tCNode> NodesLst;
    tPtrListIter<tCNode> NodesIter(NodesLst);

    cnt = 999999999;
    if (nodeList) {  // If the node list in NOT empty
        for (i = 0; i < numXsections; i++) {

            if (nodeList[i] > 0) {
                for (cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP()) {
                    if (cn->getID() == nodeList[i]) {
                        NodesLst.insertAtBack(cn);
                        cn->setChannelWidth(widthList[i]);
                        cn->DeactivateTracer(); //Deactivate tracer
                        if (cn->getID() < cnt) {
                            cnt = cn->getID();
                            cUpstream = cn; // The most Upstream node
                        }
                    }
                }
                // Check if there is a measured X-section for outlet
                if (OutletNode->getID() == nodeList[i]) {
                    OutletNode->setChannelWidth(widthList[i]);
                    NodesLst.insertAtBack(OutletNode);
                }
            }
        }
    }

    for (cn = NodesIter.FirstP(); !(NodesIter.AtEnd()); cn = NodesIter.NextP()) {
        // ------------------------------------------------------------
        // Now, do channel traverse going DOWNSTREAM:
        // 1) locate the upper node with X-section measurements
        // 2) locate the lower node with X-section measurements OR
        //    check if there is a confluence point below the upper node
        // 3) interpolate widths between the two nodes
        // ------------------------------------------------------------
        if (cn != OutletNode) {
            cprev = cmove = cn;
            flag = 0;
            cnt = 0;
            AccLgth = 0.;
            //  Go downstream until you reach a confluence/outlet node
            while (!flag) {
                cprev = cmove; // Always keeps track of the previous node
                AccLgth += cmove->getFlowEdg()->getLength();
                cmove = cmove->getDownstrmNbr(); // Downstream node...

                if (!wOption)
                    flag = IsConfluence(cmove, cprev);
                else {
                    if (cmove == OutletNode)
                        flag++;
                }

                if (cmove->NoMoreTracers())
                    flag++;
                cnt++;
            }
            // Now, interpolate widths between 'cn' and 'cmove'
            if (cnt > 1) {
                dW = cmove->getChannelWidth() - cn->getChannelWidth();
                cprev = cmove;
                cmove = cn;
                tempo = 0.;

                while (cmove != cprev) {
                    tempo += cmove->getFlowEdg()->getLength();
                    cmove = cmove->getDownstrmNbr();
                    value = dW * tempo / AccLgth + cn->getChannelWidth();
                    cmove->setChannelWidth(value);
                }
            }
        }
        // ------------------------------------------------------------
        // Now, do channel traverse going UPSTREAM:
        // 1) locate the lower node with X-section measurements
        // 2) locate the upper node with X-section measurements OR
        //    check if there is a confluence point below the upper node
        // 3) interpolate widths between the two nodes
        // ------------------------------------------------------------
        if ((!wOption || cn == cUpstream) && cn != OutletNode) {
            cprev = cmove = cn;
            cnt = 0;
            AccLgth = 0.;

            flag = IsStreamHead(cmove);
            //  Go upstream until you reach a confluence/stream head node
            while (!flag) {
                cprev = cmove; // <-- Always keeps track of the previous node
                firstedg = cmove->getFlowEdg();
                curedg = firstedg->getCCWEdg();
                while (curedg != firstedg) {
                    cmove = (tCNode *) curedg->getDestinationPtrNC();
                    // Check if it is a stream node
                    if (cmove->getBoundaryFlag() == kStream) {
                        if (cmove->getFlowEdg()->getDestinationPtrNC() == (tNode *) cprev)
                            curedg = firstedg;
                        else
                            curedg = curedg->getCCWEdg();
                    } else
                        curedg = curedg->getCCWEdg();
                }
                if (cn == OutletNode)
                    TellAboutNode(cmove);

                // Check if 'cprev' is confluence in the case
                // if the option is ON in the .in file
                // OR if node 'cn' is the most upstream
                flag = IsConfluence(cmove, cprev);

                // Go downstream one link if
                // 'cprev' is a confluence
                if (flag)
                    cmove = cprev;
                else {
                    AccLgth += cmove->getFlowEdg()->getLength();
                    // Check if it is a stream head
                    flag = IsStreamHead(cmove);

                    // Check if there is a measured X-section
                    if (cmove->NoMoreTracers())
                        flag++;
                    cnt++;
                }
            }

            if (cnt > 1) {
                dW = cmove->getChannelWidth() - cn->getChannelWidth();
                cnn = cmove;
                cprev = cmove = cn;
                tempo = 0.;

                while (cmove != cnn) {
                    cprev = cmove;
                    firstedg = cmove->getFlowEdg();
                    curedg = firstedg->getCCWEdg();
                    while (curedg != firstedg) {
                        cmove = (tCNode *) curedg->getDestinationPtrNC();
                        // Check if it is a stream node
                        if (cmove->getBoundaryFlag() == kStream) {
                            if (cmove->getFlowEdg()->getDestinationPtrNC() == (tNode *) cprev)
                                curedg = firstedg;
                            else
                                curedg = curedg->getCCWEdg();
                        } else
                            curedg = curedg->getCCWEdg();
                    }
                    tempo += cmove->getFlowEdg()->getLength();
                    value = dW * tempo / AccLgth + cn->getChannelWidth();
                    cmove->setChannelWidth(value);
                }
            }
        }
    }

    delete[] nodeList;
    delete[] widthList;
    NodesLst.Flush();
    return;
}

/*****************************************************************************
**
**  tKinemat::AllocateMemory
**
**  Allocates memory for all arrays used in the kinematic routing
**
*****************************************************************************/
void tKinemat::AllocateMemory(int NN) {
    n = NN;
    bis = new double[n];
    assert(bis != 0);
    his = new double[n];
    assert(his != 0);
    reis = new double[n];
    assert(reis != 0);
    clis = new double[n];    // ASM 2/10/2017 (2 lines)
    assert(clis != 0);
    NodeLoss = new double[n]; // ASM 2/17/17 (2 lines)
    assert(NodeLoss != 0);

    // The following three vectors contain information only for stream
    // links (n-1) in total but it is more convenient to use size 'n' instead

    ais = new double[n];    // It's n, not n-1
    assert(ais != 0);
    siis = new double[n];   // It's n, not n-1
    assert(siis != 0);
    rifis = new double[n];  // It's n, not n-1
    assert(rifis != 0);

    sumis = new double[n - 1];
    assert(sumis != 0);

    C = new double[n];
    assert(C != 0);
    Y1 = new double[n - 2];
    assert(Y1 != 0);
    Y2 = new double[n - 1];
    assert(Y2 != 0);
    Y3 = new double[n - 1];
    assert(Y3 != 0);

    return;
}

/*****************************************************************************
**
**  tKinemat::FreeMemory
**
**  Deallocates memory for all arrays used in the kinematic routing
**
*****************************************************************************/
void tKinemat::FreeMemory() {
    if (ais != NULL) delete[] ais;
    if (bis != NULL) delete[] bis;
    if (his != NULL) delete[] his;
    if (reis != NULL) delete[] reis;
    if (siis != NULL) delete[] siis;
    if (rifis != NULL) delete[] rifis;
    if (sumis != NULL) delete[] sumis;
    if (C != NULL) delete[] C;
    if (Y1 != NULL) delete[] Y1;
    if (Y2 != NULL) delete[] Y2;
    if (Y3 != NULL) delete[] Y3;
    if (clis != NULL) delete[] clis;
    if (NodeLoss != NULL) delete[] NodeLoss;

    ais = NULL;
    bis = NULL;
    his = NULL;
    reis = NULL;
    siis = NULL;
    rifis = NULL;
    sumis = NULL;
    C = NULL;
    Y1 = NULL;
    Y2 = NULL;
    Y3 = NULL;
    clis = NULL; // ASM 2/10/2017
    NodeLoss = nullptr;

    return;
}

/*****************************************************************************
**
**  tKinemat::SurfaceFlow()
**
**  Runs the Hydrologic-Kinematic Routing model
**  First, it runs the hydraulic/kinematic routing model that also contains
**  hillslope hydrologic routing model. Secondly, it runs the older hydrologic
**  routing model for the whole basin. At last, the function writes the outlet
**  streamflow value into 'theOFStream'.
**
**
*****************************************************************************/
void tKinemat::SurfaceFlow() {
    int check, it;
    it = 0; 
    Pchannel = TotChanLength = ParallelPerc = 0.0;
    RunRoutingModel(it, &check, timer->getTimeStep() * 3600.);

    tFlowNet::SurfaceFlow();
    double currentTime = timer->getCurrentTime();

#ifdef PARALLEL_TRIBS
    // If running in parallel, only partition with last reach writes
    if (tGraph::hasLastReach())
#endif
    // CJC2025: Use stream manipulators for proper floating point formatting
    theOFStream << std::fixed << std::setprecision(4) << currentTime
                << "\t" << Qout << "\t" << OutletNode->getHlevel() << "\n";
    return;
}

/*****************************************************************************
*************************** ADDED BY JECR 2015 ******************************/
void tKinemat::initialize_values(tInputFile &infile, double resTime) {
    LevelPool.setModelTimeStep(resTime * 60.0); // Converted to seconds
    int resArray = resRunTime / (resTimeStep / 60); // Converted to hours
    LevelPool.setResArraySize(resArray);
    LevelPool.SetResVariables(infile);
    LevelPool.SetResNodes(infile);
}


void tKinemat::Reservoir_Routing(int Rnode) {
    int reservoirNum = LevelPool.getNReservoirs(); // Number of Reservoirs to set the loop.

    for (int x = 0; x < reservoirNum; x++) {
        checkID = LevelPool.getNodetKinemat(x); // Read node ID from User's tables
        checkNode = Rnode; // Read node used currently by the model
        if (checkNode == checkID) {
            LevelPool.setCurrResNode(x);
            int typeID = LevelPool.getTypetKinemat(x);
            LevelPool.setCurrResType(typeID);
            LevelPool.RunLevelPoolRouting(Qin);
            Qin = LevelPool.getResDischargeOut();
            qit = 0.5 * Qin;
            H0 = pow(Qin * rifis[0] / (bis[0] * sqrt(siis[0])), 3. / 5.);
            break; // Skip to Line 834 (Prevents value of CheckID to change)
        } else {
            continue;
        }
    }
}
/*************************** END BY JECR 2015 ********************************
*****************************************************************************/


/*****************************************************************************
**
**  tKinemat::RunRoutingModel
**
**  Calls computational functions of the routing model
**  'it' is required only for an off-line test
**  'timeStep' is assumed to be in SECONDS
**
*****************************************************************************/
void tKinemat::RunRoutingModel(int it, int *check, double timeStep) {
    tCNode *cn;
    tPtrListIter<tCNode> NodesIterO(NodesLstO);
    tPtrListIter<tCNode> NodesIterH(NodesLstH);
    tListIter<int> NNodesIter(NNodes);
    int timerRes = 0;
    // Run hillslope routing model first
    RunHydrologicRouting();


    dt = timeStep;  // Computational time step

    // Update the counter for transient conditions
    if (Preach > 0.1)
        Count += 1;
    else
        Count = 0;

    // Loop through all outlets and set Q to '0' . We need to do that
    // in order to properly assign Q-s in confluence nodes

    for (cn = NodesIterO.FirstP(); !(NodesIterO.AtEnd()); cn = NodesIterO.NextP())
        cn->setQstrm(0.0);

    // Loop through all stream reaches id - Stream reach ID

    for (cn = NodesIterH.FirstP(), NodesIterO.First(), NNodesIter.First(), id = 0;
         !(NodesIterH.AtEnd());
         cn = NodesIterH.NextP(), NodesIterO.Next(), NNodesIter.Next(), id++) {

#ifdef GRAPH_TRIBS
                                                                                                                                // Process only stream reaches in local partition
    if (tGraph::inLocalPartition(id)) {
#endif

        // Initialize head and outlet for a current stream reach
        cHead = NodesIterH.DatPtr();
        cOutlet = NodesIterO.DatPtr();

        //can calculate the number of time steps to check for transient period here ASM
        CountLimit = TransientTime / (dt / 60);

        // Initialize widths, lengths, slopes, levels, C, Y1, Y2, Y3
        InitializeStreamReach(NNodesIter.DatRef(), CountLimit);

        // Initialize lateral influx array
        AssignLateralInflux();

        // Initialize upper BND condition
        AssignQin();

        // Check Reservoir Option
        if (optres == 1) {
            Reservoir_Routing(cn->getID()); // JECR 2015
        }


        // Assign reach percolation to total channel percolation ASM 2/10/2017
        Pchannel += Preach;

        // Run kinematic wave routing model
        if (NNodesIter.DatRef() == 2)
            SolveForTwoNodeReach(C, Y1, Y2, Y3, reis, his, qit, H0);
        else
            KinematWave(C, Y1, Y2, Y3, reis, his, qit, H0, check);

        // Update computed values of levels & Qs

/********************* Start of modifications by JECR 2015 **************************/
        if ((optres == 1) && (checkNode == checkID)) {
            Qout = LevelPool.getResDischargeOut(); // Skip ComputeQout();
        } else {
/********************** End of modifications by JECR 2015 ***************************/

            ComputeQout();
        }

        UpdateStreamVars();

        // Write control info
        it = TimeSteps;

        // Memory management
        FreeMemory();

#ifdef GRAPH_TRIBS
                                                                                                                                // End of processing stream reaches in local partition
    }
#endif
    }

    // Close file with reach info
    if (TimeSteps == 0) ControlOut.close();

    TimeSteps++;
    return;
}

/*****************************************************************************
**
**  tKinemat::RunHydrologicRouting()
**
**  Pseudo-Routing code
**  Loops through the active nodes: for _Hillslope_ nodes:
**     - Define current discharge in the appropriate stream node and
**       the corresponding stream velocitiy
**     - Define hillsope velocities at this moment
**     - Define the corresponding travel time
**     - Define the corresponding time step from the beginning
**       of the run at which runoff generated in a hillslope node
**       will show up in the stream node according to the set velocities
**     - Get the generated runoff volume
**     * Do the same for a _Stream_ node assuming zero travel time
**     - Store the volume in an appropriate place in stack
**       through the stream node stack
**
*****************************************************************************/
void tKinemat::RunHydrologicRouting() {
    tCNode *cn;
    tList<int> *TimeInd;  // Ptr to Time steps stack
    tList<double> *Qeff;     // Ptr to Qeff stack
    tListIter<int> IndIter;
    tListIter<double> QIter;
    tListNode<int> *indx;
    tListNode<double> *Qvalue;

    tMeshListIter<tCNode> nodIter(gridPtr->getNodeList());

    int nSStep;
    double ttime;         // travel time for a current node, SECONDS
    double vRunoff = 0.0; // runoff volume, m^3
    double Area = 0.0;    // Voronoi cell area

    if (simCtrl->Verbose_label == 'Y') {
        cout << "\nHillslope Routing scheme is running..." << endl;
    }

    int ll = 0;
    for (cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP()) {

        // Check first if node has produced runoff
        if (cn->getSrf() > 0.0) {

            Area = cn->getVArea();             // M^2
            vRunoff = cn->getSrf() * Area / 1000.; // Srf, MM --> M^3

            if (cn->getSrf() > 999999) {
                cout << "\nWarning: tKinemat: vRunoff "
                     << "> 999999: id = " << cn->getID() << "\n\t\t->Assigned to zero\n";
                vRunoff = 0.;
            }

            // 1.) Consider #Hillslope# separately from Stream
            if (cn->getBoundaryFlag() == kNonBoundary) {

                // Take an average discharge between the stream to which
                // it points to and its downstream neighbor to set the
                // travel velocity from that node
                if (cn->getStreamNode() == OutletNode)
                    setTravelVelocityKin(cn->getStreamNode()->getQstrm(),
                                         cn->getStreamNode()->getContrArea());
                else
                    setTravelVelocityKin((cn->getStreamNode()->getQstrm() +
                                          cn->getStreamNode()->getDownstrmNbr()->getQstrm()) / 2.,
                                         cn->getStreamNode()->getContrArea());

                // Update travel time if (flowexp >= 0)
                cn->setTTime(cn->getHillPath() / hillvel / 3600.);
                ttime = cn->getTTime(); //travel time for a current node, HOURS

                // Compute the corresponding # of time step for runoff shows up
                // in the node. The second argument is duration of time interval
                // to which the streamflow is related, i.e., 0.5, 1, 2 hour
                nSStep = timer->getStepForSpecifiedDT(timer->getCurrentTime() + ttime, dtReff);

                //  If the value falls in the same "box"
                if (timer->getStepForSpecifiedDT(timer->getCurrentTime(), dtReff) == nSStep) {

                    // If current time is the end of Reff interval
                    if ((int) ceil(timer->getCurrentTime() / dtReff) ==
                        (int) floor(timer->getCurrentTime() / dtReff))
                        vRunoff /= (timer->getTimeStep() * 3600.);
                        // If current time is not the end of Reff interval
                    else
                        vRunoff /= ((nSStep * dtReff - timer->getCurrentTime()
                                     + timer->getTimeStep()) * 3600.);
                } else {
                    // M^3 --> M^3/S
                    vRunoff /= (dtReff * 3600.);
                }

                // Get appropriate lists
                TimeInd = cn->getStreamNode()->getTimeIndList();
                Qeff = cn->getStreamNode()->getQeffList();
            }


                // 2.) Consider #Stream# node now
            else if (cn->getBoundaryFlag() == kStream) {
                ttime = 0.;
                nSStep = timer->getStepForSpecifiedDT(timer->getCurrentTime(), dtReff);

                // If current time is the end of Reff interval
                if ((int) ceil(timer->getCurrentTime() / dtReff) ==
                    (int) floor(timer->getCurrentTime() / dtReff))
                    vRunoff /= (timer->getTimeStep() * 3600.);
                    // If current time is not the end of Reff interval
                else
                    vRunoff /= ((nSStep * dtReff - timer->getCurrentTime()
                                 + timer->getTimeStep()) * 3600.);

                // Get appropriate lists
                TimeInd = cn->getTimeIndList();
                Qeff = cn->getQeffList();
            }

            // 3.) Now, work with the stacks
            assert(TimeInd != 0);
            assert(Qeff != 0);

            if (TimeInd->getSize() > 0) {

                // Initialize corresponding iterators
                IndIter.Reset(*TimeInd);
                QIter.Reset(*Qeff);

                // Check First if the last is less than nSStep
                IndIter.Last();
                if (IndIter.DatRef() < nSStep) {
                    TimeInd->insertAtBack(nSStep);
                    Qeff->insertAtBack(vRunoff);
                } else {
                    // Start from the first element in stack
                    // Find the one which corresponds to nSStep
                    IndIter.First();
                    QIter.First();
                    while (!(IndIter.AtEnd()) && IndIter.DatRef() < nSStep) {
                        IndIter.Next();
                        QIter.Next();
                    }

                    // If equal, just add to it
                    if (IndIter.DatRef() == nSStep) {
                        *(QIter.DatPtr()) += vRunoff;
                    }

                        // Now, it is _larger_ than nSStep put
                        // the value in the stack BEFORE
                    else {
                        indx = IndIter.NodePtr();
                        Qvalue = QIter.NodePtr();
                        TimeInd->insertAtPrev(nSStep, indx);
                        Qeff->insertAtPrev(vRunoff, Qvalue);
                    }
                }
            }

                // If there are no data in stack yet
            else {
                TimeInd->insertAtBack(nSStep);
                Qeff->insertAtBack(vRunoff);

                // For LATER versions: to start looking for
                // an appropriate place in stack...
            }
        } // (SRF > 0)

            // If runoff has not been generated update the flow field for the node
        else {

            // If this is a hillslope node
            if (cn->getBoundaryFlag() == kNonBoundary) {
                if (cn->getStreamNode() == OutletNode)
                    setTravelVelocityKin(cn->getStreamNode()->getQstrm(),
                                         cn->getStreamNode()->getContrArea());
                else
                    setTravelVelocityKin((cn->getStreamNode()->getQstrm() +
                                          cn->getStreamNode()->getDownstrmNbr()->getQstrm()) / 2.,
                                         cn->getStreamNode()->getContrArea());
            }
        }
        cn->setFlowVelocity(hillvel);

        ll++;
    }
    return;
}

/*****************************************************************************
**
**  tKinemat::setTravelVelocityKin(double)
**
**  Sets travel velocities depending on the discharge at a stream node
**  'curr_discharge'. 'CArea' is the corresponding contributing area
**  for the stream node at which streamflow is evaluated
**
*****************************************************************************/
void tKinemat::setTravelVelocityKin(double curr_discharge, double CArea) {
    // If the discharge is zero, set velocity to some value
    if (!curr_discharge) {

        if (!baseflow && flowexp) {
            cout << "\n\n\tWARNING!!! Baseflow is zero and thus the lower "
                 << "limit of stream velocity is undefined --> Set to 0.1" << endl;
            baseflow = 0.1;
        }
        // Assume linear scaling of baseflow value with the contributing area
        flowout = baseflow * CArea / OutletNode->getContrArea();
    }
        // If discharge in stream is non-zero, use its
        // value to compute the hillslope velocity
    else
        flowout = curr_discharge;

    // Non-linear model
    if (flowexp > 0.0) {
        hillvel = kincoef * pow(flowout / CArea, flowexp);
    }
        // Linear model
    else
        hillvel = velcoef / velratio;

    return;
}

/*****************************************************************************
**
**  tKinemat::InitializeStreamReach()
**
**  Assigns all appropriate information about the stream reach to be modeled
**  'n' - is the total # of stream nodes including the upper boundary
**
*****************************************************************************/
void tKinemat::InitializeStreamReach(int NN, int CountLimit) {
    int i;
    double Slope;
    tCNode *cmove, *cend;
    double ChanLength = TotWidth = 0.0; // ASM



    n = NN;        // # of nodes
    m = n - 1;
    m1 = n - 2;

    AllocateMemory(n);

    maxH = maxReff = 0.0;

    i = 0;
    cmove = cHead;  // Points to the current stream head
    cend = cOutlet; // Point to the current outlet

    // Approximate the slope for node '0' as the
    // slope of the link downstream
    // If the slope is < 0, approximate with the
    // "error" slope = 0.5ft/30m = 0.152/30
    siis[i] = cmove->getFlowEdg()->getSlope();
    if (siis[i] <= 0)
        siis[i] = 0.0050667;

    while (cmove != cend) {
        his[i] = cmove->getHlevel();
        bis[i] = cmove->getChannelWidth();
        ais[i] = cmove->getFlowEdg()->getLength();
        rifis[i] = cmove->getRoughness();
        //ASM 2/9/2017
        if (percolationOption == 1) {
            //setCoeffstest(cmove);
            //poro = soilPtr->getSoilProp(9);  // Surface hydraulic conductivity
            NodeLoss[i] = bis[i] * ais[i] * ChannelConduc * channelPorosity; // ASM testporo; w*l*poro*ksat [m3/s]
            ChanLength += ais[i]; // ASM
        } else if (percolationOption == 2) {
            // Need to get time information here
            if (Count > CountLimit - 1) {
                NodeLoss[i] = bis[i] * ais[i] * ChannelConduc * channelPorosity;
                ChanLength += ais[i];
            } else {
                NodeLoss[i] = bis[i] * ais[i] * TransientConduc * channelPorosity;
                ChanLength += ais[i];
            }
        }
        //end ASM edits

        Slope = cmove->getFlowEdg()->getSlope();

        //"Error" slope = 0.5ft/30m = 0.152/30
        if (Slope <= 0)
            siis[i + 1] = 0.0050667;
        else
            siis[i + 1] = Slope;

        if (his[i] > maxH)
            maxH = his[i];

        //PrintFlowStacks(ControlOut, cmove);

        cmove = cmove->getDownstrmNbr();
        i++;
    }

    TotChanLength += ChanLength; //ASM adds

    // Special care has to be taken regarding the outlet nodes
    // Use the separately stored outlet level from time step (t-1)
    his[i] = OutletHlev[id];           // cmove->getHlevel();
    bis[i] = cmove->getChannelWidth();
    rifis[i] = cmove->getRoughness();
    if (his[i] > maxH)
        maxH = his[i];

    // PARAMETERS
    // Slope = 0.00061140520;
    // Slope = 1.0;           // degree slope
    // Slope = tan(Slope*(atan(1)*4)/180);

    for (i = 0; i < n; i++) {
        rifis[i] = Roughness;

        // In case if uniform slope is desired
        // siis[i]  = Slope;

        // Use uniform width if desired
        if (Width > 0.)
            bis[i] = Width;
    }

    ComputeCoefficientArrays();


    // Output control
    if (TimeSteps == 0) {
        ControlOut << "## REACH ID = " << id + 1 << " ##" << "\n";
        ControlOut << "- WIDTH -\t";
        ControlPrint(ControlOut, bis, n);
        ControlOut << "- LENGTH -\t";
        ControlPrint(ControlOut, ais, m);
        ControlOut << "- ROUGHNESS -\t";
        ControlPrint(ControlOut, rifis, n);
        ControlOut << "- SLOPE -\t";
        ControlPrint(ControlOut, siis, n);
        ControlOut << "- C -\t";
        ControlPrint(ControlOut, C, n);
        ControlOut << "- Y1 -\t";
        ControlPrint(ControlOut, Y1, n - 2);
        ControlOut << "- Y2 -\t";
        ControlPrint(ControlOut, Y2, n - 1);
        ControlOut << "- Y3 -\t";
        ControlPrint(ControlOut, Y3, n - 1);
    }

    return;
}

/*****************************************************************************
**
**  tKinemat::ComputeCoefficientArrays()
**
**  Computes four coeeficient arrays required in the system of equations
**  It is assumed that the arrays ais, bis, siis, rifis containing lengths
**  widths, slopes, and roughness coefficients have been already initialized.
**  This is the only place where dt is used.
**
*****************************************************************************/
void tKinemat::ComputeCoefficientArrays() {
    int i;
    // Computation of coefficient arrays, hydraulic coef_s
    for (i = 0; i < n; i++) {

        //Computation of  Chezy friction coefficient Ci*Bi
        C[i] = 1. / rifis[i] * sqrt(siis[i]) * bis[i];

        if (i < m1) {
            Y1[i] = 1. / 6. * ais[i + 1] * bis[i + 2] / dt; // <-it starts from 2nd point
            Y2[i] = 1. / 3. * (ais[i] + ais[i + 1]) * bis[i + 1] / dt;
            Y3[i] = 1. / 6. * ais[i] * bis[i] / dt;
        } else if (i == m1) {
            Y2[m1] = 1. / 3. * ais[m1] * bis[m] / dt;
            Y3[m1] = 1. / 6. * ais[m1] * bis[m1] / dt;
        }
    }
    return;
}

/*****************************************************************************
**
**  tKinemat::AssignLateralInflux()
**
**  This option so far is only appropriate for off-line simulation, in the
**  tRIBS everything should be taken care of automatically through assigning
**  the influx to nodes/stream links
**
*****************************************************************************/
void tKinemat::AssignLateralInflux() {
    int i;
    tCNode *cmove, *cend;
    double reis1; //ASM 2/10/2017
    Preach = 0.0; //ASM

    // ---------------------------------------------------------
    // In general case, Reff is calculated
    // in the following way (Zero-th node excluded):
    // (i = 0):        Reff[i] = Rb[0]*ais[0] + Rb[1]*ais[1]
    // (0 < i < m1):   Reff[i] = Rb[i]*ais[i] + Rb[i+1]*ais[i+1]
    // (i = m1)        Reff[m1] = Rb[m1]*ais[m1]
    //   Where Rb is lateral influx on a
    //   stream reach of unit length
    // So, Reff is sort of continuous function
    // approximated by piece-wise linear polynoms
    // 'reis' does not have coefficient 1/2 in the equations
    // so, the actual influx 1/2*Reff = reis
    // ---------------------------------------------------------

    // --- Artificial variation in time ---

    for (i = 0; i < n; i++)
        reis[i] = 0.0;

    i = 0;
    cmove = cHead;  // Points to the current stream head
    cend = cOutlet; // Point to the current outlet

    while (cmove != cend) {
        reis[i] = RetrieveQeff(cmove);
        cmove = cmove->getDownstrmNbr();
        i++;
    }

    // Special care has to be taken regarding the outlet node
    // We retrieve Qeff only if this is the _Basin_ Outlet
    // Otherwise, its Qeff will be accounted when it is
    // a Stream Head for the down slope stream reach
    if (cmove == OutletNode)
        reis[i] = RetrieveQeff(cmove);

    // --------------------------------------------------
    // Now, we consider that the influx in a node
    // is actually a lateral influx to a stream
    // link lying downstream of the current node
    //             Reff(i) = Qeff(i)
    // --------------------------------------------------
    //        Qeff(0)    Qeff(1)   Qeff(2)
    //          |          |         |
    //          |          |         |
    //          |          |         |
    //          |          |         |
    //          v          v         v
    //          o----------o---------o---------o
    //

    for (i = 0; i < m; i++) {
        if (i < m1) { //ASM added {
            reis[i] = (reis[i] + reis[i + 1]) / 2.;

            //ASM Percolation:
            //ASM 2/8/2017: (Next many lines)
            if (percolationOption == 3) {        //ASM Green Ampt Method
                /*if ( (reis[i]/cmove->getVArea() ) < ChannelConduc) {
					clis[i] = reis[i];
					reis[i] = 0.0;
				}
				else {*/
                if (cmove->getFt() > 0) {    // Check if there is already an infiltration front
                    Ft = cmove->getFt();
                    Ft_previous = Ft;
                    Ft = Ft + ChannelConduc * dt;
                } else {                // if no front, estimate Ft_init
                    Ft_init = ChannelConduc * dt;
                    Ft = Ft_init;
                    Ft_previous = 0.0;
                }
                Ft_prime = 0.0;
                PsiF = ((2 * PoreInd + 3) / (2 * PoreInd + 6) * abs(PsiB)) / 1000;
                //DeltTh = channelPorosity - cmove->getRootMoisture();
                DeltTh = channelPorosity;
                test = 0;
                while (test < 100) { //Start with Ft_prime = 0
                    Ft_prime = Ft_previous + ChannelConduc * dt + PsiF * DeltTh * log(1 + Ft / (PsiF * DeltTh));
                    if (Ft_prime <= Ft * 1.1 && Ft_prime >= Ft * 0.9)
                        test = 111;
                    test += 1;
                    Ft = Ft_prime;
                }
                rate = 0.0;
                rate = ChannelConduc * ((PsiF * DeltTh + Ft) / Ft);
                NodeLoss[i] = rate * channelPorosity * cmove->getVArea();
                if (reis[i] > 0) {
                    reis1 = reis[i] - NodeLoss[i];
                    if (reis1 < 0) {
                        clis[i] = reis[i];
                        reis[i] = 0.0;
                    } else {
                        clis[i] = NodeLoss[i];
                        reis[i] = reis1;
                    }
                } else
                    clis[i] = 0.0;
                IntStormVar = cmove->getIntStormVar();
                /*if (IntStormVar < 5) //ASM originally was < IntStormMax
						cmove->setFt(Ft_prime);
					else */
                cmove->setFt(0.0);
                //}
                Preach += clis[i];
            } else if (percolationOption == 1 || percolationOption == 2) {    //ASM constant loss and transient methods
                reis1 = reis[i] - NodeLoss[i];
                if (reis1 < 0) {
                    clis[i] = reis[i];
                    reis[i] = 0.0;
                } else {
                    clis[i] = NodeLoss[i];
                    reis[i] = reis1;
                }
                Preach += clis[i]; // ASM 2/16/2017 this sums percolation in the channel
            }


        } //ASM
        else if (i == m1) {
            if (cend == OutletNode)
                reis[i] += reis[m];
            reis[i] /= 2.;

            //ASM 2/8/2017: (Next 10 lines)
            if (percolationOption == 3) {        //ASM Green Ampt Method
                /*if ( (reis[i]/cmove->getVArea() ) < ChannelConduc) {
					clis[i] = reis[i];
					reis[i] = 0.0;
				}
				else {*/
                if (cmove->getFt() > 0) {    // Check if there is already an infiltration front
                    Ft = cmove->getFt();
                    Ft_previous = Ft;
                    Ft = Ft + ChannelConduc * dt;
                } else {                // if no front, estimate Ft_init
                    Ft_init = ChannelConduc * dt;
                    Ft = Ft_init;
                    Ft_previous = 0.0;
                }
                Ft_prime = 0.0;
                PsiF = ((2 * PoreInd + 3) / (2 * PoreInd + 6) * abs(PsiB)) / 1000;
                //DeltTh = channelPorosity - cmove->getRootMoisture();
                DeltTh = channelPorosity;
                test = 0;
                while (test < 100) { //Start with Ft_prime = 0
                    Ft_prime = Ft_previous + ChannelConduc * dt + PsiF * DeltTh * log(1 + Ft / (PsiF * DeltTh));
                    if (Ft_prime <= Ft * 1.1 && Ft_prime >= Ft * 0.9)
                        test = 111;
                    test += 1;
                    Ft = Ft_prime;
                }
                rate = 0.0;
                rate = ChannelConduc * ((PsiF * DeltTh + Ft) / Ft);
                NodeLoss[i] = rate * channelPorosity * cmove->getVArea();
                if (reis[i] > 0) {
                    reis1 = reis[i] - NodeLoss[i];
                    if (reis1 < 0) {
                        clis[i] = reis[i];
                        reis[i] = 0.0;
                    } else {
                        clis[i] = NodeLoss[i];
                        reis[i] = reis1;
                    }
                } else
                    clis[i] = 0.0;
                IntStormVar = cmove->getIntStormVar();
                /*if (IntStormVar < 5) //ASM originally < IntStormMax
						cmove->setFt(Ft_prime);
					else */
                cmove->setFt(0.0);
                //}
                Preach += clis[i];
            } else if (percolationOption == 1 || percolationOption == 2) {
                if (reis[i] > 0) {
                    reis1 = reis[i] - NodeLoss[i];
                    if (reis1 < 0) {
                        clis[i] = reis[i];
                        reis[i] = 0.0;
                    } else {
                        clis[i] = NodeLoss[i];
                        reis[i] = reis1;
                    }
                } else
                    clis[i] = 0.0;
                Preach += clis[i]; // ASM 2/16/2017 this sums percolation in the channel reach
            } //end ASM
        }

    }


    return;
}

/*****************************************************************************
**
**  tKinemat::RetrieveQeff()
**
**  Based on the influx from hillsope, return Qeff for the current node cmove
**
*****************************************************************************/
double tKinemat::RetrieveQeff(tCNode *cmove) {
    int intmp;
    double dtemp, value;

    tList<int> *TimeInd;  // <-- Ptr to Time steps stack
    tList<double> *Qeff;     // <-- Ptr to Qeff stack
    tListIter<int> IndIter;
    tListIter<double> QIter;

    intmp = -999;
    dtemp = -999;

    TimeInd = cmove->getTimeIndList();

    // If the lists are NOT empty
    if (!(TimeInd->isEmpty())) {
        IndIter.Reset(*TimeInd);
        IndIter.First();

        // Average influx rate per time interval > dt
        if (IndIter.DatRef() ==
            timer->getStepForSpecifiedDT(timer->getCurrentTime(), dtReff)) {

            Qeff = cmove->getQeffList();
            QIter.Reset(*Qeff);
            QIter.First();
            value = QIter.DatRef();

            // If current time is the end of Reff interval
            if ((int) ceil(timer->getCurrentTime() / dtReff) ==
                (int) floor(timer->getCurrentTime() / dtReff)) {

                if (TimeInd->getSize() == 1) {
                    TimeInd->Flush();
                    Qeff->Flush();
                }
                    // - Or, delete the first from the stack -
                else {
                    TimeInd->removeFromFront(intmp);
                    Qeff->removeFromFront(dtemp);
                }
            }
        } else
            value = 0.;
    } else
        value = 0.;

    return value;
}

/*****************************************************************************
**
**  tKinemat::AssignQin()
**
**  Based on the influx from above stream reaches, define Qit and H0
**
*****************************************************************************/
void tKinemat::AssignQin() {
#ifdef PARALLEL_TRIBS
                                                                                                                            // If upstream reaches are on other processors, receive
  if (tGraph::hasUpstream(id) )
    tGraph::receiveUpstream(id, cHead);
#endif

    Qin = cHead->getQstrm(); // Get inflow at the upper BND node
    qit = 0.5 * Qin;           // This value is used for '0's node

    H0 = pow(Qin * rifis[0] / (bis[0] * sqrt(siis[0])), 3. / 5.); // H at the BND
    return;
}

/*****************************************************************************\
**
**  tKinemat::ComputeQout()
**
**  Computes discharge in n-th node using the Chezy formulation
**
\*****************************************************************************/
void tKinemat::ComputeQout() {
    double fk2 = 5. / 3.;
    Qout = bis[m] * pow(his[m], fk2) * sqrt(siis[m]) / rifis[m]; // Q in n-th node
    return;
}

/*****************************************************************************
**
**  tKinemat::ComputeNodeQstr()
**
**  Computes discharge in i-th node using the Chezy & Manning formulation
**
*****************************************************************************/
double tKinemat::ComputeNodeQstrm(int i) {
    return (bis[i] * pow(his[i], (5. / 3.)) * sqrt(siis[i]) / rifis[i]); // Q in i-th node
}

/*****************************************************************************
**
**  tKinemat::ComputeNodeFlowVel()
**
**  Computes flow velocity in i-th node using the Chezy & Manning formulation
**
*****************************************************************************/
double tKinemat::ComputeNodeFlowVel(int i) {
    return (pow(his[i], (2. / 3.)) * sqrt(siis[i]) / rifis[i]); //velocity in i-th node
}

/*****************************************************************************
**
**  tKinemat::UpdateStreamVars()
**
**  Updates discharge and water level in i-th node
**
*****************************************************************************/
void tKinemat::UpdateStreamVars() {
    int i;
    tCNode *cmove, *cend;

    i = 0;
    IndividualPerc = 0.0; //ASM
    cmove = cHead;  // Points to the current stream head
    cend = cOutlet; // Point to the current outlet

    while (cmove != cend) {
        cmove->setHlevel(his[i]);
        if ((optres == 1) && (checkNode == checkID) && (checkID == cmove->getID())) { //JECR 2015
            cmove->setHlevel(LevelPool.getResElevOut());
        }
        cmove->setQstrm(ComputeNodeQstrm(i));
        if ((optres == 1) && (checkNode == checkID) && (checkID == cmove->getID())) { //JECR 2015
            cmove->setQstrm(Qout);
        }
        cmove->setFlowVelocity(ComputeNodeFlowVel(i));
        if (percolationOption != 0) { // ASM ** Akram: the preceding codes are slightly different for ASM
            cmove->setChannelPerc(clis[i]); //ASM 2/10/2017
        }
        if (clis[i] > 0.0) {
            //cmove->percOccur=cmove->percOccur+floor(clis[i]*1.0E+3)+1.0E-6;
            cmove->percOccur = cmove->getPercOccur() + 1;
            cmove->avPerc = cmove->getavPerc() + clis[i];
        }
        IndividualPerc = cmove->getChannelPerc();
        ParallelPerc += IndividualPerc;
        cmove = cmove->getDownstrmNbr();
        i++;
    }

    // Special care has to be taken regarding the outlet nodes
    // We do NOT have to assign the water level for the LAST node
    // of a stream reach because it may be erroneously used as an
    // initial condition for a tributary which would also contain
    // the node - the confluence node. Instead, the level will be
    // stored in a separate array 'OutletHlev' from which it will
    // then be read. But we need to surplus the discharge Q.
    // NOTE: we must assign the water level for the basin outlet only.

    if (cOutlet == OutletNode)
        cmove->setHlevel(his[i]);
    OutletHlev[id] = his[i];

    double cOutletQstrm = ComputeNodeQstrm(i);
    cmove->addQstrm(cOutletQstrm); // 'add' not 'set'!

#ifdef PARALLEL_TRIBS
                                                                                                                            // If downstream reaches are on other processors, send
  if ( tGraph::hasDownstream(id) )
    tGraph::sendDownstream(id, cmove, cOutletQstrm);
#endif

    return;
}

/*****************************************************************************
**
**  tKinemat::ControlPrint()
**
**  Prints out an array 'a' to a destination 'Otp'
**
*****************************************************************************/
void tKinemat::ControlPrint(ofstream &Otp, double *a, int NN) {
    for (int i = 0; i < NN; i++)
        Otp << a[i] << " ";
    Otp << "\n\n";
    return;
}

/*****************************************************************************
**
**  tKinemat::PrintFlowStacks()
**
**  Prints to Otp data contained in flow stacks of node 'cmove'
**  'cmove' can only be "Stream"!
**
**
*****************************************************************************/
void tKinemat::PrintFlowStacks(ofstream &Otp, tCNode *cmove) {
    tList<int> *TimeInd;  // Ptr to Time steps stack
    tList<double> *Qeff;     // Ptr to Qeff stack
    tListIter<int> IndIter;
    tListIter<double> QIter;

    if (cmove->getBoundaryFlag() == kStream) {
        TimeInd = cmove->getTimeIndList();

        if (!(TimeInd->isEmpty())) {
            IndIter.Reset(*TimeInd);
            Qeff = cmove->getQeffList();
            QIter.Reset(*Qeff);

            Otp << "- NODE: " << cmove->getID() << " -" << endl;
            for (IndIter.FirstP(); !(IndIter.AtEnd()); IndIter.NextP())
                Otp << IndIter.DatRef() << " ";
            Otp << endl;

            for (QIter.First(); !(QIter.AtEnd()); QIter.Next())
                Otp << QIter.DatRef() << " ";
            Otp << endl << endl;
        } else
            Otp << "- NODE: " << cmove->getID() << " -\t--- ZERO STACKS ---" << endl;
    }
    return;
}

/*****************************************************************************
**
**  tKinemat::SolveForTwoNodeReach()
**
**  The function finds water stage for the outlet node when a stream reach
**  consits of just a single link containing two nodes: head and outlet
**  The unknown variable is found by finding a root of the polynomial which
**  is the first eqaution in the system
**
*****************************************************************************/
void tKinemat::SolveForTwoNodeReach(double *Cd, double *Y1d,
                                    double *Y2d, double *Y3d,
                                    double *Reff, double *HLev,
                                    double Qit, double H0d) {
    double XX;
    double c1, c2, c3;

    // Pay attention, values start from 0...

    // XX[0] is the specified by the boundary condition
    c1 = 0.5 * Cd[1];
    c2 = 0.5 * Y2d[0];  // we get 1/6*l_0*b_1/dt
    c3 = 2 * Y3d[0] * (H0d - HLev[0]) - c2 * HLev[1] - Qit - Reff[0];

    // If c3 is somehow turns out to be positive, just keep the
    // level from the preceding time step. c3 _must_ be negative
    // for the equation to have positive roots
    if (c3 > 0.0)
        XX = HLev[1];
    else
        XX = rtsafe(polynomialH, c1, c2, c3, 0., 40., 1.0E-6, H0d);

    HLev[0] = H0d;  // Comes from the BND condition for time (j+1)
    HLev[1] = XX;  // Computed levels for (j+1) time step

    return;
}

/*****************************************************************************
**
**  tKinemat::KinematWave()
**
**  Implements kinematic wave routing model for a stream reach with fixed
**  upper boundary condition Hupp for a stream reach initially having water
**  level HLev(his) and lateral influx in the nodes Reff. The geometry of the
**  channel is "contained" in arrays C, Y1, Y2, Y3 which are the coefficient
**  vectors computed as function of channel widths, lengths, slopes, roughness.
**
**  The function updates the level vector HLev and returns the status of the
**  update in 'check'
**
**  It is assumed that the following variables are "known" to the function:
**  - n, m, m1
**
*****************************************************************************/

#define MAXITS 200
#define TOLF   1.0e-4
#define TOLMIN 1.0e-6
#define TOLX   1.0e-7
#define STPMX  100.0

void tKinemat::KinematWave(double *c, double *y1,
                           double *y2, double *y3,
                           double *Reff, double *HLev,
                           double Qit, double Hupp, int *check) {
    int ITER;
    int i;
    double *X, *X1, *XR;
    double *F, *aa, *bb, *cc, *gradf;
    double den, f, fold, stpmax, sum, temp, test;

    X = new double[n];
    assert(X != 0);
    X1 = new double[n];
    assert(X1 != 0);
    XR = new double[n];
    assert(XR != 0);
    F = new double[n - 1];
    assert(F != 0);

    gradf = new double[n];
    assert(gradf != 0);

    // These will contain sparse Jacobian matrix though the actual size
    // of vectors for the problem with known upper BND condition is n-1

    aa = new double[n];
    assert(aa != 0);
    bb = new double[n];
    assert(bb != 0);
    cc = new double[n];
    assert(cc != 0);


    // Pay attention, values start from 1! Not from 0!
    for (i = 0; i < m; i++)
        X[i] = X1[i] = HLev[i + 1]; //Iterations start using levels for time (t-1)

    // The following two must be used together
    ComputeFunction(F, c, y1, y2, y3, Reff, X1, HLev, Qit, Hupp, n);

    //ControlPrint(ControlOut, F, m);
    for (sum = 0.0, i = 0; i < m; i++)
        sum += F[i];
    sum /= m;

    // Compute minimization function 'f'
    f = fmin(F, m);

    // Test the initial guess for being root
    // ...more stringent test than simply TOLF
    test = 0.0;
    for (i = 0; i < m; i++) {
        if (fabs(F[i]) > test)
            test = fabs(F[i]);
    }
    if (test < 0.01 * TOLF) {
        *check = 0;
        FREERETURN
    }

    for (sum = 0.0, i = 0; i < m; i++)
        sum += SQR(X[i]);
    stpmax = STPMX * FMAX(sqrt(sum), (double) m); // MAX step length



    for (ITER = 1; ITER <= MAXITS; ITER++) {

        // Compute tridiagonal Jacobian  matrix
        ComputeJacobian(aa, bb, cc, y1, y2, y3, c, X1, m);

        // Compute grad(f) = F*J for line search
        gradf[0] = bb[0] * F[0] + aa[1] * F[1];
        for (i = 1; i < m1; i++)
            gradf[i] = cc[i - 1] * F[i - 1] + bb[i] * F[i] + aa[i + 1] * F[i + 1];
        gradf[m1] = bb[m1] * F[m1] + cc[m1 - 1] * F[m1 - 1];

        // Update X - calculated levels of the
        // preceding iteration and also store
        for (i = 0; i < m; i++)
            X[i] = X1[i];

        fold = f; // <- ... f

        // Re-assign the right-hand side for linear eq_s: [i+1]!
        for (i = 0; i < m; i++) {
            if (F[i] == 0)
                XR[i + 1] = F[i];
            else
                XR[i + 1] = -F[i];
        }

        // Solve linear system of equations using
        // the LU decomposition for sparse matrix
        SolveLinearSystem(aa, bb, cc, XR, m);


        // LNSRCH returns new X and f; also calculates F at the new X ---
        // Take care of X1!
        lnsrch(m, X, fold, gradf, XR, X1, &f, stpmax, check,
               F, c, y1, y2, y3, Reff, HLev, Qit, Hupp);

        test = 0.0;
        for (i = 0; i < m; i++) {
            if (fabs(F[i]) > test)
                test = fabs(F[i]);  // Pick the maximum value
        }
        if (test < TOLF) { // TEST for convergence on function values
            *check = 0;
            UpdateHsShifted(X1, HLev, Hupp, m); // Update HLev for t=(j+1)

            FREERETURN
        }

        if (*check) {      // TEST for grad(f) = zero
            test = 0.0;
            den = FMAX(f, 0.5 * m);
            for (i = 0; i < m; i++) {
                temp = fabs(gradf[i]) * FMAX(fabs(X1[i]), 1.0) / den;
                if (temp > test)
                    test = temp;
            }
            *check = (test < TOLMIN ? 1 : 0);
            UpdateHsShifted(X1, HLev, Hupp, m); // <--- Update HLev for t=(j+1)

            FREERETURN
        }

        test = 0.0;
        for (i = 0; i < m; i++) {   // TEST for convergnece on dx
            temp = (fabs(X1[i] - X[i])) / FMAX(fabs(X1[i]), 1.0);
            if (temp > test)
                test = temp;
        }
        if (test < TOLX) {
            UpdateHsShifted(X1, HLev, Hupp, m); //Update HLev for t=(j+1)

            FREERETURN
        }
    }
    cerr << "MAXITS exceeded in newt" << endl << flush;
    return;
}

#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN

/*****************************************************************************
**
**  tKinemat::ComputeFunction()
**
**  Estimates vector of function values 'F' which are deviations from zero
**  in the composed non-linear system of equation
**
**  It is assumed that the following variables are "known" to the function:
**  - n, m, m1
**
*****************************************************************************/
void tKinemat::ComputeFunction(double *F, double *c, double *y1,
                               double *y2, double *y3,
                               double *Reff, double *X, double *HLev,
                               double Qit, double Hupp, int NN) {
    int i, N, M, M1;
    double *XFK2;
    double fk2 = 5. / 3.;

    N = NN;
    M = N - 1;
    M1 = N - 2;

    XFK2 = new double[N];
    assert(XFK2 != 0);

    for (i = 0; i < M; i++)
        XFK2[i] = pow(X[i], fk2); // levels in the power

    // As long we use different indexation for levels (X, X1, XFK2),
    // i.e. which starts from the 1st node (not zero-th) then the scheme
    // is slightly different from given in the book. Indices for X-s are
    // (j-1) instead of (j) as they are for HLev_s.
    // Also, the equations were multiplied by a factor of 1/2.
    // Qit is the discharge specified by the BND condition for time (t+1)
    // by the water level Hupp in the node '0':
    //          Qit = 1/2*alfa_0*b_0*(Hupp(t+1))^m
    // '- HLev' is a vector of levels for time (t) at
    //          previous time step for all 'n' nodes
    // - 'X'    is a vector of water levels being computed for
    //          time (t+1) <--- the ones we are looking for

    for (i = 0; i < M; i++) { // HLev: levels for time (t), X: time (t+1)
        if (i == 0) {
            F[0] = 0.5 * c[2] * XFK2[1] - Qit + y3[0] * (Hupp - HLev[0]) +
                   y2[0] * (X[0] - HLev[1]) +
                   y1[0] * (X[1] - HLev[2]) - Reff[0];
        } else if (i == M1) // Changed below....
            F[m1] = 0.5 * c[M] * XFK2[M1] - 0.5 * c[M1] * XFK2[M1 - 1] +
                    y3[M1] * (X[M1 - 1] - HLev[M1]) +
                    y2[M1] * (X[M1] - HLev[M]) - Reff[M1];
        else {
            F[i] = 0.5 * c[i + 2] * XFK2[i + 1] - 0.5 * c[i] * XFK2[i - 1] +
                   y3[i] * (X[i - 1] - HLev[i]) + y2[i] * (X[i] - HLev[i + 1]) +
                   y1[i] * (X[i + 1] - HLev[i + 2]) - Reff[i];
        }
    }
    delete[] XFK2;
    return;
}

/*****************************************************************************
**
**  tKinemat::fmin()
**
**  Evaluates sum of squares of deviations 'fvec' from zero
**
*****************************************************************************/
double tKinemat::fmin(double *fvec, int N) {
    int i;
    double sum;
    for (sum = 0.0, i = 0; i < N; i++)
        sum += SQR(fvec[i]);
    return 0.5 * sum;
}

/*****************************************************************************
**
**  tKinemat::ComputeJacobian
**
**  Computes Jacobian matrix for the initial non-linear system of equations
**  of the form:      dF0/dh0, dF0/dh1, ... dF0/dhN-1
**                    dF1/dh0, dF1/dh1, ... dF1/dhN-1
**                    ...............................
**                  dFN-1/dh0 dFN-1/dh1 ... dFN-1/dhN-1
**  which ends up to be a tridiagonal matrix. It is therefore much easier
**  to work with it as with a sparse matrix: vectors aa, bb, cc are used
**  to represent its three diagonals. The relevant description follows below.
**
**  It is assumed that the following variables are "known" to the function:
**
**  ------------- Jacobian matrix computation ---------------
**   'N' here is not the same as 'n' in the calling function
**   It is equal to (n-1) and defines the actual size of the
**   tridiagonal Jacobian matrix (while as we have 'n' nodes)
**   Vectors aa, bb, cc must have allocated memory of size N
**   The actual use of vector variables is:
**         - aa: from 1 to N-1
**         - bb: from 0 to N-1
**         - cc: from 0 to N-2
**   So, the Jacobian is a sparse tridiagonal matrix formed
**   by vectors aa (lower on the diagonal), bb (middle), and cc (the upper one)
**
**      | b0  c0   0      ...                       |
**      | a1  b1  c1      ...                       |
**      | 0   a2  b2  c2  ...                       |
**      |                 ...                       |
**      |                 ...   0  aN-2  bN-2  cN-2 |
**      |                 ...   0   0    aN-1  bN-1 |
**
**   Function call: ComputeJacobian(aa, bb, cc, Y1, Y2, Y3, C, X, n-1);
**
*****************************************************************************/
void tKinemat::ComputeJacobian(double *aa, double *bb, double *cc,
                               double *y1, double *y2, double *y3,
                               double *c, double *X, int N) {
    int i;
    double f1, f2;
    double fk1 = 5. / 6.;
    double fk3 = 2. / 3.;

    f1 = pow(X[0], fk3);
    for (i = 0; i < N - 1; i++) {
        f2 = pow(X[i + 1], fk3); // Here is where the index is DIFFERENT
        bb[i] = y2[i];
        cc[i] = fk1 * c[i + 2] * f2 + y1[i];
        aa[i + 1] = -fk1 * c[i + 1] * f1 + y3[i + 1]; // Y3[i+1] is CORRECTION!
        f1 = f2;
    }
    // Here is where the index is DIFFERENT -v
    bb[N - 1] = fk1 * c[N] * pow(X[N - 1], fk3) + y2[N - 1];
    return;
}

/*****************************************************************************
**
**  tKinemat::SolveLinearSystem
**
**  Using the two functions from Numerical Recipes 'bandec' and 'banbks'
**  the function solves linear system of equation A*x = XR. The matrix
**  A is sparse and is given by three vectors (diagonals) aa, bb, cc
**  They are used to write the compact form of A - matrix a which is then
**  used for LU decomposition.
**  The solution vector overwrites XR[1, N] -> XR[0,N-1]
**
*****************************************************************************/
void tKinemat::SolveLinearSystem(double *aa, double *bb,
                                 double *cc, double *XR, int N) {
    unsigned long *indx, i;
    double **a, **al, d, *tt;

    tt = new double[N];
    assert(tt != 0);
    indx = new unsigned long[N + 1];
    assert(indx != 0);
    al = new double *[N + 1];
    assert(al != 0);
    a = new double *[N + 1];
    assert(a != 0);
    for (i = 0; i < N + 1; i++) {
        a[i] = new double[4];
        assert(a[i] != 0);
        al[i] = new double[2];
        assert(al[i] != 0);
    }

    // Filling the matrix with aa, bb, cc, Zero-th elements are not used!
    // (I did this only for convenience) Indices start from '1' not form '0'!

    for (i = 0; i < 4; i++)
        a[0][i] = 999.;  // wasted space in the compact format
    a[1][0] = 999.;    // wasted space in the compact format
    a[1][1] = 999.;    // wasted space in the compact format
    a[1][2] = bb[0];
    a[1][3] = cc[0];

    for (i = 2; i < N; i++) {
        a[i][0] = 999;
        a[i][1] = aa[i - 1];
        a[i][2] = bb[i - 1];
        a[i][3] = cc[i - 1];
    }
    a[N][0] = 999;
    a[N][1] = aa[N - 1];
    a[N][2] = bb[N - 1];
    a[N][3] = 999.; // wasted space in the compact format

    bandec(a, N, 1, 1, al, indx, &d); // Solve linear equations
    //by LU decomposition for banded matrix

    banbks(a, N, 1, 1, al, indx, XR);
    for (i = 0; i < N; i++)
        XR[i] = XR[i + 1];  // <- XR[i]

    for (i = 0; i < N + 1; i++) {
        delete[] a[i];
        delete[] al[i];
    }
    delete[] a;
    delete[] al;
    delete[] indx;
    delete[] tt;
    return;
}

/*****************************************************************************
**
**  tKinemat::bandec()
**
**
*****************************************************************************/

#define SWAP(a, b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0E-20

void tKinemat::bandec(double **a, unsigned long N, int M1, int M2, double **al,
                      unsigned long indx[], double *d) {
    unsigned long i, j, k, l;
    int mm;
    double dum;

    mm = M1 + M2 + 1;
    l = M1;

    for (i = 1; i <= M1; i++) {
        for (j = M1 + 2 - i; j <= mm; j++) a[i][j - l] = a[i][j];
        l--;
        for (j = mm - l; j <= mm; j++) a[i][j] = 0.0;
    }
    *d = 1.0;
    l = M1;

    for (k = 1; k <= N; k++) {
        dum = a[k][1];
        i = k;
        if (l < N) l++;
        for (j = k + 1; j <= l; j++) {
            if (fabs(a[j][1]) > fabs(dum)) {
                dum = a[j][1];
                i = j;
            }
        }
        indx[k] = i;
        if (dum == 0.0) a[k][1] = TINY;
        if (i != k) {
            *d = -(*d);
            for (j = 1; j <= mm; j++) SWAP(a[k][j], a[i][j])
        }
        for (i = k + 1; i <= l; i++) {
            dum = a[i][1] / a[k][1];
            al[k][i - k] = dum;
            for (j = 2; j <= mm; j++) a[i][j - 1] = a[i][j] - dum * a[k][j];
            a[i][mm] = 0.0;
        }
    }
    return;
}

/*****************************************************************************
**
**  tKinemat::banbnks()
**
**
*****************************************************************************/
void tKinemat::banbks(double **a, unsigned long N, int M1, int M2, double **al,
                      unsigned long indx[], double b[]) {
    unsigned long i, k, l;
    int mm;
    double dum;

    mm = M1 + M2 + 1;
    l = M1;
    for (k = 1; k <= N; k++) {
        i = indx[k];
        if (i != k) SWAP(b[k], b[i])
        if (l < N) l++;
        for (i = k + 1; i <= l; i++) b[i] -= al[k][i - k] * b[k];
    }
    l = 1;
    for (i = N; i >= 1; i--) {
        dum = b[i];
        for (k = 2; k <= l; k++) dum -= a[i][k] * b[k + i - 1];
        b[i] = dum / a[i][1];
        if (l < mm) l++;
    }
}

#undef SWAP
#undef TINY

/*****************************************************************************
**
**  tKinemat::tridag()
**
**  Solves for a vector u[0..N-1] the tridiagonal linear set given by
**                               A*u = r
**  a[1,N-1], b[0,N-1], c[0,N-2], and r[0,N-1] are the input vectors
**  and are not modified.
**
**  Function call:      tridag(aa, bb, cc, XR, F);
**
**  We need to update levels after the function has been called:
**    for (i=0; i<m; i++)
**       X1[i] = X[i] + F[i];
**
*****************************************************************************/
void tKinemat::tridag(double a[], double b[], double c[],
                      double r[], double u[], unsigned long N) {
    long j;
    double bet;
    double *gam;

    gam = new double[N];
    assert(N != 0);

    // If this happens, we should rewrite equations as a set
    // of order (N-1), with u[1] trivially eliminated

    if (b[0] == 0.0) {
        cerr << "Error 1 in tridag: b(0) = 0, exiting..." << endl << flush;
        exit(2);
    }
    u[0] = r[0] / (bet = b[0]);
    for (j = 1; j < N; j++) {   // Decomposition and forward substitution
        gam[j] = c[j - 1] / bet;
        bet = b[j] - a[j] * gam[j];
        if (bet == 0.0) // Algorithm fails
            cerr << "Error 2 in tridag..., bet = 0.; j = " << j << endl << flush;
        u[j] = (r[j] - a[j] * u[j - 1]) / bet;
    }

    for (j = (N - 2); j >= 0; j--) { // <- Backsubstitution
        u[j] -= gam[j + 1] * u[j + 1];
    }
    delete[] gam;
}

/*****************************************************************************
**
**  tKinemat::GAUSS()
**
**  Direct Gauss method for a system of linear equations
**
*****************************************************************************/
void tKinemat::GAUSS(double *F, double **FG, int M) {
    int i, j;
    int M1 = M - 1;
    double D;

    for (i = 0; i < M1; i++) {  // Forward
        j = i + 1;
        D = FG[j][i] / FG[i][i];
        FG[j][i] = 0.;
        FG[j][j] = FG[j][j] - D * FG[i][j];
        F[j] = F[j] - D * F[i];
    }
    F[M1] = F[M1] / FG[M1][M1]; // Determination of Xm

    i = M1;  // Initial definition of K
    while (i > 0) {
        i--; // Backward
        j = i + 1;
        F[i] = (F[i] - FG[i][j] * F[j]) / FG[i][i];
    }
    return;
}

/*****************************************************************************\
**
**  tKinemat::lnsrch()
**
**  Performs Newton-Raphson procedure with backtracking if step is too
**  large
**
**  Function call: lnsrch(m, X, fold, gradf, XR, X1, &f, stpmax, check,
**                 F, C, Y1, Y2, Y3, Reff, X, HLev, Qit, H0, ComputeFunction);
**
\*****************************************************************************/
#define NRANSI
#define ALF 1.0e-4
#define TOLX 1.0e-7

void tKinemat::lnsrch(int N, double xold[], double fold,
                      double g[], double p[], double x[],
                      double *f, double stpmax, int *check,
                      double *F, double *c,
                      double *y1, double *y2, double *y3, double *Reff,
                      double *HLev, double Qit, double Hupp) {
    int i;
    double a, alam, alam2, alamin, b, disc, f2,
            rhs1, rhs2, slope, sum, temp, test, tmplam;

    f2 = alam2 = 0.0;
    *check = 0;
    for (sum = 0.0, i = 0; i < N; i++)
        sum += p[i] * p[i];
    sum = sqrt(sum);      // Square root of sum of squares of deltaX_s

    if (sum > stpmax) {   // Scale, if attempted step is too big
        for (i = 0; i < N; i++)
            p[i] *= stpmax / sum;
    }

    for (slope = 0.0, i = 0; i < N; i++) // <--- grad(f)*dx < 0!
        slope += g[i] * p[i];
    if (slope >= 0.0) {
        cout << "\nWarning: Roundoff problem in lnsrch, slope = " << slope
             << "   ... exiting to system" << endl << flush;
        exit(1);
    }

    test = 0.0;
    for (i = 0; i < N; i++) {
        temp = fabs(p[i]) / FMAX(fabs(xold[i]), 1.0);
        if (temp > test)
            test = temp;
    }

    alamin = TOLX / test; // MIN step length 10^-7/(...>1)
    alam = 1.0;         // Always try FULL Newton step first: lambda = 1


    // Infinite loop starts here ...
    for (;;) {          // Start of iteration loop...
        for (i = 0; i < N; i++) {
            x[i] = xold[i] + alam * p[i]; // <--- Updating X...
            if (x[i] < 0.0)    // Can't have hegative numbers
                x[i] = 1.0E-6;   // Bring it back close to '0'
        }

        ComputeFunction(F, c, y1, y2, y3, Reff, x, HLev, Qit, Hupp, N + 1);

        *f = fmin(F, N); // compute minimization function 'f'

        // The only two options to get out of the function
        if (alam < alamin) {  // Convergence on dx. For zero finding, the
            for (i = 0; i < N; i++)
                x[i] = xold[i];  // calling program should verify the convergence
            *check = 1;
            return;
        } else if (*f <= fold + ALF * alam * slope) //f(Xnew) <= f(Xold)+alfa*grad(f)*dx
            return;  // Sufficient function decrease --> return


        else { // Backtrack
            if (alam == 1.0)         // First time...
                tmplam = -slope / (2.0 * (*f - fold - slope));
            else {                   // Subsequent backtracks...
                rhs1 = *f - fold - alam * slope;
                rhs2 = f2 - fold - alam2 * slope;
                a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
                b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
                if (a == 0.0)
                    tmplam = -slope / (2.0 * b);
                else {
                    disc = b * b - 3.0 * a * slope;
                    if (disc < 0.0)
                        tmplam = 0.5 * alam;
                    else if (b <= 0.0)
                        tmplam = (-b + sqrt(disc)) / (3.0 * a);
                    else
                        tmplam = -slope / (b + sqrt(disc));
                }
                if (tmplam > 0.5 * alam)  //Lambda <= Lambda1
                    tmplam = 0.5 * alam;
            }
        }
        alam2 = alam;
        f2 = *f;
        alam = FMAX(tmplam, 0.1 * alam);
    }
}

#undef ALF
#undef TOLX
#undef NRANSI

/*****************************************************************************
**
**  tKinemat::UpdateHsShifted()
**
**  Updates water levels
**
*****************************************************************************/
void tKinemat::UpdateHsShifted(double *Xnew, double *Xold, double Hupp, int N) {
    Xold[0] = Hupp;          //Comes from the BND condition for time (j+1)
    for (int i = 0; i < N; i++)
        Xold[i + 1] = Xnew[i];   //Computed levels for (j+1) time step
    return;
}

/***************************************************************************
**
** tKinemat::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/

void tKinemat::writeRestart(fstream &rStr) const {
    BinaryWrite(rStr, id);
    BinaryWrite(rStr, m);
    BinaryWrite(rStr, m1);
    BinaryWrite(rStr, TimeSteps);
    BinaryWrite(rStr, dt);
    BinaryWrite(rStr, dtReff);
    BinaryWrite(rStr, qit);
    BinaryWrite(rStr, Qin);
    BinaryWrite(rStr, H0);
    BinaryWrite(rStr, Qout);
    BinaryWrite(rStr, maxH);
    BinaryWrite(rStr, maxReff);
    BinaryWrite(rStr, Roughness);
    BinaryWrite(rStr, Width);
    BinaryWrite(rStr, kincoef);

    // If this isn't dumped FlwVel, Qstrm, Hlevel are wrong
    int sz = NodesLstO.getSize();
    BinaryWrite(rStr, sz);
    for (int i = 0; i < sz; i++)
        BinaryWrite(rStr, OutletHlev[i]);

    OutletNode->writeRestart(rStr);

    tFlowNet::writeRestart(rStr);
}

/***************************************************************************
**
** tKinemat::readRestart() Function
**
***************************************************************************/

void tKinemat::readRestart(fstream &rStr) {
    BinaryRead(rStr, id);
    BinaryRead(rStr, m);
    BinaryRead(rStr, m1);
    BinaryRead(rStr, TimeSteps);
    BinaryRead(rStr, dt);
    BinaryRead(rStr, dtReff);
    BinaryRead(rStr, qit);
    BinaryRead(rStr, Qin);
    BinaryRead(rStr, H0);
    BinaryRead(rStr, Qout);
    BinaryRead(rStr, maxH);
    BinaryRead(rStr, maxReff);
    BinaryRead(rStr, Roughness);
    BinaryRead(rStr, Width);
    BinaryRead(rStr, kincoef);

    int sz;
    BinaryRead(rStr, sz);
    for (int i = 0; i < sz; i++)
        BinaryRead(rStr, OutletHlev[i]);

    OutletNode->readRestart(rStr);

    tFlowNet::readRestart(rStr);
}

#ifdef PARALLEL_TRIBS
                                                                                                                        /***************************************************************************
**
** tKinemat::openOutletFile() Function
**
***************************************************************************/

void tKinemat::openOutletFile(tInputFile &infile)
{
  // This function is used when running in parallel, since graph
  // processing is required to figure out what processor contains
  // the last reach
  if (tGraph::hasLastReach()) {
    char fullName2[kMaxNameSize+20];
    infile.ReadItem(fullName2, "OUTHYDROFILENAME" );
    strcat( fullName2, "_Outlet.qout" );
    theOFStream.open(fullName2);

    if ( !theOFStream.good() ) {
      cout<<"\nWarning: Output file not created.... "
           <<"\nExiting Program..."<<endl<<flush;
      exit(2);
    }

    if (simCtrl->Header_label=='Y')
      theOFStream<<"1-Time,hr\t "<<"2-Qstrm,m3/s\t"<<"3-Hlev,m"<<"\n";
  }

}
#endif

//=========================================================================
//
//
//                          End tKinemat.cpp
//
//
//=========================================================================
