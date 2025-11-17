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
**  tTimngs.cpp: Functions for class tTimngs (see tTimings.h)
**  This code is a modified version of the POOMA timers.
**
***************************************************************************/

// header file was missing not sure how this was previously compiled WR
#include "src/tParallel/tTimings.h"
#include <cstring>


// static data members of tTimings class
tTimings::TimerList_t tTimings::TimerList;
tTimings::TimerMap_t  tTimings::TimerMap;


//////////////////////////////////////////////////////////////////////
// default constructor
tTimings::tTimings() { }


//////////////////////////////////////////////////////////////////////
// destructor
tTimings::~tTimings() { }


//////////////////////////////////////////////////////////////////////
// create a timer, or get one that already exists
tTimings::TimerRef tTimings::getTimer(const char *nm) {
  string s(nm);
  TimerInfo *tptr = 0;
  TimerMap_t::iterator loc = TimerMap.find(s);
  if (loc == TimerMap.end()) {
    tptr = new TimerInfo;
    tptr->indx = TimerList.size();
    tptr->name = s;
    TimerMap.insert(TimerMap_t::value_type(s,tptr));
    TimerList.push_back(my_auto_ptr<TimerInfo>(tptr));
  } else {
    tptr = (*loc).second;
  }
  return tptr->indx;
}


//////////////////////////////////////////////////////////////////////
// start a timer
void tTimings::startTimer(TimerRef t) {
  if (t < 0 || t >= TimerList.size())
    return;
  TimerList[t]->start();
}


//////////////////////////////////////////////////////////////////////
// stop a timer, and accumulate it's values
void tTimings::stopTimer(TimerRef t) {
  if (t < 0 || t >= TimerList.size())
    return;
  TimerList[t]->stop();
}


//////////////////////////////////////////////////////////////////////
// clear a timer, by turning it off and throwing away its time
void tTimings::clearTimer(TimerRef t) {
  if (t < 0 || t >= TimerList.size())
    return;
  TimerList[t]->clear();
}


//////////////////////////////////////////////////////////////////////
// print out the timing results
void tTimings::print() {
  int i,j;
  if (TimerList.size() < 1)
    return;

  int nodes, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // report the average time for each timer
  if (rank == 0) {
    cout << "-----------------------------------------------------------------";
    cout << endl;
    cout << "     Timing results for " << nodes << " nodes:" << endl;
    cout << "-----------------------------------------------------------------";
    cout << endl;
  }
  for (i=0; i<1; ++i){
    TimerInfo *tptr = TimerList[i].get();
    double walltotal = 0.0, cputotal = 0.0;

    MPI_Reduce(&tptr->wallTime, &walltotal, 1, MPI_DOUBLE, MPI_MAX, 0,
        MPI_COMM_WORLD);
    MPI_Reduce(&tptr->cpuTime, &cputotal, 1, MPI_DOUBLE, MPI_MAX, 0,
        MPI_COMM_WORLD);

    if (rank == 0) {
      cout << tptr->name.c_str() << " ";
      for (j=strlen(tptr->name.c_str()); j < 20; ++j)
        cout << ".";
      cout << " Wall tot = ";
      cout.width(16);
      cout << walltotal << ", CPU tot = ";
      cout.width(16);
      cout << cputotal << endl << endl;
    }
  }

  for (i=1; i < TimerList.size(); ++i) {
    TimerInfo *tptr = TimerList[i].get();
    double wallmax = 0.0, cpumax = 0.0, wallmin = 0.0, cpumin = 0.0;
    double wallavg = 0.0, cpuavg = 0.0;

    MPI_Reduce(&tptr->wallTime, &wallmax, 1, MPI_DOUBLE, MPI_MAX, 0,
        MPI_COMM_WORLD);
    MPI_Reduce(&tptr->cpuTime, &cpumax, 1, MPI_DOUBLE, MPI_MAX, 0,
        MPI_COMM_WORLD);
    MPI_Reduce(&tptr->wallTime, &wallmin, 1, MPI_DOUBLE, MPI_MIN, 0,
        MPI_COMM_WORLD);
    MPI_Reduce(&tptr->cpuTime, &cpumin, 1, MPI_DOUBLE, MPI_MIN, 0,
        MPI_COMM_WORLD);
    MPI_Reduce(&tptr->wallTime, &wallavg, 1, MPI_DOUBLE, MPI_SUM, 0,
        MPI_COMM_WORLD);
    MPI_Reduce(&tptr->cpuTime, &cpuavg, 1, MPI_DOUBLE, MPI_SUM, 0,
        MPI_COMM_WORLD);

    if (rank == 0) {
      cout << tptr->name.c_str() << " ";
      for (j=strlen(tptr->name.c_str()); j < 20; ++j)
        cout << ".";
      cout << " Wall max = ";
      cout.width(16);
      cout << wallmax << ", CPU max = ";
      cout.width(16);
      cout << cpumax << endl;
      for (j = 0; j < 21; ++j)
        cout << " ";
      cout << " Wall avg = ";
      cout.width(16);
      cout << wallavg / nodes << ", CPU avg = ";
      cout.width(16);
      cout << cpuavg / nodes << endl;
      for (j = 0; j < 21; ++j)
        cout << " ";
      cout << " Wall min = ";
      cout.width(16);
      cout << wallmin << ", CPU min = ";
      cout.width(16);
      cout << cpumin << endl << endl;
    }
  }
  if (rank == 0) {
    cout << "-----------------------------------------------------------------";
    cout << endl;
  }
}

//=========================================================================
// 
// 
//                        End of tTimings.cpp
//    
//
//=========================================================================
