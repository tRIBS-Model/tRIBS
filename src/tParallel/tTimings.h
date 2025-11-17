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
**  tTimings.h: Header for tTimngs class and objects
**  This code is a modified version of the POOMA timers.
**
***************************************************************************/

/*************************************************************************
 * tTimings - a simple singleton class which lets the user create and
 *   timers that can be printed out at the end of the program.
 *
 * General usage
 *  1) create a timer:
 *     tTimings::TimerRef val = tTimings::getTimer("timer name");
 *  This will either create a new one, or return a ref to an existing one
 *
 *  2) start a timer:
 *     tTimings::startTimer(val);
 *  This will start the referenced timer running.  If it is already running,
 *  it will not change anything.
 *
 *  3) stop a timer:
 *     tTimings::stopTimer(val);
 *  This will stop the timer, assuming it was running, and add in the
 *  time to the accumulating time for that timer.
 *
 *  4) print out the results:
 *     tTimings::print();
 *
 *************************************************************************/

#ifndef TTIMINGS_H
#define TTIMINGS_H

// include files
#include "src/tParallel/tTimer.h"

// added by -WR to compile need to check
#include <iostream>
#include <mpi.h>
using namespace std;

//////////////////////////////////////////////////////////////////////
/*
  A simple compliant implementation of auto_ptr.
  This is from Greg Colvin's implementation posted to comp.std.c++.

  Instead of using mutable this casts away const in release.

  We have to do this because we can't build containers of these
  things otherwise.
  */
//////////////////////////////////////////////////////////////////////

template<class X>
class my_auto_ptr
{
  X* px;
public:
  my_auto_ptr() : px(0) {}
  my_auto_ptr(X* p) : px(p) {}
  my_auto_ptr(const my_auto_ptr<X>& r) : px(r.release()) {}
  my_auto_ptr& operator=(const my_auto_ptr<X>& r)
  {
    if (&r != this)
      {
	delete px;
	px = r.release();
      }
    return *this;
  }
  ~my_auto_ptr() { delete px; }
  X& operator*()  const { return *px; }
  X* operator->() const { return px; }
  X* get()        const { return px; }
  X* release()    const { X *p=px; ((my_auto_ptr<X>*)(this))->px=0; return p; }
};

#include <vector>
using std::vector;
#include <map>
using std::map;

// a simple class used to store timer values
class TimerInfo
{
public:
  // typedef for reference to a timer
  typedef int TimerRef;

  // constructor
  TimerInfo() : cpuTime(0.0), wallTime(0.0), indx(-1), name("") {
    clear();
  }

  // destructor
  ~TimerInfo() { }

  // timer operations
  void start() {
    if (!running) {
      running = true;
      t.stop();
      t.clear();
      t.start();
    }
  }

  void stop() {
    if (running) {
      t.stop();
      running = false;
      cpuTime += t.cpu_time();
      wallTime += t.clock_time();
    }
  }

  void clear() {
    t.stop();
    t.clear();
    running = false;
  }

  // the POOMA timer that this object manages
  tTimer t;

  // the name of this timer
  string name;

  // the accumulated time
  double cpuTime;
  double wallTime;

  // is the timer turned on right now?
  bool running;

  // an index value for this timer
  TimerRef indx;
};



class tTimings
{
public:
  // typedef for reference to a timer
  typedef int TimerRef;

  // a typedef for the timer information object
  typedef TimerInfo TimerInfo_t;

public:
  // Default constructor
  tTimings();

  // Destructor - clear out the existing timers
  ~tTimings();

  //
  // timer manipulation methods
  //

  // create a timer, or get one that already exists
  static TimerRef getTimer(const char *);

  // start a timer
  static void startTimer(TimerRef);

  // stop a timer, and accumulate it's values
  static void stopTimer(TimerRef);

  // clear a timer, by turning it off and throwing away its time
  static void clearTimer(TimerRef);

  // return a TimerInfo struct by asking for the name
  static TimerInfo_t *infoTimer(const char *nm) {
    return TimerMap[string(nm)];
  }

  //
  // I/O methods
  //

  // print the results to standard out
  static void print();

private:
  // type of storage for list of TimerInfo
  typedef vector<my_auto_ptr<TimerInfo> > TimerList_t;
  typedef map<string, TimerInfo *> TimerMap_t;

  // a list of timer info structs
  static TimerList_t TimerList;

  // a map of timers, keyed by string
  static TimerMap_t TimerMap;
};

//#include "tParallel/tTimings.cpp" // Not sure why it was done this way,so you only call tTimings.h, causes errors with CMAKE

#endif // TTIMINGS_H

//=========================================================================
//
//
//                          End of tTimings.h 
//
//
//=========================================================================
