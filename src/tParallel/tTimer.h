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
**  tTimer.h: Header for tTimer class and objects
**  This code is a modified version of the POOMA timers.
**
***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *	The Timer class allows for easy timing of the program.  The timer
 * tracks real (clock) time elapsed, user time, and system time.
 *
 ***************************************************************************/

#ifndef TTIMER_H
#define TTIMER_H

#ifdef __sgi
// make sure this is defined for BSD time routines
#define _BSD_TYPES
// fix a glitch in ANSI compatibility with SGI headers
#define _STAMP_T
#endif

#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <mpi.h>

#ifdef __sgi
// fix a glitch in ANSI compatibility with SGI headers
#undef _STAMP_T
#endif

class tTimer
{
public:
  tTimer();			// Constructor
  ~tTimer();                    // Destructor
  void clear();			// Set all accumulated times to 0
  void start();			// Start timer
  void stop();			// Stop timer

  double clock_time();		// Report clock time accumulated in seconds
  double user_time();		// Report user time accumlated in seconds
  double system_time();		// Report system time accumulated in seconds
  double cpu_time()
  {
    // Report total cpu_time which is just user_time + system_time
    return ( user_time() + system_time() );
  }		

  double calibration;		// Calibration time: time it takes to
                                // get in and out of timer functions
private:
  short timer_state;		// State of timer, either on or off
  long cpu_speed;		// CPU speed for times() call

  unsigned long last_secs;	// Clock seconds value when the
				// timer was last started
  long last_usecs;		// Clock useconds value when the
				// timer was last started
  unsigned long last_user_time;   // User time when timer was last started
  unsigned long last_system_time; // System time when timer was last started

  long current_secs;		// Current accumulated clock seconds
  long current_usecs;		// Current accumulated clock useconds
  long current_user_time;	// Current accumulated user time
  long current_system_time;	// Current accumulated system time

  struct tms tmsbuf;	        //  Values from call to times
  struct timeval tvbuf;	        //  Values from call to gettimeofday
  struct timezone tzbuf;        //  Timezone values from gettimeofday
	  		        //  These values aren't used for anything
};

//#include "tParallel/tTimer.cpp" //Not sure what's going on here? -WR

#endif // TTIMER_H

//=========================================================================
//
//
//                          End of tTimer.h 
//
//
//=========================================================================
