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
**  tTimer.cpp: Functions for class tTimer (see tTimer.h)
**  This code is a modified version of the POOMA timers.
**
***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *	tTimer.cpp constains the implementation of the class tTimer which
 * provides a generic way to time things.  Each Timer keeps a cumulative
 * time.  It can be started, stopped, or cleared.
 *	On the HPs, the clock time is obtained by calls to the routine
 * gettimeofday() and the user, system, and total cpu times are obtained
 * using the routine times()
 *
 * This class was originally developed by Mark Nelson in the Theoretical
 * Biophysics Group at the University of Illinois.
 ***************************************************************************/

#if defined(POOMA_T3E)
#define tick_secs(ticks, clock_ticks) \
( (double) ticks / (double) clock_ticks )          
#endif

//  TIMERON and TIMEROFF define the current state of a Timer
#define TIMEROFF	0
#define TIMERON		1

// Header file was missing WR not sure how?
#include "src/tParallel/tTimer.h"

/************************************************************************/
/*									*/
/*			FUNCTION tTimer					*/
/*									*/
/*	This is the constructor for the class tTimer.  It sets the timer */
/*  status to TIMEROFF and clears all the values.  It also makes a call */
/*  to sysconf to determine the number of clock ticks per second for    */
/*  use with the call times()						*/
/*  It also makes calibration calls.                                    */
/*									*/
/************************************************************************/

tTimer::tTimer()
{
#ifdef __MWERKS__
  // For now, stub out all Timer guts for MetroWerks
#else

#ifndef POOMA_TFLOP
  cpu_speed = sysconf(_SC_CLK_TCK);
#endif
  timer_state = TIMEROFF;
  clear();

  // Calibration:
#if defined(POOMA_T3E)
  long start_time, end_time, total_time;
  (void) rtclock();
  start_time = rtclock();
  end_time = rtclock();
  total_time = end_time - start_time;
  calibration = tick_secs(total_time, cpu_speed);
#else
  // No other machines have calibration defined yet.
  calibration = 0.0;
#endif

#endif // __MWERKS__
}
/*			END OF FUNCTION Timer				*/

tTimer::~tTimer()
{
#ifdef __MWERKS__
  // For now, stub out all Timer guts for MetroWerks
#else

  // Check to see if the timer is on
  if (timer_state == TIMERON)
    {
      //  Destroying a running Timer
      // ERRORMSG(level2 << "TRIED TO DELETE A RUNNING TIMER!\n");
      // ERRORMSG("STOPPING AND DELETING TIMER." << endl);
      timer_state=TIMEROFF;
    }

#endif // __MWERKS__
}

/************************************************************************/
/*									*/
/*			FUNCTION clear					*/
/*									*/
/*	clear sets all of the accumulated times for this timer to 0.	*/
/*  It is intended to only be used on a stopped timer.  If it is used	*/
/*  on a running timer, a warning message is printed, the timer is      */
/*  stopped and all of its values are cleared.				*/
/*									*/
/************************************************************************/

void tTimer::clear()
{
#ifdef __MWERKS__
  // For now, stub out all Timer guts for MetroWerks
  return;
#else
  // Check to see if the timer if on
  if (timer_state == TIMERON)
    {
      //  Clearing a running Timer
      // ERRORMSG(level2 << "TRIED TO CLEAR A RUNNING TIMER!\n");
      // ERRORMSG("SETTING ALL VALUES TO 0 AND STOPPING TIMER." << endl);
      timer_state = TIMEROFF;
    }

  //  Set all of the accumulated values to 0
#ifdef POOMA_TFLOP
  current_clock = 0.0;
#else
  current_secs = 0;
  current_usecs = 0;
  current_user_time = 0;
  current_system_time = 0;
#endif // POOMA_TFLOP

#endif // __MWERKS__
}
/*			END OF FUNCTION clear				*/

/************************************************************************/
/*									*/
/*			FUNCTION start					*/
/*									*/
/*	start a tTimer timing.  This will start adding time elapsed to  */
/*  the current accumulated values of the timer.  If you try to start   */
/*  a timer that is already running, a warning message is printed	*/
/*									*/
/************************************************************************/

void tTimer::start()
{
#ifdef __MWERKS__
  // For now, stub out all Timer guts for MetroWerks
  return;
#else
  //  Check to see if the timer is already running
  if (timer_state == TIMERON)
    {
      // ERRORMSG(level2 << "TRIED TO START A RUNNING TIMER!\n");
      // ERRORMSG("CONTINUING UNCHANGED." << endl);
      return;
    }

  //  Get the current time values from the system
#if defined(POOMA_T3E)
  last_secs = rtclock();
  // Omit non-real times on T3E:
  last_usecs = 0;
  last_user_time = 0;
  last_system_time = 0;
#elif defined(POOMA_TFLOP)
  last_clock = dclock();
#else
  gettimeofday(&tvbuf, &tzbuf);
  times(&tmsbuf);
  //  Set the starting values to the current time
  last_secs = tvbuf.tv_sec;
  last_usecs = tvbuf.tv_usec;
  last_user_time = tmsbuf.tms_utime;
  last_system_time = tmsbuf.tms_stime;
#endif

  //  Set the state of the Timer
  timer_state = TIMERON;
  return;
#endif // __MWERKS__
}
/*			END OF FUNCTION start				*/

/************************************************************************/
/*									*/
/*				FUNCITON stop				*/
/*									*/
/*	stop stops a tTimer from accumulating time.  If you try to stop */
/*  a stopped tTimer, a warning message is printed			*/
/*									*/
/************************************************************************/

void tTimer::stop()
{
#ifdef __MWERKS__
  // For now, stub out all Timer guts for MetroWerks
  return;
#else
  //  Check to see if the timer is already stopped
  if (timer_state == TIMEROFF)
    {
      // ERRORMSG(level2 << "TRIED TO STOP A STOPPED TIMER!\n");
      // ERRORMSG("CONTINUING UNCHANGED." << endl);
      return;
    }

  //  Get the current time values from the system and accumulate
#if defined(POOMA_T3E)
  long end_time = rtclock();

  current_secs +=  end_time - last_secs;
  current_usecs += 0;
  current_user_time += 0;
  current_system_time += 0;
#elif defined(POOMA_TFLOP)
  double end_clock = dclock();
  current_clock += end_clock - last_clock;
#else
  gettimeofday(&tvbuf, &tzbuf);
  times(&tmsbuf);

  current_secs += tvbuf.tv_sec - last_secs;
  current_usecs += tvbuf.tv_usec - last_usecs;
  current_user_time += tmsbuf.tms_utime - last_user_time;
  current_system_time += tmsbuf.tms_stime - last_system_time;
#endif

  //  Set the state of the Timer
  timer_state = TIMEROFF;
  return;
#endif // __MWERKS__
}
/*			END OF FUNCTION stop				*/

/************************************************************************/
/*									*/
/*			FUNCTION clock_time				*/
/*									*/
/*	clock_time returns the current amount of real (clock) time	*/
/*  accumulated by this timer.  If the timer is stopped, this is just	*/
/*  the total accumulated time.  If the timer is running, this is the	*/
/*  accumulated time plus the time since the timer was last started.	*/
/*									*/
/************************************************************************/

double tTimer::clock_time()
{
#ifdef __MWERKS__
  // For now, stub out all Timer guts for MetroWerks
  return 0.0;
#else

#if !defined(POOMA_TFLOP)
  long seconds;	    // seconds elapsed
  
#if !defined(POOMA_T3E)
  long useconds;    // useconds (mirco-seconds) elapsed
#endif

#endif

  double ret_val;    // time elpased

  if (timer_state == TIMEROFF)
    {
      // Timer is currently off, so just return accumulated time
#if !defined(POOMA_TFLOP)
      seconds = current_secs;
      
#if !defined(POOMA_T3E)
      useconds = current_usecs;
#endif

#else
      ret_val = current_clock;
#endif
    }
  else
    {
      // Timer is currently running, so add the elapsed
      // time since the timer was last started to the
      // accumulated time
#if defined(POOMA_T3E)
      long end_time = rtclock();
      seconds = current_secs + end_time - last_secs;
#elif defined(POOMA_TFLOP)
      double end_clock = dclock();
      ret_val = current_clock + end_clock - last_clock;
#else
      gettimeofday(&tvbuf, &tzbuf);

      seconds = current_secs + tvbuf.tv_sec - last_secs;
      useconds = current_usecs + tvbuf.tv_usec - last_usecs;
#endif
    }

  //  Convert into floating point number of seconds
#if defined(POOMA_T3E)
  ret_val = tick_secs(seconds, cpu_speed);
#elif defined(POOMA_TFLOP)
  // no need to convert
#else
  //  Adjust for the fact that the useconds may be negative.
  //  If they are, take away 1 second and add 1 million
  //  microseconds until they are positive
  while (useconds < 0)
    {
      useconds = useconds + 1000000;
      seconds = seconds - 1;
    }

  long long_ret_val = (1000000 * seconds) + useconds;
  ret_val = ( (double) long_ret_val ) / 1000000.0;
#endif

  return ret_val;
  
#endif // __MWERKS__
}
/*			END OF FUNCTION clock_time			*/

/************************************************************************/
/*									*/
/*			FUNCTION user_time				*/
/*									*/
/*	user_time reports the current amount of user cpu time           */
/*   accumulated by this tTimer.  If the timer is currently off, 	*/
/*   this is just the accumulated time.  If the tTimer is running, this  */
/*   is the accumulated time plus the time since the timer was last    */
/*   started.								*/
/*									*/
/************************************************************************/

double tTimer::user_time()
{
#ifdef __MWERKS__
// For now, stub out all Timer guts for MetroWerks
  return 0.0;
#else
  double ret_val;		//  Return value	

#if ( defined(POOMA_T3E) || defined(POOMA_TFLOP) )
  // Not defined yet on T3E or TFLOP.
  // ERRORMSG("user_time() not defined." << endl);
  ret_val = -9999.0;
#else
  if (timer_state == TIMEROFF)
    {
      //  Timer is off, just return accumulated time
      ret_val = current_user_time;
    }
  else
    {
      //  Timer is on, add current running time to accumulated time
      times(&tmsbuf);
      ret_val = current_user_time + tmsbuf.tms_utime - last_user_time;
    }

  //  Convert from clock ticks to seconds using the
  //  cpu_speed value obtained by the constructor
  ret_val = ret_val / cpu_speed;
#endif

  return ret_val;
#endif // __MWERKS__
}
/*			END OF FUNCTION user_time			*/

/************************************************************************/
/*									*/
/*			FUNCTION system_time				*/
/*									*/
/*	system_time reports the current amount of system cpu time       */
/*   accumulated by this tTimer.  If the timer is currently off, 	*/
/*   this is just the accumulated time.  If the tTimer is running, this  */
/*   is the accumulated time plus the time since the timer was last     */
/*   started.								*/
/*									*/
/************************************************************************/

double tTimer::system_time()
{
#ifdef __MWERKS__
// For now, stub out all Timer guts for MetroWerks
  return 0.0;
#else
  double ret_val;		//  Return value

#if ( defined(POOMA_T3E) || (POOMA_TFLOP) )
  // Not defined yet on T3E or TFLOP.
  // ERRORMSG("system_time() not defined." << endl);
  ret_val = -9999.0;
#else
  if (timer_state == TIMEROFF)
    {
      //  Timer is off, just return accumulated time
      ret_val = current_system_time;
    }
  else
    {
      //  Timer is on, return accumulated plus current
      times(&tmsbuf);
      ret_val = current_system_time + tmsbuf.tms_stime - last_system_time;
    }

  //  Convert from clock ticks to seconds using the
  //  cpu_speed value obtained by the constructor
  ret_val = ret_val / cpu_speed;
#endif

  return ret_val;
#endif // __MWERKS__
}

//=========================================================================
// 
// 
//                        End of tTimer.cpp
//    
//
//=========================================================================
