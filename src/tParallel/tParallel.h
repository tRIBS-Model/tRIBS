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
**  tParallel.h: Header for tParallel class and objects
**
**  tParallel Class used in tRIBS for the parallel version, providing 
**  connectivity and communications
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: tParallel Include and Define Statements
//
//
//=========================================================================

#ifndef TPARALLEL_H
#define TPARALLEL_H

#include <mpi.h>
#include <list>

//=========================================================================
//
//
//                  Section 2: tParallel Class Definitions
//
//
//=========================================================================

static const int MASTER_PROC = 0; // MPI Master node
static const int ERROR = -1;      // Returned on error

class tParallel {

public:
  /// Constructor
  tParallel();
  /// Destructor
  ~tParallel();
  /// Initialize MPI communications
  static void initialize(int& argc, char** argv);
  /// Broadcast command line args
  static void inputArgs(int& argc, char** argv);
  /// Finalize MPI communications
  static void finalize();

  /// Is this the MASTER node?
  static bool isMaster();
  /// Is this the last node?
  static bool lastProc();
  /// Return my processor
  static int getMyProc();
  /// Return number of processors
  static int getNumProcs();
  /// Barrier
  static void barrier();

  /// Send data 
  static void send(int tproc, int rtag, double* sdata, int scnt);
  /// Receive data
  static void receive(int fproc, int rtag, double* rdata, int rcnt);
  /// Delete buffers for which the send has completed
  static void freeBuffers();

  /// Global sum for a single value
  static int sum(int value);
  /// Global sum for a single value with broadcast
  static int sumBroadcast(int value);
  /// Global sum for a single value
  static double sum(double value);
  /// Global sum for a multiple values
  static double* sum(double* varray, int n);
  /// Global min for multiple values
  static double* min(double* varray, int n);
  /// Global max for multiple values
  static double* max(double* varray, int n);
  /// Global mean for multiple values
  static double* mean(double* varray, int n);
  /// Collect a value from each processor
  static int* collect(int value);

private:

  static int numProcs;              //!< # of procs
  static int myProc;                //!< # of my proc

  static std::list<double*> buffers;      //!< Buffers immediately sent
  static std::list<MPI_Request> requests; //!< Matching requests to check
};

#endif

//=========================================================================
//
//
//                          End of tParallel.h 
//
//
//=========================================================================
