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
**  tParallel.cpp: Functions for class tParallel (see tParallel.h)
**
***************************************************************************/

#include "src/tParallel/tParallel.h"

#include <iostream>
#include <cassert>
#include <cstring>

using namespace std;

int tParallel::numProcs = 0;
int tParallel::myProc = -1;

list<double*> tParallel::buffers;
list<MPI_Request> tParallel::requests;

tParallel::tParallel() {}
tParallel::~tParallel() {}

int tParallel::getMyProc()
{
  return myProc;
}

int tParallel::getNumProcs()
{
  return numProcs;
}

/*************************************************************************
**
** Initialize MPI and variables
**
*************************************************************************/

void tParallel::initialize(int& argc, char** argv) {
  // Start up MPI
  MPI_Init(&argc, &argv);

  // Get number of processors and my processor
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
  if (myProc == MASTER_PROC) { 
    cout << "tParallel: MPI initialized, " << numProcs << " processor(s)." << endl;
  }
}

/*************************************************************************
**
** Finalize MPI and variables
**
*************************************************************************/

void tParallel::finalize() {
  // Delete data and zero
  numProcs = 0;
  myProc = -1;

  // Shutdown MPI
  MPI_Finalize();
}

/*************************************************************************
**
** Broadcast command line args to all nodes
**
*************************************************************************/
                                                                                
void tParallel::inputArgs(int& argc, char** argv) {
  // If only 1 processor, done
  if (numProcs == 1) return;

  // If Master, send command line args
  int olen = 0 ;
  char* obuf;
  char* nbuf;
  int targc = 0;
  if (myProc == MASTER_PROC) {
    // Number of args
    targc = argc;
    MPI_Bcast(&targc, 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);


    // Input file name
    size_t olen = strlen(argv[1]);

    MPI_Bcast(&olen, 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);

    // Allocate a buffer that is large enough to store the input file name
    obuf = new char[olen];

    // Copy the input file name to the buffer using strncpy
    strncpy(obuf, argv[1], olen);
//    obuf[sizeof(obuf) - 1] = '\0';
//
//    // Check if the input file name was truncated
//      if (olen > sizeof(obuf) - 1) {
//          cerr << "Warning: Input file name was truncated." << endl;
//      }

    MPI_Bcast(obuf, olen, MPI_CHAR, MASTER_PROC, MPI_COMM_WORLD);

    //cout << "Master argv 1 = " << argv[1] << endl;
    delete [] obuf;
    // Input options
    if (argc == 2) return; // No-command line flags
    obuf = new char[3]; //WR debug: added additional space for good measure
    for (int i = 2; i < argc; i++) { 
      strcpy(obuf, argv[i]);
      MPI_Bcast(obuf, 2, MPI_CHAR, MASTER_PROC, MPI_COMM_WORLD); 
      //cout << "Master argv " << i << " = " << argv[i] << endl;
    }
    delete [] obuf;
  } else {
    // Number of args
    MPI_Bcast(&targc, 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
    argc = targc;
    //cout << "Non-Master argc = " << argc << endl;
    // Input file name
    MPI_Bcast(&olen, 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
    obuf = new char[olen+1];
    MPI_Bcast(obuf, olen, MPI_CHAR, MASTER_PROC, MPI_COMM_WORLD);
    obuf[olen] = '\0';
    nbuf = new char[olen+1];
    strcpy(nbuf, obuf);
    argv[1] = nbuf;
    //cout << "Non-Master argv 1 = " << argv[1] << "." << endl;
    delete [] obuf;
    // Input options
    if (argc == 2) return;
    obuf = new char[3];
    for (int i = 2; i < argc; i++) {
      MPI_Bcast(obuf, 2, MPI_CHAR, MASTER_PROC, MPI_COMM_WORLD);
      obuf[2] = '\0';
      nbuf = new char[3];
      strcpy(nbuf, obuf);
      argv[i] = nbuf;
      //cout << "Non-Master argv " << i << " = " << argv[i] << "." << endl;
    }
  }
}

/***************************************************************************
**  
** Is this the MASTER node?
**  
***************************************************************************/
bool tParallel::isMaster() 
{ 
  if (myProc == 0) return true; 
  return false;
}

/***************************************************************************
**  
** Is this the last node?
**  
***************************************************************************/
bool tParallel::lastProc() 
{ 
  if (myProc == (numProcs-1)) return true;
  return false;
}

/***************************************************************************
**  
** A barrier.
**  
***************************************************************************/
void tParallel::barrier()
{ 
  MPI_Barrier(MPI_COMM_WORLD);
}

/***************************************************************************
**
** Send data to another processor without blocking
** Store the pointer to the buffer and matching request for later cleanup
**
***************************************************************************/
                                                                                
void tParallel::send(int tproc, int rtag, double* sdata, int scnt) {
  assert((tproc >= 0) && (tproc < numProcs));

  int tag = rtag;
  MPI_Request request;
  MPI_Isend(sdata, scnt, MPI_DOUBLE, tproc, tag, MPI_COMM_WORLD, &request);

  buffers.push_back(sdata);
  requests.push_back(request);
}

/***************************************************************************
**
** Receive data from another processor blocking.
**
***************************************************************************/

void tParallel::receive(int fproc, int rtag, double* rdata, int rcnt) {
  assert((fproc >= 0) && (fproc < numProcs));

  int tag = rtag;
  MPI_Status status;
  MPI_Recv(rdata, rcnt, MPI_DOUBLE, fproc, tag, MPI_COMM_WORLD, &status);
}

/***************************************************************************
**
** Check pending send requests for completion and delete corresponding buffer
**
***************************************************************************/
                                                                                
void tParallel::freeBuffers() {

   list<double*>::iterator biter = buffers.begin();
   list<MPI_Request>::iterator riter = requests.begin();

   int complete;
   MPI_Status status;

   while (biter != buffers.end()) { 
      MPI_Test(&(*riter), &complete, &status);
      if (complete == 1) {
         double* buf = (*biter);
         biter = buffers.erase(biter);
         riter = requests.erase(riter);
         delete [] buf; //WR debug: buf should be array and not single object
      } else {
         riter++;
         biter++;
      }
   }
}

/***************************************************************************
**
** Int global summation across processors.
**
***************************************************************************/
                                                                                
int tParallel::sum(int value) {
   if (numProcs == 1) return value;

   int valueSum = 0;
   int localValue = value;
   MPI_Reduce(&localValue, &valueSum, 1, MPI_INT, MPI_SUM, MASTER_PROC,
      MPI_COMM_WORLD);

   return valueSum;
}

/***************************************************************************
**
** Int global summation across processors.
** Sum broadcast to all processors.
**
***************************************************************************/

int tParallel::sumBroadcast(int value) {
   if (numProcs == 1) return value;

   int valueSum = 0;
   int localValue = value;
   MPI_Reduce(&localValue, &valueSum, 1, MPI_INT, MPI_SUM, MASTER_PROC,
      MPI_COMM_WORLD);
   MPI_Bcast(&valueSum, 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);

   return valueSum;
}

/***************************************************************************
**
** Double global summation across processors.
** Sum broadcast to all processors.
**
***************************************************************************/
                                                                                
double tParallel::sum(double value) {
   if (numProcs == 1) return value;

   double valueSum = 0;
   double localValue = value;
   MPI_Reduce(&localValue, &valueSum, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROC,
      MPI_COMM_WORLD);
   //MPI_Bcast(&valueSum, 1, MPI_DOUBLE, MASTER_PROC, MPI_COMM_WORLD);
                                                                                
   return valueSum;
}

/***************************************************************************
**
** Double global summation of multiple values across processors.
**
***************************************************************************/
                                                                                
double* tParallel::sum(double* value, int n) {
   double* valueSum = new double[n];
   if (numProcs == 1) {
      for (int i = 0; i < n; i++) valueSum[i] = value[i];
      return valueSum;
   }

   for (int i = 0; i < n; i++) valueSum[i] = 0;
   double* localValue = new double[n];
   for (int i = 0; i < n; i++) localValue[i] = value[i];
   MPI_Reduce(localValue, valueSum, n, MPI_DOUBLE, MPI_SUM, MASTER_PROC,
      MPI_COMM_WORLD);

   delete [] localValue; //WR debug, I think  deleting local value in this context is okay as MPI_reduce is blocking, Value Sum is later deleted, also implemented in other functions.

   return valueSum;
}

/***************************************************************************
**
** Double global minimum of multiple values across processors.
**
***************************************************************************/
                                                                                
double* tParallel::min(double* value, int n) {
   double* valueSum = new double[n];
   if (numProcs == 1) {
      for (int i = 0; i < n; i++) valueSum[i] = value[i];
      return valueSum;
   }

   for (int i = 0; i < n; i++) valueSum[i] = 0;
   double* localValue = new double[n];
   for (int i = 0; i < n; i++) localValue[i] = value[i];
   MPI_Reduce(localValue, valueSum, n, MPI_DOUBLE, MPI_MIN, MASTER_PROC,
      MPI_COMM_WORLD);

   delete [] localValue;

   return valueSum;
}

/***************************************************************************
**
** Double global maximum of multiple values across processors.
**
***************************************************************************/
                                                                                
double* tParallel::max(double* value, int n) {
   double* valueSum = new double[n];
   if (numProcs == 1) {
      for (int i = 0; i < n; i++) valueSum[i] = value[i];
      return valueSum;
   } 

   for (int i = 0; i < n; i++)
       valueSum[i] = 0;

   double* localValue = new double[n];

   for (int i = 0; i < n; i++) localValue[i] = value[i];
   MPI_Reduce(localValue, valueSum, n, MPI_DOUBLE, MPI_MAX, MASTER_PROC,
      MPI_COMM_WORLD);

   delete [] localValue;

   return valueSum;
}

/***************************************************************************
**
** Double global mean of multiple values across processors.
**
***************************************************************************/
                                                                                
double* tParallel::mean(double* value, int n) {
   double* valueSum = sum(value, n);
   if (numProcs == 1) return valueSum;

   for (int i = 0; i < n; i++) valueSum[i] /= numProcs;
                                                                                
   return valueSum;
}

/***************************************************************************
**
** Integer global collection across processors.
** Sum broadcast to all processors.
**
***************************************************************************/
                                                                                
int* tParallel::collect(int value) {
   int localValue = value;
   int* collectValue = new int[numProcs];
   if (numProcs == 1) {
      collectValue[0] = value;
      return collectValue;
   }

   MPI_Gather(&localValue, 1, MPI_INT, collectValue, 1, MPI_INT, 
      MASTER_PROC, MPI_COMM_WORLD);
   MPI_Bcast(collectValue, numProcs, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
                                                                                
   return collectValue;
}

//=========================================================================
//
//
//                        End of tParallel.cpp
//
//
//=========================================================================
