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
**  tOstream.h: Header for tOstream class and objects
**
**  tOstream Class used in tRIBS for writing output for a serial or
**  parallel run.
**
***************************************************************************/
 
//=========================================================================
//
//
//                  Section 1: tOstream Include and Define Statements
//
//
//=========================================================================
 
#include "src/tInOut/tOstream.h"

#ifdef PARALLEL_TRIBS
#include "src/tParallel/tParallel.h"
#endif
 
using namespace std;

//=========================================================================
//
//
//                  Section 1: tOstream Constructor/Destructor
//
//
//=========================================================================
 
tOstream::tOstream (ostream& s) : myostream ( s ) {}
tOstream::~tOstream () {}
 
/// ostream function
tOstream& tOstream::operator<< ( ostream & (*fmanip)(ostream &)) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  fmanip( myostream ); return *this; 
}

/// bool
tOstream& tOstream::operator<< ( bool x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// short
tOstream& tOstream::operator<< ( short x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// unsigned short
tOstream& tOstream::operator<< ( unsigned short x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// int
tOstream& tOstream::operator<< ( int x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// unsigned int
tOstream& tOstream::operator<< ( unsigned int x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// long
tOstream& tOstream::operator<< ( long x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// unsigned long
tOstream& tOstream::operator<< ( unsigned long x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// float
tOstream& tOstream::operator<< ( float x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// double
tOstream& tOstream::operator<< ( double x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// long double
tOstream& tOstream::operator<< ( long double x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// char
tOstream& tOstream::operator<< ( char x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// char*
tOstream& tOstream::operator<< ( char* x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// string
tOstream& tOstream::operator<< ( string x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// const char*
tOstream& tOstream::operator<< ( const char* x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// const unsigned char*
tOstream& tOstream::operator<< ( const unsigned char* x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// const signed char*
tOstream& tOstream::operator<< ( const signed char* x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

/// setf(f1,f2)
void tOstream::setf(const ios_base::fmtflags& f1, const ios_base::fmtflags& f2) 
{
  myostream.setf(f1, f2);
}

/// T&
template <typename T>
tOstream& tOstream::operator<< ( const T & x_) 
{
#ifdef PARALLEL_TRIBS
  if (!tParallel::isMaster()) return *this;
#endif
  myostream << x_ ; return *this; 
}

//=========================================================================
//
//
//                          End of tOstream.cpp
//
//
//=========================================================================

