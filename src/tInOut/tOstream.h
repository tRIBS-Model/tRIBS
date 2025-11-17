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
 
#ifndef TOSTREAM_H
#define TOSTREAM_H
 
#include <iostream>
#include <string>
 
//=========================================================================
//
//
//                  Section 2: tOstream Class Definitions
//
//
//=========================================================================
 
class tOstream {
public:
  tOstream(std::ostream& s);
  virtual ~tOstream ();
 
  /// ostream function
  tOstream& operator<< ( std::ostream& (*fmanip)(std::ostream &));
  /// bool
  tOstream& operator<< ( bool x_);
  /// short
  tOstream& operator<< ( short x_);
  /// unsigned short
  tOstream& operator<< ( unsigned short x_);
  /// int
  tOstream& operator<< ( int x_);
  /// unsigned int
  tOstream& operator<< ( unsigned int x_);
  /// long
  tOstream& operator<< ( long x_);
  /// unsigned long
  tOstream& operator<< ( unsigned long x_);
  /// float
  tOstream& operator<< ( float x_);
  /// double
  tOstream& operator<< ( double x_);
  /// long double
  tOstream& operator<< ( long double x_);
  /// char
  tOstream& operator<< ( char x_);
  /// char*
  tOstream& operator<< ( char* x_);
  /// string
  tOstream& operator<< ( std::string x_);
  /// const char*
  tOstream& operator<< ( const char* x_);
  /// const unsigned char*
  tOstream& operator<< ( const unsigned char* x_);
  /// const signed char*
  tOstream& operator<< ( const signed char* x_);
  /// setf(f1,f2)
  void setf(const std::ios_base::fmtflags& f1, const std::ios_base::fmtflags& f2);
  /// T&
  template <typename T>
  tOstream& operator<< ( const T & x_);

private:
    std::ostream& myostream;
};

#endif
 
//=========================================================================
//
//
//                          End of tOstream.h
//
//
//=========================================================================

