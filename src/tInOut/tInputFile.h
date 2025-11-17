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
**  tInputFile.h: Header for tInputFile class and objects
**
**  tInputFile Class used in tRIBS for inputing parameters and pathnames
**  using keywords in a *.in file read within various classes.
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: tInputFile Include and Define Statements
//
//
//=========================================================================

#ifndef TINPUTFILE_H
#define TINPUTFILE_H

#include "src/Headers/Definitions.h"
#include "src/Headers/Classes.h"

#ifdef ALPHA_64
  #include <iostream>
  #include <fstream>
  #include <assert.h>
  #include <string.h>
  #include <stdlib.h>
#elif defined LINUX_32
  #include <iostream>
  #include <fstream>
  #include <cassert>
  #include <string> 
  #include <cstdlib>
#elif defined MAC
  #include <iostream>
  #include <fstream>
  #include <cassert>
  #include <string>
  #include <cstdlib>
#elif defined WIN
  #include <iostream>
  #include <fstream>
  #include <assert.h>
  #include <string.h>
  #include <stdlib.h>
#else 
  #include <iostream>
  #include <fstream>
  #include <assert.h>
  #include <string.h>
  #include <stdlib.h>
#endif

using namespace std;

//=========================================================================
//
//
//                  Section 2: tInputFile Class Definitions
//
//
//=========================================================================

class tInputFile
{
public:
  tInputFile();
  tInputFile( const char * );	
  ~tInputFile();
  int    IsItemIn( const char * );
  int    ReadItem( const int &, const char * );
  long   ReadItem( const long &, const char * );
  double ReadItem( const double &, const char * );
  void   ReadItem( char *, const char * );
  void   CloseOldAndOpenNew( const char * ); 
  char*  GetInFileName() { return InFileName; } 

private:
  ifstream infile;
  char InFileName[kMaxNameSize]; 
};

#endif

//=========================================================================
//
//
//                          End of tInputFile.h 
//
//
//=========================================================================
