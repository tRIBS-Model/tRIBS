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
**  tInputFile.cpp: Functions for class tInputFile (see tInputFile.h)
**
***************************************************************************/

#include "src/tInOut/tInputFile.h"
#include "src/Headers/globalIO.h"
#include <cstring>

//=========================================================================
//
//
//                  Section 1: tInputFile Constructors/Destructors
//
//
//=========================================================================

tInputFile::tInputFile( const char *filename )
{
	strcpy(InFileName, filename);
	infile.open( filename );
	if( !infile.good() ){
		cout << "\n\ntInputFile: Unable to open file '"
		<<InFileName<<"'."<<endl;
		cout << "Exiting Program...\n\n";
		exit(1);
	}
}

tInputFile::~tInputFile()
{
	Cout<<"tInputFile Object has been destroyed..."<<endl<<flush;
} 

void tInputFile::CloseOldAndOpenNew( const char * filename) {
	infile.close();
	strcpy(InFileName, filename); //
	infile.open( filename );
	if( !infile.good() ){
		cout << "\ntInputFile: Unable to open '"<< InFileName << "'." << endl;
		cout << "Exiting Program...\n\n";
		exit(1);
	}
}

//=========================================================================
//
//
//                  Section 2: tInputFile ReadItem Functions
//
//
//=========================================================================


/***************************************************************************
**  
**  tInputFile::ReadItem( const int &datType, const char *itemCode )
**
**  Reads one parameter from the file. The format is assumed to be a line
**  of text that begins with the keyword, followed by a line containing
**  the parameter to be read. The function is overloaded according to the
**  type of data desired. Arbitrary order of the items in the infile allowed.
**  Routine searches through list until it finds itemCode. 
**
***************************************************************************/
int tInputFile::ReadItem( const int &datType, const char *itemCode )
{
	int item;
	char headerLine[kMaxNameLength];	

	assert( infile.good() );
	
	infile.getline( headerLine, kMaxNameLength );
	while( !( infile.eof() ) &&
		   ( headerLine[0]==kCommentMark ||
			 strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
		infile.getline( headerLine, kMaxNameLength );
	
	if( !( infile.eof() ) ){
		infile >> item;
		infile.ignore( 1, '\n' );
	}
	else{
		cout<<"\nError: Expected to read the parameter '"<<itemCode
		<<"', but reached EOF first"<<endl;
		cout<<"\nMissing parameter in the input file..."<<endl;
		infile.clear(); 
		infile.seekg( 0, ios::beg );
		return -9999;
	}
	
	infile.seekg( 0, ios::beg );
	return item;
}

/***************************************************************************
**  
**  tInputFile::ReadItem( const long &datType, const char *itemCode )
** 
***************************************************************************/
long tInputFile::ReadItem( const long &datType, const char *itemCode )
{
	long item;
	char headerLine[kMaxNameLength];
	
	assert( infile.good() );
	
	infile.getline( headerLine, kMaxNameLength );
	while( !( infile.eof() ) &&
		   ( headerLine[0]==kCommentMark ||
			 strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
		infile.getline( headerLine, kMaxNameLength );
	
	if( !( infile.eof() ) ){
		infile >> item;
		infile.ignore( 1, '\n' );
	}
	else{
		cout<<"\nError: Expected to read the parameter '"<<itemCode
		<<"', but reached EOF first"<<endl;
		cout<<"\nMissing parameter in input file..."<<endl;
		infile.clear(); 
		infile.seekg( 0, ios::beg );
		item = -9999;
	}
	
	infile.seekg( 0, ios::beg );
	return item;
}

/***************************************************************************
**  
**  tInputFile::ReadItem( const double &datType, const char *itemCode )
** 
***************************************************************************/
double tInputFile::ReadItem( const double &datType, const char *itemCode )
{
	double item;
	char headerLine[kMaxNameLength];
	
	assert( infile.good() );
	
	infile.getline( headerLine, kMaxNameLength );
	while( !( infile.eof() ) &&
		   ( headerLine[0]==kCommentMark ||
			 strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
		infile.getline( headerLine, kMaxNameLength );
	
	if( !( infile.eof() ) ){
		infile >> item;
		infile.ignore( 1, '\n' );
	}
	else{
		cout<<"\nError: Expected to read the parameter '" << itemCode
		<<"', but reached EOF first" << endl;
		cout<<"\nMissing parameter in input file..."<<endl;
		infile.clear(); 
		infile.seekg( 0, ios::beg );
		item = -999999.;
	}
	infile.seekg( 0, ios::beg );
	return item;
}

/***************************************************************************
**  
**  tInputFile::ReadItem( char *theString, const char *itemCode )
** 
***************************************************************************/
void tInputFile::ReadItem( char * theString, const char *itemCode )
{
	char headerLine[kMaxNameLength];
	
	assert( infile.good() );
	
	infile.getline( headerLine, kMaxNameLength );
	while( !( infile.eof() ) &&
		   ( headerLine[0]==kCommentMark ||
			 strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
		infile.getline( headerLine, kMaxNameLength ); 
	
	if( !( infile.eof() ) ){
		infile.getline( theString, kMaxNameLength );
	}
	else {
		cout<<"\nError: Expected to read the parameter '" << itemCode
		<< "', but reached EOF first" << endl;
		cout<<"\nMissing parameter in input file..."<<endl;
		char errr[] = "-999";
		strcpy(theString, errr);
		infile.clear();
		infile.seekg( 0, ios::beg );
		return;
	}
	
	infile.seekg( 0, ios::beg );
}

/***************************************************************************
**  
**  tInputFile::IsItemIn( const char *itemCode )
** 
***************************************************************************/
int tInputFile::IsItemIn( const char *itemCode )
{
	int item;
	char headerLine[kMaxNameLength];
	
	assert( infile.good() );
	
	infile.getline( headerLine, kMaxNameLength );
	while( !( infile.eof() ) &&
		   ( headerLine[0]==kCommentMark ||
			 strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
		infile.getline( headerLine, kMaxNameLength );
	
	if ( !( infile.eof() ) ) {
		infile.seekg( 0, ios::beg );
		return 1;
	}
	else {
		infile.clear(); 
		infile.seekg( 0, ios::beg );
		return 0;
	}
}

//=========================================================================
//
//
//                        End of InputFile.cpp
//
//
//=========================================================================
