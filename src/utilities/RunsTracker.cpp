/***************************************************************************
**
**                           tRIBS Version 1.0
**
**              TIN-based Real-time Integrated Basin Simulator
**                       Ralph M. Parsons Laboratory
**                  Massachusetts Institute of Technology
**  
**
**                          Beta Release, 9/2001
**
**
**  RunTracker.cpp:  Utility Program for tRIBS used to store run information
**                   of the input *.in file and RUN_NAME_ID when using
**                   the Model Option -O
**
**  Program compiled separately from tRIBS as:
**     UNIX% CC -o <executable> RunsTracker.cpp
**     Example: CC -o ../tracker RunsTracker.cpp (inside utilities dir)
**
**  To compile the program on Tru64 (Alpha DEC Station):
**     parana:~ cxx -o <executable> RunsTracker.cpp -model ansi
** 
**  Note: RunsTracker must be run within *.in Directory
**     Example: ../tRIBS/tracker *.in runID
**     Alternatively, the tracker can be run using hill_tracker file
** 
***************************************************************************/


//=========================================================================
//
//
//                  Section 1: Include, Definitions
//                             Variable and Function Declarations
//
//
//=========================================================================

#include<iostream.h>
#include<iomanip.h>
#include<fstream.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>

#define kMaxNameLength 80
#define kMaxNameSize 80
#define kCommentMark '#'

ifstream infile;          // *.IN file for tRIBS run
ofstream theOFStream;     // Output file to store all the info

void   ExtractTable(char *);
int    ReadItem( const int &, const char * );
long   ReadItem( const long &, const char * );
double ReadItem( const double &, const char * );
void   ReadItem( char *, const char * );


//=========================================================================
//
//
//                  Section 2: Main Routine
//
//
//=========================================================================


int main(int argc, char *argv[]){
	int    BRoption;
	double tempo;
	char   FileOfInterest[kMaxNameSize];
	
	//Argument List
	
	if(argc != 3) {
		cout<<"-----------------------------------------------------------------"
		<<"-------";
		cout<<"\n\n\t\t tRIBS -- Version 1.0";
		cout<<"\n\t\t tRIBS Model: RunTracker Utility";
		cout<<"\n\t\t Ralph M. Parsons Laboratory";
		cout<<"\n\t\t Massachusetts Institute of Technology";
		cout<<"\n\n\t\t Release, 01/2002 \n\n";
		cout<<"-----------------------------------------------------------------"
			<<"-------"<<endl<<endl;
		cout<<"\nUsage: "<<argv[0]<<"  tRIBS_Input_File  Run_ID\n"<<endl;
		cout<<"RunsTracker must be run within *.in Directory"<<endl<<endl;
		cout<<"Exiting Program..."<<endl<<endl;
		exit(1);
	}
	infile.open( argv[1] );
	if (!infile) {
		cerr<<"\nFile "<<argv[1]<<" not found! \nExiting Program..."<<endl;
		exit(2);
	}
	theOFStream.open(argv[2]);
	if( !theOFStream.good() ) {
		cerr<<"\nError: Cannot create the file for output."
        <<"\nExiting Program..."<<endl<<flush;
		exit(2);
	}
	theOFStream.setf( ios::right, ios::adjustfield );
	
	//Read Items from the *.in File
	
	theOFStream<<"-----------------------------------------------------------------"
		<<"-------";
	theOFStream<<"\n\t\t tRIBS Model: RunTracker Utility";
	theOFStream<<"\n\t\t Ralph M. Parsons Laboratory";
	theOFStream<<"\n\t\t Massachusetts Institute of Technology";
	theOFStream<<"\n\t\t Release, 01/2002\n";
	theOFStream<<"-----------------------------------------------------------------"
		<<"-------"<<endl<<endl;
	
	ReadItem(FileOfInterest , "STARTINGDATETIME");
	theOFStream<<"STARTDATE:  \t"<<FileOfInterest<<endl;
	tempo = ReadItem(tempo, "RUNTIME"); 
	theOFStream<<"RUNTIME:    \t"<<tempo<<" hours"<<endl;
	tempo = ReadItem(tempo, "INTSTORMMAX");
	theOFStream<<"INTSTORMMAX:\t"<<(int)tempo<<" hours"<<endl;
	
	theOFStream<<"\n-Routing Vars-"<<endl;
	tempo = ReadItem(tempo, "BASEFLOW");
	theOFStream<<"BASEFLOW: \t"<<tempo<<endl;
	tempo = ReadItem(tempo, "KINEMVELCOEF");
	theOFStream<<"KINEMVELCOEF: \t"<<tempo;
	tempo = ReadItem(tempo, "FLOWEXP");
	theOFStream<<"\tFLOWEXP: \t"<<tempo<<endl;
	
	tempo = ReadItem(tempo, "CHANNELWIDTHCOEFF");
	if (tempo <= 0) {
		tempo = ReadItem(tempo, "CHANNELWIDTH");
		theOFStream<<"\tCHANNELWIDTH: \t"<<tempo<<endl;
	}
	else {
		theOFStream<<"WIDTHCOEFF: \t"<<tempo;
		tempo = ReadItem(tempo, "CHANNELWIDTHEXPNT");
		theOFStream<<"\tWIDTHEXPNT: \t"<<tempo<<endl;
		
		ReadItem(FileOfInterest, "CHANNELWIDTHFILE");
		ifstream InWidthFile(FileOfInterest);
		if (!InWidthFile)
			theOFStream<<"CHANNELWIDTH: \t-- VARIABLE (COMPUTED, NO MEASURED) --"<<endl;
		else {
			theOFStream<<"CHANNELWIDTH: \t-- VARIABLE (COMPUTED AND MEASURED) --";
			BRoption= ReadItem(BRoption, "WIDTHINTERPOLATION");
			if (!BRoption)
				theOFStream<<" INTERP ON"<<endl;
			else
				theOFStream<<" INTERP OFF"<<endl;
		}
		InWidthFile.close();
	}
	
	tempo = ReadItem(tempo, "CHANNELROUGHNESS");
	theOFStream<<"CHANNELROUGHNESS: \t"<<tempo<<endl;
	
	theOFStream<<"\n-GW and Bedrock info-"<<endl;
	ReadItem(FileOfInterest, "GWATERFILE" );
	theOFStream<<"GW Initial in file:\t"<<FileOfInterest<<endl;
	
	BRoption  = ReadItem(BRoption, "OPTBEDROCK");
	if ( !BRoption ) {
		tempo = ReadItem(tempo, "DEPTHTOBEDROCK");
		theOFStream<<"DEPTHTOBEDROCK: \t"<<tempo<<endl;
	}
	else if ( BRoption == 1) {
		ReadItem(FileOfInterest, "BEDROCKFILE");
		theOFStream<<"DEPTHTOBEDROCK in file\t"<<FileOfInterest<<endl;
	}
	else {
		cerr<<"Your option is "<<BRoption<<"\n...Assumed UNIFORM case..."<<endl;
		tempo = ReadItem(tempo, "DEPTHTOBEDROCK");
		theOFStream<<"DEPTHTOBEDROCK =\t"<<tempo<<endl;
	}
	theOFStream<<"\n-Soil Table-"<<endl;
	ReadItem(FileOfInterest, "SOILTABLENAME");
	ExtractTable(FileOfInterest);
	
	theOFStream<<"\n-Land Table-"<<endl;
	ReadItem(FileOfInterest, "LANDTABLENAME");
	ExtractTable(FileOfInterest);
	
	theOFStream<<"\n-Simulation Options-"<<endl;
	BRoption = ReadItem(BRoption, "RAINSOURCE");
	theOFStream<<"RAINSOURCE:\t";
	if (BRoption == 1)
		theOFStream<<"Gridded Rainfall"<<endl;
	else if (BRoption == 2)
		theOFStream<<"Rain gauge"<<endl;
	else
		theOFStream<<"-- Unknown --"<<endl;
	tempo = ReadItem(tempo, "OPTEVAPOTRANS");
	theOFStream<<"OPTEVAPOTRANS:\t"<<(int)tempo<<endl;
	tempo = ReadItem(tempo, "OPTINTERCEPT");
	theOFStream<<"OPTINTERCEPT:\t"<<(int)tempo<<endl;
	tempo = ReadItem(tempo, "GFLUXOPTION");
	theOFStream<<"GFLUXOPTION:\t"<<(int)tempo<<endl;
	
	theOFStream<<endl<<endl;
	return 0;
}

//=========================================================================
//
//
//                  Section 3: Additional Functions
//
//
//=========================================================================

/***************************************************************************
**
** ExtractTable()
**
** Used to read parameters from soil and landuse reclassification tables
**
***************************************************************************/

void ExtractTable(char *Table)
{
	int  numClass, np;
	double variable;
	
	ifstream Inp0(Table);
	if (!Inp0) {
		cerr <<"\nFile "<<Table<<" not found!"<<endl;
		cerr << "RunsTracker must be run within *.in Directory"<<endl;
		cerr <<" \nExiting Program..."<<endl;
		exit(2);
	}
	Inp0 >> numClass >> np; // Reads # of soil classes and parameters
	
	for (int i=0; i < numClass; i++) {
		for (int j=0; j < np; j++) {
			Inp0>>variable;
			if (j == 0)
				theOFStream<<dec<<variable<<" ";
			else if (j == 1)
				theOFStream<<setw(5)<<setprecision(3)<<variable<<" ";
			else if (j == np-1)
				theOFStream<<setw(8)<<setprecision(7)<<variable<<" ";
			else
				theOFStream<<setw(5)<<setprecision(4)<<variable<<" ";
		}
		theOFStream<<endl<<flush;
	}
	Inp0.close();
}

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

int ReadItem( const int &datType, const char *itemCode )
{
	int item;
	char headerLine[kMaxNameSize];
	
	assert( infile.good() );
	
	infile.getline( headerLine,  kMaxNameSize);
	while( !( infile.eof() ) &&
		   ( headerLine[0]==kCommentMark ||
			 strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
		infile.getline( headerLine, kMaxNameLength );
	
	if( !( infile.eof() ) ){
		infile >> item;
		infile.ignore( 1, '\n' );
	}
	else{
		cerr << "\nMissing Parameter '" << itemCode
		<< "' in the Input File." << endl;
		cerr << "RunsTracker must be run within *.in Directory"<<endl;
		cerr << "\nExiting Program..."<<endl<<endl;
		exit(2);
	}
	
	infile.seekg( 0, ios::beg );
	return item;
}

/***************************************************************************
**  
**  tInputFile::ReadItem( const long &datType, const char *itemCode )
** 
***************************************************************************/

long ReadItem( const long &datType, const char *itemCode )
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
		cerr << "\nMissing Parameter '" << itemCode
		<< "' in the Input File." << endl;
		cerr << "RunsTracker must be run within *.in Directory"<<endl;
		cerr << "\nExiting Program..."<<endl<<endl;
		exit(2);
	}
	infile.seekg( 0, ios::beg );
	return item;
}

/***************************************************************************
**  
**  tInputFile::ReadItem( const double &datType, const char *itemCode )
** 
***************************************************************************/

double ReadItem( const double &datType, const char *itemCode )
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
		cerr << "\nMissing Parameter '" << itemCode
		<< "' in the Input File." << endl;
		cerr << "RunsTracker must be run within *.in Directory"<<endl;
		cerr << "\nExiting Program..."<<endl<<endl;
		exit(2);
	}
	infile.seekg( 0, ios::beg );
	return item;
}

/***************************************************************************
**  
**  tInputFile::ReadItem( char *theString, const char *itemCode )
** 
***************************************************************************/

void ReadItem(  char * theString, const char *itemCode )
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
	else{
		cerr << "\nMissing Parameter '" << itemCode
		<< "' in the Input File." << endl;
		cerr << "RunsTracker must be run within *.in Directory"<<endl;
		cerr << "\nExiting Program..."<<endl<<endl;
		exit(2);
	}
	infile.seekg( 0, ios::beg );
}

//=========================================================================
//
//
//                   End of RunsTracker.cpp
//
//
//=========================================================================


