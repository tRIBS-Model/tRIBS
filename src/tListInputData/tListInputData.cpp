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
**  tListInputData.cpp:	  Functions for class tListInputData.
**
**************************************************************************/

#include "src/tListInputData/tListInputData.h"
#include "src/Headers/globalIO.h"

//=========================================================================
//
//
//                  Section 1: tListInputData Constructor and Destructor
//
//
//=========================================================================

/**************************************************************************
**
**  tListInputData constructor
**
**  The constructor takes an input file and reads from it the base name 
**  containing the triangulation. It then opens <basename>.nodes,
**  <basename>.z, <basename>.edges, and <basename>.tri. Assuming the
**  files are valid, the desired time-slice is read from infile, and
**  the start of data for that time-slice is sought in each of the four
**  triangulation files. The arrays are dimensioned as needed, and
**  GetFileEntry() is called to read the data into the arrays. Note that
**  the time in each file is identified by a space character preceding it
**  on the same line.
**
**  Modifications for tRIBS simplifies routine by not having to search
**  through multiple time slices. tRIBS has only one time slice at 
**  INPUTTIME = 0. Option for CHILD maintained. 
**
**************************************************************************/

template< class tSubNode >
tListInputData< tSubNode >::
tListInputData( tInputFile &infile )                  
{
	int righttime;                   	//flag: found the right time slice
	double time, intime;             	//current & desired time
	char basename[80],               	//base name of input files
		inname[80];                  	//full name of an input file
	char nodeHeader[kMaxNameLength]; 	//header line read from input file
	char zHeader[kMaxNameLength]; 	//header line read from input file
	char edgHeader[kMaxNameLength]; 	//header line read from input file
	char triHeader[kMaxNameLength]; 	//header line read from input file
	
	Cout<<"\nCreating tListInputData for reading in existing mesh..."<<endl;
	
	// Read base name for triangulation files from infile
	infile.ReadItem( basename, "INPUTDATAFILE" );
	
	// Open each of the four files
	strcpy( inname, basename );
	strcat( inname, ".nodes" );
	nodeinfile.open(inname);    	//Node input file pointer
	
	strcpy( inname, basename );
	strcat( inname, ".edges" );
	edgeinfile.open(inname);   	//Edge input file pointer
	
	strcpy( inname, basename );
	strcat( inname, ".tri" );
	triinfile.open( inname );   	//Triangle input file pointer
	
	strcpy( inname, basename );
	strcat( inname, ".z" );
	zinfile.open( inname );     	//Elevations input file pointer
	
	// Make sure we found them
	if( !nodeinfile.good() || !edgeinfile.good() || !triinfile.good()
		|| !zinfile.good() ){
		cerr << "Error: I can't find one or more of the following files:\n"
		<< "\t" << basename << ".nodes\n"
		<< "\t" << basename << ".edges\n"
		<< "\t" << basename << ".tri\n"
		<< "\t" << basename << ".z\n";
		cerr <<"Unable to open triangulation input file(s).";
	}
	
	// Find out which time slice we want to extract
	intime = infile.ReadItem( intime, "INPUTTIME" );
	
	Cout<<"\nReading in existing mesh files for input time slice = "<< intime<<": "<<endl;
	Cout<<"'"<< basename << ".nodes'"<<endl;
	Cout<<"'"<< basename << ".edges'"<<endl;
	Cout<<"'"<< basename << ".tri'"<<endl;
	Cout<<"'"<< basename << ".z'"<<endl;
	
	
	// For tRIBS, use INPUTTIME = 0 as flag implying only one mesh
	// list (nodes, edges, tri, z) per file opened. Will read the
	// the number of elements directly, followed by the values in
	// the GetFileEntry() Function.
	
	
	if(intime == 0){
		nodeinfile >> time;
		nodeinfile >> nnodes;
		
		zinfile >> time;
		zinfile >> nnodes;
		
		edgeinfile >> time;
		edgeinfile >> nedges;
		
		triinfile >> time;
		triinfile >> ntri;
	}
	
	
	// Keep the following for CHILD compatibility
	
	else{    // Find specified input times in input data files and read # items.
		
		//Nodes:
		righttime = 0;
		time = 0;
		while( !( nodeinfile.eof() ) && !righttime ){
			nodeinfile.getline( nodeHeader, kMaxNameLength );
			cout<<"nodeheader[0] = "<<nodeHeader[0]<<endl;
			if( nodeHeader[0] == kTimeLineMark ){
				nodeinfile.seekg( -nodeinfile.gcount(), ios::cur );
				nodeinfile >> time;
				cout << "from node file, time = " << time << endl;
				cout << "from node file, intime = " <<intime <<endl;
				if( time == intime ) righttime = 1;
			}
		}
		if( !( nodeinfile.eof() ) ) nodeinfile >> nnodes;
		else{
			cerr << "\nCouldn't find the specified input time in the node file\n";
		}
		
		//"z" values:
		righttime = 0;
		time = 0;
		while( !( zinfile.eof() ) && !righttime ){
			zinfile.getline( zHeader, kMaxNameLength );
			cout<<"zheader[0] = "<<zHeader[0]<<endl;
			if( zHeader[0] == kTimeLineMark ){
				zinfile.seekg( -zinfile.gcount(), ios::cur );
				zinfile >> time;
				cout << "from z file, time = " << time << endl;
				cout << "from z file, intime = " <<intime <<endl;
				if( time == intime ) righttime = 1;
			}
		}
		if( !( zinfile.eof() ) ) zinfile >> nnodes;
		else{
			cerr << "Couldn't find the specified input time in elevation file\n";
		}
		
		//Edges:
		righttime = 0;
		time = 0;
		while( !( edgeinfile.eof() ) && !righttime ){
			edgeinfile.getline( edgHeader, kMaxNameLength );
			cout<<"edgheader[0] = "<<edgHeader[0]<<endl;
			if( edgHeader[0] == kTimeLineMark ){
				edgeinfile.seekg( -edgeinfile.gcount(), ios::cur );
				edgeinfile >> time;
				cout << "from edg file, time = " << time << endl;
				cout << "from edg file, intime = " <<intime <<endl;
				if( time == intime ) righttime = 1;
			}
		}
		if( !( edgeinfile.eof() ) ) edgeinfile >> nedges;
		else{
			cerr << "Couldn't find the specified input time in the edge file\n";
		}
		
		//Triangles:
		righttime = 0;
		time = 0;
		while( !( triinfile.eof() ) && !righttime ){
			triinfile.getline( triHeader, kMaxNameLength );
			cout<<"triheader[0] = "<<triHeader[0]<<endl;
			if( triHeader[0] == kTimeLineMark ){
				triinfile.seekg( -triinfile.gcount(), ios::cur );
				triinfile >> time;
				cout << "from tri file, time = " << time << endl;
				cout << "from tri file, intime = " <<intime <<endl;
				if( time == intime ) righttime = 1;
			}
		}
		if( !( triinfile.eof() ) ) triinfile >> ntri;
		else{
			cerr << "Couldn't find the specified input time in the tri file\n";
		}
	}
	
	// Dimension the arrays 
	x.setSize( nnodes );
	y.setSize( nnodes );
	z.setSize( nnodes );
	edgid.setSize( nnodes );
	boundflag.setSize( nnodes );
	orgid.setSize( nedges );
	destid.setSize( nedges );
	nextid.setSize( nedges );
	p0.setSize( ntri );
	p1.setSize( ntri );
	p2.setSize( ntri );
	e0.setSize( ntri );
	e1.setSize( ntri );
	e2.setSize( ntri );
	t0.setSize( ntri );
	t1.setSize( ntri );
	t2.setSize( ntri );
	
	// Read in data from file
	GetFileEntry();
	
	// Close the files
	nodeinfile.close();
	edgeinfile.close();
	triinfile.close();
	zinfile.close();
}

template< class tSubNode >
tListInputData< tSubNode >:: ~tListInputData(){ }

//=========================================================================
//
//
//                  Section 2: tListInputData Functions
//
//
//=========================================================================

/**************************************************************************
**
**  tListInputData::GetFileEntry()
**
**  Reads node, edge, and triangle data from the four triangulation input
**  files. Assumes that each files is open and valid and that the current
**  reading point in each corresponds the start of data for the desired
**  time-slice.
**
**************************************************************************/

template< class tSubNode >
void tListInputData< tSubNode >::
GetFileEntry()                 
{
	int i;
	for( i=0; i< nnodes; i++ ){
		nodeinfile >> x[i] >> y[i] >> edgid[i] >> boundflag[i];
		zinfile >> z[i];
	}
	
	for( i=0; i<nedges; i++ )
		edgeinfile >> orgid[i] >> destid[i] >> nextid[i];
	
	for( i=0; i< ntri; i++ )
		triinfile >> p0[i] >> p1[i] >> p2[i] >> t0[i] >> t1[i] >> t2[i] 
			>> e0[i] >> e1[i] >> e2[i];
	
}


/**************************************************************************
**
**  tListInputData::GetKeyEntry()
**
**  Provides alternative keyboard entry of triangulation data for 
**  testing purposes. Not currently supported.
**
*************************************************************************/

template< class tSubNode >
void tListInputData< tSubNode >::
GetKeyEntry()                   
{
	int i;
	for( i=0; i < nnodes; i++ ){
		cout << "x y z edgid boundary:" << endl;
		cin >> x[i] >> y[i] >> z[i] >> edgid[i] >> boundflag[i];
	}
	for( i=0; i < nedges; i++ ){
		cout << "orgid destid nextid" << endl;
		cin >> orgid[i] >> destid[i] >> nextid[i];
	}
	for( i=0; i< ntri; i++ ){
		cout << "nodeids (3), edgids (3), triangleids (3)" << endl;
		cin >> p0[i] >> p1[i] >> p2[i]
			>> e0[i] >> e1[i] >> e2[i]
			>> t0[i] >> t1[i] >> t2[i];
	}
	
}

//=========================================================================
//
//
//                         End of tListInputData.cpp
//
//
//=========================================================================
