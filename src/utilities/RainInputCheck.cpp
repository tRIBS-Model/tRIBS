/***************************************************************************
**
**                           tRIBS Version 1.0
**
**              TIN-based Real-time Integrated Basin Simulator
**                       Ralph M. Parsons Laboratory
**                  Massachusetts Institute of Technology
**  
**
**                          Beta Release, 2002
**
**
**  RainInputCheck.cpp:  Utility program used for checking input
**                       precipitation files. Creates a log file and warns
**                       about values out of physically meaningful range
**
**  Program compiled separately from tRIBS as:
**
**  UNIX%  g++ -o <executable> RainInputCheck.cpp
**  Tru64% cxx -o <executable> RainInputCheck.cpp -model ansi
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
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>

#define  LOWER   0     //  <-- Lower threshold on values in grid  
#define  UPPER   25    //  <-- Upper threshold on values in grid  

static int days[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

int NR=0, MR=0;        //  <-- Dimensions of input array
double dR;             //  <-- Spatial resolution of R data
double dummy;          //  <-- NODATA value of input grid
double xllcR, yllcR;   //  <-- LL corner of input grid
double minGrid, maxGrid, sumGrid;
double  **gridIn;

int hourS, minS, dayS, monthS, yearS; 
int hourE, minE, dayE, monthE, yearE; // <-- End of specified period
int hour, day, month, year;           // <-- Current values
char mrainfileIn[30];

void printToFile(ofstream&);
void In_Mrain_Name (char *, char *, int, int, int, int);
void Destructor();
void GridStatistics();
void readInputGrid(char *);
void correctCalendarTime();
int  isTimeOver();

//=========================================================================
//
//
//                  Section 2: Main Routine
//
//
//=========================================================================

int main(int argc, char *argv[])
{
	int i, j;
	int flag = 0;
	int cflag = 0;
	double tempo;
	
	char namos[200];
	char pref[200];
	char ext[10];
	char logfile[200];
	char statfile[200];
	char tempstr[200];
	
	if(argc != 2) {
		cout <<"************************************************************\n";
		cout << "\nUsage: "<<argv[0] << "  File_containing_Info" <<endl;
		cout << "\n==> Order of information in the input file: " << endl;
		cout << "\t - Starting time:  mm/dd/yyyy/hh\n";
		cout << "\t - Finishing time: mm/dd/yyyy/hh\n";
		cout << "\t - Pathname for rainfall files and their prefix"
			<<"\n\t   (e.g., /platte/d4/Rain/Fall96/p)\n";
		cout << "\t - File extension in the *Input* files ('txt')\n";
		cout << "\t - Output log file name\n";
		cout << "\t - Output statistics file name\n"<<endl;
		cout <<"************************************************************\n";
		cout << endl << endl;
		exit(1);
	}
	
	ifstream Inp00(argv[1]);
	if (!Inp00) {
		cout <<"File "<<argv[1]<<" not found!!!"<<endl;
		exit(2);
	}
	
	Inp00 >> namos;
	sscanf(namos, "%02d/%02d/%04d/%02d", &monthS,&dayS,&yearS,&hourS);
	
	Inp00 >> namos;
	sscanf(namos, "%02d/%02d/%04d/%02d", &monthE,&dayE,&yearE,&hourE);
	
	Inp00 >> namos;
	sscanf(namos, "%200s", &pref);
	
	Inp00 >> namos;
	sscanf(namos, "%10s", &ext);
	
	Inp00 >> logfile;
	cout << "Output log file:  "<< logfile << endl;
	
	Inp00 >> statfile;
	cout << "Output statistics file:  "<< statfile << endl;
	
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//                ----- DATA CHECK -----
	
	if (monthE < monthS && yearE == yearS) {
		cout <<"\t\nERROR! Check months!\n\n";
		exit(2);
	}
	if (dayE < dayS && monthE == monthS) {
		cout << "\n\tERROR! Check days!\n\n";
		exit(2);
	}
	if (hourE < hourS && dayE == dayS && monthE == monthS) {
		cout << "\n\tERROR! Check hours!\n\n";
		exit(2);
	}
	if (hourE >= 24 || dayE > days[monthE] || dayS > days[monthS]) {
		cout << "\n\tERROR! Check hours & dates!\n\n";
		exit(2);
	}
	
	// Initialize time counter
	hour  = hourS;
	day   = dayS;
	month = monthS;
	year  = yearS;
	
	// Open program output files          
	ofstream Otp0(logfile);
	ofstream Otp1(statfile);
	Otp0<<"************************************************************\n";
	Otp0<<"Logfile for grid files located in: "<<pref<<"*"<<endl;
	Otp0<<"************************************************************\n\n";
	Otp1<<"   Time    Grid Max   Grid Min   Grid Sum"<<endl;
	
	// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	// +++++++++++++++++++++++++++  LOOP  +++++++++++++++++++++++++++
	// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	
	flag = 0;
	while ( !flag ) {
		
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
		In_Mrain_Name(pref, ext, hour, day, month, year);
		
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
		cflag = 1;
		ifstream Inp0(mrainfileIn);
		
		while (!Inp0 && !flag) {
			
			if (cflag) {
				snprintf(tempstr,sizeof(tempstr),"%02d/%02d/%04d/%02d", month, day, year, hour);
				Otp0<<"Files missing from "<<tempstr;
			}
			cflag = 0;
			
			// Mark missing data in the statistics file
			snprintf(tempstr,sizeof(tempstr),"%02d%02d%04d%02d", month, day, year, hour);
			Otp1<<tempstr<<"   -999.0\t-999.0\t-999.0"<<endl;

			// To always keep track of the last missing file
			hourS = hour;
			dayS = day;
			monthS = month;
			yearS = year;
			
			// Advance in time
			hour++;
			correctCalendarTime();
			
			In_Mrain_Name(pref, ext, hour, day, month, year);
			Inp0.open(mrainfileIn);
			
			if ( isTimeOver() )
				flag = 1;
		}
		Inp0.close();
		
		if ( flag ) {
			snprintf(tempstr,sizeof(tempstr),"%02d/%02d/%04d/%02d", month, day, year, hour);
			Otp0<<" to "<<tempstr<<endl;
			cflag = 1;
			// Mark missing data in the statistics file
			snprintf(tempstr,sizeof(tempstr),"%02d%02d%04d%02d", month, day, year, hour);
			Otp1<<tempstr<<"   -999.0\t-999.0\t-999.0"<<endl;
		}
		
		if ( !flag ) {
			
			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			
			minGrid = 1E6;
			maxGrid = -1E6; 
			sumGrid = 0.0;
			
			readInputGrid( mrainfileIn );
			
			GridStatistics();
			
			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			
			// If there was interruption in file sequence -->
			// need to output info about next file
			
			if (!cflag) {
				snprintf(tempstr,sizeof(tempstr),"%02d/%02d/%04d/%02d", monthS, dayS, yearS, hourS);
				Otp0<<" to "<<tempstr<<";  Next file Grid Max = "<<maxGrid;
				if (maxGrid > 0)
					Otp0<<"<- Correct the gap"<<endl;
				else 
					Otp0<<endl;
				cflag = 1;
			}
			
			// If values are suspicious -- warn about that int output
			if (minGrid < LOWER || maxGrid > UPPER ) {
				snprintf(tempstr,sizeof(tempstr),"%02d/%02d/%04d/%02d", month, day, year, hour);
				Otp0<<tempstr<<"\t <- Erroneous data, check the file"<<endl;
			}
			
			snprintf(tempstr,sizeof(tempstr),"%02d%02d%04d%02d", month, day, year, hour);
			Otp1<<tempstr<<"   "<<maxGrid<<"\t"<<minGrid<<"\t"<<sumGrid<<endl;
			//printToFile(Otp0); // may want to output grid sometimes
			
			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			
			Destructor();
			
			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			
			hour++;
			correctCalendarTime();
			
			if ( isTimeOver() )
				flag = 1;
		}
	}
	return 0;
}

/***************************************************************************
**
**  GridStatistics()
**
**  Computes min, max, and total value of the grid
**
***************************************************************************/
void GridStatistics()
{
	int i, j;
	for (i=0; i < NR; i++)  {
		for (j=0; j < MR; j++) {
			if (gridIn[i][j] < minGrid )
				minGrid = gridIn[i][j];
			if (gridIn[i][j] > maxGrid )
				maxGrid = gridIn[i][j];
			sumGrid += gridIn[i][j];	  
		}
	}
	return;
}

/***************************************************************************
**
**  printToFile(ofstream &Otp0)
**
**  Prints to a file current input grid
**
***************************************************************************/
void printToFile(ofstream &Otp0)
{
	int i, j;
	for (i=0; i < NR; i++)  {
		for (j=0; j < MR; j++)    {
			Otp0 <<gridIn[i][j]<<" ";
		}
		Otp0 << endl;
	}
}

/***************************************************************************
**
**  In_Mrain_Name()
**
**  Composes a string for current input grid file
**
***************************************************************************/
void In_Mrain_Name(char *pref, char *ext, int hour, int day, int month, int year)
{
	//cout<<"\n-----------------------------------------------------\n";
	snprintf(mrainfileIn, sizeof(mrainfileIn), "%s%02d%02d%04d%02d.%s", pref, month, day, year, hour, ext);
	//cout << "File IN --> "<<mrainfileIn<< endl;
	return;
}

/***************************************************************************
** 
**  Destructor()
**
**  Deallocates memory occupied by gridIN[][]
**
***************************************************************************/
void Destructor()
{
	for (int i=0; i < NR; i++)
		delete [] gridIn[i];
	delete [] gridIn;
	return;
}

/***************************************************************************
**
**  readInputGrid(char *GridIn)
**
**  The function reads a grid in the standard ArcInfo/ArcView ASCII format
**  'GridIn' is the full pathname for the input grid file which is to be
**  specified by the calling function
**
***************************************************************************/
void readInputGrid(char *GridIn)
{
	int i,j;
	char tmp[20];
	char lineIn[300];
	double tempo;
	
	ifstream Inp0(GridIn);
	if (!Inp0) {
		cout <<"File "<<GridIn<<" not found!!!"<<endl;
		exit(2);
	}
	Inp0.getline(lineIn, 256);
	//cout << "LINE IN " << 1 << ": " << lineIn << endl;
	sscanf(lineIn, "%10s %d", tmp, &MR);
	//cout <<"\n\tGRID: tmp = "<<tmp<<"; # cols  MR = "<<MR<<endl;
	
	Inp0.getline(lineIn, 256);
	//cout << "LINE IN " << 2 << ": " << lineIn << endl;
	sscanf(lineIn, "%10s %d", tmp, &NR);
	//cout <<"\tGRID: tmp = "<<tmp<<"; # rows  NR = "<<NR<<endl;
	
	Inp0 >> tmp >> xllcR;
	//cout <<"\tGRID: tmp = "<<tmp<<";   xllcR = "<<yllcR<<endl;
	Inp0 >> tmp >> yllcR;
	//cout <<"\tGRID: tmp = "<<tmp<<";   yllcR = "<<yllcR<<endl;
	Inp0 >> tmp >> dR;
	//cout <<"\tGRID: tmp = "<<tmp<<";    dR    = "<<dR<<endl;
	Inp0 >> tmp >> dummy;
	//cout <<"\tGRID: tmp = "<<tmp<<"; dummy= "<<dummy<<endl;
	
	gridIn = new double* [NR];  // <== Actual values of the grid
	assert(gridIn != 0);
	
	for (i=0; i < NR; i++)  {
		gridIn[i] = new double[MR];
		assert(gridIn[i] != 0);
		for (j=0; j < MR; j++)    {
			Inp0 >> tempo;
			gridIn[i][j] = tempo;
		}
	}
	Inp0.close();
	return;
}

/***************************************************************************
**  correctCalendarTime()
**
**  Advances time counter
**
***************************************************************************/
void correctCalendarTime()
{
	if (hour == 24) {
		day++; 
		hour = 0;
	}
	if (month == 2 && day > days[month]) { // February
		if (hour == 0) {
			if ( (year%4) != 0 || (day > 29)) {
				month++; // After/Instead February 29th
				day = 1;
			}
		}     
	}
	else if (day > days[month]) {
		month++;
		day = 1;
		if (month > 12) {
			year++;
			month = 1;
		}
	}
}

/***************************************************************************
**
**  isTimeOver()
**
**  Informs if the end time tage has been reached
**
***************************************************************************/
int isTimeOver()
{
	if (hour == hourE && day == dayE && month == monthE && year == yearE ) {
		cout << "\n\n\t <<< Time limit has been reached... >>>\n\n";
		return 1;
	}
	else 
		return 0;
}

//=========================================================================
//
//
//                   End of RainInputCheck.cpp
//
//
//=========================================================================

