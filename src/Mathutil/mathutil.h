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
**  mathutil.cpp: Math routines addition to math library. From Numerical 
**                Recipes in C by Press et. al (See mathutil.h). Includes
**                ANSI C (only) version of the Numerical Recipes utility 
**                file complex.c.
**
***************************************************************************/

#ifndef MATHUTIL_H
#define MATHUTIL_H


#ifdef ALPHA_64
  #include <math.h>
  #include<iostream>
  #include<iomanip>
  #include<assert.h>
#elif defined LINUX_32
  #include <cmath>
  #include <iostream>
  #include <cstdlib>
  #include <iomanip>
  #include <cassert>
  #include <cstdio>
#elif defined MAC
  #include <cmath>
  #include <iostream>
  #include <cstdlib>
  #include <iomanip>
  #include <cassert>
  #include <cstdio>
#elif defined WIN
  #include <math.h>
  #include<iostream>
  #include<iomanip>
  #include<assert.h>
#else 
  #include <math.h>
  #include<iostream>
  #include<iomanip>
  #include<assert.h>
#endif

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

using namespace std;

//=========================================================================
//
//
//                  Section 1: mathutil Declarations
//
//
//=========================================================================

int    rand1_00();
int    getRandomSeed();
void   setRandomSeed(int seed);

double ran3( long *idum );
double uniform();
double uniform(double, double );
double random_expon(double);
double cdfexpon(double, double);
double random_normal(double, double);
double random_gamma(double, long *);
double random_Weibull(double, double);

unsigned long random_poisson(double);
unsigned long random_binomial(unsigned long, double);
long random_hipergeom(unsigned long, unsigned long,
		      unsigned long, double);
double hipergeom_p0(unsigned long, unsigned long, unsigned long);
long double fact(unsigned long,unsigned long);
double EstimateAR1Var(double, double, double, double);

// Functions and data for creating random variables from Beta d-n 
float genbet(float aa,float bb);
float ranf();
long ignlgi();
void setall(long iseed1,long iseed2);
void initgn(long isdtyp);
void inrgcm();
long mltmod(long a,long s,long m);
void gssst(long getset,long *qset);
void gsrgs(long getset,long *qvalue);
void gscgn(long getset,long *g);

// Functions and data for creating random variables from Gamm d-n 
float gengam(float a,float r);
float sgamma(float a);
float snorm();
float sexpo();
float fsign( float num, float sign );
void  ftnstop(char* msg);

static long Xm1, Xm2, Xa1, Xa2, Xcg1[32], Xcg2[32], Xa1w, Xa2w;
static long Xig1[32], Xig2[32], Xlg1[32], Xlg2[32], Xa1vw, Xa2vw;
static long Xqanti[32];

static int rand1_00_seed;

//=========================================================================
//
//
//                  Section 2: Matrix and vector operations
//
//
//=========================================================================

void MatrixVectorProduct(double **, int, int, double *, double *);
void VectorSummation(double *, double *, int, double *);
double *Vector(long nl, long nh);
void free_Vector(double *v, long nl, long nh);
void nrerror(const char error_text[]);
void ludcmp(double a[][2], int n, int indx[], double *d);
void lubksb(double a[][2], int n, int indx[], double b[]);


//=========================================================================
//
//
//                  Section 3: complex1 Declarations
//
//
//=========================================================================

typedef struct FCOMPLEX {double r, i;} fcomplex;

void setComplex(fcomplex *a, double c, double d);
void setCompNumber(fcomplex *a, fcomplex b);
void Print(fcomplex *a);
double Cr(fcomplex *a);
double Ci(fcomplex *a);
void subtReal(fcomplex *a, double d);
void addReal(fcomplex *a, double d);
fcomplex CaddFF(fcomplex *a, double d);
fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex Complex(float re, float im);
fcomplex Conjg(fcomplex z);
fcomplex Cdiv(fcomplex a, fcomplex b);
float Cabs(fcomplex z);
fcomplex Csqrt(fcomplex z);
fcomplex RCmul(double x, fcomplex a);

double rtsafe(void (*)(double, double, double, double, double *, double *),
                       double, double, double, double, double, double, double);
void polynomialH(double, double, double, double, double *, double *);
void QuadraticFormula(double, double, double, double *, double *);

#endif


//=========================================================================
//
//
//                        End of mathutil.h
//
//
//=========================================================================
