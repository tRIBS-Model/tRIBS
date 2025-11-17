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
**  predicates.h:   Header File for predicates.cpp
**
**  Functions called from member functions of tMesh to check for line
**  segment intersection since inexact arithmetic can lead to erroneous
**  answers. (see predicates.txt for more information)
**
***************************************************************************/

#ifndef PREDICATES_H
#define PREDICATES_H

//=========================================================================
//
//
//                  Section 1: predicates Include and Define Statements
//
//
//=========================================================================


#ifdef ALPHA_64
  #include <stdio.h> 
  #include <math.h> 
  #include <sys/time.h> 
#elif defined LINUX_32
  #include <cstdio> 
  #include <cmath> 
  #include <ctime> 
#elif defined WIN
  #include <stdio.h> 
  #include <math.h> 
  #include <sys/time.h> 
#else
    #include <cstdio>    // Replaces <stdio.h>
    #include <cmath>     // Replaces <math.h>
    #include <ctime>
#endif
 
#define INEXACT                         
#define tREAL double                      // float or double 
 
// #define Absolute(a)  abs(a)              // defined in <macros.h>
#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a)) 
 
#define Fast_Two_Sum_Tail(a, b, x, y) bvirt = x - a; y = b - bvirt 
 
#define Fast_Two_Sum(a, b, x, y) x = (tREAL) (a + b); Fast_Two_Sum_Tail(a, b, x, y) 

#define Two_Sum_Tail(a, b, x, y) bvirt = (tREAL) (x - a); avirt = x - bvirt; bround = b - bvirt; around = a - avirt; y = around + bround 
 
#define Two_Sum(a, b, x, y) x = (tREAL) (a + b); Two_Sum_Tail(a, b, x, y) 
 
#define Two_Diff_Tail(a, b, x, y) bvirt = (tREAL) (a - x); avirt = x + bvirt; bround = bvirt - b; around = a - avirt; y = around + bround 
 
#define Two_Diff(a, b, x, y) x = (tREAL) (a - b); Two_Diff_Tail(a, b, x, y)
 
#define tSPLIT(a, ahi, alo) c = (tREAL) (splitter * a); abig = (tREAL) (c - a); ahi = c - abig; alo = a - ahi
 
#define Two_Product_Tail(a, b, x, y) tSPLIT(a, ahi, alo); tSPLIT(b, bhi, blo); err1 = x - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); y = (alo * blo) - err3
 
#define Two_Product(a, b, x, y) x = (tREAL) (a * b); Two_Product_Tail(a, b, x, y) 
 
#define Two_Product_Presplit(a, b, bhi, blo, x, y) x = (tREAL) (a * b); tSPLIT(a, ahi, alo); err1 = x - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); y = (alo * blo) - err3
 
#define Square_Tail(a, x, y) tSPLIT(a, ahi, alo); err1 = x - (ahi * ahi); err3 = err1 - ((ahi + ahi) * alo); y = (alo * alo) - err3
 
#define Square(a, x, y) x = (tREAL) (a * a); Square_Tail(a, x, y) 
 
#define Two_One_Sum(a1, a0, b, x2, x1, x0) Two_Sum(a0, b , _i, x0); Two_Sum(a1, _i, x2, x1) 
 
#define Two_One_Diff(a1, a0, b, x2, x1, x0) Two_Diff(a0, b , _i, x0); Two_Sum( a1, _i, x2, x1) 
 
#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) Two_One_Sum(a1, a0, b0, _j, _0, x0); Two_One_Sum(_j, _0, b1, x3, x2, x1) 
 
#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) Two_One_Diff(a1, a0, b0, _j, _0, x0); Two_One_Diff(_j, _0, b1, x3, x2, x1) 
 

//=========================================================================
//
//
//                  Section 2: predicates Class Declaration
//
//
//=========================================================================

class Predicates
{
public:
   Predicates(); // just calls exactinit()
   ~Predicates() {} // doesn't do anything
   
   // Constructor:
   void exactinit();
   
   // Basic "exact" arithmetic:
   int grow_expansion(int elen, tREAL* e, tREAL b, tREAL* h);
   int grow_expansion_zeroelim(int elen, tREAL* e, tREAL b, tREAL* h);
   int expansion_sum(int elen, tREAL* e, int flen, tREAL* f, tREAL* h);
   int expansion_sum_zeroelim1(int elen, tREAL* e, int flen, tREAL* f,
                               tREAL* h);
   int expansion_sum_zeroelim2(int elen, tREAL* e, int flen, tREAL* f,
                               tREAL* h);
   int fast_expansion_sum(int elen, tREAL* e, int flen, tREAL* f, tREAL* h);
   int fast_expansion_sum_zeroelim(int elen, tREAL* e, int flen, tREAL* f,
                                   tREAL* h);
   int linear_expansion_sum(int elen, tREAL* e, int flen, tREAL* f, tREAL* h);
   int linear_expansion_sum_zeroelim(int elen, tREAL* e, int flen, tREAL* f,
                                     tREAL* h);
   int scale_expansion(int elen, tREAL* e, tREAL b, tREAL* h);
   int scale_expansion_zeroelim(int elen, tREAL* e, tREAL b, tREAL* h);
   int compress(int elen, tREAL* e, tREAL* h);
   tREAL estimate( int elen, tREAL* e );
   
   // Lancaster functions
   double DifferenceOfProductsOfDifferences( double a, double b,
                                             double c, double d,
                                             double e, double f,
                                             double g, double h );
   double AdaptDiffOfProdsOfDiffs( double a, double b,
                                   double c, double d,
                                   double e, double f,
                                   double g, double h,
                                   double sum );

   tREAL orient2dfast(tREAL *pa, tREAL *pb, tREAL *pc);
   tREAL orient2dadapt(tREAL *pa, tREAL *pb, tREAL *pc, tREAL detsum);
   tREAL orient2d(tREAL *pa, tREAL *pb, tREAL *pc);
   tREAL incirclefast(tREAL *pa, tREAL *pb, tREAL *pc, tREAL *pd);
   tREAL incircleadapt(tREAL *pa, tREAL *pb, tREAL *pc, tREAL *pd,
                      tREAL permanent);
   tREAL incircle(tREAL *pa, tREAL *pb, tREAL *pc, tREAL *pd);
   
private:
   tREAL splitter;     /* = 2^ceiling(p / 2) + 1.  Used to split floats in half. */ 
   tREAL epsilon;                /* = 2^(-p).  Used to estimate roundoff errors. */ 
   tREAL resulterrbound; 
   tREAL ccwerrboundA, ccwerrboundB, ccwerrboundC; 
   tREAL o3derrboundA, o3derrboundB, o3derrboundC; 
   tREAL iccerrboundA, iccerrboundB, iccerrboundC; 
   tREAL isperrboundA, isperrboundB, isperrboundC; 
};

#endif

//=========================================================================
//
//
//                           End of predicates.h 
//
//
//=========================================================================
