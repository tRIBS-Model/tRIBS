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

#include "src/Mathutil/predicates.h"


Predicates::Predicates(){
	exactinit();
}

//=========================================================================
//
//
//                  Section 1: predicates Functions
//
//
//=========================================================================

/*****************************************************************************
**                                                                            
**  exactinit()   Initialize the variables used for exact arithmetic.         
**                                                                            
**  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   
**  floating-point arithmetic.  `epsilon' bounds the relative roundoff        
**  error.  It is used for floating-point error analysis.                     
**                                                                           
**  `splitter' is used to split floating-point numbers into two half-         
**  length significands for exact multiplication.                             
**                                                                            
**  I imagine that a highly optimizing compiler might be too smart for its   
**  own good, and somehow cause this routine to fail, if it pretends that     
**  floating-point arithmetic is too much like real arithmetic.              
**                                                                            
**  Don't change this routine unless you fully understand it.                 
**                                                                            
*****************************************************************************/ 

void Predicates::exactinit() 
{ 
	tREAL half; 
	tREAL check, lastcheck; 
	int every_other; 
	
	every_other = 1; 
	half = 0.5; 
	epsilon = 1.0; // Predicates private data member
	splitter = 1.0;  // Predicates private data member
	check = 1.0; 
	
	do { 
		lastcheck = check; 
		epsilon *= half; 
		if (every_other) { 
			splitter *= 2.0; 
		} 
		every_other = !every_other; 
		check = 1.0 + epsilon; 
	} while ((check != 1.0) && (check != lastcheck)); 
	splitter += 1.0; 
	
	resulterrbound = (3.0 + 8.0 * epsilon) * epsilon; 
	ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon; 
	ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon; 
	ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon; 
	o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon; 
	o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon; 
	o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon; 
	iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon; 
	iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon; 
	iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon; 
	isperrboundA = (16.0 + 224.0 * epsilon) * epsilon; 
	isperrboundB = (5.0 + 72.0 * epsilon) * epsilon; 
	isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon; 
} 

/***************************************************************************** 
**                                                                            
**  grow_expansion()   Add a scalar to an expansion.                          
**                                                                            
**  Sets h = e + b.  See the long version of my paper for details.            
**                                                                           
**  Maintains the nonoverlapping property.  If round-to-even is used (as      
**  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent     
**  properties as well.  (That is, if e has one of these properties, so       
**  will h.)                                                                  
**                                                                            
*****************************************************************************/ 

int Predicates::grow_expansion(int elen, tREAL* e, tREAL b, tREAL* h) 
{ 
	tREAL Q; 
	INEXACT tREAL Qnew; 
	int eindex; 
	tREAL enow; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	
	Q = b; 
	for (eindex = 0; eindex < elen; eindex++) { 
		enow = e[eindex]; 
		Two_Sum(Q, enow, Qnew, h[eindex]); 
		Q = Qnew; 
	} 
	h[eindex] = Q; 
	return eindex + 1; 
} 

/***************************************************************************** 
**                                                                            
**  grow_expansion_zeroelim()   Add a scalar to an expansion, eliminating     
**                              zero components from the output expansion.    
**                                                                            
**  Sets h = e + b.  See the long version of my paper for details.           
**                                                                            
**  Maintains the nonoverlapping property.  If round-to-even is used (as      
**  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent     
**  properties as well.  (That is, if e has one of these properties, so       
**  will h.)                                                                  
**                                                                            
*****************************************************************************/ 

int Predicates::grow_expansion_zeroelim(int elen, tREAL* e, tREAL b, tREAL* h)  
{ 
	tREAL Q, hh; 
	INEXACT tREAL Qnew; 
	int eindex, hindex; 
	tREAL enow; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	
	hindex = 0; 
	Q = b; 
	for (eindex = 0; eindex < elen; eindex++) { 
		enow = e[eindex]; 
		Two_Sum(Q, enow, Qnew, hh); 
		Q = Qnew; 
		if (hh != 0.0) { 
			h[hindex++] = hh; 
		} 
	} 
	if ((Q != 0.0) || (hindex == 0)) { 
		h[hindex++] = Q; 
	} 
	return hindex; 
} 

/***************************************************************************** 
**                                                                            
**  expansion_sum()   Sum two expansions.                                    
**                                                                            
**  Sets h = e + f.  See the long version of my paper for details.            
**                                                                           
**  Maintains the nonoverlapping property.  If round-to-even is used (as      
**  with IEEE 754), maintains the nonadjacent property as well.  (That is,    
**  if e has one of these properties, so will h.)  Does NOT maintain the     
**  strongly nonoverlapping property.                                         
**                                                                           
*****************************************************************************/ 

int Predicates::expansion_sum(int elen, tREAL* e, int flen, tREAL* f, tREAL* h)
{ 
	tREAL Q; 
	INEXACT tREAL Qnew; 
	int findex, hindex, hlast; 
	tREAL hnow; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	
	Q = f[0]; 
	for (hindex = 0; hindex < elen; hindex++) { 
		hnow = e[hindex]; 
		Two_Sum(Q, hnow, Qnew, h[hindex]); 
		Q = Qnew; 
	} 
	h[hindex] = Q; 
	hlast = hindex; 
	for (findex = 1; findex < flen; findex++) { 
		Q = f[findex]; 
		for (hindex = findex; hindex <= hlast; hindex++) { 
			hnow = h[hindex]; 
			Two_Sum(Q, hnow, Qnew, h[hindex]); 
			Q = Qnew; 
		} 
		h[++hlast] = Q; 
	} 
	return hlast + 1; 
} 

/***************************************************************************** 
**                                                                            
**  expansion_sum_zeroelim1()   Sum two expansions, eliminating zero          
**                              components from the output expansion.         
**                                                                            
**  Sets h = e + f.  See the long version of my paper for details.            
**                                                                            
**  Maintains the nonoverlapping property.  If round-to-even is used (as      
**  with IEEE 754), maintains the nonadjacent property as well.  (That is,    
**  if e has one of these properties, so will h.)  Does NOT maintain the     
**  strongly nonoverlapping property.                                         
**                                                                           
*****************************************************************************/ 

int Predicates::expansion_sum_zeroelim1(int elen, tREAL* e, int flen,
                                        tREAL* f, tREAL* h) 
{ 
	tREAL Q; 
	INEXACT tREAL Qnew; 
	int index, findex, hindex, hlast; 
	tREAL hnow; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	
	Q = f[0]; 
	for (hindex = 0; hindex < elen; hindex++) { 
		hnow = e[hindex]; 
		Two_Sum(Q, hnow, Qnew, h[hindex]); 
		Q = Qnew; 
	} 
	h[hindex] = Q; 
	hlast = hindex; 
	for (findex = 1; findex < flen; findex++) { 
		Q = f[findex]; 
		for (hindex = findex; hindex <= hlast; hindex++) { 
			hnow = h[hindex]; 
			Two_Sum(Q, hnow, Qnew, h[hindex]); 
			Q = Qnew; 
		} 
		h[++hlast] = Q; 
	} 
	hindex = -1; 
	for (index = 0; index <= hlast; index++) { 
		hnow = h[index]; 
		if (hnow != 0.0) { 
			h[++hindex] = hnow; 
		} 
	} 
	if (hindex == -1) { 
		return 1; 
	} else { 
		return hindex + 1; 
	} 
} 

/***************************************************************************** 
**                                                                            
**  expansion_sum_zeroelim2()   Sum two expansions, eliminating zero          
**                              components from the output expansion.         
**                                                                            
**  Sets h = e + f.  See the long version of my paper for details.            
**                                                                            
**  Maintains the nonoverlapping property.  If round-to-even is used (as      
**  with IEEE 754), maintains the nonadjacent property as well.  (That is,    
**  if e has one of these properties, so will h.)  Does NOT maintain the      
**  strongly nonoverlapping property.                                         
**                                                                            
*****************************************************************************/ 

int Predicates::expansion_sum_zeroelim2(int elen, tREAL* e,
                                        int flen, tREAL* f, tREAL* h) 
{ 
	tREAL Q, hh; 
	INEXACT tREAL Qnew; 
	int eindex, findex, hindex, hlast; 
	tREAL enow; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	
	hindex = 0; 
	Q = f[0]; 
	for (eindex = 0; eindex < elen; eindex++) { 
		enow = e[eindex]; 
		Two_Sum(Q, enow, Qnew, hh); 
		Q = Qnew; 
		if (hh != 0.0) { 
			h[hindex++] = hh; 
		} 
	} 
	h[hindex] = Q; 
	hlast = hindex; 
	for (findex = 1; findex < flen; findex++) { 
		hindex = 0; 
		Q = f[findex]; 
		for (eindex = 0; eindex <= hlast; eindex++) { 
			enow = h[eindex]; 
			Two_Sum(Q, enow, Qnew, hh); 
			Q = Qnew; 
			if (hh != 0) { 
				h[hindex++] = hh; 
			} 
		} 
		h[hindex] = Q; 
		hlast = hindex; 
	} 
	return hlast + 1; 
} 

/***************************************************************************** 
**                                                                            
**  fast_expansion_sum()   Sum two expansions.                                
**                                                                            
**  Sets h = e + f.  See the long version of my paper for details.            
**                                                                           
**  If round-to-even is used (as with IEEE 754), maintains the strongly       
**  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h    
**  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent       
**  properties.                                                               
**                                                                            
*****************************************************************************/ 

int Predicates::fast_expansion_sum(int elen, tREAL* e,
                                   int flen, tREAL* f, tREAL* h)
{ 
	tREAL Q; 
	INEXACT tREAL Qnew; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	int eindex, findex, hindex; 
	tREAL enow, fnow; 
	
	enow = e[0]; 
	fnow = f[0]; 
	eindex = findex = 0; 
	if ((fnow > enow) == (fnow > -enow)) { 
		Q = enow; 
		enow = e[++eindex]; 
	} else { 
		Q = fnow; 
		fnow = f[++findex]; 
	} 
	hindex = 0; 
	if ((eindex < elen) && (findex < flen)) { 
		if ((fnow > enow) == (fnow > -enow)) { 
			Fast_Two_Sum(enow, Q, Qnew, h[0]); 
			enow = e[++eindex]; 
		} else { 
			Fast_Two_Sum(fnow, Q, Qnew, h[0]); 
			fnow = f[++findex]; 
		} 
		Q = Qnew; 
		hindex = 1; 
		while ((eindex < elen) && (findex < flen)) { 
			if ((fnow > enow) == (fnow > -enow)) { 
				Two_Sum(Q, enow, Qnew, h[hindex]); 
				enow = e[++eindex]; 
			} else { 
				Two_Sum(Q, fnow, Qnew, h[hindex]); 
				fnow = f[++findex]; 
			} 
			Q = Qnew; 
			hindex++; 
		} 
	} 
	while (eindex < elen) { 
		Two_Sum(Q, enow, Qnew, h[hindex]); 
		enow = e[++eindex]; 
		Q = Qnew; 
		hindex++; 
	} 
	while (findex < flen) { 
		Two_Sum(Q, fnow, Qnew, h[hindex]); 
		fnow = f[++findex]; 
		Q = Qnew; 
		hindex++; 
	} 
	h[hindex] = Q; 
	return hindex + 1; 
} 

/***************************************************************************** 
**                                                                            
**  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero      
**                                  components from the output expansion.     
**                                                                            
**  Sets h = e + f.  See the long version of my paper for details.            
**                                                                            
**  If round-to-even is used (as with IEEE 754), maintains the strongly       
**  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h    
**  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      
**  properties.                                                               
**                                                                            
*****************************************************************************/ 

int Predicates::fast_expansion_sum_zeroelim(int elen, tREAL* e,
											int flen, tREAL* f, tREAL* h)
{ 
	tREAL Q; 
	INEXACT tREAL Qnew; 
	INEXACT tREAL hh; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	int eindex, findex, hindex; 
	tREAL enow, fnow; 
	
	enow = e[0]; 
	fnow = f[0]; 
	eindex = findex = 0; 
	if ((fnow > enow) == (fnow > -enow)) { 
		Q = enow; 
		enow = e[++eindex]; 
	} else { 
		Q = fnow; 
		fnow = f[++findex]; 
	} 
	hindex = 0; 
	if ((eindex < elen) && (findex < flen)) { 
		if ((fnow > enow) == (fnow > -enow)) { 
			Fast_Two_Sum(enow, Q, Qnew, hh); 
			enow = e[++eindex]; 
		} else { 
			Fast_Two_Sum(fnow, Q, Qnew, hh); 
			fnow = f[++findex]; 
		} 
		Q = Qnew; 
		if (hh != 0.0) { 
			h[hindex++] = hh; 
		} 
		while ((eindex < elen) && (findex < flen)) { 
			if ((fnow > enow) == (fnow > -enow)) { 
				Two_Sum(Q, enow, Qnew, hh); 
				enow = e[++eindex]; 
			} else { 
				Two_Sum(Q, fnow, Qnew, hh); 
				fnow = f[++findex]; 
			} 
			Q = Qnew; 
			if (hh != 0.0) { 
				h[hindex++] = hh; 
			} 
		} 
	} 
	while (eindex < elen) { 
		Two_Sum(Q, enow, Qnew, hh); 
		enow = e[++eindex]; 
		Q = Qnew; 
		if (hh != 0.0) { 
			h[hindex++] = hh; 
		} 
	} 
	while (findex < flen) { 
		Two_Sum(Q, fnow, Qnew, hh); 
		fnow = f[++findex]; 
		Q = Qnew; 
		if (hh != 0.0) { 
			h[hindex++] = hh; 
		} 
	} 
	if ((Q != 0.0) || (hindex == 0)) { 
		h[hindex++] = Q; 
	} 
	return hindex; 
} 

/***************************************************************************** 
**                                                                            
**  linear_expansion_sum()   Sum two expansions.                              
**                                                                            
**  Sets h = e + f.  See either version of my paper for details.             
**                                                                            
**  Maintains the nonoverlapping property.  (That is, if e is                 
**  nonoverlapping, h will be also.)                                          
**                                                                           
*****************************************************************************/ 

int Predicates::linear_expansion_sum(int elen, tREAL* e,
                                     int flen, tREAL* f, tREAL* h)
{ 
	tREAL Q, q; 
	INEXACT tREAL Qnew; 
	INEXACT tREAL R; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	int eindex, findex, hindex; 
	tREAL enow, fnow; 
	tREAL g0; 
	
	enow = e[0]; 
	fnow = f[0]; 
	eindex = findex = 0; 
	if ((fnow > enow) == (fnow > -enow)) { 
		g0 = enow; 
		enow = e[++eindex]; 
	} else { 
		g0 = fnow; 
		fnow = f[++findex]; 
	} 
	if ((eindex < elen) && ((findex >= flen) 
							|| ((fnow > enow) == (fnow > -enow)))) { 
		Fast_Two_Sum(enow, g0, Qnew, q); 
		enow = e[++eindex]; 
	} else { 
		Fast_Two_Sum(fnow, g0, Qnew, q); 
		fnow = f[++findex]; 
	} 
	Q = Qnew; 
	for (hindex = 0; hindex < elen + flen - 2; hindex++) { 
		if ((eindex < elen) && ((findex >= flen) 
								|| ((fnow > enow) == (fnow > -enow)))) { 
			Fast_Two_Sum(enow, q, R, h[hindex]); 
			enow = e[++eindex]; 
		} else { 
			Fast_Two_Sum(fnow, q, R, h[hindex]); 
			fnow = f[++findex]; 
		} 
		Two_Sum(Q, R, Qnew, q); 
		Q = Qnew; 
	} 
	h[hindex] = q; 
	h[hindex + 1] = Q; 
	return hindex + 2; 
} 

/***************************************************************************** 
**                                                                            
**  linear_expansion_sum_zeroelim()   Sum two expansions, eliminating zero   
**                                    components from the output expansion. 
**                                                                            
**  Sets h = e + f.  See either version of my paper for details.              
**                                                                           
**  Maintains the nonoverlapping property.  (That is, if e is                 
**  nonoverlapping, h will be also.)                                          
**                                                                            
*****************************************************************************/ 

int Predicates::linear_expansion_sum_zeroelim(int elen, tREAL* e,
											  int flen, tREAL* f, tREAL* h)
{ 
	tREAL Q, q, hh; 
	INEXACT tREAL Qnew; 
	INEXACT tREAL R; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	int eindex, findex, hindex; 
	int count; 
	tREAL enow, fnow; 
	tREAL g0; 
	
	enow = e[0]; 
	fnow = f[0]; 
	eindex = findex = 0; 
	hindex = 0; 
	if ((fnow > enow) == (fnow > -enow)) { 
		g0 = enow; 
		enow = e[++eindex]; 
	} else { 
		g0 = fnow; 
		fnow = f[++findex]; 
	} 
	if ((eindex < elen) && ((findex >= flen) 
							|| ((fnow > enow) == (fnow > -enow)))) { 
		Fast_Two_Sum(enow, g0, Qnew, q); 
		enow = e[++eindex]; 
	} else { 
		Fast_Two_Sum(fnow, g0, Qnew, q); 
		fnow = f[++findex]; 
	} 
	Q = Qnew; 
	for (count = 2; count < elen + flen; count++) { 
		if ((eindex < elen) && ((findex >= flen) 
								|| ((fnow > enow) == (fnow > -enow)))) { 
			Fast_Two_Sum(enow, q, R, hh); 
			enow = e[++eindex]; 
		} else { 
			Fast_Two_Sum(fnow, q, R, hh); 
			fnow = f[++findex]; 
		} 
		Two_Sum(Q, R, Qnew, q); 
		Q = Qnew; 
		if (hh != 0) { 
			h[hindex++] = hh; 
		} 
	} 
	if (q != 0) { 
		h[hindex++] = q; 
	} 
	if ((Q != 0.0) || (hindex == 0)) { 
		h[hindex++] = Q; 
	} 
	return hindex; 
} 

/****************************************************************************
**                                                                            
**  scale_expansion()   Multiply an expansion by a scalar.                    
**                                                                            
**  Sets h = be.  See either version of my paper for details.                 
**                                                                            
**  Maintains the nonoverlapping property.  If round-to-even is used (as      
**  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent     
**  properties as well.  (That is, if e has one of these properties, so       
**  will h.)                                                                
**                                                                            
*****************************************************************************/ 

int Predicates::scale_expansion(int elen, tREAL* e, tREAL b, tREAL* h) 
{ 
	INEXACT tREAL Q; 
	INEXACT tREAL sum; 
	INEXACT tREAL product1; 
	tREAL product0; 
	int eindex, hindex; 
	tREAL enow; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	INEXACT tREAL c; 
	INEXACT tREAL abig; 
	tREAL ahi, alo, bhi, blo; 
	tREAL err1, err2, err3;

    tSPLIT(b, bhi, blo);
	Two_Product_Presplit(e[0], b, bhi, blo, Q, h[0]); 
	hindex = 1; 
	for (eindex = 1; eindex < elen; eindex++) { 
		enow = e[eindex]; 
		Two_Product_Presplit(enow, b, bhi, blo, product1, product0); 
		Two_Sum(Q, product0, sum, h[hindex]); 
		hindex++; 
		Two_Sum(product1, sum, Q, h[hindex]); 
		hindex++; 
	} 
	h[hindex] = Q; 
	return elen + elen; 
} 

/***************************************************************************** 
**                                                                            
**  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          
**                               eliminating zero components from the         
**                               output expansion.                            
**                                                                           
**  Sets h = be.  See either version of my paper for details.                 
**                                                                            
**  Maintains the nonoverlapping property.  If round-to-even is used (as      
**  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent     
**  properties as well.  (That is, if e has one of these properties, so       
**  will h.)                                                                 
**                                                                            
*****************************************************************************/ 

int Predicates::scale_expansion_zeroelim(int elen, tREAL* e, tREAL b, tREAL* h) 
{ 
	INEXACT tREAL Q, sum; 
	tREAL hh; 
	INEXACT tREAL product1; 
	tREAL product0; 
	int eindex, hindex; 
	tREAL enow; 
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	INEXACT tREAL c; 
	INEXACT tREAL abig; 
	tREAL ahi, alo, bhi, blo; 
	tREAL err1, err2, err3;

    tSPLIT(b, bhi, blo);
	Two_Product_Presplit(e[0], b, bhi, blo, Q, hh); 
	hindex = 0; 
	if (hh != 0) { 
		h[hindex++] = hh; 
	} 
	for (eindex = 1; eindex < elen; eindex++) { 
		enow = e[eindex]; 
		Two_Product_Presplit(enow, b, bhi, blo, product1, product0); 
		Two_Sum(Q, product0, sum, hh); 
		if (hh != 0) { 
			h[hindex++] = hh; 
		} 
		Fast_Two_Sum(product1, sum, Q, hh); 
		if (hh != 0) { 
			h[hindex++] = hh; 
		} 
	} 
	if ((Q != 0.0) || (hindex == 0)) { 
		h[hindex++] = Q; 
	} 
	return hindex; 
} 

/***************************************************************************** 
**                                                                           
**  compress()   Compress an expansion.                                      
**                                                                            
**  See the long version of my paper for details.                             
**                                                                            
**  Maintains the nonoverlapping property.  If round-to-even is used (as      
**  with IEEE 754), then any nonoverlapping expansion is converted to a       
**  nonadjacent expansion.                                                    
**                                                                           
*****************************************************************************/ 

int Predicates::compress(int elen, tREAL* e, tREAL* h)
{ 
	tREAL Q, q; 
	INEXACT tREAL Qnew; 
	int eindex, hindex; 
	INEXACT tREAL bvirt; 
	tREAL enow, hnow; 
	int top, bottom; 
	
	bottom = elen - 1; 
	Q = e[bottom]; 
	for (eindex = elen - 2; eindex >= 0; eindex--) { 
		enow = e[eindex]; 
		Fast_Two_Sum(Q, enow, Qnew, q); 
		if (q != 0) { 
			h[bottom--] = Qnew; 
			Q = q; 
		} else { 
			Q = Qnew; 
		} 
	} 
	top = 0; 
	for (hindex = bottom + 1; hindex < elen; hindex++) { 
		hnow = h[hindex]; 
		Fast_Two_Sum(hnow, Q, Qnew, q); 
		if (q != 0) { 
			h[top++] = q; 
		} 
		Q = Qnew; 
	} 
	h[top] = Q; 
	return top + 1; 
} 

/***************************************************************************** 
**                                                                            
**  estimate()   Produce a one-word estimate of an expansion's value.        
**                                                                            
**  See either version of my paper for details.                              
**                                                                            
*****************************************************************************/ 

tREAL Predicates::estimate( int elen, tREAL* e ) 
{ 
	tREAL Q; 
	int eindex; 
	
	Q = e[0]; 
	for (eindex = 1; eindex < elen; eindex++) { 
		Q += e[eindex]; 
	} 
	return Q; 
} 

/***************************************************************************** 
**                                                                            
**  DifferenceOfProductsOfDifferences(...), AdaptDiffOfProdsOfDiffs(...)      
**                                                                           
**  Adaptations of orient2d and orient2dadapt                    
**                                                                           
*****************************************************************************/

double Predicates::DifferenceOfProductsOfDifferences( double a, double b,
                                                      double c, double d,
                                                      double e, double f,
                                                      double g, double h )
{
	double left;
	double right;
	double diff, errbound;
	double sum;
	
	left = ( a - b ) * ( c - d );
	right = ( e - f ) * ( g - h );
	diff = left - right;
	if( left > 0.0 )
	{
		if( right <= 0.0 ) return diff;
		else sum = left + right;
	}
	else if( left < 0.0 )
	{
		if( right >= 0.0 ) return diff;
		else sum = -left - right;
	}
	else return diff;
	errbound = ccwerrboundA * sum;
	if( diff >= errbound || -diff >= errbound ) return diff;
	return AdaptDiffOfProdsOfDiffs( a, b, c, d, e, f, g, h, sum );
}

double Predicates::AdaptDiffOfProdsOfDiffs( double terma, double termb,
                                            double termc, double termd,
                                            double terme, double termf,
                                            double termg, double termh,
                                            double sum )
{
	INEXACT tREAL diff1, diff2, diff3, diff4;
	tREAL diff1tail, diff2tail, diff3tail, diff4tail;
	INEXACT tREAL leftprod, rightprod;
	tREAL leftprodtail, rightprodtail;
	tREAL diff, errbound;
	tREAL B[4], C1[8], C2[12], D[16];
	INEXACT tREAL B3;
	int C1length, C2length, Dlength;
	tREAL u[4];
	INEXACT tREAL u3;
	INEXACT tREAL s1, t1;
	tREAL s0, t0;
	
	INEXACT tREAL bvirt;
	tREAL avirt, bround, around;
	INEXACT tREAL c;
	INEXACT tREAL abig;
	tREAL ahi, alo, bhi, blo;
	tREAL err1, err2, err3;
	INEXACT tREAL _i, _j;
	tREAL _0;
	
	diff1 = (tREAL) ( terma - termb );
	diff2 = (tREAL) ( termc - termd );
	diff3 = (tREAL) ( terme - termf );
	diff4 = (tREAL) ( termg - termh );
	
	Two_Product( diff1, diff2, leftprod, leftprodtail );
	Two_Product( diff3, diff4, rightprod, rightprodtail );
	
	Two_Two_Diff( leftprod, leftprodtail, rightprod, rightprodtail,
				  B3, B[2], B[1], B[0] );
	B[3] = B3;
	
	diff = estimate( 4, B);
	errbound = ccwerrboundB * sum;
	if( diff >= errbound || -diff >= errbound ) return diff;
	
	Two_Diff_Tail( terma, termb, diff1, diff1tail );
	Two_Diff_Tail( termc, termd, diff2, diff2tail );
	Two_Diff_Tail( terme, termf, diff3, diff3tail );
	Two_Diff_Tail( termg, termh, diff4, diff4tail );
	
	if( diff1tail == 0.0 && diff2tail == 0.0 &&
		diff3tail == 0.0 && diff4tail == 0.0 )
		return diff;
	
	errbound = ccwerrboundC * sum + resulterrbound * Absolute( diff );
	diff += ( diff1 * diff2tail + diff2 * diff1tail )
		- ( diff3 * diff4tail + diff4 * diff3tail );
	if( diff >= errbound || -diff >= errbound ) return diff;
	
	Two_Product( diff1tail, diff2, s1, s0 );
	Two_Product( diff3tail, diff4, t1, t0 );
	Two_Two_Diff( s1, s0, t1, t0, u3, u[2], u[1], u[0] );
	u[3] = u3;
	C1length = fast_expansion_sum_zeroelim( 4, B, 4, u, C1 );
	
	Two_Product( diff1, diff2tail, s1, s0 );
	Two_Product( diff3, diff4tail, t1, t0 );
	Two_Two_Diff( s1, s0, t1, t0, u3, u[2], u[1], u[0] );
	u[3] = u3;
	C2length = fast_expansion_sum_zeroelim( C1length, C1, 4, u, C2 );
	
	Two_Product( diff1tail, diff2tail, s1, s0 );
	Two_Product( diff3tail, diff4tail, t1, t0 );
	Two_Two_Diff( s1, s0, t1, t0, u3, u[2], u[1], u[0] );
	u[3] = u3;
	Dlength = fast_expansion_sum_zeroelim( C2length, C2, 4, u, D );
	
	return( D[Dlength - 1] );
}

/***************************************************************************** 
**                                                                            
**  orient2dfast()   Approximate 2D orientation test.  Nonrobust.             
**  orient2dexact()   Exact 2D orientation test.  Robust.                     
**  orient2dslow()   Another exact 2D orientation test.  Robust.              
**  orient2d()   Adaptive exact 2D orientation test.  Robust.                 
**                                                                            
**               Return a positive value if the points pa, pb, and pc occur   
**               in counterclockwise order; a negative value if they occur    
**               in clockwise order; and zero if they are collinear.  The     
**               result is also a rough approximation of twice the signed     
**               area of the triangle defined by the three points.            
**                                                                            
**  Only the first and last routine should be used; the middle two are for    
**  timings.                                                                  
**                                                                            
**  The last three use exact arithmetic to ensure a correct answer.  The      
**  result returned is the determinant of a matrix.  In orient2d() only,      
**  this determinant is computed adaptively, in the sense that exact          
**  arithmetic is used only to the degree it is needed to ensure that the     
**  returned value has the correct sign.  Hence, orient2d() is usually quite  
**  fast, but will run more slowly when the input points are collinear or     
**  nearly so.                                                                
**                                                                            
*****************************************************************************/ 

tREAL Predicates::orient2dfast(tREAL *pa, tREAL *pb, tREAL *pc) 
{ 
	tREAL acx, bcx, acy, bcy; 
	
	acx = pa[0] - pc[0]; 
	bcx = pb[0] - pc[0]; 
	acy = pa[1] - pc[1]; 
	bcy = pb[1] - pc[1]; 
	return acx * bcy - acy * bcx; 
} 


tREAL Predicates::orient2dadapt(tREAL *pa, tREAL *pb, tREAL *pc, tREAL detsum) 
{ 
	INEXACT tREAL acx, acy, bcx, bcy; 
	tREAL acxtail, acytail, bcxtail, bcytail; 
	INEXACT tREAL detleft, detright; 
	tREAL detlefttail, detrighttail; 
	tREAL det, errbound; 
	tREAL B[4], C1[8], C2[12], D[16]; 
	INEXACT tREAL B3; 
	int C1length, C2length, Dlength; 
	tREAL u[4]; 
	INEXACT tREAL u3; 
	INEXACT tREAL s1, t1; 
	tREAL s0, t0; 
	
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	INEXACT tREAL c; 
	INEXACT tREAL abig; 
	tREAL ahi, alo, bhi, blo; 
	tREAL err1, err2, err3; 
	INEXACT tREAL _i, _j; 
	tREAL _0; 
	
	acx = (tREAL) (pa[0] - pc[0]); 
	bcx = (tREAL) (pb[0] - pc[0]); 
	acy = (tREAL) (pa[1] - pc[1]); 
	bcy = (tREAL) (pb[1] - pc[1]); 
	
	Two_Product(acx, bcy, detleft, detlefttail); 
	Two_Product(acy, bcx, detright, detrighttail); 
	
	Two_Two_Diff(detleft, detlefttail, detright, detrighttail, 
				 B3, B[2], B[1], B[0]); 
	B[3] = B3; 
	
	det = estimate(4, B); 
	errbound = ccwerrboundB * detsum; 
	if ((det >= errbound) || (-det >= errbound)) { 
		return det; 
	} 
	
	Two_Diff_Tail(pa[0], pc[0], acx, acxtail); 
	Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail); 
	Two_Diff_Tail(pa[1], pc[1], acy, acytail); 
	Two_Diff_Tail(pb[1], pc[1], bcy, bcytail); 
	
	if ((acxtail == 0.0) && (acytail == 0.0) 
		&& (bcxtail == 0.0) && (bcytail == 0.0)) { 
		return det; 
	} 
	
	errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det); 
	det += (acx * bcytail + bcy * acxtail) 
		- (acy * bcxtail + bcx * acytail); 
	if ((det >= errbound) || (-det >= errbound)) { 
		return det; 
	} 
	
	Two_Product(acxtail, bcy, s1, s0); 
	Two_Product(acytail, bcx, t1, t0); 
	Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]); 
	u[3] = u3; 
	C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1); 
	
	Two_Product(acx, bcytail, s1, s0); 
	Two_Product(acy, bcxtail, t1, t0); 
	Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]); 
	u[3] = u3; 
	C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2); 
	
	Two_Product(acxtail, bcytail, s1, s0); 
	Two_Product(acytail, bcxtail, t1, t0); 
	Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]); 
	u[3] = u3; 
	Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D); 
	
	return(D[Dlength - 1]); 
} 

tREAL Predicates::orient2d(tREAL *pa, tREAL *pb, tREAL *pc)
{ 
	tREAL detleft, detright, det; 
	tREAL detsum, errbound; 
	
	detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]); 
	detright = (pa[1] - pc[1]) * (pb[0] - pc[0]); 
	det = detleft - detright; 
	
	if (detleft > 0.0)
	{ 
		if (detright <= 0.0) return det; 
		else detsum = detleft + detright;
	}
	else if (detleft < 0.0)
	{ 
		if (detright >= 0.0) return det; 
		else detsum = -detleft - detright; 
	}
	else return det;
	
	errbound = ccwerrboundA * detsum; 
	if ((det >= errbound) || (-det >= errbound))
	{ 
		return det; 
	} 
	
	return orient2dadapt(pa, pb, pc, detsum); 
} 

/***************************************************************************** 
**                                                                            
**  incirclefast()   Approximate 2D incircle test.  Nonrobust.                
**  incircleexact()   Exact 2D incircle test.  Robust.                        
**  incircleslow()   Another exact 2D incircle test.  Robust.                 
**  incircle()   Adaptive exact 2D incircle test.  Robust.                    
**                                                                            
**               Return a positive value if the point pd lies inside the      
**               circle passing through pa, pb, and pc; a negative value if   
**               it lies outside; and zero if the four points are cocircular. 
**               The points pa, pb, and pc must be in counterclockwise        
**               order, or the sign of the result will be reversed.           
**                                                                            
**  Only the first and last routine should be used; the middle two are for   
**  timings.                                                                  
**                                                                            
**  The last three use exact arithmetic to ensure a correct answer.  The      
**  result returned is the determinant of a matrix.  In incircle() only,      
**  this determinant is computed adaptively, in the sense that exact          
**  arithmetic is used only to the degree it is needed to ensure that the     
**  returned value has the correct sign.  Hence, incircle() is usually quite  
**  fast, but will run more slowly when the input points are cocircular or    
**  nearly so.                                                                
**                                                                            
*****************************************************************************/ 

tREAL Predicates::incirclefast(tREAL *pa, tREAL *pb, tREAL *pc, tREAL *pd) 
{ 
	tREAL adx, ady, bdx, bdy, cdx, cdy; 
	tREAL abdet, bcdet, cadet; 
	tREAL alift, blift, clift; 
	
	adx = pa[0] - pd[0]; 
	ady = pa[1] - pd[1]; 
	bdx = pb[0] - pd[0]; 
	bdy = pb[1] - pd[1]; 
	cdx = pc[0] - pd[0]; 
	cdy = pc[1] - pd[1]; 
	
	abdet = adx * bdy - bdx * ady; 
	bcdet = bdx * cdy - cdx * bdy; 
	cadet = cdx * ady - adx * cdy; 
	alift = adx * adx + ady * ady; 
	blift = bdx * bdx + bdy * bdy; 
	clift = cdx * cdx + cdy * cdy; 
	
	return alift * bcdet + blift * cadet + clift * abdet; 
} 


tREAL Predicates::incircleadapt(tREAL *pa, tREAL *pb, tREAL *pc, tREAL *pd,
							   tREAL permanent) 
{ 
	INEXACT tREAL adx, bdx, cdx, ady, bdy, cdy; 
	tREAL det, errbound; 
	
	INEXACT tREAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1; 
	tREAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0; 
	tREAL bc[4], ca[4], ab[4]; 
	INEXACT tREAL bc3, ca3, ab3; 
	tREAL axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32]; 
	int axbclen, axxbclen, aybclen, ayybclen, alen; 
	tREAL bxca[8], bxxca[16], byca[8], byyca[16], bdet[32]; 
	int bxcalen, bxxcalen, bycalen, byycalen, blen; 
	tREAL cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32]; 
	int cxablen, cxxablen, cyablen, cyyablen, clen; 
	tREAL abdet[64]; 
	int ablen; 
	tREAL fin1[1152], fin2[1152]; 
	tREAL *finnow, *finother, *finswap; 
	int finlength; 
	
	tREAL adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail; 
	INEXACT tREAL adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1; 
	tREAL adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0; 
	tREAL aa[4], bb[4], cc[4]; 
	INEXACT tREAL aa3, bb3, cc3; 
	INEXACT tREAL ti1, tj1; 
	tREAL ti0, tj0; 
	tREAL u[4], v[4]; 
	INEXACT tREAL u3, v3; 
	tREAL temp8[8], temp16a[16], temp16b[16], temp16c[16]; 
	tREAL temp32a[32], temp32b[32], temp48[48], temp64[64]; 
	int temp8len, temp16alen, temp16blen, temp16clen; 
	int temp32alen, temp32blen, temp48len, temp64len; 
	tREAL axtbb[8], axtcc[8], aytbb[8], aytcc[8]; 
	int axtbblen, axtcclen, aytbblen, aytcclen; 
	tREAL bxtaa[8], bxtcc[8], bytaa[8], bytcc[8]; 
	int bxtaalen, bxtcclen, bytaalen, bytcclen; 
	tREAL cxtaa[8], cxtbb[8], cytaa[8], cytbb[8]; 
	int cxtaalen, cxtbblen, cytaalen, cytbblen; 
	tREAL axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8]; 
	int axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen; 
	tREAL axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16], cytabt[16]; 
	int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen; 
	tREAL axtbctt[8], aytbctt[8], bxtcatt[8]; 
	tREAL bytcatt[8], cxtabtt[8], cytabtt[8]; 
	int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen; 
	tREAL abt[8], bct[8], cat[8]; 
	int abtlen, bctlen, catlen; 
	tREAL abtt[4], bctt[4], catt[4]; 
	int abttlen, bcttlen, cattlen; 
	INEXACT tREAL abtt3, bctt3, catt3; 
	tREAL negate; 
	
	INEXACT tREAL bvirt; 
	tREAL avirt, bround, around; 
	INEXACT tREAL c; 
	INEXACT tREAL abig; 
	tREAL ahi, alo, bhi, blo; 
	tREAL err1, err2, err3; 
	INEXACT tREAL _i, _j; 
	tREAL _0; 
	
	adx = (tREAL) (pa[0] - pd[0]); 
	bdx = (tREAL) (pb[0] - pd[0]); 
	cdx = (tREAL) (pc[0] - pd[0]); 
	ady = (tREAL) (pa[1] - pd[1]); 
	bdy = (tREAL) (pb[1] - pd[1]); 
	cdy = (tREAL) (pc[1] - pd[1]); 
	
	Two_Product(bdx, cdy, bdxcdy1, bdxcdy0); 
	Two_Product(cdx, bdy, cdxbdy1, cdxbdy0); 
	Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]); 
	bc[3] = bc3; 
	axbclen = scale_expansion_zeroelim(4, bc, adx, axbc); 
	axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, axxbc); 
	aybclen = scale_expansion_zeroelim(4, bc, ady, aybc); 
	ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, ayybc); 
	alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, adet); 
	
	Two_Product(cdx, ady, cdxady1, cdxady0); 
	Two_Product(adx, cdy, adxcdy1, adxcdy0); 
	Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]); 
	ca[3] = ca3; 
	bxcalen = scale_expansion_zeroelim(4, ca, bdx, bxca); 
	bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, bxxca); 
	bycalen = scale_expansion_zeroelim(4, ca, bdy, byca); 
	byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, byyca); 
	blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, bdet); 
	
	Two_Product(adx, bdy, adxbdy1, adxbdy0); 
	Two_Product(bdx, ady, bdxady1, bdxady0); 
	Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]); 
	ab[3] = ab3; 
	cxablen = scale_expansion_zeroelim(4, ab, cdx, cxab); 
	cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, cxxab); 
	cyablen = scale_expansion_zeroelim(4, ab, cdy, cyab); 
	cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, cyyab); 
	clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet); 
	
	ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet); 
	finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1); 
	
	det = estimate(finlength, fin1); 
	errbound = iccerrboundB * permanent; 
	if ((det >= errbound) || (-det >= errbound)) { 
		return det; 
	} 
	
	Two_Diff_Tail(pa[0], pd[0], adx, adxtail); 
	Two_Diff_Tail(pa[1], pd[1], ady, adytail); 
	Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail); 
	Two_Diff_Tail(pb[1], pd[1], bdy, bdytail); 
	Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail); 
	Two_Diff_Tail(pc[1], pd[1], cdy, cdytail); 
	if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0) 
		&& (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)) { 
		return det; 
	} 
	
	errbound = iccerrboundC * permanent + resulterrbound * Absolute(det); 
	det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) 
									   - (bdy * cdxtail + cdx * bdytail)) 
			+ 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) 
		+ ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) 
									  - (cdy * adxtail + adx * cdytail)) 
		   + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) 
		+ ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) 
									  - (ady * bdxtail + bdx * adytail)) 
		   + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx)); 
	if ((det >= errbound) || (-det >= errbound)) { 
		return det; 
	} 
	
	finnow = fin1; 
	finother = fin2; 
	
	if ((bdxtail != 0.0) || (bdytail != 0.0) 
		|| (cdxtail != 0.0) || (cdytail != 0.0)) { 
		Square(adx, adxadx1, adxadx0); 
		Square(ady, adyady1, adyady0); 
		Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]); 
		aa[3] = aa3; 
	} 
	if ((cdxtail != 0.0) || (cdytail != 0.0) 
		|| (adxtail != 0.0) || (adytail != 0.0)) { 
		Square(bdx, bdxbdx1, bdxbdx0); 
		Square(bdy, bdybdy1, bdybdy0); 
		Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]); 
		bb[3] = bb3; 
	} 
	if ((adxtail != 0.0) || (adytail != 0.0) 
		|| (bdxtail != 0.0) || (bdytail != 0.0)) { 
		Square(cdx, cdxcdx1, cdxcdx0); 
		Square(cdy, cdycdy1, cdycdy0); 
		Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]); 
		cc[3] = cc3; 
	} 
	
	if (adxtail != 0.0) { 
		axtbclen = scale_expansion_zeroelim(4, bc, adxtail, axtbc); 
		temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2.0 * adx, 
											  temp16a); 
		
		axtcclen = scale_expansion_zeroelim(4, cc, adxtail, axtcc); 
		temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, temp16b); 
		
		axtbblen = scale_expansion_zeroelim(4, bb, adxtail, axtbb); 
		temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, temp16c); 
		
		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
												 temp16blen, temp16b, temp32a); 
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
												temp32alen, temp32a, temp48); 
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
												temp48, finother); 
		finswap = finnow; finnow = finother; finother = finswap; 
	} 
	if (adytail != 0.0) { 
		aytbclen = scale_expansion_zeroelim(4, bc, adytail, aytbc); 
		temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2.0 * ady, 
											  temp16a); 
		
		aytbblen = scale_expansion_zeroelim(4, bb, adytail, aytbb); 
		temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, temp16b); 
		
		aytcclen = scale_expansion_zeroelim(4, cc, adytail, aytcc); 
		temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, temp16c); 
		
		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
												 temp16blen, temp16b, temp32a); 
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
												temp32alen, temp32a, temp48); 
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
												temp48, finother); 
		finswap = finnow; finnow = finother; finother = finswap; 
	} 
	if (bdxtail != 0.0) { 
		bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, bxtca); 
		temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2.0 * bdx, 
											  temp16a); 
		
		bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, bxtaa); 
		temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, temp16b); 
		
		bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, bxtcc); 
		temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, temp16c); 
		
		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
												 temp16blen, temp16b, temp32a); 
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
												temp32alen, temp32a, temp48); 
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
												temp48, finother); 
		finswap = finnow; finnow = finother; finother = finswap; 
	} 
	if (bdytail != 0.0) { 
		bytcalen = scale_expansion_zeroelim(4, ca, bdytail, bytca); 
		temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2.0 * bdy, 
											  temp16a); 
		
		bytcclen = scale_expansion_zeroelim(4, cc, bdytail, bytcc); 
		temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, temp16b); 
		
		bytaalen = scale_expansion_zeroelim(4, aa, bdytail, bytaa); 
		temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, temp16c); 
		
		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
												 temp16blen, temp16b, temp32a); 
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
												temp32alen, temp32a, temp48); 
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
												temp48, finother); 
		finswap = finnow; finnow = finother; finother = finswap; 
	} 
	if (cdxtail != 0.0) { 
		cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, cxtab); 
		temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2.0 * cdx, 
											  temp16a); 
		
		cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, cxtbb); 
		temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, temp16b); 
		
		cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, cxtaa); 
		temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, temp16c); 
		
		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
												 temp16blen, temp16b, temp32a); 
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
												temp32alen, temp32a, temp48); 
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
												temp48, finother); 
		finswap = finnow; finnow = finother; finother = finswap; 
	} 
	if (cdytail != 0.0) { 
		cytablen = scale_expansion_zeroelim(4, ab, cdytail, cytab); 
		temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2.0 * cdy, 
											  temp16a); 
		
		cytaalen = scale_expansion_zeroelim(4, aa, cdytail, cytaa); 
		temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, temp16b); 
		
		cytbblen = scale_expansion_zeroelim(4, bb, cdytail, cytbb); 
		temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, temp16c); 
		
		temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
												 temp16blen, temp16b, temp32a); 
		temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
												temp32alen, temp32a, temp48); 
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
												temp48, finother); 
		finswap = finnow; finnow = finother; finother = finswap; 
	} 
	
	if ((adxtail != 0.0) || (adytail != 0.0)) { 
		if ((bdxtail != 0.0) || (bdytail != 0.0) 
			|| (cdxtail != 0.0) || (cdytail != 0.0)) { 
			Two_Product(bdxtail, cdy, ti1, ti0); 
			Two_Product(bdx, cdytail, tj1, tj0); 
			Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]); 
			u[3] = u3; 
			negate = -bdy; 
			Two_Product(cdxtail, negate, ti1, ti0); 
			negate = -bdytail; 
			Two_Product(cdx, negate, tj1, tj0); 
			Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]); 
			v[3] = v3; 
			bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, bct); 
			
			Two_Product(bdxtail, cdytail, ti1, ti0); 
			Two_Product(cdxtail, bdytail, tj1, tj0); 
			Two_Two_Diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]); 
			bctt[3] = bctt3; 
			bcttlen = 4; 
		} else { 
			bct[0] = 0.0; 
			bctlen = 1; 
			bctt[0] = 0.0; 
			bcttlen = 1; 
		} 
		
		if (adxtail != 0.0) { 
			temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, temp16a); 
			axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, axtbct); 
			temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2.0 * adx, 
												  temp32a); 
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													temp32alen, temp32a, temp48); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
													temp48, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
			if (bdytail != 0.0) { 
				temp8len = scale_expansion_zeroelim(4, cc, adxtail, temp8); 
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, 
													  temp16a); 
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
														temp16a, finother); 
				finswap = finnow; finnow = finother; finother = finswap; 
			} 
			if (cdytail != 0.0) { 
				temp8len = scale_expansion_zeroelim(4, bb, -adxtail, temp8); 
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, 
													  temp16a); 
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
														temp16a, finother); 
				finswap = finnow; finnow = finother; finother = finswap; 
			} 
			
			temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail, 
												  temp32a); 
			axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, axtbctt); 
			temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2.0 * adx, 
												  temp16a); 
			temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail, 
												  temp16b); 
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													 temp16blen, temp16b, temp32b); 
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
													temp32blen, temp32b, temp64); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
													temp64, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
		} 
		if (adytail != 0.0) { 
			temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, temp16a); 
			aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, aytbct); 
			temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2.0 * ady, 
												  temp32a); 
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													temp32alen, temp32a, temp48); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
													temp48, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
			
			
			temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail, 
												  temp32a); 
			aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, aytbctt); 
			temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2.0 * ady, 
												  temp16a); 
			temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail, 
												  temp16b); 
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													 temp16blen, temp16b, temp32b); 
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
													temp32blen, temp32b, temp64); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
													temp64, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
		} 
	} 
	if ((bdxtail != 0.0) || (bdytail != 0.0)) { 
		if ((cdxtail != 0.0) || (cdytail != 0.0) 
			|| (adxtail != 0.0) || (adytail != 0.0)) { 
			Two_Product(cdxtail, ady, ti1, ti0); 
			Two_Product(cdx, adytail, tj1, tj0); 
			Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]); 
			u[3] = u3; 
			negate = -cdy; 
			Two_Product(adxtail, negate, ti1, ti0); 
			negate = -cdytail; 
			Two_Product(adx, negate, tj1, tj0); 
			Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]); 
			v[3] = v3; 
			catlen = fast_expansion_sum_zeroelim(4, u, 4, v, cat); 
			
			Two_Product(cdxtail, adytail, ti1, ti0); 
			Two_Product(adxtail, cdytail, tj1, tj0); 
			Two_Two_Diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]); 
			catt[3] = catt3; 
			cattlen = 4; 
		} else { 
			cat[0] = 0.0; 
			catlen = 1; 
			catt[0] = 0.0; 
			cattlen = 1; 
		} 
		
		if (bdxtail != 0.0) { 
			temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, temp16a); 
			bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, bxtcat); 
			temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2.0 * bdx, 
												  temp32a); 
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													temp32alen, temp32a, temp48); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
													temp48, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
			if (cdytail != 0.0) { 
				temp8len = scale_expansion_zeroelim(4, aa, bdxtail, temp8); 
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, 
													  temp16a); 
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
														temp16a, finother); 
				finswap = finnow; finnow = finother; finother = finswap; 
			} 
			if (adytail != 0.0) { 
				temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, temp8); 
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, 
													  temp16a); 
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
														temp16a, finother); 
				finswap = finnow; finnow = finother; finother = finswap; 
			} 
			
			temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail, 
												  temp32a); 
			bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, bxtcatt); 
			temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2.0 * bdx, 
												  temp16a); 
			temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail, 
												  temp16b); 
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													 temp16blen, temp16b, temp32b); 
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
													temp32blen, temp32b, temp64); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
													temp64, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
		} 
		if (bdytail != 0.0) { 
			temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, temp16a); 
			bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, bytcat); 
			temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2.0 * bdy, 
												  temp32a); 
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													temp32alen, temp32a, temp48); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
													temp48, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
			
			
			temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail, 
												  temp32a); 
			bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, bytcatt); 
			temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2.0 * bdy, 
												  temp16a); 
			temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail, 
												  temp16b); 
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													 temp16blen, temp16b, temp32b); 
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
													temp32blen, temp32b, temp64); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
													temp64, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
		} 
	} 
	if ((cdxtail != 0.0) || (cdytail != 0.0)) { 
		if ((adxtail != 0.0) || (adytail != 0.0) 
			|| (bdxtail != 0.0) || (bdytail != 0.0)) { 
			Two_Product(adxtail, bdy, ti1, ti0); 
			Two_Product(adx, bdytail, tj1, tj0); 
			Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]); 
			u[3] = u3; 
			negate = -ady; 
			Two_Product(bdxtail, negate, ti1, ti0); 
			negate = -adytail; 
			Two_Product(bdx, negate, tj1, tj0); 
			Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]); 
			v[3] = v3; 
			abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, abt); 
			
			Two_Product(adxtail, bdytail, ti1, ti0); 
			Two_Product(bdxtail, adytail, tj1, tj0); 
			Two_Two_Diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]); 
			abtt[3] = abtt3; 
			abttlen = 4; 
		} else { 
			abt[0] = 0.0; 
			abtlen = 1; 
			abtt[0] = 0.0; 
			abttlen = 1; 
		} 
		
		if (cdxtail != 0.0) { 
			temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, temp16a); 
			cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, cxtabt); 
			temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2.0 * cdx, 
												  temp32a); 
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													temp32alen, temp32a, temp48); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
													temp48, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
			if (adytail != 0.0) { 
				temp8len = scale_expansion_zeroelim(4, bb, cdxtail, temp8); 
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, 
													  temp16a); 
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
														temp16a, finother); 
				finswap = finnow; finnow = finother; finother = finswap; 
			} 
			if (bdytail != 0.0) { 
				temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, temp8); 
				temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, 
													  temp16a); 
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
														temp16a, finother); 
				finswap = finnow; finnow = finother; finother = finswap; 
			} 
			
			temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail, 
												  temp32a); 
			cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, cxtabtt); 
			temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2.0 * cdx, 
												  temp16a); 
			temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail, 
												  temp16b); 
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													 temp16blen, temp16b, temp32b); 
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
													temp32blen, temp32b, temp64); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
													temp64, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
		} 
		if (cdytail != 0.0) { 
			temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, temp16a); 
			cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, cytabt); 
			temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2.0 * cdy, 
												  temp32a); 
			temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													temp32alen, temp32a, temp48); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
													temp48, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
			
			
			temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail, 
												  temp32a); 
			cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, cytabtt); 
			temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2.0 * cdy, 
												  temp16a); 
			temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail, 
												  temp16b); 
			temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
													 temp16blen, temp16b, temp32b); 
			temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
													temp32blen, temp32b, temp64); 
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
													temp64, finother); 
			finswap = finnow; finnow = finother; finother = finswap; 
		} 
	} 
	
	return finnow[finlength - 1]; 
} 

tREAL Predicates::incircle(tREAL *pa, tREAL *pb, tREAL *pc, tREAL *pd)
{ 
	tREAL adx, bdx, cdx, ady, bdy, cdy; 
	tREAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady; 
	tREAL alift, blift, clift; 
	tREAL det; 
	tREAL permanent, errbound; 
	
	adx = pa[0] - pd[0]; 
	bdx = pb[0] - pd[0]; 
	cdx = pc[0] - pd[0]; 
	ady = pa[1] - pd[1]; 
	bdy = pb[1] - pd[1]; 
	cdy = pc[1] - pd[1]; 
	
	bdxcdy = bdx * cdy; 
	cdxbdy = cdx * bdy; 
	alift = adx * adx + ady * ady; 
	
	cdxady = cdx * ady; 
	adxcdy = adx * cdy; 
	blift = bdx * bdx + bdy * bdy; 
	
	adxbdy = adx * bdy; 
	bdxady = bdx * ady; 
	clift = cdx * cdx + cdy * cdy; 
	
	det = alift * (bdxcdy - cdxbdy) 
		+ blift * (cdxady - adxcdy) 
		+ clift * (adxbdy - bdxady); 
	
	permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift 
		+ (Absolute(cdxady) + Absolute(adxcdy)) * blift 
		+ (Absolute(adxbdy) + Absolute(bdxady)) * clift; 
	errbound = iccerrboundA * permanent; 
	if ((det > errbound) || (-det > errbound)) { 
		return det; 
	} 
	
	return incircleadapt(pa, pb, pc, pd, permanent); 
} 


//=========================================================================
//
//
//                    End of predicates.cpp
//
//
//=========================================================================
