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
**  mathutil.cpp: Functions from mathutil and complex1 (See mathutil.h). 
**
***************************************************************************/

#include "src/Mathutil/mathutil.h"

//=========================================================================
//
//
//                  Section 1: mathutil Functions
//
//
//=========================================================================

/**************************************************************************
**
**  ran3 (long *)
**
**  Random number generator from Numerical Recipes. Returns a uniform 
**  random number between 0.0 and 1.0. Set idum to any negative value to 
**  initialize or reinitialize the sequence. Parameters: idum - random seed
**
**  
**************************************************************************/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
	
	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
				inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/***************************************************************************
**
**  Random number generator
**
***************************************************************************/
int rand1_00()
{
	const int a = 16807;
	const int m = 2147483647;  
	const int q = 127773;
	const int r = 2836;
	
	int low, high, test;
	
	if (rand1_00_seed == 0) rand1_00_seed = 1;
	high = rand1_00_seed/q;
	low = rand1_00_seed%q;
	test = a*low - r*high;
	if (test > 0)
		rand1_00_seed = test;
	else
		rand1_00_seed = test + m;
	
	return rand1_00_seed;
}
int getRandomSeed() { return rand1_00_seed; }
void setRandomSeed(int seed) { rand1_00_seed = seed; }

/***************************************************************************
**
**  Generates random number between 0 and 1 using uniform distribution
**
***************************************************************************/
double uniform()
{
	const int MAX_RAND = 2147483647;
	return(static_cast<double>(rand1_00())/MAX_RAND);
}

/***************************************************************************
**
**  Generates random number between 'a' and 'b' using uniform distribution
**
***************************************************************************/
double uniform(double a, double b)
{
	const int MAX_RAND = 2147483647;
	return(a+(b-a)*static_cast<double>(rand1_00())/MAX_RAND);
}

/***************************************************************************
**
**  Generates a random value from exponential distribution
**
***************************************************************************/
double random_expon(double lambda)
{
	const int MAX_RAND = 2147483647;
	double u;
	
	u = static_cast<double>(rand1_00())/MAX_RAND;
	return (-log(u)/lambda); 
}

/***************************************************************************
**
**  Gives a cdf value of exponential distribution [ P(X <= a) ]
**
***************************************************************************/
double cdfexpon(double lambda, double a)
{
	return(1-exp(-lambda*a));
}

/***************************************************************************
**
**  Generates a random value from exponential distribution N(esp, std^2)
**  following Naylor et. al 1966 
**
***************************************************************************/
double random_normal(double esp, double std)
{
	double u1=uniform();
	double u2=uniform();
	double z=sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
	return esp+std*z;
}

/**************************************************************************
**
**  Generates a random variable drawn from a Gamma distribution with par-r m
**
**************************************************************************/
double random_gamma(double m, long * idum)
{
	double x, y,z, c,t,b,u,w,v;
	
	if (m<1){
		c = 1.0/m;
		t = 0.07 + 0.75*sqrt(1.0-m);
		b = 1.0 + exp(-t)*m/t;
		int accept = 0;
		while (accept == 0){
			u = ran3(idum);
			w = ran3(idum);
			v = b *u;
			if (v<=1.0)
			{
				x  = t * pow(v, c);
				accept = ((w<=((2.0-x)/(2.0+x))) || (w<=exp(-x)));
			}
			else
			{
				x = -log(c*t*(b-v));
				y = x/t;
				accept = (((w*(m + y - m*y)) <= 1.0) || (w<= ( pow(y, (m-1.0)))));
			}
		}
	}
	else{
		b = m-1.0;
		c = 3.0*m - 0.75;
		int accept = 0;
		while (accept == 0)
		{
			u = ran3(idum); v = ran3(idum);
			w = u* (1.0-u);
			y = sqrt(c/w) * (u - 0.5);
			x = b + y;
			if (x>= 0.0){
				z = 64*( pow(w,3.0))*v*v;
				accept = (z <= ( 1.0 - 2.0*y*y/x)) || ( log(z) <= (2.0*(b*log(x/b) - y)));
			}
		}
	}
	return x;
}

/***************************************************************************
**
**  Generates a random value from Weibull distribution w/ par-s alpha & beta
**  where 'beta' is the distribution MEAN and 'alpha' is the SHAPE FACTOR
**  See Simulation, Modelling & Analysis by Law & Kelton, pp259
**  This is the ``polar'' method.
**
***************************************************************************/
double random_Weibull(double pBeta, double pAlpha)
{
	return( pow(pBeta * ( -log(1.0 - uniform()) ), 1.0/pAlpha) );
}

/***************************************************************************
**
**  Generates a random value from Poisson distribution
**
***************************************************************************/
unsigned long random_poisson(double landa)
{
	double p=exp(-landa);
	double g=p;
	double u=uniform();
	unsigned long k=0;
	while (u>g) {
		p*=(landa/(double)(++k));
		g+=p;
	}
	return k;
}

/***************************************************************************
**
**  Generates a random value from binomial distribution
**
***************************************************************************/
unsigned long random_binomial(unsigned long n, double p)
{
	double t=p/(1.0-p);
	double u=uniform();
	double p0=pow((1.0-p),double(n));
	double g=p0;
	unsigned int k=0;
	while (u>g) {
		p0*=t*(n-k)/(k+1.0);
		g+=p0;
		k++;
	}
	return k;
}

/***************************************************************************
**
**  Generates a random value from hipergeometric distribution
**
***************************************************************************/
long random_hipergeom(unsigned long N,unsigned long D,
					  unsigned long n,double p=0.0)
{
	double u=uniform();
	long k=n-N+D;
	if (k<0) k=0;
	if (p==0.0) p=hipergeom_p0(N,D,n);
	double g=p;
	while (u>g) {
		p*=((double)(D-k)*(n-k))/(double)((k+1.0)*(N-D-n+k+1.0));
		g+=p;
		k++;
	}
	return k;
}

double hipergeom_p0(unsigned long N,unsigned long D,unsigned long n)
{
	long k=n-N+D;
	if (k<0) k=0;
	return (double)((fact(D,D-k)*fact(N-D,N-D-n+k)*fact(n,n-k))/(fact(k,1)*fact(N,N-n)));
}

long double fact(unsigned long x,unsigned long f)
{
	if (x==f) return 1.0;
	if (x<f)  return 1.0/fact(f,x);
	long double t=1.0;
	while (x>f) {
		t*=x;
		x--;
	}
	return (long double)t;
}


#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
float genbet(float aa,float bb)

/*
**********************************************************************
 float genbet(float aa,float bb)
 GeNerate BETa random deviate
 Function
 Returns a single random deviate from the beta distribution with
 parameters A and B.  The density of the beta is
 x^(a-1) * (1-x)^(b-1) / B(a,b) for 0 < x < 1
 Arguments
 aa --> First parameter of the beta distribution
 
 bb --> Second parameter of the beta distribution
 
 Method
 R. C. H. Cheng
 Generating Beta Variatew with Nonintegral Shape Parameters
 Communications of the ACM, 21:317-322  (1978)
 (Algorithms BB and BC)
 **********************************************************************
 */
{
	
#define expmax 87.49823
#define infnty 1.0E38
#define minlog 1.0E-37
	static float olda = -1.0E37;
	static float oldb = -1.0E37;
	static float genbet,a,alpha,b,beta,delta,gamma,k1,k2,r,s,t,u1,u2,v,w,y,z;
	static long qsame;
	
    qsame = olda == aa && oldb == bb;
    if(qsame) goto S20;
    if(!(aa < minlog || bb < minlog)) goto S10;
    fputs(" AA or BB < 1.0E-37 in GENBET - Abort!\n",stderr);
    fprintf(stderr," AA: %16.6E BB %16.6E\n",aa,bb);
    exit(1);
S10:
		olda = aa;
    oldb = bb;
S20:
		if(!(min(aa,bb) > 1.0)) goto S100;
	/*
		Alborithm BB
     Initialize
	 */
    if(qsame) goto S30;
    a = min(aa,bb);
    b = max(aa,bb);
    alpha = a+b;
    beta = sqrt((alpha-2.0)/(2.0*a*b-alpha));
    gamma = a+1.0/beta;
S30:
S40:
		u1 = ranf();
	/*
		Step 1
	 */
    u2 = ranf();
    v = beta*log(u1/(1.0-u1));
	/* JJV altered this */
    if(v > expmax) goto S55;
	/*
		* JJV added checker to see if a*exp(v) will overflow
	 * JJV S50 _was_ w = a*exp(v); also note here a > 1.0
	 */
    w = exp(v);
    if(w > infnty/a) goto S55;
    w *= a;
    goto S60;
S55:
		w = infnty;
S60:
		z = pow(double(u1),double(2.0))*u2;
    r = gamma*v-1.3862944;
    s = a+r-w;
	/*
		Step 2
	 */
    if(s+2.609438 >= 5.0*z) goto S70;
	/*
		Step 3
	 */
    t = log(z);
    if(s > t) goto S70;
	/*
		*   Step 4
	 *
	 *    JJV added checker to see if log(alpha/(b+w)) will 
	 *    JJV overflow.  If so, we count the log as -INF, and
	 *    JJV consequently evaluate conditional as true, i.e.
	 *    JJV the algorithm rejects the trial and starts over
	 *    JJV May not need this here since alpha > 2.0
	 */
    if(alpha/(b+w) < minlog) goto S40;
    if(r+alpha*log(alpha/(b+w)) < t) goto S40;
S70:
		/*
		Step 5
		 */
		if(!(aa == a)) goto S80;
    genbet = w/(b+w);
    goto S90;
S80:
		genbet = b/(b+w);
S90:
		goto S230;
S100:
		/*
		Algorithm BC
		 Initialize
		 */
		if(qsame) goto S110;
    a = max(aa,bb);
    b = min(aa,bb);
    alpha = a+b;
    beta = 1.0/b;
    delta = 1.0+a-b;
    k1 = delta*(1.38889E-2+4.16667E-2*b)/(a*beta-0.777778);
    k2 = 0.25+(0.5+0.25/delta)*b;
S110:
S120:
		u1 = ranf();
	/*
		Step 1
	 */
    u2 = ranf();
    if(u1 >= 0.5) goto S130;
	/*
		Step 2
	 */
    y = u1*u2;
    z = u1*y;
    if(0.25*u2+z-y >= k1) goto S120;
    goto S170;
S130:
		/*
		Step 3
		 */
		z = pow(double(u1),double(2.0))*u2; 
    if(!(z <= 0.25)) goto S160;
    v = beta*log(u1/(1.0-u1));
	/*
		*    JJV instead of checking v > expmax at top, I will check
	 *    JJV if a < 1, then check the appropriate values
	 */
    if(a > 1.0) goto S135;
	/*   JJV a < 1 so it can help out if exp(v) would overflow */
    if(v > expmax) goto S132;
    w = a*exp(v);
    goto S200;
S132:
		w = v + log(a);
    if(w > expmax) goto S140;
    w = exp(w);
    goto S200;
S135:
		/*   JJV in this case a > 1 */
		if(v > expmax) goto S140;
    w = exp(v);
    if(w > infnty/a) goto S140;
    w *= a;
    goto S200;
S140:
		w = infnty;
    goto S200;
	/*
		* JJV old code
	 *    if(!(v > expmax)) goto S140;
	 *    w = infnty;
	 *    goto S150;
	 *S140:
	 *    w = a*exp(v);
	 *S150:
	 *    goto S200;
	 */
S160:
		if(z >= k2) goto S120;
S170:
		/*
		Step 4
		 Step 5
		 */
		v = beta*log(u1/(1.0-u1));
	/*   JJV same kind of checking as above */
    if(a > 1.0) goto S175;
	/* JJV a < 1 so it can help out if exp(v) would overflow */
    if(v > expmax) goto S172;
    w = a*exp(v);
    goto S190;
S172:
		w = v + log(a);
    if(w > expmax) goto S180;
    w = exp(w);
    goto S190;
S175:
		/* JJV in this case a > 1.0 */
		if(v > expmax) goto S180;
    w = exp(v);
    if(w > infnty/a) goto S180;
    w *= a;
    goto S190;
S180:
		w = infnty;
	/*
		*   JJV old code
	 *    if(!(v > expmax)) goto S180;
	 *    w = infnty;
	 *    goto S190;
	 *S180:
	 *    w = a*exp(v);
	 */
S190:
		/*
		* JJV here we also check to see if log overlows; if so, we treat it
		 * JJV as -INF, which means condition is true, i.e. restart
		 */
		if(alpha/(b+w) < minlog) goto S120;
    if(alpha*(log(alpha/(b+w))+v)-1.3862944 < log(z)) goto S120;
S200:
		/*
		Step 6
		 */
		if(!(a == aa)) goto S210;
    genbet = w/(b+w);
    goto S220;
S210:
		genbet = b/(b+w);
S230:
S220:
		return genbet;
#undef expmax
#undef infnty
#undef minlog
}

float ranf(void)
/*
**********************************************************************
 float ranf(void)
 RANDom number generator as a Function
 Returns a random floating point number from a uniform distribution
 over 0 - 1 (endpoints of this interval are not returned) using the
 current generator
 This is a transcription from Pascal to Fortran of routine
 Uniform_01 from the paper
 L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
 with Splitting Facilities." ACM Transactions on Mathematical
 Software, 17:98-111 (1991)
 **********************************************************************
 */
{
	static float ranf;
	/*
	4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
	 and is currently 2147483563. If M1 changes, change this also.
	 */
    ranf = ignlgi()*4.656613057E-10;
    return ranf;
}

long ignlgi(void)
/*
**********************************************************************
 long ignlgi(void)
 GeNerate LarGe Integer
 Returns a random integer following a uniform distribution over
 (1, 2147483562) using the current generator.
 This is a transcription from Pascal to Fortran of routine
 Random from the paper
 L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
 with Splitting Facilities." ACM Transactions on Mathematical
 Software, 17:98-111 (1991)
 **********************************************************************
 */
{
#define numg 32L
	extern void gsrgs(long getset,long *qvalue);
	extern void gssst(long getset,long *qset);
	extern void gscgn(long getset,long *g);
	extern void inrgcm(void);
	extern long Xm1,Xm2,Xa1,Xa2,Xcg1[],Xcg2[];
	extern long Xqanti[];
	static long ignlgi,curntg,k,s1,s2,z;
	static long qqssd,qrgnin;
	/*
		IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
     THIS ROUTINE  2) A CALL TO SETALL.
	 */
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    gssst(0,&qqssd);
    if(!qqssd) setall(1234567890L,123456789L);
	/*
		Get Current Generator
	 */
    gscgn(0L,&curntg);
    s1 = *(Xcg1+curntg-1);
    s2 = *(Xcg2+curntg-1);
    k = s1/53668L;
    s1 = Xa1*(s1-k*53668L)-k*12211;
    if(s1 < 0) s1 += Xm1;
    k = s2/52774L;
    s2 = Xa2*(s2-k*52774L)-k*3791;
    if(s2 < 0) s2 += Xm2;
    *(Xcg1+curntg-1) = s1;
    *(Xcg2+curntg-1) = s2;
    z = s1-s2;
    if(z < 1) z += (Xm1-1);
    if(*(Xqanti+curntg-1)) z = Xm1-z;
    ignlgi = z;
    return ignlgi;
#undef numg
}

void setall(long iseed1,long iseed2)
/*
**********************************************************************
 void setall(long iseed1,long iseed2)
 SET ALL random number generators
 Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
 initial seeds of the other generators are set accordingly, and
 all generators states are set to these seeds.
 This is a transcription from Pascal to Fortran of routine
 Set_Initial_Seed from the paper
 L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
 with Splitting Facilities." ACM Transactions on Mathematical
 Software, 17:98-111 (1991)
 Arguments
 iseed1 -> First of two integer seeds
 iseed2 -> Second of two integer seeds
 **********************************************************************
 */
{
#define numg 32L
	extern void gsrgs(long getset,long *qvalue);
	extern void gssst(long getset,long *qset);
	extern void gscgn(long getset,long *g);
	extern long Xm1,Xm2,Xa1vw,Xa2vw,Xig1[],Xig2[];
	static long T1;
	static long g,ocgn;
	static long qrgnin;
    T1 = 1;
	/*
		TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
	 HAS BEEN CALLED.
	 */
    gssst(1,&T1);
    gscgn(0L,&ocgn);
	/*
		Initialize Common Block if Necessary
	 */
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    *Xig1 = iseed1;
    *Xig2 = iseed2;
    initgn(-1L);
    for(g=2; g<=numg; g++) {
        *(Xig1+g-1) = mltmod(Xa1vw,*(Xig1+g-2),Xm1);
        *(Xig2+g-1) = mltmod(Xa2vw,*(Xig2+g-2),Xm2);
        gscgn(1L,&g);
        initgn(-1L);
    }
    gscgn(1L,&ocgn);
#undef numg
}

void inrgcm(void)
/*
**********************************************************************
 void inrgcm(void)
 INitialize Random number Generator CoMmon
 Function
 Initializes common area  for random number  generator.  This saves
 the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
 assuring that the routine is loaded with the other routines.
 **********************************************************************
 */
{
#define numg 32L
	extern void gsrgs(long getset,long *qvalue);
	extern long Xm1,Xm2,Xa1,Xa2,Xa1w,Xa2w,Xa1vw,Xa2vw;
	extern long Xqanti[];
	static long T1;
	static long i;
	/*
		V=20;                            W=30;
     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
	 If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
	 An efficient way to precompute a**(2*j) MOD m is to start with
	 a and square it j times modulo m using the function MLTMOD.
	 */
    Xm1 = 2147483563L;
    Xm2 = 2147483399L;
    Xa1 = 40014L;
    Xa2 = 40692L;
    Xa1w = 1033780774L;
    Xa2w = 1494757890L;
    Xa1vw = 2082007225L;
    Xa2vw = 784306273L;
    for(i=0; i<numg; i++) *(Xqanti+i) = 0;
    T1 = 1;
	/*
		Tell the world that common has been initialized
	 */
    gsrgs(1L,&T1);
#undef numg
}

void initgn(long isdtyp)
/*
**********************************************************************
 void initgn(long isdtyp)
 INIT-ialize current G-e-N-erator
 Reinitializes the state of the current generator
 This is a transcription from Pascal to Fortran of routine
 Init_Generator from the paper
 L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
 with Splitting Facilities." ACM Transactions on Mathematical
 Software, 17:98-111 (1991)
 Arguments
 isdtyp -> The state to which the generator is to be set
 isdtyp = -1  => sets the seeds to their initial value
 isdtyp =  0  => sets the seeds to the first value of
 the current block
 isdtyp =  1  => sets the seeds to the first value of
 the next block
 **********************************************************************
 */
{
#define numg 32L
	extern void gsrgs(long getset,long *qvalue);
	extern void gscgn(long getset,long *g);
	extern long Xm1,Xm2,Xa1w,Xa2w,Xig1[],Xig2[],Xlg1[],Xlg2[],Xcg1[],Xcg2[];
	static long g;
	static long qrgnin;
	/*
		Abort unless random number generator initialized
	 */
    gsrgs(0L,&qrgnin);
    if(qrgnin) goto S10;
    fprintf(stderr,"%s\n",
			" INITGN called before random number generator  initialized -- abort!");
    exit(1);
S10:
		gscgn(0L,&g);
    if(-1 != isdtyp) goto S20;
    *(Xlg1+g-1) = *(Xig1+g-1);
    *(Xlg2+g-1) = *(Xig2+g-1);
    goto S50;
S20:
		if(0 != isdtyp) goto S30;
    goto S50;
S30:
		/*
		do nothing
		 */
		if(1 != isdtyp) goto S40;
    *(Xlg1+g-1) = mltmod(Xa1w,*(Xlg1+g-1),Xm1);
    *(Xlg2+g-1) = mltmod(Xa2w,*(Xlg2+g-1),Xm2);
    goto S50;
S40:
		fprintf(stderr,"%s\n","isdtyp not in range in INITGN");
    exit(1);
S50:
		*(Xcg1+g-1) = *(Xlg1+g-1);
    *(Xcg2+g-1) = *(Xlg2+g-1);
#undef numg
}

long mltmod(long a,long s,long m)
/*
**********************************************************************
 long mltmod(long a,long s,long m)
 Returns (A*S) MOD M
 This is a transcription from Pascal to Fortran of routine
 MULtMod_Decompos from the paper
 L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
 with Splitting Facilities." ACM Transactions on Mathematical
 Software, 17:98-111 (1991)
 Arguments
 a, s, m  -->
 **********************************************************************
 */
{
#define h 32768L
	static long mltmod,a0,a1,k,p,q,qh,rh;
	/*
	H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
	 machine. On a different machine recompute H
	 */
    if(!(a <= 0 || a >= m || s <= 0 || s >= m)) goto S10;
    fputs(" a, m, s out of order in mltmod - ABORT!\n",stderr);
    fprintf(stderr," a = %12ld s = %12ld m = %12ld\n",a,s,m);
    fputs(" mltmod requires: 0 < a < m; 0 < s < m\n",stderr);
    exit(1);
S10:
		if(!(a < h)) goto S20;
    a0 = a;
    p = 0;
    goto S120;
S20:
		a1 = a/h;
    a0 = a-h*a1;
    qh = m/h;
    rh = m-h*qh;
    if(!(a1 >= h)) goto S50;
    a1 -= h;
    k = s/qh;
    p = h*(s-k*qh)-k*rh;
S30:
		if(!(p < 0)) goto S40;
    p += m;
    goto S30;
S40:
		goto S60;
S50:
		p = 0;
S60:
		/*
		P = (A2*S*H)MOD M
		 */
		if(!(a1 != 0)) goto S90;
    q = m/a1;
    k = s/q;
    p -= (k*(m-a1*q));
    if(p > 0) p -= m;
    p += (a1*(s-k*q));
S70:
		if(!(p < 0)) goto S80;
    p += m;
    goto S70;
S90:
S80:
		k = p/qh;
	/*
		P = ((A2*H + A1)*S)MOD M
	 */
    p = h*(p-k*qh)-k*rh;
S100:
		if(!(p < 0)) goto S110;
    p += m;
    goto S100;
S120:
S110:
		if(!(a0 != 0)) goto S150;
	/*
		P = ((A2*H + A1)*H*S)MOD M
	 */
    q = m/a0;
    k = s/q;
    p -= (k*(m-a0*q));
    if(p > 0) p -= m;
    p += (a0*(s-k*q));
S130:
		if(!(p < 0)) goto S140;
    p += m;
    goto S130;
S150:
S140:
		mltmod = p;
    return mltmod;
#undef h
}

void gscgn(long getset,long *g)
/*
**********************************************************************
 void gscgn(long getset,long *g)
 Get/Set GeNerator
 Gets or returns in G the number of the current generator
 Arguments
 getset --> 0 Get
 1 Set
 g <-- Number of the current random number generator (1..32)
 **********************************************************************
 */
{
#define numg 32L
	static long curntg = 1;
    if(getset == 0) *g = curntg;
    else  {
        if(*g < 0 || *g > numg) {
            fputs(" Generator number out of range in GSCGN\n",stderr);
            exit(0);
        }
        curntg = *g;
    }
#undef numg
}
void gsrgs(long getset,long *qvalue)
/*
**********************************************************************
 void gsrgs(long getset,long *qvalue)
 Get/Set Random Generators Set
 Gets or sets whether random generators set (initialized).
 Initially (data statement) state is not set
 If getset is 1 state is set to qvalue
 If getset is 0 state returned in qvalue
 **********************************************************************
 */
{
	static long qinit = 0;
	
    if(getset == 0) *qvalue = qinit;
    else qinit = *qvalue;
}
void gssst(long getset,long *qset)
/*
**********************************************************************
 void gssst(long getset,long *qset)
 Get or Set whether Seed is Set
 Initialize to Seed not Set
 If getset is 1 sets state to Seed Set
 If getset is 0 returns T in qset if Seed Set
 Else returns F in qset
 **********************************************************************
 */
{
	static long qstate = 0;
    if(getset != 0) qstate = 1;
    else  *qset = qstate;
}


float gengam(float a,float r)
/*
**********************************************************************
 float gengam(float a,float r)
 GENerates random deviates from GAMma distribution
 Function
 Generates random deviates from the gamma distribution whose
 density is
 (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
 Arguments
 a --> Location parameter of Gamma distribution
 JJV   (a > 0)
 r --> Shape parameter of Gamma distribution
 JJV   (r > 0)
 Method
 Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
 instead of SUNIF.
 For details see:
 (Case R >= 1.0)
 Ahrens, J.H. and Dieter, U.
 Generating Gamma Variates by a
 Modified Rejection Technique.
 Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
 Algorithm GD
 JJV altered following to reflect argument ranges
 (Case 0.0 < R < 1.0)
 Ahrens, J.H. and Dieter, U.
 Computer Methods for Sampling from Gamma,
 Beta, Poisson and Binomial Distributions.
 Computing, 12 (1974), 223-246/
 Adapted algorithm GS.
 **********************************************************************
 */
{
	static float gengam;
	/* JJV added argument checker */
    if(a > 0.0 && r > 0.0) goto S10;
    fputs(" A or R nonpositive in GENGAM - abort!\n",stderr);
    fprintf(stderr," A value: %16.6E R value: %16.6E\n",a,r);
    exit(1);
S10:
		gengam = sgamma(r);
    gengam /= a;
    return gengam;
}



float sgamma(float a)
/*
**********************************************************************
 
 
 (STANDARD-)  G A M M A  DISTRIBUTION                             
 
 
 **********************************************************************
 **********************************************************************
 
 PARAMETER  A >= 1.0  !                                 
 
 **********************************************************************
 
 FOR DETAILS SEE:                                                 
 
 AHRENS, J.H. AND DIETER, U.                            
 GENERATING GAMMA VARIATES BY A                         
 MODIFIED REJECTION TECHNIQUE.                          
 COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  
 
 STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     
 (STRAIGHTFORWARD IMPLEMENTATION)     
 
 Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
 SUNIF.  The argument IR thus goes away.                          
 
 **********************************************************************
 
 PARAMETER  0.0 < A < 1.0  !                            
 
 **********************************************************************
 
 FOR DETAILS SEE:                                                 
 
 AHRENS, J.H. AND DIETER, U.                            
 COMPUTER METHODS FOR SAMPLING FROM GAMMA,              
 BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              
 COMPUTING, 12 (1974), 223 - 246.                       
 
 (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    
 
 **********************************************************************
 INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
 OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
 COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
 COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
 COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
 PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
 SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
 */
{
	extern float fsign( float num, float sign );
	static float q1 = 4.166669E-2;
	static float q2 = 2.083148E-2;
	static float q3 = 8.01191E-3;
	static float q4 = 1.44121E-3;
	static float q5 = -7.388E-5;
	static float q6 = 2.4511E-4;
	static float q7 = 2.424E-4;
	static float a1 = 0.3333333;
	static float a2 = -0.250003;
	static float a3 = 0.2000062;
	static float a4 = -0.1662921;
	static float a5 = 0.1423657;
	static float a6 = -0.1367177;
	static float a7 = 0.1233795;
	static float e1 = 1.0;
	static float e2 = 0.4999897;
	static float e3 = 0.166829;
	static float e4 = 4.07753E-2;
	static float e5 = 1.0293E-2;
	static float aa = 0.0;
	static float aaa = 0.0;
	static float sqrt32 = 5.656854;
	/* JJV added b0 to fix rare and subtle bug */
	static float sgamma,s2,s,d,t,x,u,r,q0,b,b0,si,c,v,q,e,w,p;
    if(a == aa) goto S10;
    if(a < 1.0) goto S120;
	/*
		STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
	 */
    aa = a;
    s2 = a-0.5;
    s = sqrt(s2);
    d = sqrt32-12.0*s;
S10:
		/*
		STEP  2:  T=STANDARD NORMAL DEVIATE,
		 X=(S,1/2)-NORMAL DEVIATE.
		 IMMEDIATE ACCEPTANCE (I)
		 */
		t = snorm();
    x = s+0.5*t;
    sgamma = x*x;
    if(t >= 0.0) return sgamma;
	/*
		STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
	 */
    u = ranf();
    if(d*u <= t*t*t) return sgamma;
	/*
		STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
	 */
    if(a == aaa) goto S40;
    aaa = a;
    r = 1.0/ a;
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
	/*
		APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
	 THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
	 C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
	 */
    if(a <= 3.686) goto S30;
    if(a <= 13.022) goto S20;
	/*
		CASE 3:  A .GT. 13.022
	 */
    b = 1.77;
    si = 0.75;
    c = 0.1515/s;
    goto S40;
S20:
		/*
		CASE 2:  3.686 .LT. A .LE. 13.022
		 */
		b = 1.654+7.6E-3*s2;
    si = 1.68/s+0.275;
    c = 6.2E-2/s+2.4E-2;
    goto S40;
S30:
		/*
		CASE 1:  A .LE. 3.686
		 */
		b = 0.463+s+0.178*s2;
    si = 1.235;
    c = 0.195/s-7.9E-2+1.6E-1*s;
S40:
		/*
		STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
		 */
		if(x <= 0.0) goto S70;
	/*
		STEP  6:  CALCULATION OF V AND QUOTIENT Q
	 */
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S50;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S60;
S50:
		q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S60:
		/*
		STEP  7:  QUOTIENT ACCEPTANCE (Q)
		 */
		if(log(1.0-u) <= q) return sgamma;
S70:
		/*
		STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
		 U= 0,1 -UNIFORM DEVIATE
		 T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
		 */
		e = sexpo();
    u = ranf();
    u += (u-1.0);
    t = b+fsign(si*e,u);
	/*
		STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
	 */
    if(t < -0.7187449) goto S70;
	/*
		STEP 10:  CALCULATION OF V AND QUOTIENT Q
	 */
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S80;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S90;
S80:
		q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S90:
		/*
		STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
		 */
		if(q <= 0.0) goto S70;
    if(q <= 0.5) goto S100;
	/*
		* JJV modified the code through line 115 to handle large Q case
	 */
    if(q < 15.0) goto S95;
	/*
		* JJV Here Q is large enough that Q = log(exp(Q) - 1.0) (for real Q)
	 * JJV so reformulate test at 110 in terms of one EXP, if not too big
	 * JJV 87.49823 is close to the largest real which can be
	 * JJV exponentiated (87.49823 = log(1.0E38))
	 */
    if((q+e-0.5*t*t) > 87.49823) goto S115;
    if(c*fabs(u) > exp(q+e-0.5*t*t)) goto S70;
    goto S115;
S95:
		w = exp(q)-1.0;
    goto S110;
S100:
		w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
S110:
		/*
		IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
		 */
		if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
S115:
		x = s+0.5*t;
    sgamma = x*x;
    return sgamma;
S120:
		/*
		ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
		 
		 JJV changed B to B0 (which was added to declarations for this)
		 JJV in 120 to END to fix rare and subtle bug.
		 JJV Line: 'aa = 0.0' was removed (unnecessary, wasteful).
		 JJV Reasons: the state of AA only serves to tell the A >= 1.0
		 JJV case if certain A-dependent constants need to be recalculated.
		 JJV The A < 1.0 case (here) no longer changes any of these, and
		 JJV the recalculation of B (which used to change with an
									 JJV A < 1.0 call) is governed by the state of AAA anyway.
		 aa = 0.0;
		 */
		b0 = 1.0+0.3678794*a;
S130:
		p = b0*ranf();
    if(p >= 1.0) goto S140;
    sgamma = exp(log(p)/ a);
    if(sexpo() < sgamma) goto S130;
    return sgamma;
S140:
		sgamma = -log((b0-p)/ a);
    if(sexpo() < (1.0-a)*log(sgamma)) goto S130;
    return sgamma;
}

float sexpo(void)
/*
**********************************************************************
 
 
 (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                
 
 
 **********************************************************************
 **********************************************************************
 
 FOR DETAILS SEE:                                                 
 
 AHRENS, J.H. AND DIETER, U.                            
 COMPUTER METHODS FOR SAMPLING FROM THE                 
 EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
 COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               
 
 ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
 'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       
 
 Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
 SUNIF.  The argument IR thus goes away.                          
 
 **********************************************************************
 Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
 (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
 */
{
	static float q[8] = {
		0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,
		.9999999
	};
	static long i;
	static float sexpo,a,u,ustar,umin;
	static float *q1 = q;
    a = 0.0;
    u = ranf();
    goto S30;
S20:
		a += *q1;
S30:
		u += u;
	/*
		* JJV changed the following to reflect the true algorithm and prevent
	 * JJV unpredictable behavior if U is initially 0.5.
	 *  if(u <= 1.0) goto S20;
	 */
    if(u < 1.0) goto S20;
    u -= 1.0;
    if(u > *q1) goto S60;
    sexpo = a+u;
    return sexpo;
S60:
		i = 1;
    ustar = ranf();
    umin = ustar;
S70:
		ustar = ranf();
    if(ustar < umin) umin = ustar;
    i += 1;
    if(u > *(q+i-1)) goto S70;
    sexpo = a+umin**q1;
    return sexpo;
}

float snorm(void)
/*
**********************************************************************
 
 
 (STANDARD-)  N O R M A L  DISTRIBUTION                           
 
 
 **********************************************************************
 **********************************************************************
 
 FOR DETAILS SEE:                                                 
 
 AHRENS, J.H. AND DIETER, U.                            
 EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             
 SAMPLING FROM THE NORMAL DISTRIBUTION.                 
 MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          
 
 ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
 (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  
 
 Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
 SUNIF.  The argument IR thus goes away.                          
 
 **********************************************************************
 THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
 H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
 */
{
	static float a[32] = {
		0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
		0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
		0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
		1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
		1.862732,2.153875
	};
	static float d[31] = {
		0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
		0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
		0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
		0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
	};
	static float t[31] = {
		7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
		1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
		2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
		4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
		9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
	};
	static float h[31] = {
		3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
		4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
		4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
		5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
		8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
	};
	static long i;
	static float snorm,u,s,ustar,aa,w,y,tt;
    u = ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long) (u);
    if(i == 32) i = 31;
    if(i == 0) goto S100;
	/*
		START CENTER
	 */
    ustar = u-(float)i;
    aa = *(a+i-1);
S40:
		if(ustar <= *(t+i-1)) goto S60;
    w = (ustar-*(t+i-1))**(h+i-1);
S50:
		/*
		EXIT   (BOTH CASES)
		 */
		y = aa+w;
    snorm = y;
    if(s == 1.0) snorm = -y;
    return snorm;
S60:
		/*
		CENTER CONTINUED
		 */
		u = ranf();
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
		tt = u;
    ustar = ranf();
S80:
		if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S70;
    ustar = ranf();
    goto S40;
S100:
		/*
		START TAIL
		 */
		i = 6;
    aa = *(a+31);
    goto S120;
S110:
		aa += *(d+i-1);
    i += 1;
S120:
		u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
		w = u**(d+i-1);
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
		tt = u;
S160:
		ustar = ranf();
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S150;
    u = ranf();
    goto S140;
}

float fsign( float num, float sign )
/* Transfers sign of argument sign to argument num */
{
	if ( ( sign>0.0f && num<0.0f ) || ( sign<0.0f && num>0.0f ) )
		return -num;
	else return num;
}

/************************************************************************
FTNSTOP:
Prints msg to standard error and then exits
************************************************************************/
void ftnstop(char* msg)
/* msg - error message */
{
	if (msg != NULL) fprintf(stderr,"%s\n",msg);
	exit(0);
}


/***************************************************************************
**
**  EstimateAR1Var
**
**  Simulates a univariate gaussian variable for time (t+1) that has 'mean' 
**  and 'std' and autocorrelation lag 1 value 'ro1' based on AR(1) model. 
**  'tt' is the current (t) value of the simulated variable 
**
***************************************************************************/
double EstimateAR1Var(double mean, double std, double ro1, double tt)
{
	return(mean + ro1*(tt-mean) + std*(sqrt(1-ro1*ro1))*random_normal(0.0,1.0));
}



//=========================================================================
//
//
//                  Section 2: Matrix and vector operations
//
//
//=========================================================================


/***************************************************************************
**
**  MatrixVectorProduct
**
**  Multiplies matrix 'a' of size n x m by vector 'b' (size m) and stores
** the results in vector 'res' (size n)
**
***************************************************************************/
void MatrixVectorProduct(double **a, int n, int m, double *b, double *res)
{
	//cerr<<"1: "<<(int)(sizeof(b))<<"  2: "<<sizeof(double)
	//    <<"  3: "<<(int)(sizeof(res))<<endl;
	//assert((int)(sizeof(b)/sizeof(double)) == m);
	//assert((int)(sizeof(res)/sizeof(double)) == n);
	
	double tempo;
	for (int i=0; i<n; i++) {
		tempo = 0.;
		for (int j=0; j<m; j++) 
			tempo += a[i][j]*b[j];
		res[i] = tempo;
	}
	return;
}

/***************************************************************************
**
**  VectorSummation
**
**  Sum two vectors 'a' and 'b' of size 'n' and stores in the resulting
**  vector 'res'
**
***************************************************************************/
void VectorSummation(double *a, double *b, int n, double *res)
{
	//assert((int)(sizeof(a)/sizeof(double)) <= n);
	//assert((int)(sizeof(b)/sizeof(double)) <= n);
	//assert((int)(sizeof(res)/sizeof(double)) <= n);
	
	for (int i = 0; i<n; i++)
		res[i] = a[i] + b[i];
	return;
}


//=========================================================================
//
//
//                  Section 3: complex1 Functions
//
//
//=========================================================================


void setComplex(fcomplex *a, double c, double d){
	a->r = c;
	a->i = d;
}

void setCompNumber(fcomplex *a, fcomplex b){
	a->r = b.r;
	a->i = b.i;
}

void Print(fcomplex *a){
	cout << setprecision(10) << setiosflags(ios::fixed | ios::showpoint)
	<< "( " << a->r << ", " << a->i << " ) " << endl;
}

double Cr(fcomplex *a){
	return a->r;
}

double Ci(fcomplex *a){
	return a->i;
}

void subtReal(fcomplex *a, double d){
	a->r = a->r - d;
}

void addReal(fcomplex *a, double d){
	a->r = a->r + d;
}

fcomplex CaddFF(fcomplex *a, double d){
	fcomplex c;
	c.r = a->r + d;
	c.i = a->i;
	return c;
}


fcomplex Cadd(fcomplex a, fcomplex b){
	fcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

fcomplex Csub(fcomplex a, fcomplex b){
	fcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}


fcomplex Cmul(fcomplex a, fcomplex b){
	fcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

fcomplex Complex(float re, float im){
	fcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

fcomplex Conjg(fcomplex z){
	fcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

fcomplex Cdiv(fcomplex a, fcomplex b){
	fcomplex c;
	float r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

float Cabs(fcomplex z){
	float x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

fcomplex Csqrt(fcomplex z){
	fcomplex c;
	float x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

fcomplex RCmul(double x, fcomplex a){
	fcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}

#define MAXITER 200

/*****************************************************************************
**  
**  rtsafe
**
**  Finds a root of the polynomial which lies in the range [x1 and x2] 
**  starting from initial guess xguess. Accuracy of estimation is xacc. 
**  c1, c2, c3 are the polynomial coefficients
**
*****************************************************************************/
double rtsafe(void (*funcd)
			  (double, double, double, double, double *, double *),
              double c1, double c2, double c3, 
              double x1, double x2, double xacc, double xguess)
{
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;
	
	// 'fl' & 'fh' below are the evaluation function values 
	// corresponding to arguments 'x1' and 'x2'
	
	(*funcd)(x1,c1,c2,c3,&fl,&df);
	
	(*funcd)(x2,c1,c2,c3,&fh,&df);
	
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		cerr<<"Root must be bracketed by negative and positive f_n values!"<<endl;
		cerr<<"fl = "<<fl<<"; fh = "<<fh<<"; xguess = "<<xguess<<endl<<flush;
	}
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {  // <-- Orient the search so that f(xl) < 0
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	//rts=0.5*(x1+x2);  // Initialize the guess for root,
	rts = xguess;       // We can input a better guess than central value
	dxold=fabs(x2-x1);  // the "stepsize before last",
	dx=dxold;           // and the last step
	
	(*funcd)(rts,c1,c2,c3,&f,&df);
	
	for (j=1; j<=MAXITER; j++) {  // Loop over allowed iterations
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) // Bisect if Newton out of range
			|| (fabs(2.0*f) > fabs(dxold*df))) {    // or not decreasing fast
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts; // Change in root is
		}                            // negligible, take it
		else {   
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts; // Convergence criterion
		
		(*funcd)(rts,c1,c2,c3,&f,&df);
		if (f < 0.0) // Maintain the bracket on the root
			xl=rts;
		else
			xh=rts;
	}
	cerr<<"Maximum number of iterations exceeded in rtsafe!"<<endl;
	return 0.0;
}
#undef MAXITER


/*****************************************************************************
**  
**  polynomialH
**
**  Evaluates the polynomial of interest as well as its derivative
**
*****************************************************************************/
void polynomialH(double x, double c1, double c2, double c3,
				 double *fn, double *df)
{
	double fk2 = 5./3.;
	double fk3 = 2./3.;
	*fn = c1*pow(x, fk2) + c2*x + c3;
	*df = fk2*c1*pow(x, fk3) + c2;
}


/*****************************************************************************
**  
**  QuadraticFormula
**
**  Finds the real roots of a quadratic equation with the arguments:
**  'a', 'b', and 'c'
**
*****************************************************************************/
void QuadraticFormula(double a, double b, double c, double *x1, double *x2)
{
	double Ds;
	// Compute discriminant
	Ds = b*b - 4*a*c;
	
	if (Ds >= 0) { 
		(*x1) = (-b + sqrt(Ds))/(2*a);
		(*x2) = (-b - sqrt(Ds))/(2*a);
	}
	else {
		cout<<"QuadraticFormula(): Negative discriminant!"
		<<"\nCheck variables!";
		//cout<<1/0.0<<endl;
		exit(2);
	}
	return;
}

#define NRANSI
#define TINY 1.0e-20
#define NR_END 1
#define FREE_ARG char*

void ludcmp(double a[][2], int n, int indx[], double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;
	
	vv=Vector(1,n);
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != (n-1)) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	free_Vector(vv,1,n);
}
#undef TINY
#undef NRANSI

void lubksb(double a[][2], int n, int indx[], double b[])
{
	int i,ii=0,ip,j;
	double sum;
	
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii+1)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=(n-1);i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

double *Vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
	
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in Vector()");
	return v-nl+NR_END;
}

void free_Vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void nrerror(const char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

#undef NR_END
#undef FREE_ARG

//=========================================================================
//
//
//                        End of mathutil.cpp
//
//
//=========================================================================
