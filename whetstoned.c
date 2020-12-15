/*
 * C Converted Whetstone Double Precision Benchmark
 *		Version 1.2	22 March 1998
 *
 *	(c) Copyright 1998 Painter Engineering, Inc.
 *		All Rights Reserved.
 *
 *		Permission is granted to use, duplicate, and
 *		publish this text and program as long as it
 *		includes this entire comment block and limited
 *		rights reference.
 *
 * Converted by Rich Painter, Painter Engineering, Inc. based on the
 * www.netlib.org benchmark/whetstoned version obtained 16 March 1998.
 *
 * A novel approach was used here to keep the look and feel of the
 * FORTRAN version.  Altering the FORTRAN-based array indices,
 * starting at element 1, to start at element 0 for C, would require
 * numerous changes, including decrementing the variable indices by 1.
 * Instead, the array E1[] was declared 1 element larger in C.  This
 * allows the FORTRAN index range to function without any literal or
 * variable indices changes.  The array element E1[0] is simply never
 * used and does not alter the benchmark results.
 *
 * The major FORTRAN comment blocks were retained to minimize
 * differences between versions.  Modules N5 and N12, like in the
 * FORTRAN version, have been eliminated here.
 *
 * An optional command-line argument has been provided [-c] to
 * offer continuous repetition of the entire benchmark.
 * An optional argument for setting an alternate LOOP count is also
 * provided.  Define PRINTOUT to cause the POUT() function to print
 * outputs at various stages.  Final timing measurements should be
 * made with the PRINTOUT undefined.
 *
 * Questions and comments may be directed to the author at
 *			r.painter@ieee.org
 */
/*
C**********************************************************************
C     Benchmark #2 -- Double  Precision Whetstone (A001)
C
C     o	This is a REAL*8 version of
C	the Whetstone benchmark program.
C
C     o	DO-loop semantics are ANSI-66 compatible.
C
C     o	Final measurements are to be made with all
C	WRITE statements and FORMAT sttements removed.
C
C**********************************************************************   
*/

/* standard C library headers required */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* the following is optional depending on the timing function used */
#include <sys/time.h>
#include <time.h>

/* map the FORTRAN math functions, etc. to the C versions */
#define DSIN	sin
#define DCOS	cos
#define DATAN	atan
#define DLOG	log
#define DEXP	exp
#define DSQRT	sqrt
#define IF		if

/*
	COMMON T,T1,T2,E1(4),J,K,L
*/

long double_precision(void)
{
	/* used in the FORTRAN version */
	long I;
	long N1, N2, N3, N4, N6, N7, N8, N9, N10, N11;
	double X1,X2,X3,X4,X,Y,Z;
	long LOOP;
	int II, JJ;
	double T,T1,T2,E1[5];
	int J,K,L;


	long loopstart;
	long startsec, finisec;
	float KIPS;
	int continuous;
	struct timeval s_tv, e_tv;

	loopstart = 5000000;		/* see the note about LOOP below */
	continuous = 0;

	II = 1;		/* start at the first arg (temp use of II here) */


LCONT:
/*
C
C	Start benchmark timing at this point.
C
*/
	gettimeofday(&s_tv, NULL);
	startsec = time(0);

/*
C
C	The actual benchmark starts here.
C
*/
	T  = .499975;
	T1 = 0.50025;
	T2 = 2.0;
/*
C
C	With loopcount LOOP=10, one million Whetstone instructions
C	will be executed in EACH MAJOR LOOP..A MAJOR LOOP IS EXECUTED
C	'II' TIMES TO INCREASE WALL-CLOCK TIMING ACCURACY.
C
	LOOP = 1000;
*/
	LOOP = loopstart;
	II   = 1;

	JJ = 1;

IILOOP:
	N1  = 0;
	N2  = 12 * LOOP;
	N3  = 14 * LOOP;
	N4  = 345 * LOOP;
	N6  = 210 * LOOP;
	N7  = 32 * LOOP;
	N8  = 899 * LOOP;
	N9  = 616 * LOOP;
	N10 = 0;
	N11 = 93 * LOOP;
/*
C
C	Module 1: Simple identifiers
C
*/
	X1  =  1.0;
	X2  = -1.0;
	X3  = -1.0;
	X4  = -1.0;

	for (I = 1; I <= N1; I++) {
	    X1 = (X1 + X2 + X3 - X4) * T;
	    X2 = (X1 + X2 - X3 + X4) * T;
	    X3 = (X1 - X2 + X3 + X4) * T;
	    X4 = (-X1+ X2 + X3 + X4) * T;
	}

/*
C
C	Module 2: Array elements
C
*/
	E1[1] =  1.0;
	E1[2] = -1.0;
	E1[3] = -1.0;
	E1[4] = -1.0;

	for (I = 1; I <= N2; I++) {
	    E1[1] = ( E1[1] + E1[2] + E1[3] - E1[4]) * T;
	    E1[2] = ( E1[1] + E1[2] - E1[3] + E1[4]) * T;
	    E1[3] = ( E1[1] - E1[2] + E1[3] + E1[4]) * T;
	    E1[4] = (-E1[1] + E1[2] + E1[3] + E1[4]) * T;
	}


/*
C
C	Module 3: Array as parameter
C
*/
	for (I = 1; I <= N3; I++)
	{
			J = 0;

		L10:
			E1[1] = ( E1[1] + E1[2] + E1[3] - E1[4]) * T;
			E1[2] = ( E1[1] + E1[2] - E1[3] + E1[4]) * T;
			E1[3] = ( E1[1] - E1[2] + E1[3] + E1[4]) * T;
			E1[4] = (-E1[1] + E1[2] + E1[3] + E1[4]) / T2;
			J += 1;

			if (J < 6)
				goto L10;			
	}

/*
C
C	Module 4: Conditional jumps
C
*/
	J = 1;
	for (I = 1; I <= N4; I++) {
		if (J == 1)
			J = 2;
		else
			J = 3;

		if (J > 2)
			J = 0;
		else
			J = 1;

		if (J < 1)
			J = 1;
		else
			J = 0;
	}


/*
C
C	Module 5: Omitted
C 	Module 6: Integer arithmetic
C
*/

	J = 1;
	K = 2;
	L = 3;

	for (I = 1; I <= N6; I++) {
	    J = J * (K-J) * (L-K);
	    K = L * K - (L-J) * K;
	    L = (L-K) * (K+J);
	    E1[L-1] = J + K + L;
	    E1[K-1] = J * K * L;
	}


/*
C
C	Module 7: Trigonometric functions
C
*/
	X = 0.5;
	Y = 0.5;

	for (I = 1; I <= N7; I++) {
		X = T * DATAN(T2*DSIN(X)*DCOS(X)/(DCOS(X+Y)+DCOS(X-Y)-1.0));
		Y = T * DATAN(T2*DSIN(Y)*DCOS(Y)/(DCOS(X+Y)+DCOS(X-Y)-1.0));
	}

/*
C
C	Module 8: Procedure calls
C
*/
	X = 1.0;
	Y = 1.0;

	for (I = 1; I <= N8; I++)
	{
		double X1, Y1;

		X1 = X;
		Y1 = Y;
		X1 = T * (X1 + Y1);
		Y1 = T * (X1 + Y1);
		Z = (X1 + Y1) / T2;		
	}	

/*
C
C	Module 9: Array references
C
*/
	J = 1;
	K = 2;
	L = 3;
	E1[1] = 1.0;
	E1[2] = 2.0;
	E1[3] = 3.0;

	for (I = 1; I <= N9; I++)
	{
		E1[J] = E1[K];
		E1[K] = E1[L];
		E1[L] = E1[J];
	}


/*
C
C	Module 10: Integer arithmetic
C
*/
	J = 2;
	K = 3;

	for (I = 1; I <= N10; I++) {
	    J = J + K;
	    K = J + K;
	    J = K - J;
	    K = K - J - J;
	}


/*
C
C	Module 11: Standard functions
C
*/
	X = 0.75;

	for (I = 1; I <= N11; I++)
		X = DSQRT(DEXP(DLOG(X)/T1));


/*
C
C      THIS IS THE END OF THE MAJOR LOOP.
C
*/
	if (++JJ <= II)
		goto IILOOP;

/*
C
C      Stop benchmark timing at this point.
C
*/
	gettimeofday(&e_tv, NULL);
	finisec = time(0);

/*
C----------------------------------------------------------------
C      Performance in Whetstone KIP's per second is given by
C
C	(100*LOOP*II)/TIME
C
C      where TIME is in seconds.
C--------------------------------------------------------------------
*/
	printf("#2 Double Precision Whetstone in [ %4.4lf ] seconds @ Loops: %ld, Iterations: %d.\n",
			(float) ((e_tv.tv_sec - s_tv.tv_sec)*1000000L + e_tv.tv_usec - s_tv.tv_usec) / 1000000L, LOOP, II);

	KIPS = (100.0*LOOP*II)/(float)(finisec-startsec);
	if (KIPS >= 1000.0)
		printf("  -- C Converted Double Precision Whetstones: %.1f MIPS\n", KIPS/1000.0);
	else
		printf("  -- C Converted Double Precision Whetstones: %.1f KIPS\n", KIPS);

	if (continuous)
		goto LCONT;

	return loopstart;
}

