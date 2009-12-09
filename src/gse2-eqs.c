/*
 * The GSE equations, implemented in C
 *
 * Usage outside of the package:
 *    $ gcc -c -I/usr/share/R/include -fPIC gse2-eqs.c
 *    $ gcc -shared -o gse2-eqs.so gse2-eqs.o
 *    creates gse2-eqs.so
 *    > dyn.load('gse2-eqs.so')   # load in R
 */

#include <R.h>

static double parms[7];

#define sA parms[0]     /* speciation in region A */
#define sB parms[1]     /* speciation in region B */
#define sAB parms[2]    /* allopatric speciation */
#define xA parms[3]     /* extinction in region A */
#define xB parms[4]     /* extinction in region B */
#define dA parms[5]     /* dispersal from A to B */
#define dB parms[6]     /* dispersal from B to A */

void gse2_initmod(void (* odeparms)(int *, double *))
{
	int N = 7;
	odeparms(&N, parms);
}

/* states:
 * 1 = in both regions
 * 2 = region A endemic
 * 3 = region B endemic
 */

void gse2_derivs(int *neq, double *t, double *y, double *ydot, double *yout,
                 int *ip)
{
	double E_1 = y[0];
	double E_2 = y[1];
	double E_3 = y[2];

	double D_N1 = y[3];
	double D_N2 = y[4];
	double D_N3 = y[5];
	    		 
	/*  dE_1 / dt  */
	ydot[0] = -(sA + sB + xA + xB + sAB) * E_1 
	          + xA * E_3 + xB * E_2 
			+ sA * E_1 * E_2 + sB * E_1 * E_3 + sAB * E_2 * E_3;

	/*  dE_2 / dt  */
	ydot[1] = -(sA + dA + xA) * E_2 
	          + xA + dA * E_1 + sA * E_2 * E_2;

	/*  E_3 / dt  */
	ydot[2] = -(sB + dB + xB) * E_3 
	          + xB + dB * E_1 + sB * E_3 * E_3;

	/*  dD_N1 / dt  */
	ydot[3] = -(sA + sB + sAB + xA + xB) * D_N1 
	          + xA * D_N3 + xB * D_N2 
			+ sA * (E_2 * D_N1 + E_1 * D_N2) 
			+ sB * (E_3 * D_N1 + E_1 * D_N3)
			+ sAB * (E_2 * D_N3 + E_3 * D_N2);

	/*  dD_N2 / dt  */
	ydot[4] = -(sA + dA + xA) * D_N2 
	          + dA * D_N1 + 2 * sA * D_N2 * E_2;

	/*  dD_N3 / dt  */
	ydot[5] = -(sB + dB + xB) * D_N3 
	          + dB * D_N1 + 2 * sB * D_N3 * E_3;
}
