/*
 * STATUS: I think this will work for gradual and punctuated change 
 * with any number of states.
 *
 * For debugging:
 *	 $ gcc -c -I/usr/share/R/include -fPIC gp-bisse-eqs.c
 *	 $ gcc -shared -o gp-bisse-eqs.so gp-bisse-eqs.o
 *	 or instead:
 *	 $ R CMD SHLIB gp-bisse-eqs.c
 *	 That creates gp-bisse-eqs.so.
 *	 In R, dyn.load("gp-bisse-eqs.so").
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static double *lambda;
static double *mu;
static double *q;
static int *nstates;

/* get the list element named str, or return NULL; from R-exts.c */
SEXP getListElement(SEXP list, const char *str)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;

	for (i = 0; i < length(list); i++)
	{
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	 }
	return elmt;
}

void initmod_gp(void (* odeparms)(int *, double *))
{
	DL_FUNC get_deSolve_gparms;
	SEXP gparms;
	SEXP temp;
	const char *varnames[4] = {"nstates", "lambda", "mu", "q"};

	/* acquire R's list of parameters */
	get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
	gparms = get_deSolve_gparms();

	/* should have 4 elements: nstates, lambda, mu, q */
	if (LENGTH(gparms) != 4)
	{
		PROBLEM "Confusion over the length of parms." ERROR;
	}
	else
	{
		/* the number of states; an integer */
		temp = getListElement(gparms, varnames[0]);
		if (temp == R_NilValue)
			PROBLEM "nstates not found" ERROR;
		nstates = INTEGER(temp);
		/* (ERROR terminates, so not bothering with else's and {}'s */

		/* speciation rates; indexed by parent, daught1, daught2, 
 		 * so a nstates x nstates x nstates array */
		temp = getListElement(gparms, varnames[1]);
		if (temp == R_NilValue)
			PROBLEM "lambda not found" ERROR;
		lambda = REAL(temp);

		/* extinction rates; indexed by lineage state, 
 		 * so a length-nstates vector */
		temp = getListElement(gparms, varnames[2]);
		if (temp == R_NilValue)
			PROBLEM "mu not found" ERROR;
		mu = REAL(temp);

		/* transition rates; indexed by oldstate, newstate, 
 		 * so a nstates x nstates array */
		temp = getListElement(gparms, varnames[3]);
		if (temp == R_NilValue)
			PROBLEM "q not found" ERROR;
		q = REAL(temp);
	}
} 

void derivs_gp(int *neq, double *t, double *y, double *ydot, 
               double *yout, int *ip)
{
	double Di, Dj, Dk, Ei, Ej, Ek, q_ij, lam_ijk;
	int i, j, k;
	const int ns = *(nstates);  /* just to keep notation shorter */

	/* 
 	 * the E's are elements 0 ... ns-1  of y
 	 * the D's are elements ns ... 2ns-1 of y
	 *
	 * the [i,j] element of q is *(q + i + j*ns)
	 *    (transition from i to j)
	 * the [i, j, k] element of lambda is *(lambda + i + j*ns + k*ns*ns)
	 *   (parent = i, daughters = j and k)
	 */

	for (i=0; i<ns; i++)
	{
		Di = y[ns+i];
		Ei = y[i];

		/*** extinction ***/

		/* start dDi / dt */
		ydot[ns+i] = -mu[i] * Di;

		/* start dEi / dt */
		ydot[i] = mu[i] - mu[i]  * Ei;

		/*** transition ***/

		for (j=0; j<ns; j++)
		{
			if (j != i)
			{
				q_ij = *(q + i+j*ns);
				Dj = y[ns+j];
				Ej = y[j];

				/* continue dDi / dt */
				ydot[ns+i] += q_ij * Dj - q_ij * Di;

				/* continue dEi / dt */
				ydot[i] += q_ij * Ej - q_ij * Ei;
			}
		}

		/*** speciation ***/

		for (j=0; j<ns; j++)
		{
			Dj = y[ns+j];
			Ej = y[j];

			/* order of daughters doesn't matter, so assume that j <= k; 
			 * elements with j > k should be zero/ignored */
			for (k=0; k<=j; k++)
			{
				Dk = y[ns+k];
				Ek = y[k];
				lam_ijk = *(lambda + i + j*ns + k*ns*ns);

				/* continue dDi / dt */
				ydot[ns+i] += -lam_ijk * Di + lam_ijk * (Dj*Ek + Dk*Ej);

				/* continue dEi / dt */
				ydot[i] += -lam_ijk * Ei + lam_ijk * Ej * Ek;
			}
		}
	}
}
