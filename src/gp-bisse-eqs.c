/*
 * IN PROGRESS
 *
 * For debugging:
 *	 $ gcc -c -I/usr/share/R/include -fPIC gp-bisse-eqs.c
 *	 $ gcc -shared -o gp-bisse-eqs.so gp-bisse-eqs.o
 *	 creates gp-bisse-eqs.so
 *	 in R, dyn.load('gp-bisse-eqs.so')
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static double *lambda;
static double *mu;
static double *q;

/* will eventually need to index these dynamically */
#define lambda0 lambda[0]
#define lambda1 lambda[7]
#define mu0 mu[0]
#define mu1 mu[1]
#define q01 q[2]
#define q10 q[1]

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
	const char *varnames[3] = {"lambda", "mu", "q"};

	/* acquire R's list of parameters */
	get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
	gparms = get_deSolve_gparms();

	/* should have 3 elements: lambda, mu, q */
	if (LENGTH(gparms) != 3)
	{
		PROBLEM "Confusion over the length of parms" ERROR;
	}
	else
	{
		/* a 3D array; parent, daught1, daught2 */
		lambda = REAL(getListElement(gparms, varnames[0]));
		//printf("%f, %f, %f, %f, %f, %f, %f, %f\n", lambda[0], lambda[1], lambda[2], lambda[3], lambda[4], lambda[5], lambda[6], lambda[7]);

		/* a vector */
		mu = REAL(getListElement(gparms, varnames[1]));

		/* a 2D array; read col by col */
		q = REAL(getListElement(gparms, varnames[2]));
	}
} 

void derivs_gp(int *neq, double *t, double *y, double *ydot, 
               double *yout, int *ip)
{
	double E0 = y[0];
	double E1 = y[1];
	double D0 = y[2];
	double D1 = y[3];

	ydot[0] = -(mu0 + q01 + lambda0) * E0 + lambda0 * E0 * E0 + mu0 + q01 * E1;
	ydot[1] = -(mu1 + q10 + lambda1) * E1 + lambda1 * E1 * E1 + mu1 + q10 * E0;
	ydot[2] = -(mu0 + q01 + lambda0) * D0 + 2 * lambda0 * E0 * D0 + q01 * D1;
	ydot[3] = -(mu1 + q10 + lambda1) * D1 + 2 * lambda1 * E1 * D1 + q10 * D0;
}
