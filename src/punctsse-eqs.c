/* PunctSSE compiled code */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include "util.h"

void do_derivs_punctsse(int n, double *pars, double *y, double *ydot, 
		              int jk_array[][2]);

static double *parms_punctsse;
void initmod_punctsse(void (* odeparms)(int *, double *))
{
  /* parameter length checking will be done in R */
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  parms_punctsse = REAL(get_deSolve_gparms());
} 

void derivs_punctsse(int *neq, double *t, double *y, double *ydot,
                     double *yout, int *ip)
{
  /* saves a little time to pre-compute indices outside of do_derivs_punctsse
   * would be better to do it outside of here, too */
  int n = *neq / 2;
  int len_lam_i = n * (n + 1) / 2;
  int jk_array[len_lam_i][2];
  int j, k, m;
  m = 0;
  for (j=0; j<n; j++)
  {
    for (k=j; k<n; k++)
    {
	 jk_array[m][0] = j;
	 jk_array[m][1] = k;
	 m++;
    }
  }

  do_derivs_punctsse(n, parms_punctsse, y, ydot, jk_array);
}

void do_derivs_punctsse(int n, double *pars, double *y, double *ydot, int jk_array[][2])
{
/* note: n = num states is called k elsewhere, but k is used as an index here */
  double *E = y, *D = y + n;
  double *dEdt = ydot, *dDdt = ydot + n;
  int len_lam_i = n * (n + 1) / 2;
  int len_lam = n * len_lam_i;
  double *lambda = pars, *lambda_i;
  double *mu = pars + len_lam, *Q = pars + len_lam + n;
  double tmp, Ei, Di;
  int i, j, k, m;

  for (i=0; i<n; i++)
  {
    Ei = E[i];
    Di = D[i];

    /* extinction (start ydot) */
    dEdt[i] = mu[i] * (1 - Ei);
    dDdt[i] = -mu[i] * Di;

    /* speciation (continue ydot) */
    lambda_i = lambda + i * len_lam_i;
    for (m=0; m<len_lam_i; m++)
    {
      j = jk_array[m][0];
      k = jk_array[m][1];
      dEdt[i] += lambda_i[m] * (-Ei + E[j] * E[k]);
      dDdt[i] += lambda_i[m] * (-Di + D[j] * E[k] + D[k] * E[j]);
    }
  }

  /* transitions (complete ydot by adding Q y to it) */
  do_gemm2(Q, n, n, y, n, 2, ydot);
}

/* no time-dependence yet */
