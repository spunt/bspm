#include "mex.h"
#include "math.h"

/* $Id: modwtj.c 300 2004-03-18 02:15:10Z ccornish $ */

/* To compile manually, type:
     mex -v -f $MATLAB/bin/mexopts.sh modwtj.c
 */

void modwtj(double *Vin, int N, int j, int L, double *ht, double *gt, 
	    double *Wout, double *Vout)
{

  int k, n, t;

  for (t = 0; t < N; t++) {
    k = t;
    Wout[t] = ht[0] * Vin[k];
    Vout[t] = gt[0] * Vin[k];
    for (n = 1; n < L; n++) {
      k -= (int) pow(2.0, (double) j - 1.0);
      if (k < 0) {
	k += N;
      }
      Wout[t] += ht[n] * Vin[k];
      Vout[t] += gt[n] * Vin[k];
    }
  }

}

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
  int N, j, L;
  int m, n;
  double *Wout, *Vout;
  double *Vin, *ht, *gt;
  
  /* Check for proper number of arguments */
  
  if (nrhs != 4) {
    mexErrMsgTxt("modwtj requires four input arguments.");
  } else if (nlhs > 2) {
    mexErrMsgTxt("modwtj requires two output arguments.");
  }

  /* Read input arguments */
  
  Vin = mxGetPr(prhs[0]);
  ht = mxGetPr(prhs[1]);
  gt = mxGetPr(prhs[2]);

  N = mxGetNumberOfElements(prhs[0]);
  L = mxGetNumberOfElements(prhs[1]);

  j = (int) mxGetScalar(prhs[3]);
  #ifdef DEBUG
    mexPrintf("j = %d\n", j);
  #endif

  m = mxGetM(prhs[0]);
  #ifdef DEBUG
    mexPrintf("m = %d\n", m);
  #endif

  n = mxGetN(prhs[0]);
  #ifdef DEBUG
    mexPrintf("n = %d\n", n);
  #endif
  
  /* Create arrays for the output arguments */
  
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
  
  /* Assign pointers to output arrays. */
  
  Wout = mxGetPr(plhs[0]);
  Vout = mxGetPr(plhs[1]);
  
  /* Compute the modwt at the jth level. */
  
  modwtj(Vin, N, j, L, ht, gt, Wout, Vout);

  return;
}
