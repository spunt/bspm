#include "mex.h"
#include "math.h"

/* mex -v -f ./mexopts.sh imodwtj.c */

void imodwtj(double *Win, double *Vin, int *N, int *j, int *L, 
            double *ht, double *gt, double *Vout)
{

  int k, n, t;

  for(t = 0; t < *N; t++) {
    k = t;
    Vout[t] = (ht[0] * Win[k]) + (gt[0] * Vin[k]);
    for(n = 1; n < *L; n++) {
      k += (int) pow(2.0, (double) *j - 1.0);
      if(k >= *N) k -= *N;
      Vout[t] += (ht[n] * Win[k]) + (gt[n] * Vin[k]);
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
  int N, J, L;
  int m, n;
  double *Vout;
  double *Win, *Vin, *ht, *gt;
  
  /* Check for proper number of arguments */
  
  if (nrhs != 5) {
    mexErrMsgTxt("DWT requires five input arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("DWT requires one output argument.");
  }
  
  Win = mxGetPr(prhs[0]);
  Vin = mxGetPr(prhs[1]);
  ht = mxGetPr(prhs[2]);
  gt = mxGetPr(prhs[3]);

  J = (int) mxGetScalar(prhs[4]);
  #ifdef DEBUG
    mexPrintf("J = %d\n", J);
  #endif

  m = mxGetM(prhs[0]);
  #ifdef DEBUG
    mexPrintf("m = %d\n", m);
  #endif

  n = mxGetN(prhs[0]);
  #ifdef DEBUG
    mexPrintf("n = %d\n", n);
  #endif
  
  /* Create matrices for the return arguments */
  
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  
  /* Assign pointers to the various parameters */
  
  Vout = mxGetPr(plhs[0]);
  
  N = mxGetNumberOfElements(prhs[0]);
  L = mxGetNumberOfElements(prhs[2]);
  
  /* Do the actual computations in a subroutine */
  
  imodwtj(Win, Vin, &N, &J, &L, ht, gt, Vout);
  return;
}





