#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <float.h>
#include "viterbi517.h"
#define N 3

/************************************************************************
 * Main function
 ************************************************************************/
void mexFunction(
        int nlhs,               // number of outputs
        mxArray *plhs[],        // outputs vector
        int nrhs,               // number of inputs
        const mxArray *prhs[]   // inputs vector
        )
{
    // useful variables
    double *r, *u_hat, sigma_w, mode;
    int n;

    // check input length
    if (nrhs != 3) // we need the received vector and the noise variance
        mexErrMsgTxt("Three input arguments required");
    // check output length
    if (nlhs != 1) // return the decoded vector
        mexErrMsgTxt("One output argument required");
    
    //read input and prepare output
    // received values
    r = mxGetPr(prhs[0]);
    // vector length
    n = mxGetN(prhs[0]);
    // noise std deviation
    sigma_w = mxGetScalar(prhs[1]);
    // SD or HD
    mode = mxGetScalar(prhs[2]);    
    // out size
    int u_hat_size = n/N;
    plhs[0] = mxCreateDoubleMatrix(1, u_hat_size, mxREAL);
    u_hat = mxGetPr(plhs[0]);
    
    // decode
    viterbi517(r,sigma_w,n,u_hat,mode);

    return; 
}