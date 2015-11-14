#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <float.h>
#include "viterbi517_windowed.h"
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
    int windowSize;


    //

    /* 1. Check validity of expressions */
    
    // check input length
    if (nrhs != 4) // we need the received vector and the noise variance
        mexErrMsgTxt("Two input arguments required");
    // check output length
    if (nlhs != 1) // return the decoded vector
        mexErrMsgTxt("One output argument required");
    
    
    /* 2. Read inputs */
    
    // received values
    r = mxGetPr(prhs[0]);
    // vector length
    n = mxGetN(prhs[0]);
    // noise std deviation
    sigma_w = mxGetScalar(prhs[1]);
    // SD or HD
    mode = mxGetScalar(prhs[2]);
    // window size
    windowSize = mxGetScalar(prhs[3]);
    
    //printf("n = %d, sigw = %f\n", n, sigma_w);
    
    
    /* 3. Prepare output vectors */
    int u_hat_size = n/N; // TODO sanitize length
    plhs[0] = mxCreateDoubleMatrix(1, u_hat_size, mxREAL);
    u_hat = mxGetPr(plhs[0]);
    //mexPrintf("%d\n", u_hat_size);
    
    
    /* 4. Run the algorithm */
    viterbi517_windowed(r,sigma_w,n,u_hat,mode,windowSize);
    
    
    /* 5. Exit */
    
    return;
    
}