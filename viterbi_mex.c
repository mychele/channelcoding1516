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
    double *r, *u_hat, sigma_w;
    int n;

    //

    /* 1. Check validity of expressions */
    
    // check input length
    if (nrhs != 2) // we need the received vector and the noise variance
        mexErrMsgTxt("Two input arguments required");
    // check output length
    if (nlhs != 1) // return the decoded vector
        mexErrMsgTxt("One output argument required");
    
    
    /* 2. Read inputs */
    
    // received values
    r = mxGetPr(prhs[0]);
    // vector length
    n = mxGetN(prhs[0]);
    // noise variance
    sigma_w = mxGetScalar(mxGetField(prhs[1], 0, "sigw"));
    
    printf("n = %d, sigw = %f\n",n,sigma_w);
    
    
    /* 3. Prepare output vectors */
    int u_hat_size = n/N; // TODO sanitize length
    plhs[0] = mxCreateDoubleMatrix(1, u_hat_size, mxREAL);
    u_hat = mxGetPr(plhs[0]);
    
    
    /* 4. Run the algorithm */
    
    viterbi517(r,sigma_w,n,u_hat);
    
    
    /* 5. Exit */
    
    return;
    
}