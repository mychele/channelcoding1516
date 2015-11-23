#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <float.h>
#include "encoder577.h"
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
    double *y, *u;
    int n;

    //

    /* 1. Check validity of expressions */
    
    // check input length
    if (nrhs != 1) // we need the received vector and the noise variance
        mexErrMsgTxt("Two input arguments required");
    // check output length
    if (nlhs != 1) // return the decoded vector
        mexErrMsgTxt("One output argument required");
    
    
    /* 2. Read inputs */
    
    // input vector
    u = mxGetPr(prhs[0]);
    // vector length
    n = mxGetN(prhs[0]);
    
    /* 3. Prepare output vectors */
    int y_size = n*N; 
    plhs[0] = mxCreateDoubleMatrix(1, y_size, mxREAL);
    y = mxGetPr(plhs[0]);    
    
    /* 4. Encode */
    encoder577(u, y, n);
    
    
    /* 5. Exit */
    
    return;
    
}