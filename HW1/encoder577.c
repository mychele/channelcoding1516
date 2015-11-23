#include "math.h"
#include <stdlib.h>
#include <float.h>
#include <mex.h>

static int const N = 3; // codeword length
static int const mem = 2; // memory length
static int const M = 2; // binary constellation
static int const numStates = 4;

void encoder577(double* u, double* y, int n)
{
	//int y_size = n*N;
	int u_index = 0;
	int state[2] = {0,0};
	for(; u_index < n; u_index++)
	{
		*y++ = ((int)u[u_index] + *(state+1))%mem;
		*y++ = ((int)u[u_index] + *(state) + *(state+1))%mem;
		*y++ = ((int)u[u_index] + *(state) + *(state+1))%mem;
		*(state+1) = *state;
		*state = (int)u[u_index];
	}
}