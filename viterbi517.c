#include "math.h"
#include <stdlib.h>
#include <float.h>
#include <mex.h>

// LUT for the output, hardcoded for a 517 code
static int const y_lut[4][6] = {	  //u=0		//u=1
									{ 0, 0, 0, 1, 0, 1 }, // state 0
									{ 1, 1, 1, 0, 1, 0 }, // state 1
									{ 0, 0, 1, 1, 0, 0 }, // state 2
									{ 1, 1, 0, 0, 1, 1 }  // state 3
								};

// LUT for the neighbors, hardcoded for a 517 code
static int const neighbors[4][2] = {
										{0, 1},
										{2, 3},
										{0, 1},
										{2, 3}
									};

static int const N = 3; // codeword length
static int const mem = 2; // memory length
static int const M = 2; // binary constellation
static int const numStates = 4;
// static char const verb = 0;


/**
 * Compute cost of the transition from neighID to stateID with symbol u
 * @param the input symbol
 * @param the prev state
 * @param the codeword associated with the transition
 * @return the cost
 */
double getCost(int symbol, int neighID, double *codeword)
{
	double transCost = 0;
	int l = 0;
	for(; l < 3; l++)
	{
		double r_i = *(codeword + l);
		transCost += -2*r_i*y_lut[neighID][symbol*N + l];
		// if (verb)
		// {
		// 	printf("y%d=%d, r%d=%f\n", l, y_lut[neighID][symbol*N + l], l, -2*r_i);
		// }
	}
	// if(verb)
	// {
	// 	printf("cost=%f\n", transCost);
	// }
	return transCost;
}

void viterbi517(double *r, double sigma_w, int n, double *u_hat)
{
	int const outSize = n/N;							

	// initialize useful quantities
	double gammaPrev[numStates] = {0, DBL_MAX, DBL_MAX, DBL_MAX};
	double gamma[numStates] = {DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX};

	char prevState[outSize][numStates];

	int l;
	// cycle on received vector
	for(l = 0; l < n; l = l + 3)
	{
		// TODO consider if it is a valuable improvement to compute LLR just once
		// if(verb)
		// {
		// 	printf("%d\n", l);
		// 	for(int in = 0; in <= 3; in++)
		// 	{
		// 		printf("gamma(%d) = %f\n", in, gammaPrev[in]);
		// 	}
		// }
		// cycle on the states
		int stateID;
		double minCost = DBL_MAX; // this will always be updated because of 
		// the paths from 0 to 0 and 2
		for (stateID = 0; stateID < numStates; stateID++)
		{
			// compute the input symbol that can bring to this state
			// for 0 and 1 the input symbol that bring to the state is always
			// 0, for 2 and 3 the input symbol is 1
			int u_poss = stateID/mem;
			// if (verb)
			// {
			// 	printf("stateID = %d, possible input=%d\n", stateID, u_poss);
			// }
			double cost = DBL_MAX;
			// each node has mem neighbors, cycle on them
			int maxNeighID = neighbors[stateID][mem - 1];
			int neighID = neighbors[stateID][0];
			for(; neighID <= maxNeighID; neighID++)
			{
				// if (verb)
				// {
				// 	printf("neighID=%d\n", neighID);
				// }
				double newCost = gammaPrev[neighID] + getCost(u_poss, neighID, r + l)/pow(sigma_w,2);
				if (newCost < cost)
				{
					cost = newCost;
					prevState[l/N][stateID] = neighID;
					gamma[stateID] = cost;
				}
			}
			if (cost == DBL_MAX) // no paths going to this state, this happens 
				// for state 1 and 3 at the first step
			{
				prevState[l/N][stateID] = 0;
				gamma[stateID] = cost;
			}
			if (cost < minCost) // to find the minimum cost
			{
				minCost = cost;
			}

		}

		// consider using a windowed version

		// normalize gamma to avoid overflow
		int gammaIndex;
		for(gammaIndex = 0; gammaIndex < numStates; gammaIndex++)
		{
			gammaPrev[gammaIndex] = gamma[gammaIndex] - minCost;
			// if (verb)
			// {
			// 	printf("prevState=%d for stateID=%d\n", prevState[l/N][gammaIndex], gammaIndex);
			// }
		}
	}

	// backtrack from state 0
	int stateID = 0;
	l = outSize - 1;
	for(; l >= 0; l--)
	{
		// if (verb)
		// {
		// 	printf("stateID=%d\n", stateID);
		// }
		u_hat[l] = stateID/mem;
		stateID = prevState[l][stateID];
	}
}



