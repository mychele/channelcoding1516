#include "math.h"
#include <stdlib.h>
#include <float.h>
#include <mex.h>
#include "viterbi517.h"

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
static char const verb = 0;


void viterbi517_windowed(double *r, double sigma_w, int n, double *u_hat, double mode, int windowSize)
{
	int const outSize = n/N;	
	// printf("outSize=%d\n", outSize);						

	// initialize useful quantities
	double gammaPrev[numStates] = {0, DBL_MAX/2, DBL_MAX/2, DBL_MAX/2};
	double gamma[numStates] = {DBL_MAX/2, DBL_MAX/2, DBL_MAX/2, DBL_MAX/2};

	char circularPrevStateMatrix[windowSize][numStates];
	int currentIndex = 0;

	int l;
	// cycle on received vector
	for(l = 0; l < n; l = l + 3)
	{
		if(verb)
			printf("l=%d, l/N=%d\n", l, l/N);
		
		// cycle on the states
		int stateID;
		int minCostState;
		double minCost = DBL_MAX/2; // this will always be updated because of 
		// the paths from 0 to 0 and 2
		for (stateID = 0; stateID < numStates; stateID++)
		{
			// compute the input symbol that can bring to this state
			// for 0 and 1 the input symbol that bring to the state is always
			// 0, for 2 and 3 the input symbol is 1
			int u_poss = stateID/mem;
			if (verb)
			{
				printf("stateID = %d, possible input=%d\n", stateID, u_poss);
			}
			double cost = DBL_MAX/2;
			// each node has mem neighbors, cycle on them
			int maxNeighID = neighbors[stateID][mem - 1];
			int neighID = neighbors[stateID][0];
			for(; neighID <= maxNeighID; neighID++)
			{
				if (verb)
				{
					printf("neighID=%d\n", neighID);
				}
				double newCost = gammaPrev[neighID] + getCost(u_poss, neighID, r + l, sigma_w, mode);
				if (newCost < cost)
				{
					cost = newCost;
					circularPrevStateMatrix[currentIndex][stateID] = neighID;
					gamma[stateID] = cost;
				}
			}
			if (cost == DBL_MAX/2) // no paths going to this state, this happens 
				// for state 1 and 3 at the first step
			{
				circularPrevStateMatrix[l/N][stateID] = 0;
				gamma[stateID] = cost;
			}
			if (cost < minCost) // to find the minimum cost
			{
				minCost = cost;
				minCostState = circularPrevStateMatrix[currentIndex][stateID];
			}

		}

		// normalize gamma to avoid overflow
		int gammaIndex;
		for(gammaIndex = 0; gammaIndex < numStates; gammaIndex++)
		{
			gammaPrev[gammaIndex] = gamma[gammaIndex] - minCost;
			if (gammaPrev[gammaIndex] == 0)
			{
				minCostState = gammaIndex;
			}
			if (verb)
			{
				printf("in buffer prevState=%d for stateID=%d\n", circularPrevStateMatrix[l/N][gammaIndex], gammaIndex);
			}
		}

		if(verb)
		{
			for(int in = 0; in <= 3; in++)
			{
				printf("gamma(%d) = %g\n", in, gammaPrev[in]);
			}
		}

		if(l >= (windowSize-1)*3) // the circular matrix is full
		{
			// backtrack
			int bckIndex = currentIndex;
			int bckIter = 0;
			for(; bckIter < windowSize - 1; bckIter++)
			{
				if(verb)
					printf("bckIndex=%d, minCostState=%d\n", bckIndex, minCostState);
				minCostState = circularPrevStateMatrix[bckIndex][minCostState];
				bckIndex = (bckIndex - 1) >= 0 ? bckIndex - 1 : windowSize-1;
			}
			if(verb)
			{
				printf("minCostState=%d\n", minCostState);
				printf("decide on l=%d\n", l/N - windowSize + 1);
			}
			u_hat[l/N - windowSize + 1] = minCostState/mem;
		}
		currentIndex = (currentIndex+1)%windowSize;
	}
	
	//backtrack from state 0 to get the last windowSize symbols
	//backtrack
	int bckIndex = (currentIndex - 1) > 0 ? currentIndex - 1 : windowSize-1;
	int bckIter = 0;
	int minCostState = 0;
	for(; bckIter < windowSize - 1; bckIter++)
	{
		u_hat[l/N - bckIter - 1] = minCostState/mem;
		if(verb)
			printf("decide on l=%d, bckIter=%d, bckIndex=%d, minCostState=%d\n", l/N - bckIter - 1, bckIter, bckIndex, minCostState);
		minCostState = circularPrevStateMatrix[bckIndex][minCostState];
		bckIndex = (bckIndex - 1) >= 0 ? bckIndex - 1 : windowSize - 1;
	}
	if(verb)
		printf("minCostState=%d\n", minCostState);
	u_hat[l/N - bckIter - 1] = minCostState/mem;
	
}



