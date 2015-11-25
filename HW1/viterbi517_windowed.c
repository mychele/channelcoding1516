#include "math.h"
#include <stdlib.h>
#include <float.h>
#include <mex.h>
#include "viterbi517.h"

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

void viterbi517_windowed(double *r, double sigma_w, int n, double *u_hat, double mode, int windowSize)
{
	int const outSize = n/N;				

	// initialize useful quantities
	double gammaPrev[numStates] = {0, DBL_MAX/2, DBL_MAX/2, DBL_MAX/2};
	double gamma[numStates] = {DBL_MAX/2, DBL_MAX/2, DBL_MAX/2, DBL_MAX/2};

	// use a circular matrix
	char circularPrevStateMatrix[windowSize][numStates];
	int currentIndex = 0;

	int l;
	// cycle on received vector
	for(l = 0; l < n; l = l + 3)
	{
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
			double cost = DBL_MAX/2;
			// each node has mem neighbors, cycle on them
			int maxNeighID = neighbors[stateID][mem - 1];
			int neighID = neighbors[stateID][0];
			for(; neighID <= maxNeighID; neighID++)
			{
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
		}

		if(l >= (windowSize-1)*3) // the circular matrix is full
		{
			// backtrack
			int bckIndex = currentIndex;
			int bckIter = 0;
			for(; bckIter < windowSize - 1; bckIter++)
			{
				minCostState = circularPrevStateMatrix[bckIndex][minCostState];
				bckIndex = (bckIndex - 1) >= 0 ? bckIndex - 1 : windowSize-1;
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
		minCostState = circularPrevStateMatrix[bckIndex][minCostState];
		bckIndex = (bckIndex - 1) >= 0 ? bckIndex - 1 : windowSize - 1;
	}
	u_hat[l/N - bckIter - 1] = minCostState/mem;
	
}



