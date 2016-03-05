#include <bitset>
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>

#ifndef LDPC_COMMON
#define LDPC_COMMON

#define INFO_ROWS 105
#define PC_ROWS 7
#define ALL_ROWS 112
#define ALL_COLUMNS 293
#define ALL_BIT 32816
#define INFO_BIT 30592
#define INIT_ZERO_BIT 173
#define ALL_INFO_BIT 30765
#define CODE_WORD 32640
#define IND_EQ 2045
#define ALL_EQ 2051
#define MAX_ATTEMPTS 99

typedef std::vector< std::bitset<ALL_INFO_BIT> > EncodingMatrix;
typedef std::bitset< CODE_WORD > CodeWord;
typedef std::bitset< ALL_INFO_BIT > InfoWord;
typedef std::pair< int, int > line; // line.first = slope_index, line.second = c
typedef std::pair< int, int > coordinates; // (a, b)

static const int slopes[7] = {1,2,3,5,7,11,13};

// map j to the k of the scrambled input vector
inline static int mapJtoK(int j) {
   	int col_r = ALL_COLUMNS;
	int q = j + (int)INIT_ZERO_BIT; // j + 172, with j from 1 to 30592
	int r = q/col_r; // floor(q/293) since both q and 293 are positive int, the result is the floor
	int col_index = col_r*r + col_r - 1 - q; // 293*r + 292 - q
	return col_r*r + col_index;
}

// compute x in (sum_term+x)%div = remainder, with x >= 0
inline static int reverseModulo(int div, int sum_term, int remainder) {
   	if(sum_term < remainder) {
      	return remainder - sum_term;
   	}
  	else {
  		int possible = div + remainder - sum_term;
  		while (possible < 0) {
  			possible += div;
  		}
  		if(possible == (int)ALL_COLUMNS) {possible = 0;}
  		return possible;
  	}
   
}

// phiTilde function
inline static double phiTilde(double x) {
	if (x > 38) return 0;
	if (x < 1.0e-300) return std::numeric_limits<double>::infinity();
	return -std::log(std::tanh(x/2));
}

#endif
// LDPC_COMMON