#include <stdio.h>
#include <stdlib.h>
#include "math.h"

/**
 * This function performs encoding of a given binary vector
 * @param input the binary input vector
 * @param input_size the length of the input vector
 * @return a pointer to a float vector of codewords that can be passed to the conform map
 */
float* encode(char *input, int input_size);

/**
 * This function conforms the codeword (1->1, 0->-1) and transmits it through the channel
 * @param codeword the codeword
 * @param codeword_size the length of the codeword vector
 */
void conformAndChannel(float *codeword, int codeword_size);

/**
 * This function performs encoding of a given binary vector
 * @param codeword the codeword vector to be decoded
 * @param codeword_size the length of the codeword vector
 * @return a pointer to a char vector of decoded output symbols
 */
char* decode(float *codeword, int codeword_size);

int main()
{
	char verb = 0;
	int inv_rate = 3;
	int input_size;
	input_size = 5*pow(10, 8);
	// execute
	char* vec;
	vec = (char *) malloc(input_size*sizeof(char));
	char *cursor_vec;
	for(cursor_vec = vec; cursor_vec < vec + input_size; cursor_vec++)
		*cursor_vec = rand()%2;
	char *input = vec;
	float *codeword = encode(input, input_size);
	if (verb)
	{ 
		for(int i = 0; i < input_size*inv_rate; i++)
			printf("%f\n", codeword[i]);
	}
	conformAndChannel(codeword, inv_rate*input_size);

	free(codeword);
	free(input);
}

float* encode(char *input, int input_size)
{
	// allocate memory for the codeword vector
	float *codeword = (float *) malloc(3*input_size*sizeof(float));

	// the state will be updated on the go
	char s1 = 0;
	char s2 = 0;

	// to iterate through the vector
	char *cursor_input = input;
	float *cursor_cw = codeword;
	while(cursor_cw < codeword + 3*input_size && cursor_input < input + input_size)
	{
		/* 
		 * y1 = u + s2
		 * y2 = s2
		 * y3 = u + s1 + s2
		 */
		*cursor_cw++ = (*cursor_input + s2)%2;
		*cursor_cw++ = s2;
		*cursor_cw++ = (*cursor_input + s1 + s2)%2;

		s2 = s1;
		s1 = *cursor_input++;
	}

	return codeword;
}

void conformAndChannel(float *codeword, int codeword_size)
{
	float *cursor_cw = codeword;
	while(cursor_cw < codeword + codeword_size)
	{
		if (*cursor_cw == 0) {
				*cursor_cw = -1;
		}
		// add a random noise

		cursor_cw++;
	}
}



/*char* decode(float *codeword, int codeword_size)
{

}*/

float* generateGaussianNoise(float sigma)
{
	int still_one_available = 1;
	float u1 = rand()/RAND_MAX;
	float u2 = rand()/RAND_MAX;
	double x1 = sigma*sqrt(-2*log(u1))*cos(2*M_PI*u2);
	double x2 = sigma*sqrt(-2*log(u1))*sin(2*M_PI*u2);
}

