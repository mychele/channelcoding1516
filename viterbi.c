#include <stdio.h>
#include <stdlib.h>
#include "math.h"

/**
 * This function performs encoding of a given binary vector
 */
char* encode(char *input, int input_size);

int main()
{
	int inv_rate = 3;
	int input_size;
	input_size = 5.5*pow(10, 8);
	// execute
	char* vec;
	vec = (char *) malloc(input_size*sizeof(char));
	char *cursor_vec;
	for(cursor_vec = vec; cursor_vec < vec + input_size; cursor_vec++)
		*cursor_vec = rand()%2;
	char *input = vec;
	char *codeword = encode(input, input_size);
	//for(int i = 0; i < input_size*inv_rate; i++)
	//	printf("%d", codeword[i]);
	free(codeword);
	free(input);
}

char* encode(char *input, int input_size)
{
	// allocate memory for the codeword vector
	char *codeword = (char *) malloc(3*input_size*sizeof(char));

	// the state will be updated on the go
	char s1 = 0;
	char s2 = 0;

	// to iterate through the vector
	char *cursor_input = input;
	char *cursor_cw = codeword;
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

