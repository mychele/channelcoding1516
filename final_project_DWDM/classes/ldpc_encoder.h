#include "ldpc_common.h"

#ifndef LDPC_ENCODER
#define LDPC_ENCODER

class LdpcEncoder
{
public:
	LdpcEncoder();
	//~LdpcEncoder();

	/**
	 * Read matrix K from file
	 * @return -1 if problems in opening the file, 0 otherwise
	 */
	int setup();

	// TODO condider if it is a valuable improvement encoding without the all 0 bit
	/**
	 * Encode the given information word
	 * @ InfoWord, a bitset with the information bit and the 173 0 bit in (292-d, 292) 0<=d<=172
	 */
	CodeWord encode(InfoWord info_word);


private:
	EncodingMatrix m_encodingMatrix;

	
};

#endif
// LDPC_ENCODER