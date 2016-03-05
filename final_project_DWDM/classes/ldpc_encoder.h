#ifndef LDPC_ENCODER
#define LDPC_ENCODER

#include "ldpc_common.h"

/**
 * An LdpcEncoder object represents the encoder, once setup its method encode
 * accepts and InfoWord and returns a CodeWord
 */
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