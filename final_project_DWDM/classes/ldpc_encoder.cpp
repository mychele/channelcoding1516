#include <fstream>
#include <iostream>
#include "ldpc_encoder.h"
#include "bit_io.h"

LdpcEncoder::LdpcEncoder() {
	// init private variables
	EncodingMatrix m_encodingMatrix = EncodingMatrix();
}

// LdpcEncoder::~LdpcEncoder() {
// 	delete m_encodingMatrix;
// }

int
LdpcEncoder::setup() {
	std::cout << "Setup\n";
	// read K from K.bin
	std::ifstream bin_file_in("bin_files/K.bin", std::ios::binary);

	// create BitIo object
	BitIo<ALL_INFO_BIT> bit_io;

	// check if file is open
	if(bin_file_in.is_open()) {
		// read the whole file into bit_io object
		bin_file_in >> bit_io;
		bin_file_in.close();
		// save the bit_io in m_encodingMatrix, splitted into rows
		m_encodingMatrix.resize(2045);
		for(int row_index = 0; row_index < (int)IND_EQ; row_index++) { 
			m_encodingMatrix[row_index] = bit_io.pop_front();
		}
		return 0;
	} else {
		return -1; // error in opening the file
	}
}

CodeWord 
LdpcEncoder::encode(InfoWord info_word) {
	// info_word has 30765 entries. Apply the encoding matrix to get 2045 parity check bit, then create the codeword

	// encode using bitset
	std::bitset<IND_EQ> parity_check_bitset;
	for(int ar_index = 0; ar_index < (int)IND_EQ; ar_index++) {
		parity_check_bitset[ar_index] = ( (m_encodingMatrix[ar_index]&info_word).count() & 1);
												// multiply and sum					// if sum%2 = 1 -> 1, else -> 0
	}

	// create codeword
	CodeWord codeword = CodeWord();
	// copy information bit
	int bit_index = 0;
	for(; bit_index < (int)ALL_COLUMNS-(int)INIT_ZERO_BIT; bit_index++) { // cycle up to bit 119
		codeword[bit_index] = info_word[bit_index];
	}
	int offset_info = (int)INIT_ZERO_BIT;
	for(; bit_index < (int)ALL_INFO_BIT-(int)INIT_ZERO_BIT; bit_index++) { // cycle in the codeword up to bit 30765 - 173 - 1, in the info_word up to bit 30765 - 1
		codeword[bit_index] = info_word[bit_index + offset_info];
	}

	// copy parity check bit
	int pcb_index = 0;
	int offset_cw = 3; // 3 unused bits set to 0
	bit_index += offset_cw;
	for(; bit_index < (int)CODE_WORD; bit_index++) { // in the codeword from 30765 - 173, copy the 2045 parity check bit
		codeword[bit_index] = parity_check_bitset[pcb_index++];
	}

	return codeword;
}