#include "ldpc_encoder.h"
#include "ldpc_decoder.h"
#include <bitset>
#include <iostream>
#include <random>
#include <chrono>


int main(int argc, char const *argv[])
{
	// rng stuff
	std::random_device rd; // device entropy
    std::mt19937 m_rng(rd()); // initialize our mersenne twister with a random seed
    std::uniform_int_distribution<int> bit_generator(0,1);
    std::normal_distribution<double> noise_generator(0,1);

	//create information word
	InfoWord infoword;
	for(int bit_index = 0; bit_index < (int)INFO_BIT; ++bit_index) { 
		infoword[mapJtoK(bit_index)] = (bit_generator(m_rng)==1);
	}

	// encode
	LdpcEncoder encoder = LdpcEncoder();
	std::cout << "LdpcEncoder created\n";
	// import matrix
	if(encoder.setup() == 0) {std::cout << "File opened and read successfully\n";}
	else {std::cout << "Error in opening matrix file, abort\n"; return 2;}
	CodeWord codeword = encoder.encode(infoword);

	// add noise
	double ebn0 = 7;
	double sigma_w = 1/std::sqrt(2*std::pow(10, (ebn0/10)) * INFO_BIT/CODE_WORD);
	std::cout << "sigma_w " << sigma_w << "\n";
	std::vector<double> *received_signal = new std::vector<double>;
	received_signal->resize(ALL_COLUMNS*ALL_ROWS, -1);

	// copy the first 120 bit of the codeword, then skip 173 zeros
	int cw_bit_index = 0;
	for(; cw_bit_index < ALL_COLUMNS - INIT_ZERO_BIT; cw_bit_index++) {
		received_signal->at(cw_bit_index) = ((codeword[cw_bit_index] == 0) ? -1:1)  + sigma_w*noise_generator(m_rng);
	}
	int all_offset = (int)INIT_ZERO_BIT;
	// copy up to						30765 - 1, which is bit 30592 - 1 in the codeword
	for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT; cw_bit_index++) { 
		received_signal->at(cw_bit_index + all_offset) = ((codeword[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_generator(m_rng);
	}
	// skip 3 zeros in the codeword, but not in all_bit
	cw_bit_index += 3;
	all_offset -= 3;
	for(int i = 1; i < 7; ++i) {
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + (i)*(int)ALL_COLUMNS - 1; cw_bit_index++) { 
		// skip bit 31057, 31350, 31643, 31936, 32229, 32522 in all_bits
			received_signal->at(cw_bit_index + all_offset) = ((codeword[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_generator(m_rng);
		}
		all_offset += 1; 
	}
	// copy the last 293 bit
	for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + 7*(int)ALL_COLUMNS; cw_bit_index++) { // up to the last bit 
		received_signal->at(cw_bit_index + all_offset) = ((codeword[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_generator(m_rng);
	}

	// decode
	LdpcDecoder decoder = LdpcDecoder();
	std::chrono::time_point<std::chrono::system_clock> begin = std::chrono::system_clock::now();
	std::vector<bool> *decoded_symbols = decoder.decode(received_signal, std::pow(sigma_w, 2));
	std::chrono::microseconds duration = std::chrono::system_clock::now() - begin;
	std::cout << "Decoding time " << (double)duration.count()/1000 << " ms\n";

	// check
	int num_error = 0;
	for(int i = 0; i < 120; i++) {
		if(!(decoded_symbols->at(i) == infoword[i])) {
			//std::cout << "Index " << i << " " << decoded_symbols[i] << " while info_bit " << infoword[i] << " codeword " << codeword[i] << " received_signal " << received_signal->at(i) << "\n";
			num_error++;
		}
	}
	for(int i = 293; i < 30592; i++) {
		if(!(decoded_symbols->at(i) == (bool)infoword[i+173])) {
			//std::cout << "Index " << i << " " << decoded_symbols[i] << " while info_bit " << infoword[i+173] << " codeword " << codeword[i] << " received_signal " << received_signal->at(i + 173) << "\n";
			num_error++;
		}
	}
	std::cout << "BER " << (double)num_error/INFO_BIT << "\n";

	return 0;
}