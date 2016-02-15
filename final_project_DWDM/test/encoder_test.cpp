#include "ldpc_encoder.h"
#include <bitset>
#include <iostream>
#include <random>
#include <chrono>


int main(int argc, char const *argv[])
{
	// test the LdpcEncoder class

	// rng stuff
	std::random_device rd; // device entropy
    std::mt19937 m_rng(rd()); // initialize our mersenne twister with a random seed
    std::uniform_int_distribution<int> int_uni_gen(0,1);


	// create LdpcEncoder object
	std::chrono::time_point<std::chrono::system_clock> tic = std::chrono::system_clock::now(); 
	LdpcEncoder encoder = LdpcEncoder();
	std::cout << "LdpcEncoder created\n";
	// import matrix
	if(encoder.setup() == 0) {std::cout << "File opened and read successfully\n";}
	else {std::cout << "Error in opening matrix file, abort\n"; return 2;}
	std::cout << "Elapsed time to read K matrix " << (std::chrono::system_clock::now() - tic).count()/1000 << " ms\n";

	// create a random information word, encode it and test it against the standard specifications
	// repeat the procedure N times
	int N = 1000;
	int tot_line_ok = 0;
	std::vector< std::chrono::microseconds > encoding_time;
	std::chrono::time_point<std::chrono::system_clock> begin_encoding;
	std::vector< std::chrono::microseconds > unif_generation_time;
	std::chrono::time_point<std::chrono::system_clock> begin_info_generation;
	for (int attempt = 0; attempt < N; ++attempt) {
		// generate random uncoded word
		begin_info_generation = std::chrono::system_clock::now();
		InfoWord infoword;
		for(int bit_index = 0; bit_index < (int)INFO_BIT; ++bit_index) { 
			infoword[mapJtoK(bit_index)] = (int_uni_gen(m_rng)==1);
		}
		unif_generation_time.push_back(std::chrono::system_clock::now() - begin_info_generation);

		// encode it
		// std::cout << "Encode\n";
		begin_encoding = std::chrono::system_clock::now();
		CodeWord codeword = encoder.encode(infoword);
		encoding_time.push_back(std::chrono::system_clock::now() - begin_encoding);

		// std::cout << "Test\n";
		// test the encoded word by filling the standard matrix and checking every line
		// for practical reasons, create a 112*293 bitset and fill it with the correct bits
		std::bitset<(int)ALL_COLUMNS*((int)INFO_ROWS + (int)PC_ROWS)> all_bits;

		// copy the first 120 bit of the codeword, then skip 173 zeros
		int cw_bit_index = 0;
		for(; cw_bit_index < (int)ALL_COLUMNS-(int)INIT_ZERO_BIT; cw_bit_index++) {
			all_bits[cw_bit_index] = codeword[cw_bit_index];
		}
		int all_offset = (int)INIT_ZERO_BIT;
		// copy up to						30765 - 1, which is bit 30592 - 1 in the codeword
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT; cw_bit_index++) { 
			all_bits[cw_bit_index + all_offset] = codeword[cw_bit_index];
		}
		// skip 3 zeros in the codeword, but not in all_bit
		cw_bit_index += 3;
		all_offset -= 3;
		// copy 292 bit
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + (int)ALL_COLUMNS - 1; cw_bit_index++) { // skip bit 31057 in all_bits
			all_bits[cw_bit_index + all_offset] = codeword[cw_bit_index];
		}

		all_offset += 1; 
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + 2*(int)ALL_COLUMNS - 1; cw_bit_index++) { // skip bit 31350 in all_bits
			all_bits[cw_bit_index + all_offset] = codeword[cw_bit_index];
		}

		all_offset += 1;
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + 3*(int)ALL_COLUMNS - 1; cw_bit_index++) { // skip bit 31643 in all_bits
			all_bits[cw_bit_index + all_offset] = codeword[cw_bit_index];
		}

		all_offset += 1;
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + 4*(int)ALL_COLUMNS - 1; cw_bit_index++) { // skip bit 31936 in all_bits
			all_bits[cw_bit_index + all_offset] = codeword[cw_bit_index];
		}

		all_offset += 1; 
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + 5*(int)ALL_COLUMNS - 1; cw_bit_index++) { // skip bit 32229 in all_bits
			all_bits[cw_bit_index + all_offset] = codeword[cw_bit_index];
		}

		all_offset += 1; 
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + 6*(int)ALL_COLUMNS - 1; cw_bit_index++) { // skip bit 32522 in all_bits
			all_bits[cw_bit_index + all_offset] = codeword[cw_bit_index];
		}

		all_offset += 1; 
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + 7*(int)ALL_COLUMNS; cw_bit_index++) { // up to the last bit 
			all_bits[cw_bit_index + all_offset] = codeword[cw_bit_index];
		}

		// fill the matrix
		int std_M[INFO_ROWS + PC_ROWS][ALL_COLUMNS];
		for(int bit_index = 0; bit_index < all_bits.size(); bit_index++) {
			int row_index = bit_index/((int)ALL_COLUMNS); // floor(bit_index/293) since both q and 293 are positive int, the result is the floor
			int col_index = bit_index % ((int)ALL_COLUMNS); // (bit_index%293)
			std_M[row_index][col_index] = all_bits[bit_index];
		}

		// check if the 2051 lines sum to 0
		int col_r = (int)ALL_COLUMNS;
		bool line_ok = 1;
		for(int slope_index = 0; slope_index < (int)PC_ROWS; slope_index++) {
			// cycle on the offsets
			for(int c = 0; c < (int)col_r; c++) {
				// cycle on a from 0 to 111
				int lineValue = 0;
				for(int a = 0; a < (int)PC_ROWS + (int)INFO_ROWS; a++) {
					int col_index = (a*slopes[slope_index] + c)%col_r;
					if(std_M[a][col_index] == 1) {lineValue++;}
				}
				if(lineValue%2 != 0) {std::cout << "s " << slopes[slope_index] << " c " << c << "\n"; line_ok = 0;}
			}
		}
		tot_line_ok+=line_ok;
		// std::cout << "Encoding " << ((line_ok) ? "is":"is not") << " successfull!\n";

	}

	// print stats
	std::cout << "Correct encoding in " << (double)tot_line_ok/N*100 << "% of the cases\n";
    double total_duration_generation;
    for(auto iter : unif_generation_time) {
    	total_duration_generation += iter.count();
    }
    std::cout << "Average time to generate an infoword (30592 uniform_int (0,1) samples) " << total_duration_generation/(N * 1000) << " ms\n";

    double total_duration_encoding;
    for(auto iter : encoding_time) {
    	total_duration_encoding += iter.count();
    }
    std::cout << "Average time to encode an infoword " << total_duration_encoding/(N * 1000) << " ms\n";



	return 0;
}