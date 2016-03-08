#include "ldpc_decoder.h"
#include "ldpc_encoder.h"
#include <sstream>
#include <string>
#include <bitset>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <thread>
#include <mutex>

std::mutex lock_variable;


void simulate(int snr_ind, double sigma_w, LdpcDecoder *decoder, 
	std::vector< std::vector<double> > *num_error_matrix, std::vector< std::vector<double> > *decoding_time,
	std::vector<double> *noise_vector, CodeWord *codeword, InfoWord *infoword) 
{
	std::vector<double> *received_signal = new std::vector<double>(ALL_COLUMNS*ALL_ROWS, -1);

	// create the received signal

	// copy the first 120 bit of the codeword, then skip 173 zeros
	int cw_bit_index = 0;
	for(; cw_bit_index < ALL_COLUMNS - INIT_ZERO_BIT; cw_bit_index++) {
		received_signal->at(cw_bit_index) = (((*codeword)[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_vector->at(cw_bit_index);
	}
	int all_offset = (int)INIT_ZERO_BIT;
	// copy up to						30765 - 1, which is bit 30592 - 1 in the codeword
	for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT; cw_bit_index++) { 
		received_signal->at(cw_bit_index + all_offset) = (((*codeword)[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_vector->at(cw_bit_index + all_offset);
	}
	// skip 3 zeros in the codeword, but not in all_bit
	cw_bit_index += 3;
	all_offset -= 3;
	for(int i = 1; i < 7; ++i) {
		for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + (i)*(int)ALL_COLUMNS - 1; cw_bit_index++) { 
		// skip bit 31057, 31350, 31643, 31936, 32229, 32522 in all_bits
			received_signal->at(cw_bit_index + all_offset) = (((*codeword)[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_vector->at(cw_bit_index + all_offset);
		}
		all_offset += 1; 
	}
	// copy the last 293 bit
	for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + 7*(int)ALL_COLUMNS; cw_bit_index++) { // up to the last bit 
		received_signal->at(cw_bit_index + all_offset) = (((*codeword)[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_vector->at(cw_bit_index + all_offset);
	}

	// decode
	std::chrono::time_point<std::chrono::system_clock> begin = std::chrono::system_clock::now();
	std::vector<bool> *decoded_symbols = decoder->decode(received_signal, std::pow(sigma_w,2));
	std::chrono::nanoseconds duration = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - begin);

	// check
	int num_error = 0;
	for(int i = 0; i < 120; i++) {
		if(!(decoded_symbols->at(i) == (*infoword)[i])) {
			//std::cout << "Index " << i << " " << decoded_symbols[i] << " while info_bit " << infoword[i] << " codeword " << codeword[i] << " received_signal " << received_signal->at(i) << "\n";
			num_error++;
		}
	}
	for(int i = 293; i < 30765; i++) {
		if(!(decoded_symbols->at(i) == (bool)(*infoword)[i])) {
			//std::cout << "Index " << i << " " << decoded_symbols[i] << " while info_bit " << infoword[i+173] << " codeword " << codeword[i] << " received_signal " << received_signal->at(i + 173) << "\n";
			num_error++;
		}
	}

	delete received_signal;

	// write operation: lock and unlock
	lock_variable.lock();
	num_error_matrix->at(snr_ind).push_back(num_error);
	decoding_time->at(snr_ind).push_back(duration.count());
	//std::cout << "snr = " << sigma_w << " num_error_matrix = " << num_error_matrix->at(snr_ind).back() << "\n";

	lock_variable.unlock();
}

/**
 * Script to launch a simulation with encoding
 */

int main(int argc, char const *argv[])
{
	// rng stuff
	std::random_device rd; // device entropy
    std::mt19937 m_rng(rd()); // initialize our mersenne twister with a random seed
    std::uniform_int_distribution<int> bit_generator(0,1);
    std::normal_distribution<double> noise_generator(0,1);

    const int num_SNR = 6;
    double ebn0_vec[num_SNR] = {2,3,4,5,6,6.5};//,6.7, 6.75, 6.8, 6.85, 6.9,  6.95,  7,    7.05, 	7.1};
    const int num_attempt_per_snr[num_SNR] = {100,100,100,100,100,100};//,100, 1000, 1000, 1000, 1000, 10000, 10000, 10000, 10000};

    const int N = 100; // attempts
    std::vector< std::vector<double> > *num_error_matrix = new std::vector< std::vector<double> >(num_SNR);
    std::vector< std::vector<double> > *decoding_time = new std::vector< std::vector<double> >(num_SNR);
    std::vector< std::thread > snr_thread_vec;

	LdpcEncoder encoder = LdpcEncoder();
	std::cout << "LdpcEncoder created\n";
	// import matrix
	if(encoder.setup() == 0) {std::cout << "File opened and read successfully\n";}
	else {std::cout << "Error in opening matrix file, abort\n"; return 2;}

    std::vector< LdpcDecoder > *decoder_vector = new std::vector<LdpcDecoder>(num_SNR); //one decoder per SNR
    std::cout << "Decoders initialized\n";

    for(int attempt = 0; attempt < N; ++attempt) {
    	std::cout << attempt << "\n";

    	// Create information word
		InfoWord infoword;
		for(int bit_index = 0; bit_index < (int)INFO_BIT; ++bit_index) { 
			infoword[mapJtoK(bit_index)] = (bit_generator(m_rng)==1);
		}

		// Encode
		CodeWord codeword = encoder.encode(infoword);

		// Generate the same noise for all the SNR, with random N(0,1) variables
		std::vector<double> *noise_vector = new std::vector<double>(ALL_BIT);
		for(int noise_index = 0; noise_index < noise_vector->size(); noise_index++) {
			noise_vector->at(noise_index) = noise_generator(m_rng);
		}

		for(int snr_ind = 0; snr_ind < num_SNR; ++snr_ind) {
			if(num_attempt_per_snr[snr_ind] >= attempt) { 

				// ---------------------------------------- add noise and decode --------------------------------------
				double sigma_w = 1/std::sqrt(2*std::pow(10, (ebn0_vec[snr_ind]/10)) * INFO_BIT/CODE_WORD);
				snr_thread_vec.push_back(std::thread(simulate, snr_ind, sigma_w, &(decoder_vector->at(snr_ind)), num_error_matrix, decoding_time, noise_vector, &codeword, &infoword));
			}
		}

		for(std::vector<std::thread>::iterator th_iter = snr_thread_vec.begin(); th_iter != snr_thread_vec.end(); ++th_iter) {
			th_iter->join();
		}
		// remove the threads that just joined
		snr_thread_vec.clear();

		if((attempt+1)%100 == 0) {
			for(int snr_write = 0; snr_write < num_SNR; snr_write++) {
				if(num_attempt_per_snr[snr_write] >= attempt + 1) {
					std::stringstream conc_filename;
					conc_filename << "simulation_results_multi_3_" << ebn0_vec[snr_write] << ".txt";
					std::ofstream output_file (conc_filename.str().c_str(), std::ios::out | std::ios::app);
					for(int attempt_index = attempt + 1 - 100; attempt_index < attempt + 1; attempt_index++) {
						output_file << num_error_matrix->at(snr_write).at(attempt_index) << "\n";
					}
					output_file.close();
				}
			}
		}
		delete noise_vector;
    }

    // elaborate num_error_matrix
    double mean_BER[num_SNR] = {0};
    double mean_decoding[num_SNR] = {0};
    for (int snr_ind = 0; snr_ind < num_SNR; snr_ind++) {
    	double BER_sum = 0;
    	double time_sum = 0;
    	for(int i = 0; i < num_attempt_per_snr[snr_ind]; ++i) {
    		BER_sum += num_error_matrix->at(snr_ind).at(i)/INFO_BIT;
    		time_sum += decoding_time->at(snr_ind).at(i);
    	}
    	mean_BER[snr_ind] = BER_sum/num_attempt_per_snr[snr_ind];
    	mean_decoding[snr_ind] = time_sum/num_attempt_per_snr[snr_ind];
    }

	
	std::ofstream output_file ("BER_and_time_multi_3.txt", std::ios::out | std::ios::app);
    for (int snr_ind = 0; snr_ind < num_SNR; ++snr_ind)
    {
    	output_file << ebn0_vec[snr_ind] << " " << mean_BER[snr_ind] << "\n";
    }
    for (int snr_ind = 0; snr_ind < num_SNR; ++snr_ind)
    {
    	output_file << ebn0_vec[snr_ind] << " " << mean_decoding[snr_ind] << "\n";
    }
    output_file.close();
	return 0;
}