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
    std::mt19937 m_rng(10); // initialize our mersenne twister with a random seed
    std::uniform_int_distribution<int> bit_generator(0,1);
    std::normal_distribution<double> noise_generator(0,1);

    const int N = 50; // attempts
    const int num_SNR = 2;
    double ebn0_vec[num_SNR] = {6.6, 7.2};
    double BER[N][num_SNR] = {{0}};
    double decoding_time[N][num_SNR] = {{0}};
    std::vector<double> decoding_time_nocw;

	LdpcEncoder encoder = LdpcEncoder();
	std::cout << "LdpcEncoder created\n";
	// import matrix
	if(encoder.setup() == 0) {std::cout << "File opened and read successfully\n";}
	else {std::cout << "Error in opening matrix file, abort\n"; return 2;}

    LdpcDecoder decoder = LdpcDecoder();

    for(int attempt = 0; attempt < N; ++attempt) {
    	std::cout << attempt << "\n";
    	//create information word
		InfoWord infoword;
		for(int bit_index = 0; bit_index < (int)INFO_BIT; ++bit_index) { 
			infoword[mapJtoK(bit_index)] = (bit_generator(m_rng)==1);
		}

		// encode
		CodeWord codeword = encoder.encode(infoword);

		// generate noise N(0,1)
		std::vector<double> noise_vector;
		noise_vector.resize(ALL_COLUMNS*ALL_ROWS);
		for(int noise_index = 0; noise_index < noise_vector.size(); noise_index++) {
			noise_vector[noise_index] = noise_generator(m_rng);
		}

		for(int snr_ind = 0; snr_ind < num_SNR; ++snr_ind) {
			double sigma_w = 1/std::sqrt(2*std::pow(10, (ebn0_vec[snr_ind]/10)) * INFO_BIT/CODE_WORD);

			std::vector<double> *received_signal = new std::vector<double>;
			received_signal->resize(ALL_COLUMNS*ALL_ROWS, -1);

			// create the received signal

			// copy the first 120 bit of the codeword, then skip 173 zeros
			int cw_bit_index = 0;
			for(; cw_bit_index < ALL_COLUMNS - INIT_ZERO_BIT; cw_bit_index++) {
				received_signal->at(cw_bit_index) = ((codeword[cw_bit_index] == 0) ? -1:1)  + sigma_w*noise_vector[cw_bit_index];
			}
			int all_offset = (int)INIT_ZERO_BIT;
			// copy up to						30765 - 1, which is bit 30592 - 1 in the codeword
			for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT; cw_bit_index++) { 
				received_signal->at(cw_bit_index + all_offset) = ((codeword[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_vector[cw_bit_index + all_offset];
			}
			// skip 3 zeros in the codeword, but not in all_bit
			cw_bit_index += 3;
			all_offset -= 3;
			for(int i = 1; i < 7; ++i) {
				for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + (i)*(int)ALL_COLUMNS - 1; cw_bit_index++) { 
				// skip bit 31057, 31350, 31643, 31936, 32229, 32522 in all_bits
					received_signal->at(cw_bit_index + all_offset) = ((codeword[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_vector[cw_bit_index + all_offset];
				}
				all_offset += 1; 
			}
			// copy the last 293 bit
			for(; cw_bit_index + all_offset < (int)ALL_INFO_BIT + 7*(int)ALL_COLUMNS; cw_bit_index++) { // up to the last bit 
				received_signal->at(cw_bit_index + all_offset) = ((codeword[cw_bit_index] == 0) ? -1:1) + sigma_w*noise_vector[cw_bit_index + all_offset];
			}

			// decode
	
			std::chrono::time_point<std::chrono::system_clock> begin = std::chrono::system_clock::now();
			std::vector<bool> *decoded_symbols = decoder.decode(received_signal, std::pow(sigma_w,2));
			std::chrono::nanoseconds duration = std::chrono::system_clock::now() - begin;
			decoding_time[attempt][snr_ind] = duration.count();

			// check
			int num_error = 0;
			for(int i = 0; i < 120; i++) {
				if(!(decoded_symbols->at(i) == infoword[i])) {
					//std::cout << "Index " << i << " " << decoded_symbols[i] << " while info_bit " << infoword[i] << " codeword " << codeword[i] << " received_signal " << received_signal->at(i) << "\n";
					num_error++;
				}
			}
			for(int i = 293; i < 30765; i++) {
				if(!(decoded_symbols->at(i) == (bool)infoword[i])) {
					//std::cout << "Index " << i << " " << decoded_symbols[i] << " while info_bit " << infoword[i+173] << " codeword " << codeword[i] << " received_signal " << received_signal->at(i + 173) << "\n";
					num_error++;
				}
			}

			BER[attempt][snr_ind] = (double)num_error/INFO_BIT;
			std::cout << "snr = " << ebn0_vec[snr_ind] << " BER = " << BER[attempt][snr_ind] << "\n";
			if(num_error > 0) {
				decoding_time_nocw.push_back(duration.count());
			}
		}
    }

    // elaborate BER
    double mean_BER[num_SNR] = {0};
    double mean_decoding[num_SNR] = {0};
    for (int snr_ind = 0; snr_ind < num_SNR; snr_ind++) {
    	double BER_sum = 0;
    	double time_sum = 0;
    	for(int i = 0; i < N; ++i) {
    		BER_sum += BER[i][snr_ind];
    		time_sum += decoding_time[i][snr_ind];
    		std::cout << BER[i][snr_ind] << " ";
    	}
    	std::cout << "\n";
    	mean_BER[snr_ind] = BER_sum/N;
    	mean_decoding[snr_ind] = time_sum/N;
    }

    double mean_nocw = 0;
    for(int i = 0; i < decoding_time_nocw.size(); ++i) {
    	mean_nocw += decoding_time_nocw[i];
    }
    mean_nocw = mean_nocw/decoding_time_nocw.size();

    for (int snr_ind = 0; snr_ind < num_SNR; ++snr_ind)
    {
    	std::cout << "eb/N0 " << ebn0_vec[snr_ind] << " mean BER " << mean_BER[snr_ind] << "\n";
    }
    for (int snr_ind = 0; snr_ind < num_SNR; ++snr_ind)
    {
    	std::cout << "decoding time (avg) for eb/N0 " << ebn0_vec[snr_ind] << " is " << mean_decoding[snr_ind]/1e6 << "ms\n";
    }

    std::cout << "in case cw is not found " << mean_nocw/1e6 << "ms\n";
	return 0;
}