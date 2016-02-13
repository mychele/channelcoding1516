#include <iostream>
#include <random>
#include <chrono>


int main(int argc, char const *argv[])
{
	int num_samples = 32816;
	int repetitions = 1000;

	// rng stuff
	std::random_device rd; // device entropy
    std::mt19937 m_rng(rd()); // initialize our mersenne twister with a random seed
    std::uniform_int_distribution<int> int_uni_gen(0,1);
           
    std::vector< std::chrono::microseconds > execution_times;

	std::chrono::time_point<std::chrono::system_clock> beginning;
    // generate num_samples uniform for 100 times
    for (int i = 0; i < repetitions; i++) {
    	beginning = std::chrono::system_clock::now();
    	int array[num_samples];
    	for(int s = 0; s < num_samples; s++) {
    		array[s] = int_uni_gen(m_rng);
    	}
    	execution_times.push_back(std::chrono::system_clock::now() - beginning);
    }

    // stats
    double total_duration;
    for(std::vector< std::chrono::microseconds >::iterator iter = execution_times.begin(); iter != execution_times.end(); iter++) {
    	total_duration += iter->count();
    }

    std::cout << "Average time to generate 32816 uniform_int (0,1) samples " << total_duration/(repetitions * 1000) << " ms\n";

    std::normal_distribution<double> norm_gen(0,2);
           
    execution_times.clear();

    // generate num_samples uniform for 100 times
    for (int i = 0; i < repetitions; i++) {
    	beginning = std::chrono::system_clock::now();
    	double array[num_samples];
    	for(int s = 0; s < num_samples; s++) {
    		array[s] = norm_gen(m_rng);
    	}
    	execution_times.push_back(std::chrono::system_clock::now() - beginning);
    }

    // stats
    total_duration = 0;
    for(std::vector< std::chrono::microseconds >::iterator iter = execution_times.begin(); iter != execution_times.end(); iter++) {
    	total_duration += iter->count();
    }

    std::cout << "Average time to generate 32816 normal_double samples " << total_duration/(repetitions * 1000) << " ms\n";

	/* code */
	return 0;
}