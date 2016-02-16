#include "variable_node.h"
#include "check_node.h"
#include <iostream>
#include <chrono>
#include <random>
#include <array>


int main(int argc, char const *argv[])
{
	// rng stuff
	std::random_device rd; // device entropy
    std::mt19937 m_rng(rd()); // initialize our mersenne twister with a random seed
    std::uniform_real_distribution<double> double_gen(1e-300,38);

	// ---------------------------------- check correctness of checkNodes computations in variableNodes --------------------------
	for(int a = 0; a < 111; ++a) {
		for(int b = 0; b < 293; ++b) {
			VariableNode node(a, b);
			for(int i = 0; i < 7; ++i) { // a*slope +    c  		
				if( ((a*slopes[node.getCheckNodes()[i].first]  + node.getCheckNodes()[i].second)%293 != b) || node.getCheckNodes()[i].second < 0) {
					if (!(node.getCheckNodes()[i].second < 0 && reverseModulo((int)ALL_COLUMNS, a*slopes[i], b) == 292)) { // c = 292 -> c = -1
						std::cout << "a = " << a << " b = " << b << " s_index = " << node.getCheckNodes()[i].first << " c = " << node.getCheckNodes()[i].second << "\n";
					}
				}

			}
		}
	}

	// ------------------------------------------- test phiTilde function for large values ----------------------------------------
	for(int i = 30; i < 40; ++i) {
		std::cout << "i = " << i << " phiTilde = " << phiTilde(i) << " isNaN? " << isnan(phiTilde(i)) << "\n";
	}

	// ----------------------------------------- test phiTilde function for very small values -------------------------------------
	for(int i = 280; i < 320; i++) {
		std::cout << "10^{-" << i << "} = " << std::pow(10, -i) << " phiTilde = " << phiTilde(std::pow(10, -i)) << "\n";
	}

	// ------------------------------------------------------ test infinity sum ---------------------------------------------------
	std::cout << "Expecting inf -> " << phiTilde(std::pow(10, -400)) + 2 << "\n";
	std::cout << "Expecting -inf -> " << -phiTilde(std::pow(10, -400)) + 2 << "\n";
	std::cout << "Expecting 2 -> " <<  phiTilde(phiTilde(std::pow(10, -400))) + 2 << "\n";

	// ---------------------------------------------------- phiTilde performances -------------------------------------------------
	const int N = 1e6;
	// create N random values
	std::array<double, N> random_array;
	for(int i = 0; i < N; ++i) {
		random_array[i] = double_gen(m_rng);
	}

	std::chrono::time_point<std::chrono::system_clock> begin = std::chrono::system_clock::now();
	double val;
	for(int i = 0; i < N; ++i) {
		val = phiTilde(random_array[i]);
	}
	std::chrono::nanoseconds duration = std::chrono::system_clock::now() - begin;
	std::cout << "Average time to compute phiTilde(x) " << (double)duration.count()/N << " ns\n";

	// ----------------------------------------- test the update LLR function on variableNode --------------------------------------
	int a = 10;
	int b = 20;
	VariableNode node(a, b);
	node.setLLR(3.4);
	std::cout << "LLR = " << node.getLLR() << "\n";
	for(int i = 0; i < 7; ++i) {
		std::cout << slopes[node.getCheckNodes()[i].first] << " c " << node.getCheckNodes()[i].second << "\n";
	}
	// create checkNodesVector
	std::vector<CheckNode> *checkNodesVector = new std::vector<CheckNode>;
	checkNodesVector->resize(2051);
	for(int slope_index = 0; slope_index < (int)PC_ROWS - 1; ++slope_index) { // fill the first six blocks
		for(int c = 0; c < (int)ALL_COLUMNS - 1; ++c) { // with 292 valid nodes each
			checkNodesVector->at(slope_index*(int)ALL_COLUMNS + c).setLine(line(slope_index, c));
		}
	}
	// fill the last block
	for(int c = 0; c < (int)ALL_COLUMNS; ++c) { // with 293 valid nodes 
		checkNodesVector->at(((int)PC_ROWS - 1)*(int)ALL_COLUMNS + c).setLine(line((int)PC_ROWS - 1, c));
	}

	// init the LLR of the blocks for (10,20) -> 0+10, 293+0, 293*2+283, 293*3+273, 293*4+263, 293*5 + 253, 293*6+243
	checkNodesVector->at(10).setLLR(3);
	checkNodesVector->at(293).setLLR(3);
	checkNodesVector->at(293*2+283).setLLR(4);
	checkNodesVector->at(293*3+273).setLLR(3);
	checkNodesVector->at(293*4+263).setLLR(3);
	checkNodesVector->at(293*5+253).setLLR(3);
	checkNodesVector->at(293*6+243).setLLR(0);

	node.updateLLR(checkNodesVector);
	std::cout << "LLR = " << node.getLLR() << " expecting 19 " << "\n";

	delete checkNodesVector;

	// ----------------------------------------- test the update LLR function on checkNode ----------------------------------------
	// TODO


	return 0;
}