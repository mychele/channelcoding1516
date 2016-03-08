#include "variable_node.h"
#include "check_node.h"
#include <iostream>
#include <chrono>
#include <random>
#include <array>

/**
 * Script to test operations on variable and check nodes
 */

int main(int argc, char const *argv[])
{
	// rng stuff
	std::random_device rd; // device entropy
    std::mt19937 m_rng(rd()); // initialize our mersenne twister with a random seed
    std::uniform_real_distribution<double> double_gen(1e-300,38);

	// ---------------------------------- check correctness of checkNodes computations in VariableNode --------------------------
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
	for(int i = 300; i < 320; i++) {
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
	int a = 1;
	int b = 0;
	VariableNode node(a, b);
	for(int i = 0; i < 7; ++i) {
		std::cout << slopes[node.getCheckNodes()[i].first] << " c " << node.getCheckNodes()[i].second << "\n";
	}
	// set channel LLR and propagate it
	node.setChannelLLR(6);
	for(int i = 0; i < PC_ROWS + 1; ++i) {
		std::cout << "LLR = " << node.getLLRat(i) << "\n";
	}
	
	// create checkNodesVector
	std::vector<CheckNode> *checkNodesVector = new std::vector<CheckNode>(2051);
	for(int slope_index = 0; slope_index < (int)PC_ROWS - 1; ++slope_index) { // init the first six blocks
		for(int c = 0; c < (int)ALL_COLUMNS - 1; ++c) { // with 292 valid nodes each
			checkNodesVector->at(slope_index*(int)ALL_COLUMNS + c).setLine(line(slope_index, c));
		}
	}
	// init the last block
	for(int c = 0; c < (int)ALL_COLUMNS; ++c) { // with 293 valid nodes 
		checkNodesVector->at(((int)PC_ROWS - 1)*(int)ALL_COLUMNS + c).setLine(line((int)PC_ROWS - 1, c));
	}

	// fill the LLR of the blocks for (1,0) -> 0+293, 293+291, 293*2+290, 293*3+289, 293*4+288, 293*5 + 287, 293*6+286, each at row 1
	checkNodesVector->at(293).setLLRat(3, 1); // this will not be summed, since this check node is not valid
	checkNodesVector->at(293+291).setLLRat(3, 1);
	checkNodesVector->at(293*2+290).setLLRat(2, 1);
	checkNodesVector->at(293*3+289).setLLRat(-4, 1);
	checkNodesVector->at(293*4+288).setLLRat(-1, 1);
	checkNodesVector->at(293*5+287).setLLRat(3, 1);
	checkNodesVector->at(293*6+286).setLLRat(0, 1);

	node.updateLLR(checkNodesVector);
	std::cout << "LLR = " << node.getLLRat(0) << " expecting 9\n";
	std::cout << "LLR = " << node.getLLRat(1) << " expecting 6\n";
	std::cout << "LLR = " << node.getLLRat(2) << " expecting 7\n";
	std::cout << "LLR = " << node.getLLRat(3) << " expecting 13\n";
	std::cout << "LLR = " << node.getLLRat(4) << " expecting 10\n";
	std::cout << "LLR = " << node.getLLRat(5) << " expecting 6\n";
	std::cout << "LLR = " << node.getLLRat(6) << " expecting 9\n";
	std::cout << "LLR = " << node.getLLRat(7) << " expecting 3\n";

	delete checkNodesVector;

	// ----------------------------------------- test the update LLR function on checkNode ----------------------------------------
	// Create VariableNode vector
	std::vector<VariableNode> *variableNodeVector = new std::vector<VariableNode>;
	variableNodeVector->resize(ALL_COLUMNS*ALL_ROWS);
	// Init its coordinates
	for(int node_index = 0; node_index < variableNodeVector->size(); ++node_index) {
															// row 						// column
		variableNodeVector->at(node_index).setCoordinates(node_index/(int)ALL_COLUMNS, node_index%(int)ALL_COLUMNS);
	}

	int line_index = 2;
	int node_c = 3;
	CheckNode check_node(line(line_index,node_c));
	double all_llr = -11;
	// init LLR on the corresponding nodes
	for(int a = 0; a < (int)ALL_ROWS; ++a) {
		variableNodeVector->at(a*ALL_COLUMNS + (slopes[line_index]*a + node_c)%ALL_COLUMNS).setLLRat(all_llr, line_index);
	}
	// set an LLR on a different branch to 0 and check that it does not change the results
	variableNodeVector->at(106*ALL_COLUMNS + (slopes[line_index]*106 + node_c)%ALL_COLUMNS).setLLRat(0, line_index+1);

	double phiTildeNode = phiTilde(std::abs(all_llr));
	double exp_llr = (all_llr > 0) ? phiTilde(111*phiTildeNode) : -phiTilde(111*phiTildeNode);
	begin = std::chrono::system_clock::now();
	check_node.updateLLR(variableNodeVector);
	std::cout << "updateLLR time " << (double)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - begin)).count()/1000 << " ms\n";
	std::cout << "Expecting " << exp_llr << " computed " << check_node.getLLRat(1) << " " <<  check_node.getLLRat(2) <<"\n";

	variableNodeVector->at(110*ALL_COLUMNS + (slopes[line_index]*110 + node_c)%ALL_COLUMNS).setLLRat(-all_llr, line_index);
	phiTildeNode = phiTilde(std::abs(all_llr));
	exp_llr = (all_llr > 0) ? -phiTilde(111*phiTildeNode) : phiTilde(111*phiTildeNode);
	begin = std::chrono::system_clock::now();
	check_node.updateLLR(variableNodeVector);
	std::cout << "updateLLR time " << (double)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - begin)).count()/1000 << " ms\n";
	std::cout << "Expecting " << exp_llr << " computed " << check_node.getLLRat(1) << " " <<  check_node.getLLRat(2) <<"\n";
	std::cout << "For node 110, expecting " << -exp_llr << " computed " << check_node.getLLRat(110) <<"\n";

	variableNodeVector->at(107*ALL_COLUMNS + (slopes[line_index]*107 + node_c)%ALL_COLUMNS).setLLRat(-all_llr, line_index);
	phiTildeNode = phiTilde(std::abs(all_llr));
	exp_llr = (all_llr > 0) ? phiTilde(111*phiTildeNode) : -phiTilde(111*phiTildeNode);
	check_node.updateLLR(variableNodeVector);
	std::cout << "Expecting " << exp_llr << " computed " << check_node.getLLRat(1) << " " <<  check_node.getLLRat(2) <<"\n";
	std::cout << "For node 110, expecting " << -exp_llr << " computed " << check_node.getLLRat(110) <<"\n";
	std::cout << "For node 107, expecting " << -exp_llr << " computed " << check_node.getLLRat(107) <<"\n";



	// set one LLR to 0 -> also the resulting llr should be 0
	variableNodeVector->at(106*ALL_COLUMNS + (slopes[line_index]*106 + node_c)%ALL_COLUMNS).setLLRat(0, line_index);
	phiTildeNode = phiTilde(std::abs(all_llr));
	exp_llr = (all_llr > 0) ? phiTilde(111*phiTildeNode) : -phiTilde(111*phiTildeNode);
	check_node.updateLLR(variableNodeVector);
	std::cout << "Expecting " << 0 << " computed " << check_node.getLLRat(1) << " " <<  check_node.getLLRat(2) <<"\n";
	std::cout << "For node 110, expecting " << 0 << " computed " << check_node.getLLRat(110) <<"\n";
	std::cout << "For node 107, expecting " << 0 << " computed " << check_node.getLLRat(107) <<"\n";
	std::cout << "For node 106, expecting " << exp_llr << " computed " << check_node.getLLRat(106) <<"\n";


	// set one LLR to infinity -> the result is given only by the other ones
	variableNodeVector->at(106*ALL_COLUMNS + (slopes[line_index]*106 + node_c)%ALL_COLUMNS).setLLRat(std::numeric_limits<double>::infinity(), line_index);
	phiTildeNode = phiTilde(std::abs(all_llr));
	exp_llr = phiTilde(110*phiTildeNode);
	check_node.updateLLR(variableNodeVector);
	std::cout << "Expecting " << exp_llr << " computed " << check_node.getLLRat(1) << " " <<  check_node.getLLRat(2) <<"\n";
	std::cout << "For node 110, expecting " << -exp_llr << " computed " << check_node.getLLRat(110) <<"\n";
	std::cout << "For node 107, expecting " << -exp_llr << " computed " << check_node.getLLRat(107) <<"\n";
	std::cout << "For node 106, expecting " << ((all_llr > 0) ? phiTilde(111*phiTildeNode) : -phiTilde(111*phiTildeNode)) << " computed " << check_node.getLLRat(106) <<"\n";


	for(int a = 0; a < (int)ALL_ROWS; ++a) {
		variableNodeVector->at(a*ALL_COLUMNS + (slopes[line_index]*a + node_c)%ALL_COLUMNS).setLLRat(std::numeric_limits<double>::infinity(), line_index);
	}
	phiTildeNode = phiTilde(std::abs(std::numeric_limits<double>::infinity()));
	exp_llr = (std::numeric_limits<double>::infinity() > 0) ? phiTilde(111*phiTildeNode) : -phiTilde(111*phiTildeNode);
	check_node.updateLLR(variableNodeVector);
	std::cout << "Expecting " << exp_llr << " computed " << check_node.getLLRat(1) << " " <<  check_node.getLLRat(2) <<"\n";

	// ----------------------------------------- test the min sum update LLR function on checkNode ----------------------------------------

	std::cout << "Min Sum test\n";
	line_index = 2;
	node_c = 3;
	CheckNode check_node_ms = CheckNode(line(line_index,node_c));
	all_llr = -11;
	// init LLR on the corresponding nodes
	for(int a = 0; a < (int)ALL_ROWS; ++a) {
		variableNodeVector->at(a*ALL_COLUMNS + (slopes[line_index]*a + node_c)%ALL_COLUMNS).setLLRat(all_llr, line_index);
	}
	// set an LLR on a different branch to 0 and check that it does not change the results
	variableNodeVector->at(106*ALL_COLUMNS + (slopes[line_index]*106 + node_c)%ALL_COLUMNS).setLLRat(0, line_index+1);

	begin = std::chrono::system_clock::now();
	check_node_ms.updateLLRminSum(variableNodeVector);
	std::cout << "updateLLR time " << (double)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - begin)).count()/1000 << " ms\n";
	std::cout << "Expecting " << all_llr << " computed " << check_node_ms.getLLRat(1) << " " <<  check_node_ms.getLLRat(2) <<"\n";

	variableNodeVector->at(110*ALL_COLUMNS + (slopes[line_index]*110 + node_c)%ALL_COLUMNS).setLLRat(-all_llr, line_index);
	exp_llr = -all_llr;
	begin = std::chrono::system_clock::now();
	check_node_ms.updateLLRminSum(variableNodeVector);
	std::cout << "updateLLR time " << (double)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - begin)).count()/1000 << " ms\n";
	std::cout << "Expecting " << exp_llr << " computed " << check_node_ms.getLLRat(1) << " " <<  check_node_ms.getLLRat(2) <<"\n";
	std::cout << "For node 110, expecting " << -exp_llr << " computed " << check_node_ms.getLLRat(110) <<"\n";


	variableNodeVector->at(107*ALL_COLUMNS + (slopes[line_index]*107 + node_c)%ALL_COLUMNS).setLLRat((all_llr > 0 ? -(all_llr - 1) : -(all_llr + 1)), line_index);
	exp_llr = (all_llr > 0) ? all_llr-1 : -(all_llr-1);
	check_node_ms.updateLLRminSum(variableNodeVector);
	std::cout << "Expecting " << (all_llr > 0 ? (all_llr - 1) : (all_llr + 1)) << " computed " << check_node_ms.getLLRat(1) << " " <<  check_node_ms.getLLRat(2) <<"\n";
	std::cout << "For node 110, expecting " << (all_llr < 0 ? -(all_llr + 1) : -(all_llr - 1))   << " computed " << check_node_ms.getLLRat(110) <<"\n";
	std::cout << "For node 107, expecting " << -all_llr << " computed " << check_node_ms.getLLRat(107) <<"\n";



	// set one LLR to 0 -> also the resulting llr should be 0
	variableNodeVector->at(106*ALL_COLUMNS + (slopes[line_index]*106 + node_c)%ALL_COLUMNS).setLLRat(0, line_index);
	check_node_ms.updateLLRminSum(variableNodeVector);
	std::cout << "Expecting " << 0 << " computed " << check_node_ms.getLLRat(1) << " " <<  check_node_ms.getLLRat(2) <<"\n";
	std::cout << "For node 110, expecting " << 0 << " computed " << check_node_ms.getLLRat(110) <<"\n";
	std::cout << "For node 107, expecting " << 0 << " computed " << check_node_ms.getLLRat(107) <<"\n";
	std::cout << "For node 106, expecting " << (all_llr > 0 ? (all_llr - 1) : (all_llr + 1)) << " computed " << check_node_ms.getLLRat(106) <<"\n";


	// set one LLR to infinity -> the result is given only by the other ones
	variableNodeVector->at(106*ALL_COLUMNS + (slopes[line_index]*106 + node_c)%ALL_COLUMNS).setLLRat(std::numeric_limits<double>::infinity(), line_index);
	check_node_ms.updateLLRminSum(variableNodeVector);
	std::cout << "Expecting " << (all_llr > 0 ? (all_llr - 1) : -(all_llr + 1)) << " computed " << check_node_ms.getLLRat(1) << " " <<  check_node_ms.getLLRat(2) <<"\n";
	std::cout << "For node 110, expecting " << -(all_llr > 0 ? (all_llr - 1) : -(all_llr + 1)) << " computed " << check_node_ms.getLLRat(110) <<"\n";
	std::cout << "For node 107, expecting " << (all_llr > 0 ? -all_llr : all_llr) << " computed " << check_node_ms.getLLRat(107) <<"\n";
	std::cout << "For node 106, expecting " << (all_llr > 0 ? (all_llr - 1) : (all_llr + 1)) << " computed " << check_node_ms.getLLRat(106) <<"\n";


	for(int a = 0; a < (int)ALL_ROWS; ++a) {
		variableNodeVector->at(a*ALL_COLUMNS + (slopes[line_index]*a + node_c)%ALL_COLUMNS).setLLRat(std::numeric_limits<double>::infinity(), line_index);
	}
	check_node_ms.updateLLRminSum(variableNodeVector);
	std::cout << "Expecting " << std::numeric_limits<double>::infinity() << " computed " << check_node_ms.getLLRat(1) << " " <<  check_node_ms.getLLRat(2) <<"\n";
	
	return 0;
}