#include "ldpc_decoder.h"
#include "variable_node.h"
#include "check_node.h"
#include <vector>
#include <random>


int main(int argc, char const *argv[])
{
	// rng stuff
	std::random_device rd; // device entropy
    std::mt19937 m_rng(rd()); // initialize our mersenne twister with a random seed
    std::uniform_real_distribution<double> double_gen(-10,10);

	// create a LdpcDecoder object
	LdpcDecoder decoder = LdpcDecoder();
	std::vector<double> *received_signal = new std::vector<double>;
	received_signal->resize(ALL_COLUMNS*ALL_ROWS);
	std::cout << "Fill vector\n";
	for(std::vector<double>::iterator iter = received_signal->begin(); iter != received_signal->end(); ++iter) {
		(*iter) = double_gen(m_rng);
	}
	std::cout << "Start decoder\n";
	decoder.decode(received_signal, 1);

	// const std::vector<VariableNode> *variableNodeVector = decoder.getVariableNodeVector();

	// //------------------------------------------- check if nodes that must have LLR = inf have LLR = inf--------------------------
	// for(int v_index = 120; v_index < 293; v_index++) {
	// 	if(!isinf(variableNodeVector->at(v_index).getLLR())) 
	// 		{std::cout << "Index " << v_index << " value " << variableNodeVector->at(v_index).getLLR();}
	// }
	// for(int r_index = 105; r_index < 111; r_index++) {
	// 	int v_index = r_index*293 + 292;
	// 	if(!isinf(variableNodeVector->at(v_index).getLLR())) 
	// 		{std::cout << "Index " << v_index << " value " << variableNodeVector->at(v_index).getLLR();}
	// }
	// std::cout << variableNodeVector->at(6536).getLLR() << " " << -received_signal->at(6536)*2 << "\n";
	// std::cout << variableNodeVector->at(variableNodeVector->size()-1).getLLR() << " " << -received_signal->at(32815)*2 << "\n";

	return 0;
	
}