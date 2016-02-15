#include "variable_node.h"
#include <iostream>


int main(int argc, char const *argv[])
{
	// check correctness of checkNodes computations in variableNodes
	for(int a = 0; a < 111; ++a) {
		for(int b = 0; b < 293; ++b) {
			VariableNode node(a, b);
			for(int i = 0; i < 7; ++i) { // a*slope +    c  		
				if( ((a*node.getCheckNodes()[i].first  + node.getCheckNodes()[i].second)%293 != b) || node.getCheckNodes()[i].second < 0) {
					std::cout << "a = " << a << " b = " << b << " s = " << node.getCheckNodes()[i].first << " c = " << node.getCheckNodes()[i].second << "\n";
				}

			}
		}
	}

	return 0;
}