#include "ldpc_common.h"
#include "check_node.h"
#include <iostream>

#ifndef VARIABLE_NODE
#define VARIABLE_NODE

class VariableNode
{
public:
	VariableNode(int a, int b);

	void updateLLR(std::vector<CheckNode> *checkNodeVector);

	void setLLR(double llr);

	void setCoordinates(int a, int b);

	double getLLR() const;

	std::pair<int, int> getCoordinates() const;

	std::vector<line> getCheckNodes() const;

private:
	double m_llr;
	int m_j;
	int m_a;
	int m_b;
	std::vector<line> m_checkNodes;
	
};

#endif
// VARIABLE_NODE