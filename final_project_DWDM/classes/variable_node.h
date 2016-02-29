#ifndef VARIABLE_NODE
#define VARIABLE_NODE

#include "ldpc_common.h"
#include "check_node.h"
#include <iostream>

class CheckNode;

class VariableNode
{
public:
	VariableNode();
	VariableNode(int a, int b);
	~VariableNode();

	void updateLLR(std::vector<CheckNode> *checkNodeVector);

	void setAllLLR(std::vector<double> *llrVector);

	void setLLRat(double llr, int branch_index);

	void setChannelLLR(double llr);

	void setCoordinates(int a, int b);

	double getLLRat(int branch_index) const;

	/**
	 * Get the complete vector of LLRs
	 */
	const std::vector<double>* getAllLLR() const;

	std::pair<int, int> getCoordinates() const;

	std::vector<line> getCheckNodes() const;

private:
	std::vector<double> *m_llrVector;
	double m_channelLLR;
	int m_j;
	int m_a;
	int m_b;
	std::vector<line> m_checkNodes;

	void updateLLRat(int branch_index, std::vector<CheckNode> *checkNodeVector);
	
};

#endif
// VARIABLE_NODE