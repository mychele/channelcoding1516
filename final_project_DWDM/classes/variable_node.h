#ifndef VARIABLE_NODE
#define VARIABLE_NODE

#include "ldpc_common.h"
#include "check_node.h"
#include <iostream>

class CheckNode;

class VariableNode
{
public:
	/**
	 * Default public constructor
	 */
	VariableNode();

	/**
	 * Public constructor
	 * @param the row in the standard matrix
	 * @param the column in the standard matrix
	 */
	VariableNode(int a, int b);

	/**
	 * Public destructor
	 */
	~VariableNode();

	/**
	 * Update all the LLR, given the CheckNode vector with ALL_EQ entries
	 * @param a vector of CheckNode, with ALL_EQ entries
	 */
	void updateLLR(std::vector<CheckNode> *checkNodeVector);

	/**
	 * Set all the LLR with a given vector
	 * @param a vector of LLR
	 */
	void setAllLLR(std::vector<double> *llrVector);

	/**
	 * Set the LLR with of a given branch
	 * @param an LLR
	 * @param the index of the branch to set
	 */
	void setLLRat(double llr, int branch_index);

	/**
	 * Set all the LLR with the same value
	 * @param an LLR
	 */
	void setChannelLLR(double llr);

	/**
	 * Set the coordinates 
	 * @param the row in the standard matrix
	 * @param the column in the standard matrix
	 */
	void setCoordinates(int a, int b);

	/**
	 * Get the LLR with of a given branch
	 * @param the index of the branch to get
 	 * @return an LLR
	 */
	double getLLRat(int branch_index) const;

	/**
	 * Get the complete vector of LLRs
	 * @return a pointer to the LLR vector
	 */
	const std::vector<double>* getAllLLR() const;

	/**
	 * Get the coordinates
	 * @return a pair (a,b)
	 */
	std::pair<int, int> getCoordinates() const;

	/**
	 * Get the lines of connected check nodes
	 * @return a vector of lines
	 */
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