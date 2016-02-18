#ifndef CHECK_NODE
#define CHECK_NODE

#include "ldpc_common.h"
#include "variable_node.h"
#include <vector>

class VariableNode;

class CheckNode
{
public:
	CheckNode();
	CheckNode(line lineEq);

	void updateLLR(std::vector<VariableNode> *variableNodeVector);

	void setLine(line lineEq);

	/**
	 * Set the LLR of the row_index outgoing branch
	 */
	void setLLRat(double llr, int row_index);

	void setAllLLR(std::vector<double> *llrVector);

	/**
	 * Get the complete vector of LLRs
	 */
	const std::vector<double>* getAllLLR() const;

	/**
	 * Get the LLR of row_index outgoing branch
	 */
	double getLLRat(int row_index) const;

	line getLine() const;
	
private:
	line m_line; // if m_line has both entries negative, the check node is one of the invalid ones
	// TODO consider if it makes any difference to precompute the variable nodes and store them in an array or compute them on the fly
	/**
	 * Vector that stores all the LLR, one entry per outgoing branch
	 */
	std::vector<double> *m_llrVector;

	void updateLLRat(int row_index, std::vector<VariableNode> *variableNodeVector);


};

#endif
//CHECK_NODE