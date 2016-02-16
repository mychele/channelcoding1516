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

	void setLLR(double llr);

	double getLLR() const;

	line getLine() const;
	
private:
	line m_line; // if m_line has both entries negative, the check node is one of the invalid ones
	// TODO consider if it makes any difference to precompute the variable nodes and store them in an array or compute them on the fly
	double m_llr;


};

#endif
//CHECK_NODE