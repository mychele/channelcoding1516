#ifndef CHECK_NODE
#define CHECK_NODE

#include "ldpc_common.h"
#include "variable_node.h"
#include <vector>
#include <array>


class VariableNode;

class CheckNode
{
public:
	CheckNode();
	CheckNode(line lineEq);
	~CheckNode();

	void updateLLR(std::vector<VariableNode> *variableNodeVector);
	void updateLLRminSum(std::vector<VariableNode> *variableNodeVector);

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

	/**
	 * Return 1 if there is an odd number of 1 in the incoming branches
	 */
	bool areOnesOdd(std::vector<bool> *decisionVector);
	
private:
	line m_line; // if m_line has both entries negative, the check node is one of the invalid ones
	std::vector<int> m_variableNodeColumnIndex;
	// store the incoming LLR
	std::vector<double> *m_incomingLlrArray;
	// store phiTilde of the abs of incoming LLR
	std::vector<double> *m_phiArray;
	// store the sign of the incoming LLR
	std::vector<int> *m_signArray;
	// sum of the elements in m_phiArray
	double m_phiSum;
	// product of the sign in m_signArray
	int m_signProd;

	/**
	 * Vector that stores all the LLR, one entry per outgoing branch
	 */
	std::vector<double> *m_llrVector;

	inline void updateLLRat(int row_index);
	void updateLLRatMinSum(int row_index, std::vector<VariableNode> *variableNodeVector);

};

#endif
//CHECK_NODE