#ifndef CHECK_NODE
#define CHECK_NODE

#include <vector>
#include <array>
#include "ldpc_common.h"
#include "variable_node.h"

class VariableNode;

/**
 * CheckNode class
 */
class CheckNode
{
public:
	/**
	 * Default constructor, it creates a check node with (-1, -1) line
	 */
	CheckNode();

	/**
	 * Public constructor, it creates a check node with line lineEq
	 * @param a line
	 */
	CheckNode(line lineEq);

	/**
	 * Public destructor
	 */
	~CheckNode();

	/**
	 * Update each branch's LLR with Sum Product
	 * @param a pointer to a vector of ALL_BITS variable nodes
	 */
	void updateLLR(std::vector<VariableNode> *variableNodeVector);
	/**
	 * Update each branch's LLR with Min Sum
	 * @param a pointer to a vector of ALL_BITS variable nodes
	 */
	void updateLLRminSum(std::vector<VariableNode> *variableNodeVector);

	/**
	 * Set the line of the node
	 * @param a line
	 */
	void setLine(line lineEq);

	/**
	 * Set the LLR of the row_index outgoing branch
	 * @param an LLR
	 * @param the index of the branch to update
	 */
	void setLLRat(double llr, int row_index);

	/**
	 * Set the LLR of each outgoing branch
	 * @param a vector of 112 LLR
	 */
	void setAllLLR(std::vector<double> *llrVector);

	/**
	 * Get the complete vector of LLRs
	 * @return a vector of 112 LLR
	 */
	const std::vector<double>* getAllLLR() const;

	/**
	 * Get the LLR of row_index outgoing branch
	 * @param the index of the branch
	 * @return the LLR of the given branch
	 */
	double getLLRat(int row_index) const;

	/**
	 * Returns the line of this check node
	 * @return a line
	 */
	line getLine() const;

	/**
	 * Return 1 if there is an odd number of 1 in the incoming branches
	 * @param a vector of bool (marginalized values)
	 * @return 1 if the number of 1 is odd
	 */
	bool areOnesOdd(std::vector<bool> *decisionVector);
	
private:
	line m_line; // if m_line has both entries negative, the check node is one of the invalid ones
	std::vector<int> m_variableNodeColumnIndex;
	// store the incoming LLR
	std::vector<double> *m_incomingLlrArray;
	// store phiTilde of the abs of incoming LLR, or abs of incoming LLR if minsum is used
	std::vector<double> *m_phiArray;
	// store the sign of the incoming LLR
	std::vector<int> *m_signArray;
	// sum of the elements in m_phiArray
	double m_phiSum;
	// product of the sign in m_signArray
	int m_signProd;
	// minimum LLR, used with min_sum
	double m_minLLR;
	// second minimum LLR, used with min_sum
	double m_secondMinLLR;

	/**
	 * Vector that stores all the LLR, one entry per outgoing branch
	 */
	std::vector<double> *m_llrVector;

	inline void updateLLRat(int row_index);
	void updateLLRatMinSum(int row_index);

};

#endif
//CHECK_NODE