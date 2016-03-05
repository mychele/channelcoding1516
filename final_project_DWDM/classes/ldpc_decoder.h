#ifndef LDPC_DECODER
#define LDPC_DECODER

#include "variable_node.h"
#include "check_node.h"
#include "ldpc_common.h"
#include <vector>


/**
 * This class represents the decoder, either Sum Product or Min Sum. It assumes equally likely input symbols.
 */
class LdpcDecoder
{
public:
	/**
	 * Public constructor
	 */ 
	LdpcDecoder();
	/**
	 * Public destructor
	 */ 
	~LdpcDecoder();

	/**
	 * Decode receivedData using Sum Product decoder, given the estimate of sigma_w^2
	 * @param a pointer to the vector of received data, whose length is ALL_BITS
	 * @param the value of sigma_w^2
	 * @return a pointer to a vector of bool with the decisions
	 */
	std::vector<bool>* decode(std::vector<double> *receivedData, double sigma_w2);

	/**
	 * Decode receivedData using Min Sum decoder, given the estimate of sigma_w^2
	 * @param a pointer to the vector of received data, whose length is ALL_BITS
	 * @param the value of sigma_w^2
	 * @return a pointer to a vector of bool with the decisions
	 */
	std::vector<bool>* decodeMinSum(std::vector<double> *receivedData, double sigma_w2);

	/**
	 * Get a pointer to VariableNode vector
	 */
	const std::vector<VariableNode>* getVariableNodeVector() const; 
	/**
	 * Get a pointer to CheckNode vector
	 */
	const std::vector<CheckNode>* getCheckNodeVector() const; 

private:
	std::vector<VariableNode> *m_variableNodeVector;
	std::vector<CheckNode> *m_checkNodeVector;
	std::vector<double> *m_receivedLLR;
	std::vector<bool> *m_decisionVector;
	double m_sigmaw2;
	double m_alpha; //-2/sigma_w2

	inline void initializeVariableNodes();
	void updateVariableNodes();
	void updateCheckNodes();
	void updateCheckNodesMinSum();
	bool isCodewordFound();
	
};

#endif
//LDPC_DECODER