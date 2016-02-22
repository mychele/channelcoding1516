#ifndef LDPC_DECODER
#define LDPC_DECODER

#include "variable_node.h"
#include "check_node.h"
#include "ldpc_common.h"
#include <vector>

class LdpcDecoder
{
public:
	LdpcDecoder();
	~LdpcDecoder();

	std::vector<bool>* decode(std::vector<double> *receivedData, double sigma_w);
	std::vector<bool>* decodeMinSum(std::vector<double> *receivedData, double sigma_w);

	// this is not a good practice, but it is useful for testing purposes
	// TODO remove once testing is done
	const std::vector<VariableNode>* getVariableNodeVector() const; 
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