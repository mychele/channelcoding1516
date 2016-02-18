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

	std::vector<bool> decode(std::vector<double> *receivedData, double sigma_w);

	// this is not a good practice, but it is useful for testing purposes
	// TODO remove once testing is done
	const std::vector<VariableNode>* getVariableNodeVector() const; 
	const std::vector<CheckNode>* getCheckNodeVector() const; 

private:
	std::vector<VariableNode> *m_variableNodeVector;
	std::vector<CheckNode> *m_checkNodeVector;

	void initializeVariableNodes(std::vector<double> *receivedData, double sigma_w);
	inline void updateVariableNodes();
	inline void updateCheckNodes();
	bool marginalizeCheckNodes();
	std::vector<bool> marginalizeVariableNodes(std::vector<double> *receivedData, double sigma_w);
	
};

#endif
//LDPC_DECODER