#include "variable_node.h"
#include "check_node.h"
#include <limits>

VariableNode::VariableNode() {
	// Create LLR vector
	m_llrVector = new std::vector<double>(PC_ROWS+1, 0); // size 8
}

VariableNode::VariableNode(int a, int b) {
	m_a = a;
	m_b = b;
	m_j = a*(int)ALL_COLUMNS + b;	
	// Fill check nodes vector
	m_checkNodes = std::vector<line>((int)PC_ROWS);	// 7
	for(int slope_index = 0; slope_index < (int)PC_ROWS - 1; ++slope_index) {
		int c = reverseModulo((int)ALL_COLUMNS, m_a*slopes[slope_index], m_b);
		m_checkNodes[slope_index] = (c != 292) ? line(slope_index, c) : line(slope_index, -1); // for slope_index from 0 to 5 remove the eq with c=292
	}
	m_checkNodes[(int)PC_ROWS - 1] = line((int)PC_ROWS - 1, reverseModulo((int)ALL_COLUMNS, m_a*slopes[(int)PC_ROWS - 1], m_b));
	// Create LLR vector
	m_llrVector = new std::vector<double>(PC_ROWS+1, 0); // size 8
}

VariableNode::~VariableNode() {
	delete m_llrVector;
}

void 
VariableNode::updateLLR(std::vector<CheckNode> *checkNodeVector) {
	// Update each LLR 
	for(int branch_index = 0; branch_index < m_llrVector->size(); ++branch_index) {
		updateLLRat(branch_index, checkNodeVector);
	}
}

void
VariableNode::updateLLRat(int branch_index, std::vector<CheckNode> *checkNodeVector) {
	
	double llr = 0;
	int slope_index = 0;
	int block_index = 0;
	int branch_block_index = branch_index*ALL_COLUMNS;
	for(; block_index < branch_block_index; block_index += ALL_COLUMNS) {
		if(m_checkNodes[slope_index].second > -1) {
			llr += checkNodeVector->at(block_index + m_checkNodes[slope_index].second).getLLRat(m_a);
		}
		slope_index++;
	}
	slope_index++;
	block_index += ALL_COLUMNS;
	for(; block_index < ALL_EQ; block_index += ALL_COLUMNS) {
		if(m_checkNodes[slope_index].second > -1) {
			llr += checkNodeVector->at(block_index + m_checkNodes[slope_index].second).getLLRat(m_a);
		}
		slope_index++;
	}

	if(branch_index < PC_ROWS) { // branch PC_ROWS is the one connected to leaf node
		llr += m_channelLLR;
	}

	m_llrVector->at(branch_index) = llr;
}

void 
VariableNode::setCoordinates(int a, int b) {
	m_a = a;
	m_b = b;
	// when a and b are changed, m_checkNodes needs to change too
	m_checkNodes = std::vector<line>((int)PC_ROWS);	// 7
	for(int slope_index = 0; slope_index < (int)PC_ROWS - 1; ++slope_index) {
		int c = reverseModulo((int)ALL_COLUMNS, m_a*slopes[slope_index], m_b);
		m_checkNodes[slope_index] = (c < 292) ? line(slope_index, c) : line(slope_index, -1); // for slope_index from 0 to 5 remove the eq with c=292
	}
	m_checkNodes[(int)PC_ROWS - 1] = line((int)PC_ROWS - 1, reverseModulo((int)ALL_COLUMNS, m_a*slopes[(int)PC_ROWS - 1], m_b));
}

void 
VariableNode::setAllLLR(std::vector<double> *llrVector) {
	m_llrVector = llrVector;
}

void
VariableNode::setLLRat(double llr, int branch_index) {
	m_llrVector->at(branch_index) = llr;
}

void 
VariableNode::setChannelLLR(double llr) {
	// set all the LLR to the channel LLR, for the first step of the iterative procedure
	for(int branch_index = 0; branch_index < m_llrVector->size(); branch_index++) {
		m_llrVector->at(branch_index) = llr;
	}
	m_channelLLR = llr; // this won't be changed
}

double 
VariableNode::getLLRat(int branch_index) const {
	return m_llrVector->at(branch_index);
}

const std::vector<double>*
VariableNode::getAllLLR() const {
	return m_llrVector;
}

std::pair<int, int>
VariableNode::getCoordinates() const {
	return coordinates(m_a, m_b);
}

std::vector<line> 
VariableNode::getCheckNodes() const {
	return m_checkNodes;
}