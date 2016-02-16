#include "variable_node.h"
#include "check_node.h"

VariableNode::VariableNode(int a, int b) {
	//TODO check intervals!
	m_a = a;
	m_b = b;
	m_j = a*(int)ALL_COLUMNS + b;	
	m_checkNodes = std::vector<line>((int)PC_ROWS);	// 7
	for(int slope_index = 0; slope_index < (int)PC_ROWS; ++slope_index) {
		int c = reverseModulo((int)ALL_COLUMNS, m_a*slopes[slope_index], m_b);
		m_checkNodes[slope_index] = line(slope_index, c);
	}
}

void
VariableNode::updateLLR(std::vector<CheckNode> *checkNodeVector) {
	// sum the LLR of the corresponding check nodes

}

void 
VariableNode::setLLR(double llr) {
	m_llr = llr;
}

void 
VariableNode::setCoordinates(int a, int b) {
	m_a = a;
	m_b = b;
	// when a and b are changed, m_checkNodes needs to change too
	m_checkNodes = std::vector<line>((int)PC_ROWS);	// 7
	for(int slope_index = 0; slope_index < (int)PC_ROWS; ++slope_index) {
		int c = reverseModulo((int)ALL_COLUMNS, m_a*slopes[slope_index], m_b);
		m_checkNodes[slope_index] = line(slope_index, c);
	}
}

double 
VariableNode::getLLR() const {
	return m_llr;
}

std::pair<int, int>
VariableNode::getCoordinates() const {
	return coordinates(m_a, m_b);
}

std::vector<line> 
VariableNode::getCheckNodes() const {
	return m_checkNodes;
}