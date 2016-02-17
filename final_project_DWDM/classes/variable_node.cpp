#include "variable_node.h"
#include "check_node.h"

VariableNode::VariableNode() {

}

VariableNode::VariableNode(int a, int b) {
	//TODO check intervals!
	m_a = a;
	m_b = b;
	m_j = a*(int)ALL_COLUMNS + b;	
	m_checkNodes = std::vector<line>((int)PC_ROWS);	// 7
	for(int slope_index = 0; slope_index < (int)PC_ROWS - 1; ++slope_index) {
		int c = reverseModulo((int)ALL_COLUMNS, m_a*slopes[slope_index], m_b);
		m_checkNodes[slope_index] = (c != 292) ? line(slope_index, c) : line(slope_index, -1); // for slope_index from 0 to 5 remove the eq with c=292
	}
	m_checkNodes[(int)PC_ROWS - 1] = line((int)PC_ROWS - 1, reverseModulo((int)ALL_COLUMNS, m_a*slopes[(int)PC_ROWS - 1], m_b));
}

void
VariableNode::updateLLR(std::vector<CheckNode> *checkNodeVector) {
	// sum the LLR of the corresponding check nodes
	double llr = 0;
	for(int slope_index = 0; slope_index < (int)PC_ROWS; ++slope_index) {
		if(m_checkNodes[slope_index].second > -1) { // check if it is a valid check node
			llr += checkNodeVector->at(slope_index*(int)ALL_COLUMNS + m_checkNodes[slope_index].second).getLLR(); // since all check nodes (even redundant) are in checkNodeVector,
																								// then the indices will be slope_index*293 + c
		}
	}
	m_llr = llr;
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
	for(int slope_index = 0; slope_index < (int)PC_ROWS - 1; ++slope_index) {
		int c = reverseModulo((int)ALL_COLUMNS, m_a*slopes[slope_index], m_b);
		m_checkNodes[slope_index] = (c < 292) ? line(slope_index, c) : line(slope_index, -1); // for slope_index from 0 to 5 remove the eq with c=292
	}
	m_checkNodes[(int)PC_ROWS - 1] = line((int)PC_ROWS - 1, reverseModulo((int)ALL_COLUMNS, m_a*slopes[(int)PC_ROWS - 1], m_b));
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