#include "check_node.h"

CheckNode::CheckNode(line lineEq) {
	m_line = lineEq;
}

void 
CheckNode::updateLLR() {

}

void 
CheckNode::setLine(line lineEq) {
	m_line = lineEq;
}

double 
CheckNode::getLLR() const {
	return m_llr;
}

line 
CheckNode::getLine() const {
	return m_line;
}
