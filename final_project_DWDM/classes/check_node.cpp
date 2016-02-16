#include "check_node.h"
#include "variable_node.h"
#include <cmath>

CheckNode::CheckNode() {
	m_line = line(-1, -1);
}

CheckNode::CheckNode(line lineEq) {
	m_line = lineEq;
}

void 
CheckNode::updateLLR(std::vector<VariableNode> *variableNodeVector) {
	// cycle on the variableNodes of the line, compute sum(phiTilde(|x|)) and prod(sgn(x))
	double sumPhiTilde = 0;
	int prodSgn = 1;
	int slope = slopes[m_line.first];
	int b;
	int variable_node_index;
	double llr_var;
	for (int a = 0; a < (int)ALL_ROWS; ++a) {
					//(a*s + c)%293
		b = (a*slope + m_line.second)%ALL_COLUMNS;
		variable_node_index = ALL_COLUMNS*a + b;
		llr_var = variableNodeVector->at(variable_node_index).getLLR();
		if(llr_var > 0) {
			sumPhiTilde += phiTilde(llr_var);
		} else if(llr_var < 0) {
			sumPhiTilde += phiTilde(-llr_var);
			prodSgn = -prodSgn;
		} else { // llr_var is equal to 0. Then phiTilde(llr_var) -> inf, and also the sum is inf. Then phiTilde(inf) = 0
			sumPhiTilde = std::numeric_limits<double>::infinity(); //just a flag now
			break; // no need to compute other variable nodes
		}
	}
	if(isinf(sumPhiTilde)) {
		m_llr = 0;
	} else {
		if(prodSgn > 0) {
			m_llr = phiTilde(sumPhiTilde);
		} else {
			m_llr = - phiTilde(sumPhiTilde);
		}
	}
}

void 
CheckNode::setLine(line lineEq) {
	m_line = lineEq;
}

void 
CheckNode::setLLR(double llr) {
	m_llr = llr;
}

double 
CheckNode::getLLR() const {
	return m_llr;
}

line 
CheckNode::getLine() const {
	return m_line;
}
