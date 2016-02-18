#include "check_node.h"
#include "variable_node.h"
#include <cmath>

CheckNode::CheckNode() {
	m_line = line(-1, -1);
	m_llrVector = new std::vector<double>();
	m_llrVector->resize(ALL_ROWS);
}

CheckNode::CheckNode(line lineEq) {
	m_line = lineEq;
	m_llrVector = new std::vector<double>();
	m_llrVector->resize(ALL_ROWS);
}

void 
CheckNode::updateLLR(std::vector<VariableNode> *variableNodeVector) {
	// cycle on all the outgoing branch, one per line, and update the LLR of each one
	for(int row_index = 0; row_index < ALL_ROWS; ++row_index) {
		//std::cout << "Updating LLR for outgoing branch " << row_index << "\n";
		updateLLRat(row_index, variableNodeVector);
	}
}

void
CheckNode::updateLLRat(int row_index, std::vector<VariableNode> *variableNodeVector) {
	// cycle on the other incoming branches, compute phiTilde(sum(phiTilde(|x|)) and prod(sgn(x)))

	// verbose mode for debugging
	bool verb = 0; //(m_line.first==0 && m_line.second==29);

	// variables
	double sumPhiTilde = 0;
	int prodSgn = 1;
	int slope = slopes[m_line.first];
	int b;
	int variable_node_index;
	double llr_var;

	// cycle on all the rows, with the exception of row_index
	for (int a = 0; a < row_index; ++a) {
					//(a*s + c)%293
		b = (a*slope + m_line.second)%ALL_COLUMNS;
		variable_node_index = ALL_COLUMNS*a + b;
		llr_var = variableNodeVector->at(variable_node_index).getLLR();
		
		if(verb) {std::cout << "LLR of node " << a << "," << b << " = " << llr_var << " for check_node " << m_line.first << " " << m_line.second <<"\n";}
		
		if(llr_var > 0) {
			sumPhiTilde += phiTilde(llr_var);
		} else if(llr_var < 0) {
			sumPhiTilde += phiTilde(-llr_var);
			prodSgn = - prodSgn;
		} else { // llr_var is equal to 0. Then phiTilde(llr_var) -> inf, and also the sum is inf. Then phiTilde(inf) = 0
			sumPhiTilde = std::numeric_limits<double>::infinity(); //just a flag now
			break; // no need to compute other variable nodes
		}
		if(verb) {std::cout << "sumPhiTilde = " << sumPhiTilde <<"\n";}
		if(verb) {std::cout << "prodSgn = " << prodSgn <<"\n";}
	}

	if(isinf(sumPhiTilde)) { // don't need to cycle also on the other rows, the sum will be surely +inf, and phiTilde(+inf)=0
		// save the outgoing llr
		m_llrVector->at(row_index) = 0;
	} else { // continue to cycle
		for (int a = row_index + 1; a < ALL_ROWS; ++a) {

			b = (a*slope + m_line.second)%ALL_COLUMNS;
			variable_node_index = ALL_COLUMNS*a + b;
			llr_var = variableNodeVector->at(variable_node_index).getLLR();
			
			if(verb) {std::cout << "LLR of node " << a << "," << b << " = " << llr_var << " for check_node " << m_line.first << " " << m_line.second <<"\n";}
			
			if(llr_var > 0) {
				sumPhiTilde += phiTilde(llr_var);
			} else if(llr_var < 0) {
				sumPhiTilde += phiTilde(-llr_var);
				prodSgn = - prodSgn;
			} else { // llr_var is equal to 0. Then phiTilde(llr_var) -> inf, and also the sum is inf. Then phiTilde(inf) = 0
				sumPhiTilde = std::numeric_limits<double>::infinity(); //just a flag now
				break; // no need to compute other variable nodes
			}
			if(verb) {std::cout << "sumPhiTilde = " << sumPhiTilde <<"\n";}
			if(verb) {std::cout << "prodSgn = " << prodSgn <<"\n";}
		}

		// compute the outgoing llr
		if(isinf(sumPhiTilde)) {
			m_llrVector->at(row_index) = 0;
		} else {
			if(prodSgn > 0) {
				m_llrVector->at(row_index) = phiTilde(sumPhiTilde);
			} else {
				m_llrVector->at(row_index) = - phiTilde(sumPhiTilde);
			}
		}
	}
	if(verb) {std::cout << "for outgoing branch " << row_index << " of check_node " << m_line.first << " " << m_line.second << " llr " << m_llrVector->at(row_index) <<"\n";}
}

void 
CheckNode::setLine(line lineEq) {
	m_line = lineEq;
}

void 
CheckNode::setAllLLR(std::vector<double> *llrVector) {
	m_llrVector = llrVector;
}

void
CheckNode::setLLRat(double llr, int row_index) {
	m_llrVector->at(row_index) = llr;
}

double 
CheckNode::getLLRat(int row_index) const {
	return m_llrVector->at(row_index);
}

const std::vector<double>*
CheckNode::getAllLLR() const {
	return m_llrVector;
}

line 
CheckNode::getLine() const {
	return m_line;
}
