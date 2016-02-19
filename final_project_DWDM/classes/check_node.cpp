#include "check_node.h"
#include "variable_node.h"
#include <cmath>

CheckNode::CheckNode() {
	m_line = line(-1, -1);
	m_llrVector = new std::vector<double>(ALL_ROWS, 0);
	m_variableNodeColumnIndex = std::vector<int>(ALL_ROWS, 0);

}

CheckNode::CheckNode(line lineEq) {
	m_line = lineEq;
	m_llrVector = new std::vector<double>(ALL_ROWS, 0);
	m_variableNodeColumnIndex = std::vector<int>(ALL_ROWS, 0);
	int slope = slopes[m_line.first];
	for(int a = 0; a < ALL_ROWS; ++a) {
		m_variableNodeColumnIndex[a] = (a*slope + m_line.second)%ALL_COLUMNS;
	}
}

CheckNode::~CheckNode() {
	delete m_llrVector;
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

	// variables
	double sumPhiTilde = 0;
	int prodSgn = 1;
	int slope = slopes[m_line.first];
	double llr_var;

	// cycle on all the rows, with the exception of row_index
	int block_index = 0;
	int a = 0;
	int row_block_index = row_index*ALL_COLUMNS;
	for (; block_index < row_block_index; block_index += ALL_COLUMNS) {
															//(a*s + c)%293
		llr_var = variableNodeVector->at(block_index + m_variableNodeColumnIndex[a++]).getLLR();
				
		if(llr_var > 0) {
			sumPhiTilde += phiTilde(llr_var);
		} else if(llr_var < 0) {
			sumPhiTilde += phiTilde(-llr_var);
			prodSgn = - prodSgn;
		} else { // llr_var is equal to 0. Then phiTilde(llr_var) -> inf, and also the sum is inf. Then phiTilde(inf) = 0
			sumPhiTilde = L_INFINITY; //just a flag now
			break; // no need to compute other variable nodes
		}
	}
	block_index += ALL_COLUMNS; // skip row_index row
	a++;
	if(sumPhiTilde == L_INFINITY) { // don't need to cycle also on the other rows, the sum will be surely +inf, and phiTilde(+inf)=0
		// save the outgoing llr
		m_llrVector->at(row_index) = 0;
	} else { // continue to cycle
		for (; block_index < ALL_BIT; block_index += ALL_COLUMNS) {
															//(a*s + c)%293
		llr_var = variableNodeVector->at(block_index + m_variableNodeColumnIndex[a++]).getLLR();
					
			if(llr_var > 0) {
				sumPhiTilde += phiTilde(llr_var);
			} else if(llr_var < 0) {
				sumPhiTilde += phiTilde(-llr_var);
				prodSgn = - prodSgn;
			} else { // llr_var is equal to 0. Then phiTilde(llr_var) -> inf, and also the sum is inf. Then phiTilde(inf) = 0
				sumPhiTilde = L_INFINITY; //just a flag now
				break; // no need to compute other variable nodes
			}
		}

		// compute the outgoing llr
		if(sumPhiTilde == L_INFINITY) {
			m_llrVector->at(row_index) = 0;
		} else {
			if(prodSgn > 0) {
				m_llrVector->at(row_index) = phiTilde(sumPhiTilde);
			} else {
				m_llrVector->at(row_index) = - phiTilde(sumPhiTilde);
			}
		}
	}
}

void 
CheckNode::setLine(line lineEq) {
	m_line = lineEq;
	int slope = slopes[m_line.first];
	for(int a = 0; a < ALL_ROWS; ++a) {
		m_variableNodeColumnIndex[a] = (a*slope + m_line.second)%ALL_COLUMNS;
	}
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
