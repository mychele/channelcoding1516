#include "check_node.h"
#include "variable_node.h"
#include <cmath>
#include <limits>
//#include <chrono>

CheckNode::CheckNode() {
	m_line = line(-1, -1);
	m_llrVector = new std::vector<double>(ALL_ROWS, 0);
	m_incomingLlrArray = new std::vector<double>(ALL_ROWS, 0);
	m_phiArray = new std::vector<double>(ALL_ROWS, 0);
	m_signArray = new std::vector<int>(ALL_ROWS, 1);
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
	m_incomingLlrArray = new std::vector<double>(ALL_ROWS, 0);
	m_phiArray = new std::vector<double>(ALL_ROWS, 0);
	m_signArray = new std::vector<int>(ALL_ROWS, 1);
}

CheckNode::~CheckNode() {
	delete m_signArray;
	delete m_phiArray;
	delete m_llrVector;
}

void 
CheckNode::updateLLR(std::vector<VariableNode> *variableNodeVector) {
	// get all the LLR that will be used in the computation in an array, so that reading is performed just once
	int block_index = 0;
	int a = 0;
	// we can compute a vector of phiTilde of |llr| and a vector of sign
	// then we can sum the vector of phiTilde(|llr|) once, and subtract the corresponding llr when computing each branch
	// we can compute the product of the sign. Then if the outgoing LLR has a positive sign the sign doesn't change, otherwise it becomes - the prev sign
	double llr_temp;
	double phi_temp;
	m_phiSum = 0;
	m_signProd = 1;
	for(; block_index < ALL_BIT; block_index += ALL_COLUMNS) {
														//(a*s + c)%293
		llr_temp = variableNodeVector->at(block_index + m_variableNodeColumnIndex[a]).getLLRat(m_line.first);
		// store the llr to handle inf cases
		m_incomingLlrArray->at(a) = llr_temp;
		if(llr_temp > 0) {
			phi_temp = phiTilde(llr_temp);
			m_phiArray->at(a) = phi_temp;
			m_signArray->at(a) = 1;
			m_phiSum += phi_temp;
			// m_signProd doesn't change
		} else if(llr_temp < 0) {
			phi_temp = phiTilde(-llr_temp);
			m_phiArray->at(a) = phi_temp;
			m_signArray->at(a) = -1;
			m_phiSum += phi_temp;
			m_signProd = -m_signProd;
		} else { // this LLR is 0 -> the other LLR will be 0, but this one maybe not. When there are LLR=0 use m_incomingLlrArray
			m_phiArray->at(a) = std::numeric_limits<double>::infinity();
			m_signArray->at(a) = 0;
			m_phiSum = std::numeric_limits<double>::infinity();
			m_signProd = 0;
		}
		a++;
	}
	// cycle on all the outgoing branch, one per line, and update the LLR of each one
	for(int row_index = 0; row_index < ALL_ROWS; ++row_index) {
		updateLLRat(row_index);
	}
}

void
CheckNode::updateLLRat(int row_index) {
	// check if m_signProd is diff from 0 (i.e. no LLR equal to 0)
	if(m_signProd != 0) {
		if((*m_signArray)[row_index] > 0) {
			// then the outgoing llr will have sign m_signProd
			m_llrVector->at(row_index) = m_signProd*phiTilde(m_phiSum - (*m_phiArray)[row_index]);
		} else { // it will be surely < 0 -> the outgoing llr will have sign -m_signProd
			m_llrVector->at(row_index) = - m_signProd*phiTilde(m_phiSum - (*m_phiArray)[row_index]);
		}
	} else { // there is one LLR that is 0
		if((*m_signArray)[row_index] != 0) { // the LLR=inf is not the one this outgoing branch -> this outgoing branch will have LLR = 0
			m_llrVector->at(row_index) = 0;
		} else { // we need to compute the LLR step by step
			// cycle on the other incoming branches, compute phiTilde(sum(phiTilde(|x|)) and prod(sgn(x)))
			// variables
			double sumPhiTilde = 0;
			int prodSgn = 1;
			double llr_var;

			// cycle on all the rows, with the exception of row_index
			int a = 0;
			//std::chrono::time_point<std::chrono::system_clock> begin = std::chrono::system_clock::now();

			for (; a < row_index; a++) {
				llr_var = (*m_incomingLlrArray)[a];
			
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
			a++;

			if(!isinf(sumPhiTilde)) { // continue to cycle
				for (; a < ALL_ROWS; ++a) {
					llr_var = (*m_incomingLlrArray)[a]; 

					if(llr_var > 0) {
						sumPhiTilde += phiTilde(llr_var);
					} else if(llr_var < 0) {
						sumPhiTilde += phiTilde(-llr_var);
						prodSgn = - prodSgn;
					} else { // llr_var is equal to 0. Then phiTilde(llr_var) -> inf, and also the sum is inf. Then phiTilde(inf) = 0
						sumPhiTilde = std::numeric_limits<double>::infinity(); //just a flag now
						break; // no need to compute other variable nodes
					}
				}
				// compute the outgoing llr
				if(sumPhiTilde < std::numeric_limits<double>::infinity()) {
					if(prodSgn > 0) {
						m_llrVector->at(row_index) = phiTilde(sumPhiTilde);
					} else {
						m_llrVector->at(row_index) = - phiTilde(sumPhiTilde);
					}
				} else {
					m_llrVector->at(row_index) = 0;
				}
			} else { // don't need to cycle also on the other rows, the sum will be surely +inf, and phiTilde(+inf)=0
				// save the outgoing llr
				m_llrVector->at(row_index) = 0;
			}

		}

	}

}

void 
CheckNode::updateLLRminSum(std::vector<VariableNode> *variableNodeVector) {
	// cycle on all the outgoing branch, one per line, and update the LLR of each one
	for(int row_index = 0; row_index < ALL_ROWS; ++row_index) {
		updateLLRatMinSum(row_index, variableNodeVector);
	}
}

void
CheckNode::updateLLRatMinSum(int row_index, std::vector<VariableNode> *variableNodeVector) {
	// cycle on the other incoming branches, compute phiTilde(sum(phiTilde(|x|)) and prod(sgn(x)))

	// variables
	double minLLR = std::numeric_limits<double>::infinity();
	int prodSgn = 1;
	int slope = slopes[m_line.first];
	double llr_var;

	// cycle on all the rows, with the exception of row_index
	int block_index = 0;
	int a = 0;
	int row_block_index = row_index*ALL_COLUMNS;
	for (; block_index < row_block_index; block_index += ALL_COLUMNS) {
															//(a*s + c)%293
		llr_var = variableNodeVector->at(block_index + m_variableNodeColumnIndex[a++]).getLLRat(m_line.first);
				
		if(llr_var > 0) {
			if(llr_var < minLLR) {
				minLLR = llr_var;
			}
		} else if(llr_var < 0) {
			if(-llr_var < minLLR) {
				minLLR = -llr_var;
			}
			prodSgn = - prodSgn;
		} else { // llr_var is equal to 0. Then min(|llr|)=0
			minLLR = 0; //just a flag now
			break; // no need to compute other variable nodes
		}
	}
	
	block_index += ALL_COLUMNS; // skip row_index row
	a++;

	if(minLLR > 0) { // continue to cycle
		for (; block_index < ALL_BIT; block_index += ALL_COLUMNS) {
														//(a*s + c)%293
			llr_var = variableNodeVector->at(block_index + m_variableNodeColumnIndex[a++]).getLLRat(m_line.first);
				
			if(llr_var > 0) {
				if(llr_var < minLLR) {
					minLLR = llr_var;
				}
			} else if(llr_var < 0) {
				if(-llr_var < minLLR) {
					minLLR = -llr_var;
				}
				prodSgn = - prodSgn;
			} else { // llr_var is equal to 0. Then min(|llr|)=0
				minLLR = 0; //just a flag now
				break; // no need to compute other variable nodes
			}
		}
		// compute the outgoing llr
		if(minLLR > 0) {
			if(prodSgn > 0) {
				m_llrVector->at(row_index) = minLLR;
			} else {
				m_llrVector->at(row_index) = - minLLR;
			}
		} else {
			m_llrVector->at(row_index) = 0;
		}
	} else { // don't need to cycle also on the other rows, the LLR will surely be 0
		// save the outgoing llr
		m_llrVector->at(row_index) = 0;
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

bool
CheckNode::areOnesOdd(std::vector<bool> *decisionVector) {
	// for rows a from 0 to 111 count the number of 1 in decisionVector
	int numOnes = 0;
	int block_index = 0;
	int a = 0;
	for (; block_index < ALL_BIT; block_index += ALL_COLUMNS) {
															//(a*s + c)%293
		numOnes += decisionVector->at(block_index + m_variableNodeColumnIndex[a++]);
	}
	
	return numOnes%2; // if numOnes is even, return 0
}
