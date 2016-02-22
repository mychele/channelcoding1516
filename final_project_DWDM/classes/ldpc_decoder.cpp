#include "ldpc_decoder.h"
#include <vector>
#include <algorithm>
#include <limits>
#include <cstring>

LdpcDecoder::LdpcDecoder() {
	// initialize the m_variableNodeVector
	m_variableNodeVector = new std::vector<VariableNode>(ALL_COLUMNS*ALL_ROWS);
	// set coordinates for each node
	for(int node_index = 0; node_index < m_variableNodeVector->size(); ++node_index) {
															// row 						// column
		m_variableNodeVector->at(node_index).setCoordinates(node_index/(int)ALL_COLUMNS, node_index%(int)ALL_COLUMNS);
	}
	// set the LLR to +inf for nodes [120, 292] and 293*105-1, ... 293*110-1. This will never be updated!
	for(int node_index = ALL_COLUMNS - INIT_ZERO_BIT; node_index < ALL_COLUMNS; ++node_index) {
		m_variableNodeVector->at(node_index).setLLR(std::numeric_limits<double>::infinity());
	}					//105						//111
	for(int row_index = INFO_ROWS; row_index < INFO_ROWS + PC_ROWS - 1; ++row_index) {
		m_variableNodeVector->at(row_index*ALL_COLUMNS + ALL_COLUMNS - 1).setLLR(std::numeric_limits<double>::infinity());
	}

	//initialize the m_checkNodeVector
	m_checkNodeVector = new std::vector<CheckNode>(ALL_EQ);
	// set lines in each valid node
	for(int slope_index = 0; slope_index < (int)PC_ROWS - 1; ++slope_index) { // fill the first six blocks
		for(int c = 0; c < (int)ALL_COLUMNS - 1; ++c) { // with 292 valid nodes each
			m_checkNodeVector->at(slope_index*(int)ALL_COLUMNS + c).setLine(line(slope_index, c));
		}
	}
	// fill the last block
	for(int c = 0; c < (int)ALL_COLUMNS; ++c) { // with 293 valid nodes 
		m_checkNodeVector->at(((int)PC_ROWS - 1)*(int)ALL_COLUMNS + c).setLine(line((int)PC_ROWS - 1, c));
	}

	m_decisionVector = new std::vector<bool>(ALL_BIT, 0);
}

LdpcDecoder::~LdpcDecoder() {
	// clean up
	delete m_variableNodeVector;
	delete m_checkNodeVector;
	//delete m_receivedLLR;
	delete m_decisionVector;
}

// received data must be ALL_BIT bits
std::vector<bool>*
LdpcDecoder::decode(std::vector<double> *receivedData, double sigma_w2) {
	// assign m_sigmaw2
	m_sigmaw2 = sigma_w2;
	m_alpha = -2/m_sigmaw2;
	// initialize variable nodes LLR with the channel LLR
	m_receivedLLR = new std::vector<double>();
	m_receivedLLR->insert(m_receivedLLR->begin(), receivedData->begin(), receivedData->end());
	std::transform(m_receivedLLR->begin(), m_receivedLLR->end(), m_receivedLLR->begin(),
               std::bind1st(std::multiplies<double>(),m_alpha));
	initializeVariableNodes();
	// until the maximum number of attempt is reached
	int attempt_index = 0;
	do {
		updateCheckNodes();
		updateVariableNodes();
		// int zero_llr = 0;
		// for (int i = 0; i < m_variableNodeVector->size(); i++) {
		// 	if(m_variableNodeVector->at(i).getLLR()==0) {zero_llr++;}
		// }
		// std::cout << zero_llr << "\n";
	} while(attempt_index++ < MAX_ATTEMPTS && !isCodewordFound());
	delete m_receivedLLR;
	return m_decisionVector;
}

// received data must be ALL_BIT bits
std::vector<bool>*
LdpcDecoder::decodeMinSum(std::vector<double> *receivedData, double sigma_w2) {
	// assign m_sigmaw2
	m_sigmaw2 = sigma_w2;
	m_alpha = -2/m_sigmaw2;
	// initialize variable nodes LLR with the channel LLR
	m_receivedLLR = new std::vector<double>();
	m_receivedLLR->insert(m_receivedLLR->begin(), receivedData->begin(), receivedData->end());
	std::transform(m_receivedLLR->begin(), m_receivedLLR->end(), m_receivedLLR->begin(),
               std::bind1st(std::multiplies<double>(),m_alpha));
	initializeVariableNodes();
	// until the maximum number of attempt is reached
	int attempt_index = 0;
	do {
		updateCheckNodesMinSum();
		updateVariableNodes();
		// int zero_llr = 0;
		// for (int i = 0; i < m_variableNodeVector->size(); i++) {
		// 	if(m_variableNodeVector->at(i).getLLR()==0) {zero_llr++;}
		// }
		// std::cout << zero_llr << "\n";
	} while(attempt_index++ < MAX_ATTEMPTS && !isCodewordFound());
	delete m_receivedLLR;
	return m_decisionVector;
}

const std::vector<VariableNode>* 
LdpcDecoder::getVariableNodeVector() const {
	return m_variableNodeVector;
}
const std::vector<CheckNode>* 
LdpcDecoder::getCheckNodeVector() const {
	return m_checkNodeVector;
}

void 
LdpcDecoder::initializeVariableNodes() {
	// hypothesis of equally likely input symbols -> ln(p(u=0)/p(u=1)) = 0
	// LLRg_>c = -2r_l/sigma_w^2
	
	// LLR for [0, 120-1]
	int rx_index = 0;
	for(; rx_index < ALL_COLUMNS - INIT_ZERO_BIT; ++rx_index) {
		m_variableNodeVector->at(rx_index).setLLR(m_receivedLLR->at(rx_index));
	}

	rx_index = ALL_COLUMNS;
	for(; rx_index < ALL_INFO_BIT; ++rx_index) {
		m_variableNodeVector->at(rx_index).setLLR(m_receivedLLR->at(rx_index));
	} // up to 30765 - 1

	for(int pc_row_index = 1; pc_row_index < PC_ROWS; ++pc_row_index) {
		for(; rx_index < ALL_INFO_BIT + pc_row_index*ALL_COLUMNS - 1; ++rx_index) { 
		// skip bit 31057, 31350, 31643, 31936, 32229, 32522, since they must be left at +inf
			m_variableNodeVector->at(rx_index).setLLR(m_receivedLLR->at(rx_index));
		}
		rx_index++; 
	}
	// copy the last 293 bit
	for(; rx_index < ALL_INFO_BIT + ALL_EQ; rx_index++) { // up to the last bit 
		m_variableNodeVector->at(rx_index).setLLR(m_receivedLLR->at(rx_index));
	}
}

void 
LdpcDecoder::updateVariableNodes() {
	// update the first 120 nodes
	int rx_index = 0;
	for(; rx_index < ALL_COLUMNS - INIT_ZERO_BIT; ++rx_index) {
		m_variableNodeVector->at(rx_index).updateLLR(m_checkNodeVector);
	}
	rx_index = ALL_COLUMNS;
	for(; rx_index < ALL_INFO_BIT; ++rx_index) {
		m_variableNodeVector->at(rx_index).updateLLR(m_checkNodeVector);
	} // up to 30765 - 1

	for(int pc_row_index = 1; pc_row_index < PC_ROWS; ++pc_row_index) {
		for(; rx_index < ALL_INFO_BIT + pc_row_index*ALL_COLUMNS - 1; ++rx_index) { 
		// skip bit 31057, 31350, 31643, 31936, 32229, 32522, since they must be left at +inf
			m_variableNodeVector->at(rx_index).updateLLR(m_checkNodeVector);
		}
		rx_index++; 
	}
	// update the last 293 bit
	for(; rx_index < ALL_INFO_BIT + ALL_EQ; rx_index++) { // up to the last bit 
		m_variableNodeVector->at(rx_index).updateLLR(m_checkNodeVector);
	}
}

void 
LdpcDecoder::updateCheckNodes() {
	// update 292 nodes for block [0, 292], [293, 585] ... and 293 nodes for the last block. Redudant check nodes are skipped
	int block_index = 0;
	int index = 0;
	int six_rows = (PC_ROWS - 1)*ALL_COLUMNS;
	for(; block_index < six_rows; block_index+=ALL_COLUMNS) {
		for(int c = 0; c < ALL_COLUMNS - 1; ++c) {
			m_checkNodeVector->at(block_index + c).updateLLR(m_variableNodeVector);
		}
	}
	for(int c = 0; c < ALL_COLUMNS; ++c) {
		m_checkNodeVector->at(block_index + c).updateLLR(m_variableNodeVector);
	}
}

void 
LdpcDecoder::updateCheckNodesMinSum() {
	// update 292 nodes for block [0, 292], [293, 585] ... and 293 nodes for the last block. Redudant check nodes are skipped
	int block_index = 0;
	int index = 0;
	int six_rows = (PC_ROWS - 1)*ALL_COLUMNS;
	for(; block_index < six_rows; block_index+=ALL_COLUMNS) {
		for(int c = 0; c < ALL_COLUMNS - 1; ++c) {
			m_checkNodeVector->at(block_index + c).updateLLRminSum(m_variableNodeVector);
		}
	}
	for(int c = 0; c < ALL_COLUMNS; ++c) {
		m_checkNodeVector->at(block_index + c).updateLLRminSum(m_variableNodeVector);
	}
}

bool 
LdpcDecoder::isCodewordFound() {
	// marginalize, then cycle on the 2045 checkNodes and check if each has an even number of 1
	int rx_index = 0;
	for(; rx_index < ALL_COLUMNS - INIT_ZERO_BIT; ++rx_index) {
		m_decisionVector->at(rx_index) = (m_variableNodeVector->at(rx_index).getLLR() + m_receivedLLR->at(rx_index) >= 0) ? 0 : 1;
	}
	rx_index = ALL_COLUMNS;
	for(; rx_index < ALL_INFO_BIT; ++rx_index) {
		m_decisionVector->at(rx_index) = (m_variableNodeVector->at(rx_index).getLLR() + m_receivedLLR->at(rx_index) >= 0) ? 0 : 1;
	} // up to 30765 - 1

	for(int pc_row_index = 1; pc_row_index < PC_ROWS; ++pc_row_index) {
		for(; rx_index < ALL_INFO_BIT + pc_row_index*ALL_COLUMNS - 1; ++rx_index) { 
		// skip bit 31057, 31350, 31643, 31936, 32229, 32522, since they must be left at +inf
			m_decisionVector->at(rx_index) = (m_variableNodeVector->at(rx_index).getLLR() + m_receivedLLR->at(rx_index) >= 0) ? 0 : 1;
		}
		rx_index++; 
	}
	// update the last 293 bit
	for(; rx_index < ALL_INFO_BIT + ALL_EQ; rx_index++) { // up to the last bit 
		m_decisionVector->at(rx_index) = (m_variableNodeVector->at(rx_index).getLLR() + m_receivedLLR->at(rx_index) >= 0) ? 0 : 1;
	}

	// cycle on the 2045 independet check nodes
	int block_index = 0;
	int six_rows = (PC_ROWS - 1)*ALL_COLUMNS;
	int numNonEvenCheckNodes = 0;
	for(; block_index < six_rows; block_index+=ALL_COLUMNS) {
		for(int c = 0; c < ALL_COLUMNS - 1; ++c) {
			if(m_checkNodeVector->at(block_index+c).areOnesOdd(m_decisionVector)) {
				++numNonEvenCheckNodes;
			}
		}
	}
	for(int c = 0; c < ALL_COLUMNS; ++c) {
		if(m_checkNodeVector->at(block_index+c).areOnesOdd(m_decisionVector)) {
				++numNonEvenCheckNodes;
		}
	}
	return ((numNonEvenCheckNodes > 0) ? 0 : 1);
}

