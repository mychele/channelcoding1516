#include "ldpc_decoder.h"
#include <vector>

LdpcDecoder::LdpcDecoder() {
	// initialize the m_variableNodeVector
	m_variableNodeVector = new std::vector<VariableNode>;
	m_variableNodeVector->resize(ALL_COLUMNS*ALL_ROWS);
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
	m_checkNodeVector = new std::vector<CheckNode>;
	m_checkNodeVector->resize(PC_ROWS*ALL_COLUMNS);
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
}

LdpcDecoder::~LdpcDecoder() {
	// clean up
	delete m_variableNodeVector;
	delete m_checkNodeVector;
}

InfoWord
LdpcDecoder::decode(std::vector<double> *receivedData, double sigma_w) {
	// initialize variable nodes LLR with the channel LLR
	int attempt_index = 0;
	do {
	initializeVariableNodes(receivedData, sigma_w);
	updateCheckNodes();
	updateVariableNodes();
	} while(!marginalizeCheckNodes() && attempt_index++ < MAX_ATTEMPTS);

	std::cout << "Decoding completed" << "\n";

	return 0;
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
LdpcDecoder::initializeVariableNodes(std::vector<double> *receivedData, double sigma_w) {
	// hypothesis of equally likely input symbols -> ln(p(u=0)/p(u=1)) = 0
	// LLRg_>c = -2r_l/sigma_w^2
	double sigma_w2 = std::pow(sigma_w, 2);
	double alpha = -2/sigma_w2;
	// LLR for [0, 120-1]
	int rx_index = 0;
	for(; rx_index < ALL_COLUMNS - INIT_ZERO_BIT; ++rx_index) {
		m_variableNodeVector->at(rx_index).setLLR((receivedData->at(rx_index))*alpha);
	}
	rx_index = ALL_COLUMNS;
	for(; rx_index < ALL_INFO_BIT; ++rx_index) {
		m_variableNodeVector->at(rx_index).setLLR((receivedData->at(rx_index))*alpha);
	} // up to 30765 - 1

	for(int pc_row_index = 1; pc_row_index < PC_ROWS; ++pc_row_index) {
		for(; rx_index < ALL_INFO_BIT + pc_row_index*ALL_COLUMNS - 1; ++rx_index) { 
		// skip bit 31057, 31350, 31643, 31936, 32229, 32522, since they must be left at +inf
			m_variableNodeVector->at(rx_index).setLLR((receivedData->at(rx_index))*alpha);
		}
		rx_index++; 
	}
	// copy the last 293 bit
	for(; rx_index < ALL_INFO_BIT + PC_ROWS*ALL_COLUMNS; rx_index++) { // up to the last bit 
		m_variableNodeVector->at(rx_index).setLLR((receivedData->at(rx_index))*alpha);
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
	// copy the last 293 bit
	for(; rx_index < ALL_INFO_BIT + PC_ROWS*ALL_COLUMNS; rx_index++) { // up to the last bit 
		m_variableNodeVector->at(rx_index).updateLLR(m_checkNodeVector);
	}
}

void 
LdpcDecoder::updateCheckNodes() {
	// update 292 nodes for block [0, 292], [293, 585] ... and 293 nodes for the last block. Redudant check nodes are skipped
	int block_index = 0;
	for(int slope_index = 0; slope_index < PC_ROWS - 1; ++slope_index) {
		block_index = slope_index*ALL_COLUMNS;
		for(int c = 0; c < ALL_COLUMNS - 1; ++c) {
			m_checkNodeVector->at(block_index + c).updateLLR(m_variableNodeVector);
		}
	}
	block_index = (PC_ROWS - 1)*ALL_COLUMNS;
	for(int c = 0; c < ALL_COLUMNS; ++c) {
		m_checkNodeVector->at(block_index + c).updateLLR(m_variableNodeVector);
	}
}

bool 
LdpcDecoder::marginalizeCheckNodes() {
	int block_index = 0;
	bool all_zero = 1;
	for(int slope_index = 0; slope_index < PC_ROWS - 1; ++slope_index) {
		block_index = slope_index*ALL_COLUMNS;
		for(int c = 0; c < ALL_COLUMNS - 1; ++c) {
			if(m_checkNodeVector->at(block_index + c).getLLR() < 0) {all_zero = 0; return 0;}
		}
	}
	block_index = (PC_ROWS - 1)*ALL_COLUMNS;
	for(int c = 0; c < ALL_COLUMNS; ++c) {
		if(m_checkNodeVector->at(block_index + c).getLLR() < 0) {all_zero = 0; return 0;}
	}
	return all_zero;
}

std::vector<bool> 
LdpcDecoder::marginalizeVariableNodes(std::vector<double> *receivedData, double sigma_w) {
	std::vector<bool> decisionVector;
	decisionVector.resize(INFO_BIT);
	double alpha = -2/std::pow(sigma_w, 2);
	int rx_index = 0;
	bool u = 0;
	for(; rx_index < ALL_COLUMNS - INIT_ZERO_BIT; ++rx_index) {
		u = (m_variableNodeVector->at(rx_index).getLLR() + receivedData->at(rx_index)*alpha > 0) ? 0 : 1;
		decisionVector.push_back(u);
	}
	rx_index = ALL_COLUMNS;
	for(; rx_index < ALL_INFO_BIT; ++rx_index) {
		u = (m_variableNodeVector->at(rx_index).getLLR() + receivedData->at(rx_index)*alpha > 0) ? 0 : 1;
		decisionVector.push_back(u);
	} // up to 30765 - 1
	return decisionVector;
}
