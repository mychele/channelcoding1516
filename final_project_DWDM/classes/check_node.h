#include "ldpc_common.h"

#ifndef CHECK_NODE
#define CHECK_NODE

class CheckNode
{
public:
	CheckNode(line lineEq);

	void updateLLR();

	void setLine(line lineEq);

	double getLLR() const;

	line getLine() const;
	
private:
	line m_line;
	// TODO consider if it makes any difference to precompute the variable nodes and store them in an array or compute them on the fly
	double m_llr;
};

#endif
//CHECK_NODE