/*
 * This function performs ML decoding for a 517 binary convolutional code
 * @param a pointer to the received vector
 * @param the noise std_deviation of the channel
 * @param the length of the received vector
 * @param the pointer to the array where received values will be stored
 */

void viterbi517(double *r, double sigma_w, int n, double *u_hat, double mode);

/**
 * Compute cost of the transition from neighID to stateID with symbol u
 * @param the input symbol
 * @param the prev state
 * @param the codeword associated with the transition
 * @param sigma_w
 * @param use SD (0) or HD (1)
 * @return the cost
 */
double getCost(int symbol, int neighID, double *codeword, double sigma_w, double mode);

double getSign(double value);