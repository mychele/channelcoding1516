/*
 * This function performs ML decoding for a 517 binary convolutional code
 * @param a pointer to the received vector
 * @param the noise std_deviation of the channel
 * @param the length of the received vector
 * @param the pointer to the array where received values will be stored
 */

void viterbi517(double *r, double sigma_w, int n, double *u_hat);