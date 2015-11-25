/*
 * This function performs ML decoding for a 517 binary convolutional code
 * @param a pointer to the received vector
 * @param the noise std_deviation of the channel
 * @param the length of the received vector
 * @param the pointer to the array where received values will be stored
 * @param use SD or HD
 * @param the size of the window
 */

void viterbi517_windowed(double *r, double sigma_w, int n, double *u_hat, double mode, int windowSize);