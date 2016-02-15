###################################################
# Makefile for Channel Coding 1516 final project - code only
###################################################

test: invert_H encoder_test rng

invert_H: ./test/invert_H.cpp
	@clang++ -std=c++11 -c -I/Users/Michele/Dev/lib/ntl/include -I./classes $?
	@clang++ -std=c++11 -o invert_H invert_H.o -L/Users/Michele/Dev/lib/ntl/lib -lntl -lm

rng: ./test/rng_test.cpp
	@clang++ -std=c++11 -o rng_test $?

encoder: ./classes/ldpc_encoder.cpp
	@clang++ -std=c++11 -c -I. $?

encoder_test: ./test/encoder_test.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/ldpc_encoder.cpp
	@clang++ -std=c++11 -c -I./classes -g $?
	@clang++ -std=c++11 -o encoder_test encoder_test.o ldpc_encoder.o


clean:
	rm invert_H encoder_test rng_test *.o
