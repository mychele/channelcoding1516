###################################################
# Makefile for Channel Coding 1516 final project - code only
###################################################

test: invert_H encoder_test rng node_test

invert_H: ./test/invert_H.cpp
	@clang++ -std=c++11 -c -I/Users/Michele/Dev/lib/ntl/include -I./classes $?
	@clang++ -std=c++11 -o invert_H invert_H.o -L/Users/Michele/Dev/lib/ntl/lib -lntl -lm

rng: ./test/rng_test.cpp
	@clang++ -std=c++11 -o rng_test $?

encoder: ./classes/ldpc_encoder.cpp
	@clang++ -std=c++11 -c -I./classes $?

variable_node: ./classes/variable_node.cpp
	@clang++ -std=c++11 -c -I./classes $?

encoder_test: ./test/encoder_test.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/ldpc_encoder.cpp
	@clang++ -std=c++11 -c -I./classes -g $?
	@clang++ -std=c++11 -o encoder_test encoder_test.o ldpc_encoder.o

node_test: ./test/node_test.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/variable_node.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/check_node.cpp
	@clang++ -std=c++11 -c -I./classes -g $?
	@clang++ -std=c++11 -o node_test node_test.o variable_node.o check_node.o
 

clean:
	rm invert_H encoder_test rng_test *.o
