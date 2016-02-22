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

decoder: ./classes/ldpc_decoder.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/variable_node.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/check_node.cpp
	@clang++ -std=c++11 -c -I./classes -g $?

variable_node: ./classes/variable_node.cpp
	@clang++ -std=c++11 -c -I./classes $?

encoder_test: ./test/encoder_test.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/ldpc_encoder.cpp
	@clang++ -std=c++11 -c -I./classes -g $?
	@clang++ -std=c++11 -o encoder_test encoder_test.o ldpc_encoder.o

node_test: ./test/node_test.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/variable_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/check_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g $?
	@clang++ -std=c++11 -O3 -o node_test variable_node.o check_node.o node_test.o

decoder_test: ./test/decoder_test.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/variable_node.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/check_node.cpp
	@clang++ -std=c++11 -c -I./classes -g ./classes/ldpc_decoder.cpp
	@clang++ -std=c++11 -c -I./classes -g $?
	@clang++ -std=c++11 -o decoder_test decoder_test.o ldpc_decoder.o variable_node.o check_node.o


dwdm: ./test/dwdm_encode_decode.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/variable_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/check_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/ldpc_decoder.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/ldpc_encoder.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g $?
	@clang++ -std=c++11 -O3 -o dwdm_encode_decode dwdm_encode_decode.o ldpc_decoder.o variable_node.o check_node.o ldpc_encoder.o

sim: ./test/simulation.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/variable_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/check_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/ldpc_decoder.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/ldpc_encoder.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g $?
	@clang++ -std=c++11 -O3 -o simulation simulation.o ldpc_decoder.o variable_node.o check_node.o ldpc_encoder.o

sim_noenc: ./test/simulation_noenc.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/variable_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/check_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/ldpc_decoder.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g $?
	@clang++ -std=c++11 -O3 -o simulation_noenc simulation_noenc.o ldpc_decoder.o variable_node.o check_node.o

sim_multi: ./test/simulation_multi.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/variable_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/check_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/ldpc_decoder.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g $?
	@clang++ -std=c++11 -O3 -o simulation_multi simulation_multi.o ldpc_decoder.o variable_node.o check_node.o

sim_minsum: ./test/simulation_minsum.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/variable_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/check_node.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g ./classes/ldpc_decoder.cpp
	@clang++ -std=c++11 -O3 -c -I./classes -g $?
	@clang++ -std=c++11 -O3 -o simulation_minsum simulation_minsum.o ldpc_decoder.o variable_node.o check_node.o
 

clean:
	rm invert_H encoder_test rng_test dwdm_encode_decode simulation decoder_test node_test simulation_multi simulation_noenc *.o
