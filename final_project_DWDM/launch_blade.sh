#!/bin/bash

g++ -std=c++11 -O3 -c -I./classes -g ./classes/variable_node.cpp
g++ -std=c++11 -O3 -c -I./classes -g ./classes/check_node.cpp
g++ -std=c++11 -O3 -c -I./classes -g ./classes/ldpc_decoder.cpp
g++ -std=c++11 -O3 -c -I./classes -g ./test/simulation_multi.cpp
g++ -std=c++11 -O3 -o simulation_multi simulation_multi.o ldpc_decoder.o variable_node.o check_node.o

./simulation_multi