###################################################
# Makefile for Channel Coding 1516 final project - code only
###################################################

all: invert_H

invert_H: invert_H.cpp
	@clang++ -std=c++11 -c -I/Users/Michele/Dev/lib/ntl/include $?
	@clang++ -std=c++11 -o invert_H invert_H.o -L/Users/Michele/Dev/lib/ntl/lib -lntl -lm

rng: normal_rng_test.cpp
	@clang++ -std=c++11 -o rng_test $?

debug: invert_H.cpp
	@clang++ -std=c++11 -c -I/Users/Michele/Dev/lib/ntl/include -g $?
	@clang++ -std=c++11 -o invert_H invert_H.o -L/Users/Michele/Dev/lib/ntl/lib -lntl -lm

clear:
	rm invert_H *.o
