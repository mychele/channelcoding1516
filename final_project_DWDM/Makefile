###################################################
# Makefile for Channel Coding 1516 final project - code only
###################################################

all: invert_H

invert_H: invert_H.cpp
	@g++ -std=c++11 -c -I/Users/Michele/Dev/lib/ntl/include $?
	@g++ -std=c++11 -o invert_H invert_H.o -L/Users/Michele/Dev/lib/ntl/lib -lntl -lm

clear:
	rm invert_H *.o
