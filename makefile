make: Code/encoders.cpp
	g++ -g -O2 -std=c++11 -pthread -march=native Code/encoders.cpp -o encoders -lntl -lgmp -lm -lgsl -lgslcblas