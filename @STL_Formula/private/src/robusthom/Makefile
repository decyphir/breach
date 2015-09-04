OBJ= robustness.o signal.o

all: test_rob testSignal testRobustness

robustness.o: robustness.cpp robustness.h
	g++ -c -O2 robustness.cpp

signal.o: signal.cpp signal.h
	g++ -c -O2 signal.cpp

test_rob: robustness.o signal.o test.cpp
	g++ -o test test.cpp robustness.o signal.o

testSignal: $(OBJ) testSignal.cpp
	g++ -o testSignal testSignal.cpp $(OBJ)

testRobustness: $(OBJ) testRobustness.cpp
	g++ -o testRobustness testRobustness.cpp $(OBJ)

