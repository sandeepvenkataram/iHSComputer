CC=g++


iHS:
	$(CC) -O3 -o bin/iHS src/iHSComputer_singlethreaded.cpp src/iHSComputer.h src/iHSComputer.cpp
iHS_multithreaded:
	$(CC) -std=c++0x -pthread -o bin/iHS_multithreaded src/iHSComputer_multithreaded.cpp src/iHSComputer.h src/iHSComputer.cpp
normalizer:
	$(CC) -O3 -o bin/normalizer src/normalization.cpp
clean:
	rm -f bin/iHS bin/iHS_multithreaded bin/normalizer

 
