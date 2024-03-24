all: seism.bin

seism.bin: ktilt
	./ktilt 
ktilt: ktilt.cpp
	g++ -fopenmp $^ -o $@
