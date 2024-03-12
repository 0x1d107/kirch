OPENMP=-fopenmp -DOMP
all: data.png

data.png: data.txt plot.sh
	bash plot.sh data.txt data.png
data.txt: ktilt
	./ktilt > $@
ktilt: ktilt.cpp
	g++ $(OPENMP) $^ -o $@
