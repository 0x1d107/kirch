all: data.png

data.png: data.txt plot.sh
	bash plot.sh data.txt data.png
data.txt: kflat
	./kflat > $@
kflat: kflat.cpp
	g++ -fopenmp $^ -o $@
