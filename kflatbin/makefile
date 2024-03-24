all: slice.png xz.png

%.png: %.txt
	bash plot.sh $< $@
xz.txt: seism.bin makefile xzslice.py
	python xzslice.py > $@ 
slice.txt: seism.bin makefile binslice.py
	python binslice.py > $@ 
seism.bin: kflat
	./kflat
kflat: kflat.cpp
	g++ -fopenmp $^ -o $@
