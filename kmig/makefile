all:xz.png mig.png

xz.png:xz.txt
	bash plot.sh $< $@
mig.png:mig.txt
	bash plotidx.sh $< $@
xz.txt: mig.txt makefile
	awk 'NR%101==50' $< > $@
mig.txt: kmig seism.bin
	./kmig >mig.txt
kmig: kmig.cpp
	$(CXX) -o $@ -fopenmp $^ 
