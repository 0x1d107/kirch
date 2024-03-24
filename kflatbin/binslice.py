#!/bin/env python
import sys,struct
Nx = 50
Ny = 50
Nt = 100
it = Nt//2
if len(sys.argv) >1:
	it = int(sys.argv[1])
double_size=8
fmt = f'{Nx*Ny}d'
with open('seism.bin','rb') as f:
	f.seek(double_size*Nx*Ny*it)
	buf = f.read(struct.calcsize(fmt))
	data= struct.unpack(fmt,buf)
	for j in range(Ny):
		for i in range(Nx):
			print(data[Nx*j+i],end=' ')
		print()
