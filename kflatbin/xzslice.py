#!/bin/env python
import sys,struct
Nx = 50
Ny = 50
Nt = 100
it = 25
iy = 20
if len(sys.argv) >1:
	it = int(sys.argv[1])
double_size=8
fmt = f'{Nx}d'
with open('seism.bin','rb') as f:
	for it in range(Nt):
		f.seek(double_size*(Nx*Ny*it+iy*Nx))
		buf = f.read(struct.calcsize(fmt))
		data= struct.unpack(fmt,buf)
		print(*data)
