#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
double C=1.0;
double C2=0.1;

double rick(double t,double f){
	return (1-2*M_PI*M_PI*f*f*t*t)*exp(-M_PI*M_PI*f*f*t*t);
}
double Refl(double ix,double iy,double dx,double dy){
	return ( C2-C)/(C+C2);
}

double K(std::vector<double> &seism,double xp,double yp,double zp,double dx,double dy,double dt,int Nx,int Ny,int Nt){
	double sum = 0;
	#pragma omp parallel for reduction(+:sum)
	for(int iys =0;iys<Ny;iys++){
		for(int ixs =0 ;ixs<Nx;ixs++){
			double xs = ixs * dx;
			double ys = iys * dy;
			double zs = 0;
			double L  = sqrt(pow(xp - xs,2)+pow(yp - ys,2)+pow(zp - zs,2));
			double st = L/C;
			int it = min(max((int)(st/dt),0),Nt-1);
			double P = seism[Nx*Ny*it+Nx*iys+ixs];
			sum+=P/L;

		}
	} 
	return sum;
	
}
int main(){
	int Nx = 50;
	int Ny = 50;
	int Nt = 100;
	double dt = 2.0/Nt;
	double dz = 2.0/Nt;
	double dx = 1.0/Nx;
	double dy = 1.0/Ny;
	double z1 = 0.5;
	double z0 = 0;
	double sx = 0.5;
	double sy = 0.5;
	double ry = 0.5;
	double h=1e-5;
	std::ifstream inpf("seism.bin",std::ios::binary|std::ios::in);
	if(!inpf.is_open()){
		std::cerr<<"Can't open file seism.bin"<<std::endl;
		return 1;
	}
	std::vector<double> seismvec(Nx*Ny*Nt,0); 
	for(int i=0;i<Nx*Ny*Nt;i++){
		inpf.read((char*)&seismvec[i],sizeof(double));
	}
	for(int iz = 0; iz < Nt;iz++){
		for(int iy = 0; iy<Ny;iy++){
			for(int ix = 0;ix < Nx;ix++){
				double zp = iz*dz;
				double zp0 = zp - h;
				double zp2 = zp + h;
				double xp = ix*dx ;
				double yp = iy*dy;
				double sum0 = K(seismvec,xp,yp,zp0,dx,dy,dt,Nx,Ny,Nt);
				double sum2 = K(seismvec,xp,yp,zp2,dx,dy,dt,Nx,Ny,Nt);
				double sum = (sum2-sum0)/2/h;
				std::cout << sum<<' ';

			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}
