#include <cmath>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>


using namespace std;
double C=1.0;
double C2=0.5;

double rick(double t,double f){
	return (1-2*M_PI*M_PI*f*f*t*t)*exp(-M_PI*M_PI*f*f*t*t);
}
double K(double z0,double z1,double sx,double sy,double rx,double ry,double rt,double x,double y){
	double as = 1.0/M_PI/4.0/sqrt(pow(x-sx,2)+pow(y-sy,2)+pow(z1-z0,2));
	double ar = 1.0/M_PI/4.0/sqrt(pow(x-rx,2)+pow(y-ry,2)+pow(z1-z0,2));
	double tau_s = sqrt(pow(x-sx,2)+pow(y-sy,2)+pow(z1-z0,2))/C;
	double tau_r = sqrt(pow(x-rx,2)+pow(y-ry,2)+pow(z1-z0,2))/C;
	
	return ar*as*rick(rt - tau_r - tau_s,10.0);
}
double Refl(double ix,double iy,double dx,double dy){
	return ( C-C2)/(C+C2);
}
struct vec2d{
	double x;
	double y;
};

int main(){
	int Nx=50;
	int Ny=50;
	int Nt=100;
	double dt = 2.0/Nt;
	double dx = 1.0/Nx;
	double dy = 1.0/Ny;
	double z1 = 0.5;
	double z0 = 0;
	std::vector<vec2d> ss = {{0.5,0.5}};
	/*
	for(int j=0;j<5;j++)
		for(int i=0;i<5;i++)
			ss.push_back({i/5.0+1.0/10,j/5.0+1.0/5});
	*/	
	double sx = 0.5;
	double sy = 0.5;
	double h=1e-5;
	ofstream outf("seism.bin",std::ios::out |std::ios::binary);
	for(int irt = 0; irt < Nt;irt++){
		for(int iry = 0; iry<Ny;iry++){
			for(int irx = 0;irx < Nx;irx++){
				double rt = irt*dt;
				double rx = irx*dx ;
				double ry = iry*dy;
				double sum = 0.0;
				for(const vec2d &vs: ss){
					double sx = vs.x;
					double sy = vs.y;
					#pragma omp parallel for reduction(+:sum)
					for(int iy = 0; iy <Ny; iy++){
						for(int ix = 0;ix < Nx;ix++){
							double x = dx*ix;
							double y = dy*iy;
							sum+=Refl(ix,iy,dx,dy)*(K(z0,z1+h,sx,sy,rx,ry,rt,x,y) - K(z0,z1-h,sx,sy,rx,ry,rt,x,y))/2.0/h;
						}
					}
				}
				sum *= dx*dy/2/M_PI;
				outf.write(reinterpret_cast<char*>(&sum),sizeof(double));
			} 
		}
		std::cout << "irt = "<<irt<<std::endl;
	}
	outf.close();
}
