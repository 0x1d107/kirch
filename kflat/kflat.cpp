#include <cmath>
#include <iostream>
#include <math.h>
#include <omp.h>

using namespace std;
double C=1.0;
double C2=0.1;

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
	return ( C2-C)/(C+C2);
}

int main(){
	double dt = 0.01;
	double dx = 0.01;
	double dy = 0.01;
	double z1 = 0.5;
	double z0 = 0;
	double sx = 0.5;
	double sy = 0.5;
	double ry = 0.5;
	double h=1e-5;
	for(int irt = 0; irt < 200;irt++){
		for(int irx = 0;irx < 100;irx++){
			double rt = irt*dt;
			double rx = irx*dx ;
			double sum = 0.0;
			#pragma omp parallel for reduction(+:sum)
			for(int iy = 0; iy <100; iy++){
				for(int ix = 0;ix < 100;ix++){
					double x = dx*ix;
					double y = dy*iy;
					sum+=Refl(ix,iy,dx,dy)*(K(z0,z1+h,sx,sy,rx,ry,rt,x,y) - K(z0,z1-h,sx,sy,rx,ry,rt,x,y))/2.0/h;
				}
			}
			cout << sum*dx*dy << ' ';
		} 
		cout << endl;
	}
}
