#include <cmath>
#include <iostream>
#include <math.h>
#ifdef OMP
#include <omp.h>
#endif

using namespace std;
double C=1.0;
double C2=0.1;
double EPS = 1e-5;
// BA*x+BB*y+BC*z=BD
double BA = 1;
double BB = 1;
double BC = -1;
double BD = 2;
typedef struct {double x[3];} vec3_t;
double vec3_rn2(vec3_t rp,vec3_t r){
	
	double rn2 = 0;
	for(int i=0;i<3;i++){
		rn2+=(rp.x[i]-r.x[i])*(rp.x[i]-r.x[i]);
	}
	return rn2;
}
double rick(double t,double f){
	return (1-2*M_PI*M_PI*f*f*t*t)*exp(-M_PI*M_PI*f*f*t*t);
}
double rickdel(double t,double f,vec3_t rp, vec3_t r, vec3_t rs){
	return rick(t - sqrt(vec3_rn2(r,rs))/C - sqrt(vec3_rn2(rp,r))/C,f);
}
double dnrickdel(double t,double f,vec3_t n,vec3_t rp, vec3_t r, vec3_t rs){
	double s=0;
	for(int i=0;i<3;i++){
		vec3_t rp2,rp0;
		rp0=rp2=rp;
		rp0.x[i]-=EPS;
		rp2.x[i]+=EPS;
		double dkrick = (rickdel(t,f,rp2,r,rs) - rickdel(t,f,rp0,r,rs))/2/EPS;
		s+=dkrick*n.x[i];
	}
	return s;
}
double dtrickdel(double t,double f,vec3_t rp, vec3_t r, vec3_t rs){
	double t2 = t+EPS;
	double t0 = t-EPS;
	return (rickdel(t2,f,rp,r,rs) - rickdel(t2,f,rp,r,rs))/2/EPS;
}
double K(double z0,double z1,double sx,double sy,double rx,double ry,double rt,double x,double y){
	double as = 1.0/M_PI/4.0/sqrt(pow(x-sx,2)+pow(y-sy,2)+pow(z1-z0,2));
	double ar = 1.0/M_PI/4.0/sqrt(pow(x-rx,2)+pow(y-ry,2)+pow(z1-z0,2));
	double tau_s = sqrt(pow(x-sx,2)+pow(y-sy,2)+pow(z1-z0,2))/C;
	double tau_r = sqrt(pow(x-rx,2)+pow(y-ry,2)+pow(z1-z0,2))/C;
	
	return ar*as*rick(rt - tau_r - tau_s,10.0);
}
double Refl(double ix,double iy,double dx,double dy){
	return (C - C2)/(C+C2);
}
double dndr(vec3_t rp,vec3_t r,vec3_t n){
    double s = 0;
	double rn2 = 0;
	for(int i=0;i<3;i++){
		rn2+=(rp.x[i]-r.x[i])*(rp.x[i]-r.x[i]);
	}
    for(int k=0;k<3;k++){
        s+= (rp.x[k] - r.x[k])/sqrt(rn2)*n.x[k]; //TODO
    }
	return s;
    
}
double dnidr(vec3_t rp,vec3_t r,vec3_t n){
    double s = 0;
	double rn2 = 0;
	for(int i=0;i<3;i++){
		rn2+=(rp.x[i]-r.x[i])*(rp.x[i]-r.x[i]);
	}
    for(int k=0;k<3;k++){
        s+= -(rp.x[k] - r.x[k])/sqrt(rn2)/rn2*n.x[k]; //TODO
    }
	return s;
    
}
int main(){
	double dt = 0.01;
	double dx = 0.01;
	double dy = 0.01;
	double z1 = 0.5;
	double z0 = 0;
	double sx = 0.5;
	double sy = 0.5;
	vec3_t sr = {{sx,sy,0}};
	double ry = 0.0;
	double h=1e-5;
	vec3_t Bn = {{BA,BB,BC}};
	double Bnl = sqrt(BA*BA+BB*BB+BC*BC);
	for(int k=0;k<3;k++){
		Bn.x[k]/=Bnl;
	}
	for(int irt = 0; irt < 200;irt++){
		for(int irx = 0;irx < 100;irx++){
			double rt = irt*dt;
			double rx = irx*dx ;
			vec3_t rp = {{rx,ry,0}};
			double sum = 0.0;
			#pragma omp parallel for reduction(+:sum)
			for(int iy = 0; iy <100; iy++){
				for(int ix = 0;ix < 100;ix++){
					double x = dx*ix;
					double y = dy*iy;
					double z = (BD-BA*x-BB*y)/BC;
					vec3_t r = {{x,y,z}};
					double K1 = dndr(rp,r,Bn)*dtrickdel(rt,15,rp,r,sr)/C/sqrt(vec3_rn2(rp,r));
					double K2 = -dnidr(rp,r,Bn)*rickdel(rt,15,rp,r,sr);
					double K3 = dnrickdel(rt,15,Bn,rp,r,sr)/sqrt(vec3_rn2(r,rp));
					sum+=K1+K2+K3;
				}
			}
			cout << sum*dx*dy/4/M_PI << ' ';
		} 
		cout << endl;
	}
}
