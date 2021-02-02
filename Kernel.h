#include<bits/stdc++.h>
#define h 0.05
using namespace std;


double ker(double ri, double rj){
	double alpha=(1.0/(h*sqrt(PI)));
	double R=(ri-rj)/h;
	return alpha*exp(-R*R);
}
double dker(double ri, double rj){
	double alpha=(1.0/(h*sqrt(PI))), beta=(2.0/h);
	double R=(ri-rj)/h;
	return alpha*exp(-R*R)*(-R*beta);
}
double ddker(double ri, double rj){
	double alpha=(1.0/(h*sqrt(PI))), beta=(2.0/h), delta=2.0/(h*h);
	double R=(ri-rj)/h;
	return alpha*exp(-R*R)*(R*R*beta*beta-delta);
}
