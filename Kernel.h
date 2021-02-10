#include<bits/stdc++.h>
#define h0 100.0
using namespace std;


double ker(int N, double ri, double rj){
	double h=h0/(double)N;
	double alpha=(1.0/(h*sqrt(PI)));
	double R=(ri-rj)/h;
	return alpha*exp(-R*R);
}
double dker(int N, double ri, double rj){
	double h=h0/(double)N;
	double alpha=(1.0/(h*sqrt(PI))), beta=(2.0/h);
	double R=(ri-rj)/h;
	return alpha*exp(-R*R)*(-R*beta);
}
double ddker(int N, double ri, double rj){
	double h=h0/(double)N;
	double alpha=(1.0/(h*sqrt(PI))), beta=(2.0/h), delta=2.0/(h*h);
	double R=(ri-rj)/h;
	return alpha*exp(-R*R)*(R*R*beta*beta-delta);
}
