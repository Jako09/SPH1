#include<bits/stdc++.h>

using namespace std;
double DA(int Ptype, double Ri){
        double z;
	if(Ptype==1){
		z=(2.0/a)*sin(PI*Ri/a)*sin(PI*Ri/a);
	}
	if(Ptype==2){
		z=(1.0/sqrt(PI))*exp(-Ri*Ri);
	}
	if(Ptype==3){
		z=1.0;
	}
	return z;
}
double L2(int Ptype, int N, double R[], double D[]){
	double z;
	for(int i=0; i<N; ++i){
		double e=D[i]-DA(Ptype, R[i]);
		z+=e*e;
	}	
	return sqrt(z);
}
void error(int Ptype, int NError, int NL2[], double L2A[]){
	int N=0;
	for(int i=0; i<NError; ++i){
		NL2[i]=0.0;
		L2A[i]=0.0;
	}
	for(int i=0; i<NError; ++i){
		N+=100;
		NL2[i]=N;
		double m[N], R[N], D[N];
		double Na, Nb;
		bool ERes=false, NPV=false;
		Na=(double)N;
		if(NPV==false){
			Nb=(double)N;
		}else{
		}
		InitialD(Ptype, Na, Nb, R);
		Densidad0(ERes, N, m, R, D);
		L2A[i]=L2(Ptype, N, R, D);
	}
}




