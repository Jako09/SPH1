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
		if(i>0.25*N && i<0.75*N){
		double e=D[i]-DA(Ptype, R[i]);
		z+=e*e;
		}
	}	
	return sqrt(z);
}
void error(int Ptype, int NError){
	int N=200000;
	ofstream file("Error.dat");
	for(int i=0; i<NError; ++i){
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
		file << N << " " << L2(Ptype, N, R, D) << '\n';
		cout << N << '\n';
		N+=10000;
	}
	file.close();
}




