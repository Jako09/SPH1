#include<bits/stdc++.h>
#include"Densidad.h"

using namespace std;

double Press(bool ERes, int N, double m[], double R[], double D[], double Dx[], double Dxx[], double P[]){
	int Max;
	double z;
	if(ERes==true){
	
	}else{
		Max=N;
	}
	for(int i=0; i<Max; ++i){
		P[i]=0.0;
		for(int j=0; j<Max; ++j){
			double Q;
			Q=Dx[j]*Dx[j]/D[j]-Dxx[j];
			P[i]+=(m[j]/(D[j]*4.0))*ker(R[i], R[j])*(Q);
		}
	}
	return z;
}
void PrintPress(int N, double R[], double P[]){
	ofstream file("PvsR.dat");
	for(int i=0; i<N; ++i){
		file << R[i] << " " << P[i] << '\n';
	}
	file.close();
}




