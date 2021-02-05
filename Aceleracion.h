#include<bits/stdc++.h>
#include"Presion.h"

using namespace std;

void Ace(bool ERes, int Ptype, int N, double m[], double R[], double D[], double P[], double A[]){
	int Max;
	double z;
	if(ERes==true){
	
	}else{
		Max=N;
	}
	for(int i=0; i<Max; i++){
		A[i]=0.0;
		for(int j=0; j<Max; j++){
			double Vj, Wi;
			if(Ptype==1){
				Vj=0;
			}
			if(Ptype==2){
				Vj=(1.0/2.0)*(R[j]*R[j]);
			}
			if(Ptype==3){
				Vj=0;
			}
				Wi=ker(R[i], R[j]);
				A[i]+=-m[j]*(P[i]/(D[i]*D[i])+P[j]/(D[j]*D[j])+Vj/D[j])*Wi;
		}
	}
}
void PrintAce(int N, double R[], double A[]){
	ofstream file("AvsR.dat");
	for(int i=0; i<N; ++i){
		file << R[i] << " " << A[i] << '\n';
	}
}







