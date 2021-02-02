#include<bits/stdc++.h>
#include"Presion.h"

using namespace std;
double baseAce(int Ptype, double R[], double D[], int i, int j, int N, bool ERes){
	double m=1.0/(double)N;
	double Di, Dj, Vj, Pi, Pj, Wi;
	Di=D[i];
	Dj=D[j];
	if(Ptype==1){
		Vj=0;
	}
	if(Ptype==2){
		Vj=(1.0/2.0)*(R[j]*R[j]);
	}
	if(Ptype==3){
		Vj=0;
	}
	Pi=Press(i, R, D, N, ERes);
	Pj=Press(j, R, D, N, ERes);
	Wi=ker(R[i], R[j]);
	return m*(Pi/(Di*Di)+Pj/(Dj*Dj)+Vj/Dj)*Wi;
	
}
double Ace(int Ptype, double R[], double D[], int i, int N, bool ERes){
	int Max;
	double z;
	if(ERes==true){
	
	}else{
		Max=N;
	}
	for(int j=0; j<Max; j++){
			z+=-baseAce(Ptype, R, D, i, j, N, ERes);
	}
	return z;
}
void PrintAce(int Ptype, double R[], double D[], int N, bool ERes){
	ofstream file("AvsR.dat");
	for(int i=0; i<N; ++i){
		file << R[i] << " " << Ace(Ptype, R, D, i, N, ERes) << '\n';
	}
}







