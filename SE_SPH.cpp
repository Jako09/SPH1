#include<bits/stdc++.h>
#include"InitialD.h"
#include"Aceleracion.h"
#include"MCuadrados.h"

using namespace std;

void printAll(int N, double R[], double D[], double Dx[], double Dxx[], double P[], double A[]){
	PrintDensidades(N, R, D, Dx, Dxx);
	PrintPress(N, R, P);
	PrintAce(N, R, A);
}

void SPH(int Ptype, int N){
	double m[N], R[N], D[N], Dx[N], Dxx[N], P[N], A[N];
	double Na, Nb;
	bool ERes=false, NPV=false;
	Na=(double)N;
	if(NPV==false){
		Nb=(double)N;
	}else{
	}
	InitialD(Ptype, Na, Nb, R);
	Densidad0(ERes, N, m, R, D);
	Densidad1(ERes, N, m, R, D, Dx);
	Densidad2(ERes, N, m, R, D, Dxx);
	Press(ERes, N, m, R, D, Dx, Dxx, P);
	Ace(ERes, Ptype, N, m, R, D, P, A);
	printAll(N, R, D, Dx, Dxx, P, A);
}

int main(){
	int Ptype=2, N=100, NError=20;
	MCuadrados(Ptype, NError);
	return 0;
}


