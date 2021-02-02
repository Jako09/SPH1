#include<bits/stdc++.h>
#include"InitialD.h"
#include"Aceleracion.h"

using namespace std;

void printAll(int Ptype, double D[], double R[], int N, bool ERes){
	PrintDensidades(R, D, N, ERes);
	PrintPress(R, D, N, ERes);
	PrintAce(Ptype, R, D, N, ERes);
}



int main(){
	int N=100, Ptype=3;
	double R[N], Na;
	double D[N];
	bool ERes=false;
	Na=(double)N;
	InitialD(Ptype, Na, 100.0, R);
	Densidad(R, N, D, ERes);
//	for(int i=0; i<N; ++i){
//		cout << D[i] << '\n';
//	}
	printAll(Ptype, D, R, N, ERes);

	return 0;
}


