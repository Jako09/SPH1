#include<bits/stdc++.h>
#include"InitialD.h"
#include"Aceleracion.h"

using namespace std;

void printAll(int Ptype, double R[], int N, bool ERes){
	PrintDensidades(R, N, ERes);
	PrintPress(R, N, ERes);
	PrintAce(Ptype, R, N, ERes);
}

int main(){
	int N=100, Ptype=3;
	double R[N], Na;
	Na=(double)N;
	InitialD(Ptype, Na, 100.0, R);
	printAll(Ptype, R, N, false);

	return 0;
}


