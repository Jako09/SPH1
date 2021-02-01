#include<bits/stdc++.h>
#include"Densidad.h"

using namespace std;
double basePress(int i, double R[], int j, int N, bool ERes){
	double m=1.0/(double)N;
	double Dx, Dxx, Q, Rhok;
	Dx=Densidad(1, ERes, R, j, N);
	Dxx=Densidad(2, ERes, R, j, N);
	Rhok=Densidad(0, ERes, R, j, N);
	Q=Dx*Dx/Rhok-Dxx;
		
	return	(m/(Rhok*4.0))*ker(R[i], R[j])*(Q);
}
double Press(int i, double R[], int N, bool ERes){
	int Max;
	double z;
	if(ERes==true){
	
	}else{
		Max=N;
	}
	for(int j=0; j<Max; ++j){
		if(i!=j){
			z+=basePress(i, R, j, N, ERes);
		}
	}
	return z;
}
void PrintPress(double R[], int N, bool ERes){
	ofstream file("PvsR.dat");
	for(int i=0; i<N; ++i){
		file << R[i] << " " << Press(i, R, N, ERes) << '\n';
	}
	file.close();
}




