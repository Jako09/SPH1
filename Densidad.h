#include<bits/stdc++.h>
#include"Kernel.h"

using namespace std;


double Densidad0(bool ERes, double R[], int i, int N){
	double z;
	int Max;
	if(ERes==true){
		
	}else{
		Max=N;
	}
	double m=1.0/(double)N;
	for(int j=0; j<Max; ++j){
			z+=m*ker(R[i], R[j]);
	}
	return z;
}

void Densidad(double R[], int N, double D[], bool ERes){
	double z;
	for(int i=0; i<N; i++){
		z=Densidad0(ERes, R, i, N);
		D[i]=z;
		cout << D[i] << '\n'; 
	}
	for(int i=0; i<N; i++){
		cout << D[i] << '\n';
	}
}

double Densidad1(bool ERes, double D[], double R[], int i, int N){
	double z;
	int Max;
	if(ERes==true){
		
	}else{
		Max=N;
	}
	double m=1.0/(double)N;
	for(int j=0; j<Max; ++j){
		double q=(D[j]-D[i])/D[j];
		z+=m*q*dker(R[i], R[j]);
	}
	return z;
}
double Densidad2(bool ERes, double D[], double R[], int i, int N){
	double z;
	int Max;
	if(ERes==true){
	
	}else{
		Max=N;
	}
	double m=1.0/(double)N;
	for(int j=0; j<Max; ++j){
		double q=(D[j]-D[i])/D[j];
		z+=m*q*ddker(R[i], R[j]);
	}
	return z;
}	



void PrintDensidades(double R[], double D[], int N, bool ERes){
	ofstream file0("DvsR.dat");
	ofstream file1("DxvsR.dat");
	ofstream file2("DxxvsR.dat");
	for(int i=0; i<N; ++i){
		file0 << R[i] << " " << D[i] << '\n';
		file1 << R[i] << " " << Densidad1(ERes, D, R, i, N) << '\n';
		file2 << R[i] << " " << Densidad2(ERes, D, R, i, N) << '\n';
	}
	file2.close();
	file1.close();
	file0.close();
}
/*
double Den(vector<double> W){
	double m=1.0/NPV, b2;
	for(int j=0; j<W.size();++j){
		b2+=m*W.at(j);
	}
	return b2;
}
double Denx(vector<double> Wx){
	double m=1.0/NPV, b2;
	for(int j=0; j<Wx.size(); ++j){
		b2+=m*Wx.at(j);
	}
	return b2;
}
double Denxx(vector<double> Wxx){
	double m=1.0/NPV, b2;
	for(int j=0; j<Wxx.size(); ++j){
		b2+=m*Wxx.at(j);
	}
	return b2;
}
*/
