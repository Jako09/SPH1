#include<bits/stdc++.h>
#include"Kernel.h"

using namespace std;


void Densidad0(bool ERes, int N, double m[], double R[], double D[]){
	int Max;

//	if(ERes==true){
		
//	}else{
		Max=N;
//	}
	for(int i=0; i<Max; i++){
		m[i]=1.0/(double)N;
	}
	for(int i=0; i<Max; i++){
		D[i]=0.0;
		for(int j=0; j<Max; j++){
				D[i]+=m[j]*ker(N, R[i], R[j]);
		}
		
	}
}

void Densidad1(bool ERes, int N, double m[], double R[], double D[], double Dx[]){
	int Max;
//	if(ERes==true){
		
//	}else{
		Max=N;
//	}
	for(int i=0; i<Max; i++){
		Dx[i]=0.0;
		for(int j=0; j<Max; j++){
//			double q=(D[j]-D[i])/D[j];
//			Dx[i]+=m[j]*q*dker(N, R[i], R[j]);
			Dx[i]+=m[j]*dker(N,R[i],R[j]);
		}
	}
}
void Densidad2(bool ERes, int N , double m[], double R[], double D[], double Dxx[]){
	int Max;
//	if(ERes==true){
	
//	}else{
		Max=N;
//	}
	for(int i=0; i<N; i++){
		Dxx[i]=0.0;
		for(int j=0; j<Max; j++){
			double q=(D[j]-D[i])/D[j];
			Dxx[i]+=m[j]*q*ddker(N, R[i], R[j]);
//			Dxx[i]+=m[j]*ddker(N, R[i], R[j]);
		}
	}
}	



void PrintDensidades( int N, double R[], double D[], double Dx[], double Dxx[]){
	ofstream file0("DvsR.dat");
	ofstream file1("DxvsR.dat");
	ofstream file2("DxxvsR.dat");
	for(int i=0; i<N; i++){
		file0 << R[i] << " " << D[i] << '\n';
		file1 << R[i] << " " << Dx[i] << '\n';
		file2 << R[i] << " " << Dxx[i] << '\n';
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
