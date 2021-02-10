#include<bits/stdc++.h>


using namespace std;

double NEficiente(int N, double R[]){
	double NEffA[N];
	double NEff;
	for(int i=0; i<N; ++i){
		for(int j=0; j<N; ++j){
			double o=pow(1.0,-14.0);
			double k=ker(N, R[i], R[j]);
			if(k>o){
				NEffA[i]+=1;
			}
		}
	}
	for(int i=0; i<N; ++i){
		NEff+=NEffA[i];
	}
	return NEff/(double)N;
}

void PrintNEff(int Ptype, int NError){
	int N=1000;
	ofstream file("NEff.dat");
	for(int i=0; i<NError; ++i){
		double m[N], R[N];
		double Na, Nb, h1;
		bool ERes=false, NPV=false;
		Na=(double)N;
		h1=100.0/(double)N;
		if(NPV==false){
			Nb=(double)N;
		}else{
		}
		InitialD(Ptype, Na, Nb, R);
		file << N << " " << NEficiente(N, R) << '\n';
		N=N*2;
	}
	file.close();
}
