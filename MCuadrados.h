#include<bits/stdc++.h>
#include"L2.h"

using namespace std;



double MCPowb(int NError, int x[], double y[]){
	double b, c, d, e;
	for(int i=0; i < NError; i++){
		b+=log(x[i])*log(y[i]);
		c+=log(x[i]);
		d+=log(y[i]);
		e+=log(x[i])*log(x[i]);
	}	
	return ((double)NError*b-c*d)/((double)NError*e-c*c);
}

double MCPowa(int NError, int x[], double y[]){
	double b, c;
	for(int i=0; i<NError; ++i){
		b+=log(y[i]);
		c+=log(x[i]);
	}
	return (b-MCPowb(NError, x, y)*c)/(double)NError;
}

void MCuadrados(int Ptype, int NError){
	double L2A[NError];
	int NL2[NError];
	error(Ptype, NError, NL2, L2A);
	ofstream file("Error.dat");
	for(int i=0; i<NError; ++i){
		file << NL2[i] << " " << L2A[i] << '\n';
	}
	file.close();
	ofstream file0("Error.gnu");
	file0 << "p " << exp(MCPowa(NError, NL2, L2A)) << "*x**("<< MCPowb(NError, NL2, L2A) << "), \"Error.dat\" w lp" <<  '\n';
	file0 << "set ylabel \"L2(Error de la densidad)\"" << '\n';
	file0 << "set xlabel \"N-Numero de particulas\"" << '\n';
	file0 << "set title \"Errores\"" << '\n';
	file0 << "set logscale y" << '\n';
	file0 << "set logscale x" << '\n';
	file0 << "replot" << '\n';
	file0 << "pause -1" << '\n'; 
	file0.close();
}

