#include<bits/stdc++.h>
#include<vector>
#define PI 3.1415926535897932 
#define EPSILON 0.0001
#define a 5.0

using namespace std;
/*Para esta función será necesario pedir el tipo de problema (typeP), el número de partículas (N), con una boolean (htype) que determine el tipo de parámetro de suavizado a usar, puede ser estático o variable*/
double NR1(double x, double NPV, int typeP);
double NR2(double x, int typeP);
double aux1(double x, double NPV, int typeP);
double aux2(double x1, double x2, double NPV, int typeP);
double NR3(int i, double x1, double x2, double NPV, int typeP);
double dis0(int i, double r1, double r2, double N, double NPV, int typeP);

void InitialD(int typeP, double N, double NPV, double R[]){
	double IDR[(int)N];
	double IDR2[(int)N];
	double z;
	if(typeP==1){
		z=2.5;
	}if(typeP==2){
		z=0.0;
	}
	for(int i=0; i<N;i++){
		if(typeP==3){
				IDR[i]=1.0/(2.0*N)+(double)i/N;
				R[i]=IDR[i];
		}else{
			if(i==0){
				IDR[i]=NR3(i, z, 0, NPV, typeP);
				IDR2[i]=NR3(i, z, 0, NPV, typeP);
			}else{
				IDR[i]=NR3(i, z, NR1(IDR[i-1], NPV, typeP), NPV, typeP);
				IDR2[i]=NR3(i, z, NR1(IDR[i-1], NPV, typeP), NPV, typeP);
			}
		}
	}

	for(int i=0; i<N; i++){
		if(typeP==3){	
		}else{
			if(i==0){
				IDR[i]=dis0(i, IDR2[i], 0, N, NPV, typeP);
				R[i]=IDR[i];
			}else{
				IDR[i]=dis0(i, IDR2[i], IDR2[i-1], N, NPV, typeP);
				R[i]=IDR[i];
			}
		}
	}
}

void PrintPosicion(double R[], int N){
	ofstream file1("RvsIh01t0.dat");
	for(int i=0; i<N; i++){
		file1 << i << " " << R[i] << '\n';
	}
	file1.close();
}


double NR1(double x, double NPV, int typeP){ 
	double z;
	if(typeP == 1){
		z =(1.0/a)*(x)-(1.0/(2.0*PI))*sin(2.0*PI*x/a)-1.0/NPV;  /*N-NPV*/
	}
	if(typeP == 2){
		z =(0.5)*erf(x)-1.0/NPV;
	}
	if(typeP == 3){
		z =0.0;
	}
	return z;
}

double NR2(double x, int typeP){
	double z;
	if(typeP == 1){
		z= (1.0/a)*(1.0-cos(2.0*PI*x/a));
	}
	if(typeP == 2){
		z= (1.0/sqrt(PI))*exp(-x*x);		
	}
	if(typeP == 3){
		z= 0.0;
	}
	return z;
}
double aux1(double r, double NPV, int typeP){
	double z;
	if(typeP == 1){
		z=NR1(r, NPV, typeP)/NR2(r, typeP);
	}
	if(typeP == 2){
		z=(NR1(r,NPV, typeP)-(0.5)*erf(-10.0))/NR2(r, typeP);
	}
	if(typeP == 3){
		z= 0.0;
	}
	return z;
}
double aux2(double r, double fi, double NPV,int typeP){
	return (NR1(r, NPV, typeP)-(fi+1.0/NPV))/(NR2(r, typeP));
}
//-------------------------Newton-Rapson Method---------------------------------------------------------------------------
double NR3(int i, double r, double fi, double NPV, int typeP){
	double z;
		if(i==0){
			if(typeP == 1){
				double e= aux1(r, NPV, typeP); /*Para la primera posición*/
				while(abs(e) >= EPSILON){
					e=aux1(r, NPV, typeP);
					r=r-e;
				}
			z= r;
			}
			if(typeP == 2){
				double e= aux1(r, NPV, typeP); /*Para la primera posición*/
				while(abs(e) >= EPSILON){
					e=aux1(r, NPV, typeP);
					r=r-e;
				}
				z= r;
			}
			if(typeP == 3){
				z= 0.0;
			}
		}else{	/*Para las posiciones siguientes, es diferente debido a los límites de integración*/
			double e= aux2(r, fi, NPV, typeP);
			while(abs(e)>=EPSILON){
				e=aux2(r, fi, NPV, typeP); /*funct(p[i-1].r)--fi1*/
				r=r-e;
			}	
			z= r;
		}
	return z;
}

double dis0(int i, double r, double rim, double N, double NPV, int typeP){
	double b1;
	if(N==NPV){
			if(i==0){
				if(typeP==1){
					b1=r/2.0;  /*Para la primera posición*/
				}
				if(typeP==2){
					b1=r;
				}
				if(typeP==3){
					return 0.0;		
				}
			//	file << i << " " << p[i].r << '\n';
			}else{				/*Para las demás posiciones*/
				b1=(r+rim)/2.0; /*surge de tomar en cuenta (b[i].r-b[i-1].r)/2+b[i-1].r*/
			//	file << i << " " << p[i].r << '\n';
			}
	}else{
		if(i<(int)NPV){			/*N-NPV*/
			if(i==0){
				b1=r/2.0-a;  /*Para la primera posición*/
			//	file << i << " " << p[i].r << '\n';
			}else{				/*Para las demás posiciones*/
				b1=(r+rim)/2.0-a; /*surge de tomar en cuenta (b[i].r-b[i-1].r)/2+b[i-1].r*/
			//	file << i << " " << p[i].r << '\n';
			}	
		}
		/*----Partículas estaticas-------*/
		if(i >= (int)NPV  && i<(int)(2.0*NPV)){			/*N-NPV*/
			b1=rim+a;
		}
		
		if(i >= (int)(2.0*NPV) && i<(int)(3.0*NPV)){			/*N-NPV*/
			b1=rim+2.0*a;
		}
		
	}
	return b1;
}


