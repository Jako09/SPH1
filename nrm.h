#include<bits/stdc++.h>
#define flag 2
#define NPV 100.0
#define PI 3.1415926535897932 
#define EPSILON 0.0001

using namespace std;

double NR1(double x){ /*Definimos la ecuación que se va a resolver en el método Newton-Raphson*/
//	return x*x*x-x*x+2;
	if(flag == 1){
		return (1.0/a)*(x)-(1.0/(2.0*PI))*sin(2.0*PI*x/a)-1.0/NPV;  /*N-NPV*/
	}
	if(flag == 2){
		return (0.5)*erf(x)-1.0/NPV;
	}
	if(flag == 3){
		return 0;
	}
}
double NR2(double x){
//	return 3*x*x-2*x; /*La derivada de la ecuación a encontrar para el método Newton-Raphson*/
	if(flag == 1){
		return (1.0/a)*(1.0-cos(2.0*PI*x/a));
	}
	if(flag == 2){
		return (1.0/sqrt(PI))*exp(-x*x);		
	}
	if(flag == 3){
		return 0;
	}
}
double aux1(double r){
	if(flag == 2){
		return (NR1(r)-(0.5)*erf(-4.0))/NR2(r);
	}
	if(flag == 1){
		return NR1(r)/NR2(r);
	}
		if(flag == 3){
		return 0;
	}
}
double aux2(double r, double fi){
	return (NR1(r)-(fi+1.0/NPV))/(NR2(r));
}
//-------------------------Newton-Rapson Method---------------------------------------------------------------------------
double NR3(int i, double r, double fi){
		if(i==0){
			if(flag == 1){
				double e= aux1(r); /*Para la primera posición*/
				while(abs(e) >= EPSILON){
					e=aux1(r);
					r=r-e;
				}
			return r;
			}
			if(flag == 2){
				double e= aux1(r); /*Para la primera posición*/
				while(abs(e) >= EPSILON){
					e=aux1(r);
					r=r-e;
				}
				return r;
					}
			if(flag == 3){
				return 0;
			}
		}else{	/*Para las posiciones siguientes, es diferente debido a los límites de integración*/
			double e= aux2(r, fi);
			while(abs(e)>=EPSILON){
				e=aux2(r, fi); /*funct(p[i-1].r)--fi1*/
				r=r-e;
			}	
			return r;
		}
}

