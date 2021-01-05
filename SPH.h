#include<bits/stdc++.h>
#include<vector>
#define EPSILON 0.0001
#define a 5.0
#define PI 3.1415926535897932 
#define NPV 100.0  /*Al tomar NPV y N con el mismo valor, se indica que no hay partículas estáticas virtuales cerca de las fronteras*/
#define N 100.0
#define tmax 1
#define dt 0.01
#define flag 2 /*flag1/flag2/flag3*/
#define h 0.4 /*0.5/0.4/0.1*/
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

//---------------------------------------------Distribución inicial----------------------------------------------------
double dis0(int i, double r, double rim){
	double b1;
	if(N==NPV){
			if(i==0){
				if(flag==1){
					b1=r/2.0;  /*Para la primera posición*/
				}
				if(flag==2){
					b1=r;
				}
				if(flag == 3){
					return 0;		
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
//-------------------------Busqueda de vecinos----------------------------------------------------------------------------
int cell(double ri){
	if(flag == 1){
		return (int)(ri/h);
	}
	if(flag == 2){
		return (int)((ri+4.0)/h);
	}
	if(flag == 3){
		return 0;
	}
}
int key(int n, int celli){	
		return n+1+(celli)*(int)N;
}
//--------------------Método sort radix----------------------------------------------------------------------------------
int getmax(int arr[]){
	int max=arr[0];
	for(int i=0; i< (int)N; ++i){
		if(arr[i]>max){
			max=arr[i];
		}
	}
	return max;
}
void countSort(int arr[], int pot){
	/*primero se definen los arreglos*/

	int output[(int)N], count[(int)((a/h)*(N))]={0};
	/*Después se define el número de repeticiones que hay de cada número*/
	for(int i=0; i<N;++i){
		count[(arr[i]/pot)%10]++;
	}
	for(int i=0; i<(int)((a/h)*(N));++i){
		count[i]+=count[i-1];
	}
	for(int i=N-1; i>=0; --i){
		output[count[(arr[i]/pot)%10]-1]=arr[i];
		count[(arr[i]/pot)%10]--;
	}
	for(int i=0; i<N; ++i){
		arr[i]=output[i];
	}
}
void radixsort(int arr[]){
	int m=getmax(arr);
	for(int pot=1; (m/pot)>0; pot*=10 ){
		countSort(arr, pot);
	}
}


//-------------------------Funcion de suavizado y derivadas o Kernel a h constante tipo Gaussiana--------------------------
double ker(double ri, double rj){
	double alpha=(1.0/(h*sqrt(PI)));
	double R=(ri-rj)/h;
	return alpha*exp(-R*R);
}
double dker(double ri, double rj){
	double alpha=(1.0/(h*sqrt(PI))), beta=(2.0/h);
	double R=(ri-rj)/h;
	return alpha*exp(-R*R)*(-R*beta);
}
double ddker(double ri, double rj){
	double alpha=(1.0/(h*sqrt(PI))), beta=(2.0/h), delta=2.0/(h*h);
	double R=(ri-rj)/h;
	return alpha*exp(-R*R)*(R*R*beta*beta-delta);
}
//------------------------------------------Densidad y Derivadas-----------------------------------------------------------
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
//------------------------------------------Presión-----------------------------------------------------------------------

double Press(int j, vector<double> W, vector<double> Wk, vector<double> Wxk, vector<double> Wxxk){
	double m=1.0/NPV;
		double Dx, Dxx, Q;
		Dx=Denx(Wxk);
		Dxx=Denx(Wxxk);
		double Rhok=Den(Wk);
		Q=Dx*Dx/Rhok-Dxx;
		
		return	(m/(Rhok*4.0))*W.at(j)*(Q);
}

//------------------------------------------Aceleración---------------------------------------------------------------------

double Ace(int j, double Pi, double Pj, double rj, vector<double> Wi, vector<double> Wj, vector<double> Wx){
	double m=1.0/NPV;
	double Di, Dj;
	Di=Den(Wi);
	Dj=Den(Wj);
	double Vj=(1.0/2.0)*(rj*rj);
	
	return m*(Pi/(Di*Di)+Pj/(Dj*Dj)+Vj/Dj)*Wx.at(j);
}

double error(double numerico,double analitic){
	return analitic-numerico;
}
