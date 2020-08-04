#include<bits/stdc++.h>
#include<vector>
#define EPSILON 0.001
#define a 5.0
#define PI 3.1415926535897932 
#define NPV 100.0  /*Al tomar NPV y N con el mismo valor, se indica que no hay partículas estáticas virtuales cerca de las fronteras*/
#define N 300.0

using namespace std;

struct particle{
	double r, Rho; /*r es igual a la posicion de la particula virtual*/ /*Rho es la densidad de cada partícula virtual*/
	vector<double> W, Wx, Wxx; /*Arreglo para gusradar el valor de la función de suavizado y sus derivadas*/
};

double NR1(double x){ /*Definimos la ecuación que se va a resolver en el método Newton-Raphson*/
//	return x*x*x-x*x+2;
	return (1.0/a)*(x)-(1.0/(2.0*PI))*sin(2.0*PI*x/a)-1.0/NPV;  /*N-NPV*/
}

double NR2(double x){
//	return 3*x*x-2*x; /*La derivada de la ecuación a encontrar para el método Newton-Raphson*/
	return (1.0/a)*(1.0-cos(2.0*PI*x/a));
}
//-------------------------Newton-Rapson Method--------------------------
double NR3(int i, double r, double fi1){
		if(i==0){
			double e= NR1(r)/NR2(r); /*Para la primera posición*/
			while(abs(e) >= EPSILON){
				e=NR1(r)/NR2(r);
				r=r-e;
			}
			return r;
		}else{				/*Para las posiciones siguientes, es diferente debido a los límites de integración*/
			double e= (NR1(r)-fi1-1.0/NPV)/(NR2(r));
			while(abs(e)>=EPSILON){
				e= (NR1(r)-fi1-1.0/NPV)/(NR2(r)); /*funct(p[i-1].r)--fi1*/
				r=r-e;
			}	
			return r;
		}
}
//--------------------Distribución inicial, Particula cuántica en una caja------------------
double dis0(int i, double r, double rim){
	double b1;
	if(N==NPV){
			if(i==0){
				b1=r/2.0;  /*Para la primera posición*/
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
//------------------------------Funcion de suavizado o Kernel a h constante tipo Gaussiana-------------------------------------
double ker(double ri, double rj){
	double h=0.5, alpha=(1.0/(h*sqrt(PI)));
	double R=(ri-rj)/h;
	return alpha*exp(-R*R);
}
double dker(double ri, double rj){
	double h=0.5, alpha=(1.0/(h*sqrt(PI))), beta=(2.0/h);
	double R=(ri-rj)/h;
	return alpha*exp(-R*R)*(-R*beta);
}
double ddker(double ri, double rj){
	double h=0.5, alpha=(1.0/(h*sqrt(PI))), beta=(2.0/h), delta=2.0/(h*h);
	double R=(ri-rj)/h;
	return alpha*exp(-R*R)*(R*R*beta*beta-delta);
}
//------------------------------------------Densidad------------------------------------------------

double Den(vector<double> W){
	double m=1/NPV, b2;
	for(int j=0; j<W.size();++j){
		b2+=m*W.at(j);
	}
	return b2;
}

//------------------------------------------SPH Method ----------------------------------------------
void run(){
	particle p[(int)N];
//	ofstream file0("r0iN100.dat");
	for(int i=0; i<(int)NPV; ++i){
		if(i==0){
			p[i].r=NR3(i, 2.5, 0);
		}else{
			p[i].r=NR3(i, 2.5, NR1(p[i-1].r));
		}
//		file0 << i << " " << p[i].r << '\n';
	}
//	file0.close();
//	ofstream file("riN100.dat");
	double d[(int)N];
	for(int i=0; i<(int)N; ++i){
		d[i]=p[i].r;
		if(N==NPV){
			if(i==0){
				p[i].r=dis0(i, d[i], 0);
			}else{
				p[i].r=dis0(i, d[i], d[i-1]);
			}
			
		}else{
			if(i<(int)NPV){
				if(i==0){
					p[i].r=dis0(i, d[i], 0);
				}else{
					p[i].r=dis0(i, d[i], d[i-1]);
				}
			}else{
				if(i >= (int)NPV  && i<(int)(2.0*NPV)){
					p[i].r=dis0(i, 0, p[i-(int)NPV].r);
				}
				if(i >= (int)(2.0*NPV) && i<(int)(3.0*NPV)){
					p[i].r=dis0(i, 0, p[i-(int)(2.0*NPV)].r);
				}
			}
		}
//		file << i << " " << p[i].r << '\n';
	}
//	file.close();
	for(int i=0; i<(int)N;++i){
		for(int j=0; j<(int)N;++j){
			if(i!=j){
				double W=ker(p[i].r, p[j].r), Wx=dker(p[i].r, p[j].r), Wxx=ddker(p[i].r, p[j].r);
				p[i].W.push_back(W);
				p[i].Wx.push_back(Wx);
				p[i].Wxx.push_back(Wxx);
			}
		}
	}
	ofstream file("DPN100.dat");
	for(int i=0; i<(int)N; ++i){
		if(N==NPV){
			p[i].Rho=Den(p[i].W);
			
		}else{
			if(i < (int)NPV){			/*___N-NPV*/
					p[i].Rho=-Den(p[i].W);
			}
			if(i >= (int)(NPV) && i<(int)(2.0*NPV)){
				p[i].Rho=Den(p[i].W);
			}
			if(i >= (int)(2.0*NPV) && i<(int)(3.0*NPV)){
				p[i].Rho=-Den(p[i].W);
			}
		}
		file << p[i].r << " " << p[i].Rho << '\n';
	}
	file.close();
}




int main(){
	run();
	return 0;
}


