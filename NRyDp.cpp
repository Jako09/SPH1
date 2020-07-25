#include<bits/stdc++.h>
#include<vector>
#define EPSILON 0.001
#define a 5.0
#define PI 3.1415926535897932 


double N=1600.0;/*El # de particulas virtuales*/
using namespace std;

/*Se desea resolver la ecuación */

struct particle{
	double r, Rho; /*r es igual a la posicion de la particula virtual*/ /*Rho es la densidad de cada partícula virtual*/
	vector<int> in; /*Arreglo para determinar todas las partículas virtuales que estarán dentro de la longitud de suavizado*/
	vector<double> W, Wx, Wxx; /*Arreglo para gusradar el valor de la función de suavizado y sus derivadas*/
};

double funct(double x){ /*Definimos la ecuación que se va a resolver en el método Newton-Raphson*/
//	return x*x*x-x*x+2;
	return (1.0/a)*(x)-(1.0/(2.0*PI))*sin(2.0*PI*x/a)-1.0/N;
}

double dfunct(double x){
//	return 3*x*x-2*x; /*La derivada de la ecuación a encontrar para el método Newton-Raphson*/
	return (1.0/a)*(1.0-cos(2.0*PI*x/a));
}


int main(){
	particle p[(int)N]; /*Definimos el número de partículas*/
	for(int i=0; i<(int)N; ++i){
		double r=2.5; /*corresponde a un valor inicial, para hacer el calculo de Newton-Raphson*/
		if(i==0){
			double e= funct(r)/dfunct(r); /*Para la primera posición*/
			while(abs(e) >= EPSILON){
				e=funct(r)/dfunct(r);
				r=r-e;
			}
			p[i].r=r;
		}else{					/*Para las posiciones siguientes, es diferente debido a los límites de integración*/
			double e= (funct(r)-funct(p[i-1].r)-1.0/N)/(dfunct(r));
			while(abs(e)>=EPSILON){
				e= (funct(r)-funct(p[i-1].r)-1.0/N)/(dfunct(r));
				r=r-e;
			}	
			p[i].r=r;
		}
	}
	
	ofstream file("r.dat");
	double b[(int)N];  /*Un arreglo solo para hacer cáculos*/
	for(int i=0; i<(int)N;++i){
		b[i]=p[i].r;
	}
	for(int i=0; i<(int)N; ++i){
		if(i==0){
			p[i].r=p[i].r/2.0;  /*Para la primera posición*/
			file << i << " " << p[i].r << '\n';
		}else{				/*Para las demás posiciones*/
			p[i].r=(p[i].r+b[i-1])/2.0; /*surge de tomar en cuenta (b[i].r-b[i-1].r)/2+b[i-1].r*/
			file << i << " " << p[i].r << '\n';
		}
	}
	file.close();
	double  h[(int)N] /*200/N*/, alpha[(int)N] /*(1.0/(h*sqrt(PI)))*/, beta[(int)N]/*(2.0/h)*/, delta[(int)N] /*2.0/(h*h)*/,  n, m=1.0/N; /*h[N] es el valor de la distancia de suavizado para cada partícula, alpha, beta y delta solo son valores que ayudan a reducir la expresión*/
	
	for(int i=0; i<(int)N;++i){
		double s=sin(PI*p[i].r*0.2);
		h[i]=1.1*m/(0.4*s*s);
		alpha[i]=(1.0/(h[i]*sqrt(PI)));
		beta[i]=(2.0/h[i]);
		delta[i]=2.0/(h[i]*h[i]);
	}
	/*Se determinan los valores de la función de suavizado o Kernel*/
	for(int i=0; i<(int)N;i++){
			p[i].in.clear();
			p[i].W.clear();
			p[i].Wx.clear();
			p[i].Wxx.clear();
			for(int j=0; j<(int)N; j++){
				if(i!=j){
					double R, W, Wx, Wxx, eta=1e-7;
					R=(p[i].r-p[j].r)/h[i];
//					file3 << i << " " << j << " " << R << '\n';
					if(R<1){
						p[i].in.push_back(j);
						if(R!=0){/*Aquí se toman los valores donde no hay superposición de partículas*/
							W=alpha[i]*exp(-R*R);
							Wx=alpha[i]*exp(-R*R)*(-R*beta[i]);
							Wxx=alpha[i]*exp(-R*R)*(R*R*beta[i]*beta[i]-delta[i]);
							p[i].W.push_back(W);
							p[i].Wx.push_back(Wx);
							p[i].Wxx.push_back(Wxx);
//							file4 << i << " " << j << " " << W << " " << Wx << " " << Wxx << '\n';
						}else{/*Aquí se toman los valores donde hay superposición*/
							W=alpha[i]*exp(eta*eta);
							Wx=alpha[i]*exp(-eta*eta)*(-eta*beta[i]);
							Wxx=alpha[i]*exp(-eta*eta)*(eta*eta*beta[i]*beta[i]-delta[i]);
							p[i].W.push_back(W);
							p[i].Wx.push_back(Wx);
							p[i].Wxx.push_back(Wxx);
//							file4 << i << " " << j << " " << W << " " << Wx << " " << Wxx << '\n';
						}
					}	
				}
			}
	}
	/*Se calcula el valor de la densidad*/
		ofstream file2("dvsrN1600.dat");
		for(int i=0; i<(int)N; ++i){
			p[i].Rho=0;
			for(int j=0; j<p[i].in.size(); ++j){
				p[i].Rho+=m*p[i].W.at(j);
			}
			file2 << p[i].r << " " << p[i].Rho << '\n'; 
		}
		file2.close();
	/*Se cálcula el error*/
		ofstream file3("errorN1600.dat");
		for(int i=0; i<(int)N;++i){
			double s=sin(PI*p[i].r*0.2);
			double den=0.4*s*s;
			file3 << p[i].r << " " << abs(den-p[i].Rho) << '\n'; 
		}
		file3.close();
	return 0;
}













