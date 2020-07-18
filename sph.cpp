#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>

using namespace std;

int N=100, tmax=10; /*N es la cantidad de partículas, tmax, es la cantidad de dt que se realizará la simulación*/
double m=1/(double)N; /*m es la masa definida en el artículo m=1/N*/
double PI= 3.1415926535897932;
double h=2;/*Parametro de suavizado*/ 
double g=0;
double alpha=(1/(h*sqrt(PI))), beta=(2/h), delta=2/(h*h);

struct part{
	double r, v, a, Rho, P, V; /*r es la posicion, v es la velocidad, a es la aceleración, Rho es la densidad, P la presión y V el potencial*/
	vector<int> in;
	vector<double> W, dxW, dxxW; /*in son los indices de las partículas dentro del rango de interacción de la función de suavizado, W es la función de suavizado, dxW es el valor de la primera derivada de W, dxxW es el valor de la segunda derivada de W*/
};


/*Es de notar que primero se dan las posiciones y velocidades, además de asegurarnos que no se superpongan las partículas*/
/*Después se determinan los valores de*/
int main(){
	part sph[N];
	double dt=0.004;
	ofstream file9("vtmax.dat");
	ofstream file10("rtmax.dat");
	for(int t=0; t<tmax;++t){
		/*Primer paso, veloccidades y posciones*/
		if(t==0){
		ofstream file1 ("r0.dat");
		ofstream file2 ("v0.dat");
			for(int i=0; i<N;++i){
				sph[i].v=0.0;    /*VELOCIDADES*/
				file2 << i << " " << sph[i].v << '\n';
				sph[i].r=-5+0.1*i; /*POSICIONES*/			
				file1 << i << " " << sph[i].r << '\n';
			}
		file1.close();
		file2.close();
		}
			
		
		/*Primero se calcula la distancia relativa*/
		/*Se determinan que partículas se encuentran dentro de la región de interacción*/
		/*Se calcula el valor de la función de suavizado y de sus derivadas*/
//		ofstream file3 ("R.dat");
//		ofstream file4 ("W.dat");
		for(int i=0; i<N;i++){
			sph[i].in.clear();
			sph[i].W.clear();
			sph[i].dxW.clear();
			sph[i].dxxW.clear();
			for(int j=0; j<N; j++){
				if(i!=j){
					double R, W, Wx, Wxx, eta=1e-7;
					R=(sph[i].r-sph[j].r)/h;
//					file3 << i << " " << j << " " << R << '\n';
					if(R<1){
						sph[i].in.push_back(j);
						if(R!=0){/*Aquí se toman los valores donde no hay superposición de partículas*/
							W=alpha*exp(-R*R);
							Wx=alpha*exp(-R*R)*(-R*beta);
							Wxx=alpha*exp(-R*R)*(R*R*beta*beta-delta);
							sph[i].W.push_back(W);
							sph[i].dxW.push_back(Wx);
							sph[i].dxxW.push_back(Wxx);
//							file4 << i << " " << j << " " << W << " " << Wx << " " << Wxx << '\n';
						}else{/*Aquí se toman los valores donde hay superposición*/
							W=alpha*exp(eta*eta);
							Wx=alpha*exp(-eta*eta)*(-eta*beta);
							Wxx=alpha*exp(-eta*eta)*(eta*eta*beta*beta-delta);
							sph[i].W.push_back(W);
							sph[i].dxW.push_back(Wx);
							sph[i].dxxW.push_back(Wxx);
//							file4 << i << " " << j << " " << W << " " << Wx << " " << Wxx << '\n';
						}
					}	
				}
			}
		}
//		file3.close();
//		file4.close();
		/*Se calculan los valores de la densidad*/
//		ofstream file5("rho.dat");
		for(int i=0; i<N; ++i){
			sph[i].Rho=0;
			for(int j=0; j<sph[i].in.size(); ++j){
				sph[i].Rho+=m*sph[i].W.at(j);
			}
//			file5 << i << " " << sph[i].Rho << '\n'; 
		}
//		file5.close();
		/*Se calculan los valores de la presión*/
//		ofstream file6("p.dat");
		for(int i=0; i<N; ++i){
			sph[i].P=0;
			for(int j=0; j<sph[i].in.size(); ++j){
				double dxrho=0, dxxrho=0;
				int l=sph[i].in.at(j);/*Indice de la partícula en el arreglo in de la partícula i*/
				double rhoj=sph[l].Rho;
				/*Se calculan las derivadas de la densidad*/
				for(int k=0; k<sph[l].in.size();++k){
					dxrho+=m*sph[l].dxW.at(k);
					dxxrho+=m*sph[l].dxxW.at(k);
				}
				sph[i].P+=(m/rhoj)*(0.25)*(dxrho*dxrho/rhoj-dxxrho)*sph[i].W.at(j);
			}
//			file6 << i << " " << sph[i].P << '\n'; 
		}
//		file6.close();
		/*Se calcula el potencial*/
//		ofstream file7("V.dat");
		for(int i=0; i<N; ++i){
			sph[i].V=0;
			for(int j=0; j<sph[i].in.size(); ++j){
				int l=sph[i].in.at(j);
				double rhoj=sph[l].Rho;
				sph[i].V+=(m/rhoj)*(0.5)*sph[l].r*sph[l].r*sph[i].W.at(j);
			}
//			file7 << i << " " << sph[i].V << '\n';
		}
//		file7.close();
		/*Se calculan las aceleraciones*/
//		ofstream file8("a.dat");
		for(int i=0; i<N; ++i){
			sph[i].a=0;
			double rhoi=sph[i].Rho;
			for(int j=0; j<sph[i].in.size(); ++j){
				int l=sph[i].in.at(j);
				double rhoj=sph[l].Rho;
				double Vj=sph[l].V;
				sph[i].a+=-m*sph[i].dxW.at(j)*(g+(1/rhoj)*Vj+(sph[i].P/(rhoi*rhoi)+sph[l].P/(rhoj*rhoj)));
			}
//			file8 << i << " " << sph[i].a << '\n';
		}
//		file8.close();


	for(int i=0; i<N; ++i){
		sph[i].v+=sph[i].a*dt;
		cout << sph[i].v << '\n';
		file9 << t << " " << i << " " << sph[i].v << '\n';
		sph[i].r+=sph[i].v*dt;
		cout << sph[i].r << '\n';
		file10 << t << " " << i << " " << sph[i].r << '\n';
	}

	}
	file9.close();
	file10.close();
}





