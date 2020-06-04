#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>

using namespace std;

int N=100;
double h=1; /*El parametro de suavizado es de 5, constante*/
double m=1/(double)N;
double PI = 3.1415926535897932;

struct particula{
	double r;
	double v;
	vector<int> in;
	vector<double> R;
	vector<double> W;
	double Rho;
	double Pres;
	double ac;
};

double R_relativa(double x1, double x2);
double densidad(vector<double> W);

int main(){
	particula sph[N];
	/*Se dan las posiciones iniciales--------------------------*/
	ofstream file1 ("posiciones.dat");
	for(int i=0; i<N; ++i){
		sph[i].r=((double)(rand()%10001)-5000)/1000;
		file1 << i << " " << sph[i].r << '\n';
	}
	file1.close();
	/*Se calcula la R_relativa y W-suavizado---------------------------------*/
	for(int i=0; i<N;++i){
		for(int j=0; j<N;++j){
			if(i!=j){
				double R;
				double Wij;
				R=R_relativa(sph[i].r,sph[j].r);
				if(R<=h){
					sph[i].in.push_back(j);
					sph[i].R.push_back(R);
					Wij=((1/(h*sqrt(PI)))*(1/(h*sqrt(PI)))*(1/(h*sqrt(PI))))*exp(-(R*R)/(h*h));
					sph[i].W.push_back(Wij);
				}
			}
		}
	}
	/*calculo de la densidad para cada partícula*/
	for(int i=0; i<N;++i){
		sph[i].Rho=densidad(sph[i].W);
	}
	/*Calculo de la presión, solo se presenta en una dimensión*/
	double dxWij=0, dxxWij=0;
	for(int i=0; i<N; ++i){
		for(int j=0; j<sph[i].in.size();++j){
			double rij=(sph[i].r-sph[sph[i].in.at(j)].r);
			dxWij+=((1/(h*sqrt(PI)))*(1/(h*sqrt(PI)))*(1/(h*sqrt(PI))))*exp(-(rij*rij)/(h*h))*(-2*rij/(h*h));
			dxxWij+=((1/(h*sqrt(PI)))*(1/(h*sqrt(PI)))*(1/(h*sqrt(PI))))*exp(-(rij*rij)/(h*h))*(4*(rij/(h*h))*(rij/(h*h))-2/(h*h));
		}
	}
	
	/*Impresion de resultados*/
	ofstream file2 ("R_relativa.dat");
	ofstream file3 ("F_suavizado.dat");
	ofstream file4 ("Rho_densidad.dat");
	for(int i=0; i<N;++i){
		for(int j=0;j<sph[i].in.size() ;++j){
			if(i!=j){
				file2 << i << " " << sph[i].in.at(j) << " " << sph[i].R.at(j) << '\n'; 
				file3 << i << " " << sph[i].in.at(j) << " " << sph[i].W.at(j) << '\n'; 
			}
		}
		file4 << i << " " << sph[i].Rho << '\n';
	}
	file2.close();
	file3.close();
	file4.close();
	

}


double R_relativa(double x1, double x2){
	return sqrt((x1-x2)*(x1-x2));
}
double densidad(vector<double> W){
	double rho;
	for(int i=0; i<W.size();++i){
		rho+=(m)*(W.at(i));
	}
	return rho;
}



