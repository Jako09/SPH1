#include<bits/stdc++.h>
#include<vector>
#include "SPH.h"
#define EPSILON 0.0001
#define a 5.0
#define PI 3.1415926535897932 
#define NPV 100.0  /*Al tomar NPV y N con el mismo valor, se indica que no hay partículas estáticas virtuales cerca de las fronteras*/
#define N 100.0
#define tmax 1
#define dt 0.01
#define flag 2
using namespace std;

struct particle{
	double r, v, Rho, P, Phi, A; /*r es igual a la posicion de la particula virtual*/ /*Rho es la densidad de cada partícula virtual*/
	vector<int> in;
	vector<double> W, Wx, Wxx; /*Arreglo para gusradar el valor de la función de suavizado y sus derivadas*/
};

int main(){
particle p[(int)N];
bool search = false;
for(int t=0; t<tmax; ++t){
	if(t==0){
//	ofstream file3("f2t0i100.dat");
	int z1=0;
	for(int i=0; i<(int)NPV; ++i){
//		cout << i << '\n';	
		if(i==0){
			if(flag == 1){
				p[i].r=NR3(i, 2.5, 0);
			}
			if(flag == 2){
				p[i].r=NR3(i,0.1, 0);	
//				file3 << i << " " << p[i].r << '\n';
			}
			if(flag == 3){
				p[i].r=1.0/(2.0*N);
			}
		}else{
			if(flag == 1){
				p[i].r=NR3(i, 2.5, NR1(p[i-1].r));
			}
			if(flag == 2){
				if(i<(int)(NPV/2)){
					p[i].r=NR3(i, 0.1, NR1(p[i-1].r));
				}else{
					p[(int)(N/2.0)+z1].r=-p[(int)(N/2.0-1)-z1].r;
					z1++;
				}
			}
			if(flag == 3){
				p[i].r=1/(2.0*N)+((double)i)/(N);
			}
		}
	}
	if(flag != 3){
//	file3.close();
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
			if(i<(int)(NPV)/2){
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
	}
	
	if(flag == 2){
		for(int i=0; i<(int)(N/2.0); ++i){
			p[(int)(N/2.0)+i].r=-p[(int)(N/2.0-1)-i].r;
		}	
	}
	
	for(int i=0; i<(int)N; ++i){
		p[i].v=0;
	}
	}
	}else{
		for(int i=0; i<(int)N; ++i){
			p[i].v+=p[i].A*dt;
			p[i].r+=p[i].v*dt;
		}
	}
	/*--Se introduce la función cell y key--*/
	if(search == true){
		int keyo[(int)N];
		ofstream filek("f2keyw.dat");
		for(int i=0; i<(int)N; ++i){
			keyo[i]=key(i, cell(p[i].r));
			filek << i << " " << keyo[i] << '\n';
		}
		filek.close();
	/*------Se ordenan los valores del Key, mediante sort radix------*/
		radixsort(keyo);
		ofstream fileko("f2keyo.dat");
		for(int i=0; i< (int)N; ++i){
			fileko << i << " " << keyo[i] << '\n';
		}
		fileko.close();
	/*------------------------------Se determina que intervalos están ocupados--------------------------*/
		int lap;
		if(flag==1){
			lap=10*1;
		}
		if(flag ==2){
			lap=10*2;
		}
			int on[lap]={0};	
	/*-------------------------------cambiar el for por el while, los indices*/

/*	for(int i=0; i<N;++i){
		int k1=keyo[i];
		bool prueba=false;
		int j=0;	
		while(prueba == false){
			if(k1 >= j*(int)N && k1 < (j+1)*(int)N){
				on[j]=1;
				prueba=true;
			}
			j++;
		}	
	}
*/
		for(int i=0; i<lap; ++i){
			bool prueba = true;
			int j=0;
			while(j< (int)N){
				if(keyo[j] > i*(int)N && keyo[j] <= (i+1)*(int)N){
					on[i]++;
				}
				j++;
			}
		}
	/*-----------Se muestran que celdas están ocupadas-----------------*/
		ofstream fileon("f2on.dat");
		for(int i=0; i<lap; ++i){
			fileon << i << " " << on[i] << '\n';
		}
		fileon.close();
	}
	
	/*------------------------------------------------------------------------------------------------------*/
	/*-------Implementación del buscador eficiente de vecinos------*/
	
	if(search == true){
		for(int i=0; i<(int)N; ++i){
			p[i].in.clear();
			p[i].W.clear();
			p[i].Wx.clear();
			p[i].Wxx.clear();
			int j=0;
			
		}	
	}else{
		for(int i=0; i<(int)N;++i){
			p[i].in.clear();
			p[i].W.clear();
			p[i].Wx.clear();
			p[i].Wxx.clear();
			for(int j=0; j<(int)N;++j){
				if(i!=j){
					double W=ker(p[i].r, p[j].r), Wx=dker(p[i].r, p[j].r), Wxx=ddker(p[i].r, p[j].r);
					p[i].in.push_back(j);
					p[i].W.push_back(W);
					p[i].Wx.push_back(Wx);
					p[i].Wxx.push_back(Wxx);
				}
			}
		}
	}
		
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
		cout << p[i].Rho << '\n';	
	}

	
	for(int i=0; i<(int)N; ++i){
		p[i].P=0;
		for(int j=0; j<p[i].in.size(); ++j){			
			p[i].P+=Press(j, p[i].W, p[p[i].in.at(j)].W, p[p[i].in.at(j)].Wx, p[p[i].in.at(j)].Wxx);
		}
	}
	for(int i=0; i<(int)N; ++i){
		double m=1.0/NPV;
		p[i].Phi=0;
		for(int j=0; j<p[i].in.size();++j){
			double Phi0=(1.0/2.0)*(p[p[i].in.at(j)].r)*(p[p[i].in.at(j)].r);
			p[i].Phi+=(m/Den(p[p[i].in.at(j)].W))*Phi0*p[i].Wx.at(j);
		}
	}
	for(int i=0; i<(int)N; ++i){
		p[i].A=0;
		if( i<10 || i>=90 ){
			p[i].A=0;
		}else{
		for(int j=0; j<p[i].in.size(); ++j){
			p[i].A+=-Ace(j, p[i].P, p[p[i].in.at(j)].P, p[p[i].in.at(j)].r,  p[i].W, p[p[i].in.at(j)].W, p[i].Wx);
		}
		}
	}
	/*Obtención del histograma de frecuencia
	double error[N];
	for(int i=0; i<(int)N; ++i){
		error[i]=error(Den(p[i].W),(1.0/sqrt(PI))*exp(-p[i].r*p[i].r));
	}
	double max=error[0], min=error[0];
	for(int i=0; i<(int)N; ++i){
		if(error[i]>max){
			max=error[i];
		}
		if(error[i]<min){
			min=error[i];
		}
	}
	int range= (int)(max-min);
	int classes= (int)sqrt(N);
	int amp=(int) range/clasees;
	int histogram[classes];
	for(int i=0; i<(int)N; ++i){
		bool proof2=false;
		int j=0;
		while(proof2==false){
			double left=min+j*amp, right=min+(j+1)*amp;
			if(error[i]>left || error[i]<right){
				proof2=true;
				histogram[j]++;
			}else{
				++j;
			}
		}
	}
	*/
	
	
	
	
	if(t==0){
		ofstream file0("f3t0DRN100h01.dat");
		ofstream file1("f3t0PRN100h01.dat");
		ofstream file2("f3t0ARN100h01.dat");
		ofstream file3("f3t0PhiRN100h01.dat");
		ofstream file7("f3RhovsPosicionN100h01.dat");
		ofstream file8("f3RhoxvsPosicionN100h01.dat");
		ofstream file9("f3RhoxxvsPosicionN100h01.dat");
		ofstream file10("f3RhoErrorN100h01.dat");
		ofstream file11("f3histogramN100h01.dat");
		for(int i=0; i<(int)N; ++i){
			file0 << p[i].r << " " << p[i].Rho << '\n';
			file1 << p[i].r << " " << p[i].P << '\n';
			file2 << p[i].r << " " << abs(p[i].A) << '\n';
			file3 << p[i].r << " " << (1.0/2.0)*(p[i].r)*(p[i].r) << '\n';
			file7 << p[i].r << " " << p[i].Rho << '\n';
			file8 << p[i].r << " " << Denx(p[i].Wx) << '\n';
			file9 << p[i].r << " " << Denxx(p[i].Wxx) << '\n';
			if(flag == 1){
				file10 << p[i].r << " " << abs(error(Den(p[i].W),(2.0/a)*sin(PI*p[i].r/a)*sin(PI*p[i].r/a)));
				file11 << abs(error(Den(p[i].W),(2.0/a)*sin(PI*p[i].r/a)*sin(PI*p[i].r/a))) << '\n';
			}
			if(flag == 2){
				file10 << p[i].r << " " << abs(error(Den(p[i].W),(1.0/sqrt(PI))*exp(-p[i].r*p[i].r))) << '\n';
				file11 << abs(error(Den(p[i].W),(1.0/sqrt(PI))*exp(-p[i].r*p[i].r))) << '\n';
			}
			if(flag == 3){
				file10 << p[i].r << " " << abs(error(Den(p[i].W),1.0)) << '\n';
				if(p[i].r < 0.8 || p[i].r > 0.2 ){
				file11 << error(Den(p[i].W),1.0) << '\n';
				}
			}

		}
		file11.close();
		file10.close();
		file9.close();
		file8.close();
		file7.close();
		file3.close();
		file2.close();
		file1.close();
		file0.close();
		
//		ofstream file11("f2HistogramerrorRho.dat");
//		for(int i=0; i<classes; i++){
//			double left=min+j*amp, right=min+(j+1)*amp;
//			file11 <<  << 
//		}
//		file11.close();
	}
/*	if(t==tmax-1){
		ofstream file6("f2t100DRN100.dat");
		ofstream file4("f2t100PRN100.dat");
		ofstream file5("f2t100ARN100.dat");
		for(int i=0; i<(int)N; ++i){
			file6 << p[i].r << " " << p[i].Rho << '\n';
			file4 << p[i].r << " " << p[i].P << '\n';
			file5 << p[i].r << " " << abs(p[i].A) << '\n';
		}
		file5.close();
		file4.close();
		file6.close();
		
	}
*/

}

}
