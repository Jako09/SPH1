#include<bits/stdc++.h>
#include<vector>
#define EPSILON 0.0001
#define a 5.0
#define PI 3.1415926535897932 
#define NPV 100.0  /*Al tomar NPV y N con el mismo valor, se indica que no hay partículas estáticas virtuales cerca de las fronteras*/
#define N 100.0
#define tmax 1
#define dt 0.001
#define flag 2
using namespace std;

struct particle{
	double r, v, Rho, P, A; /*r es igual a la posicion de la particula virtual*/ /*Rho es la densidad de cada partícula virtual*/
	vector<int> in;
	vector<double> W, Wx, Wxx; /*Arreglo para gusradar el valor de la función de suavizado y sus derivadas*/
};

double NR1(double x){ /*Definimos la ecuación que se va a resolver en el método Newton-Raphson*/
//	return x*x*x-x*x+2;
	if(flag == 1){
		return (1.0/a)*(x)-(1.0/(2.0*PI))*sin(2.0*PI*x/a)-1.0/NPV;  /*N-NPV*/
	}
	if(flag == 2){
		return (0.5)*erf(x)-1.0/NPV;
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
}
double aux1(double r){
	if(flag == 2){
		return (NR1(r)-(0.5)*erf(-4.0))/NR2(r);
	}
	if(flag == 1){
		return NR1(r)/NR2(r);
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
		double h=0.5;
		return (int)(ri/h);
	}
	if(flag == 2){
		double h=0.4;
		return (int)((ri+4.0)/h);
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
	double h=0.5;
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

double Ace(int j, double Pi, double Pj, vector<double> Wi, vector<double> Wj, vector<double> Wx){
	double m=1.0/NPV;
	double Di, Dj;
	Di=Den(Wi);
	Dj=Den(Wj);
	
	return m*(Pi/(Di*Di)+Pj/(Dj*Dj))*Wx.at(j);
}

//------------------------------------------SPH Method ----------------------------------------------------------------------
void run(){
particle p[(int)N];
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
		}else{
			if(flag == 1){
				p[i].r=NR3(i, 2.5, NR1(p[i-1].r));
			}
			if(flag == 2){
				if(i<(int)(NPV/2)){
					p[i].r=NR3(i, 0.1, NR1(p[i-1].r));
				}else{
					p[50+z1].r=-p[49-z1].r;
					z1++;
				}
			}
		}
	}
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
			p[50+i].r=-p[49-i].r;
		}	
	}
	
	for(int i=0; i<(int)N; ++i){
		p[i].v=0;
	}
	}else{
		for(int i=0; i<(int)N; ++i){
			p[i].v+=p[i].A*dt;
			p[i].r+=p[i].v*dt;
		}
	}
	/*--Se introduce la función cell y key--*/
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
	/*------------------------------------------------------------------------------------------------------*/
	/*-----------Es necesario encontrar cual es el promedio de partículas dentro de cada celda--------------*/
	/*-------------------------------------------------------------------------------------------------------*/
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
		p[i].A=0;
		if( i<10 || i>=90 ){
			p[i].A=0;
		}else{
		for(int j=0; j<p[i].in.size(); ++j){
			p[i].A+=Ace(j, p[i].P, p[p[i].in.at(j)].P, p[i].W, p[p[i].in.at(j)].W, p[i].Wx);
		}
		}
	}
	
//	if(t==0){
		ofstream file0("f2t0DRN100.dat");
		ofstream file1("f2t0PRN100.dat");
		ofstream file2("f2t0ARN100.dat");
//		ofstream file3("f2t0i100.dat");
		for(int i=0; i<(int)N; ++i){
			file0 << p[i].r << " " << p[i].Rho << '\n';
			file1 << p[i].r << " " << p[i].P << '\n';
			file2 << p[i].r << " " << abs(p[i].A) << '\n';
//			file3 << i << " " << p[i].r << '\n';
		}
//		file3.close();
		file2.close();
		file1.close();
		file0.close();

//	}
/*	if(t==tmax-1){
		ofstream file3("t100DRN100.dat");
		ofstream file4("t100PRN100.dat");
		ofstream file5("t100ARN100.dat");
		for(int i=0; i<(int)N; ++i){
			file3 << p[i].r << " " << p[i].Rho << '\n';
			file4 << p[i].r << " " << p[i].P << '\n';
			file5 << p[i].r << " " << abs(p[i].A) << '\n';
		}
		file5.close();
		file4.close();
		file3.close();
		
	}
*/

}
}




int main(){
	run();
	return 0;
}


