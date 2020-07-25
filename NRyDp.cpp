#include<bits/stdc++.h>
#include<vector>
#define EPSILON 0.001
#define a 5.0
#define PI 3.1415926535897932 /*Sería para 100 partículas virtuales*/


double N=1600.0;
using namespace std;

/*Se desea resolver la ecuación */

struct particle{
	double r, Rho;
	vector<int> in;
	vector<double> W, Wx, Wxx;
};

double funct(double x){
//	return x*x*x-x*x+2;
	return (1.0/a)*(x)-(1.0/(2.0*PI))*sin(2.0*PI*x/a)-1.0/N;
}

double dfunct(double x){
//	return 3*x*x-2*x;
	return (1.0/a)*(1.0-cos(2.0*PI*x/a));
}


int main(){
	particle p[(int)N];
	for(int i=0; i<(int)N; ++i){
		double r=2.5;
		if(i==0){
			double e= funct(r)/dfunct(r);
			while(abs(e) >= EPSILON){
				e=funct(r)/dfunct(r);
				r=r-e;
			}
			p[i].r=r;
		}else{
			double e= (funct(r)-funct(p[i-1].r)-1.0/N)/(dfunct(r));
			while(abs(e)>=EPSILON){
				e= (funct(r)-funct(p[i-1].r)-1.0/N)/(dfunct(r));
				r=r-e;
			}	
			p[i].r=r;
		}
	}
	
	ofstream file("r.dat");
	double b[(int)N];
	for(int i=0; i<(int)N;++i){
		b[i]=p[i].r;
	}
	for(int i=0; i<(int)N; ++i){
		if(i==0){
			p[i].r=p[i].r/2.0;
			file << i << " " << p[i].r << '\n';
		}else{
			p[i].r=(p[i].r+b[i-1])/2.0;
			file << i << " " << p[i].r << '\n';
		}
	}
	file.close();
	double  h[(int)N] /*200/N*/, alpha[(int)N] /*(1.0/(h*sqrt(PI)))*/, beta[(int)N]/*(2.0/h)*/, delta[(int)N] /*2.0/(h*h)*/,  n, m=1.0/N;
	
	for(int i=0; i<(int)N;++i){
		double s=sin(PI*p[i].r*0.2);
		h[i]=1.5*m/(0.4*s*s);
		alpha[i]=(1.0/(h[i]*sqrt(PI)));
		beta[i]=(2.0/h[i]);
		delta[i]=2.0/(h[i]*h[i]);
	}
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
	
		ofstream file2("dvsrN1600.dat");
		for(int i=0; i<(int)N; ++i){
			p[i].Rho=0;
			for(int j=0; j<p[i].in.size(); ++j){
				p[i].Rho+=m*p[i].W.at(j);
			}
			file2 << p[i].r << " " << p[i].Rho << '\n'; 
		}
		file2.close();
		
		ofstream file3("errorN1600.dat");
		for(int i=0; i<(int)N;++i){
			double s=sin(PI*p[i].r*0.2);
			double den=0.4*s*s;
			file3 << p[i].r << " " << abs(den-p[i].Rho) << '\n'; 
		}
		file3.close();
	return 0;
}













