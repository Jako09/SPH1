#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>

using namespace std;
int l=50, N=1000;
double dpi[50]; 
double dp[50];
double h=1, PI= 3.1415926535897932, alpha=(1/(h*sqrt(PI))), beta=(2/h), delta=2/(h*h), eps=0.1, a=5.0, b=0.0, n, m=1/(double)N;
struct part{
	double r, Rho;
	vector<int> in;
	vector<double> W, Wx, Wxx;
};

int main(){
	part sq[N];
	ofstream filedp("dp.dat");	
	for(int j=0; j<l; ++j){
		double S=sin(PI*j*eps*0.2);
		dp[j]=0.4*S*S;
		dpi[j]=dp[j]+b;
		b+=dp[j];
		filedp << j << " " << dp[j] << " " << dpi[j] << '\n'; 
	}
	filedp.close();
	ofstream filer("r.dat");
	for(int z=0; z<N; ++z){
		bool proof=false;
		int k=0;
		n=(double)(1+rand()%(100))*0.1;
			while(proof==false){
				if(n>dpi[k] && n<dpi[k+1]){
					if(z<N/2){
						proof=true;
						sq[z].r=(a/PI)*asin(sqrt(dp[k]*2.5));
						filer << z << " " << sq[z].r << '\n';
					}else{
						proof=true;
						sq[z].r=5-(a/PI)*asin(sqrt(dp[k]*2.5));
						filer << z << " " << sq[z].r << '\n';
					}
				}else{
					k+=1;
				}
			}
	}
	filer.close();
	for(int i=0; i<N;i++){
			sq[i].in.clear();
			sq[i].W.clear();
			sq[i].Wx.clear();
			sq[i].Wxx.clear();
			for(int j=0; j<N; j++){
				if(i!=j){
					double R, W, Wx, Wxx, eta=1e-7;
					R=(sq[i].r-sq[j].r)/h;
//					file3 << i << " " << j << " " << R << '\n';
					if(R<1){
						sq[i].in.push_back(j);
						if(R!=0){/*Aquí se toman los valores donde no hay superposición de partículas*/
							W=alpha*exp(-R*R);
							Wx=alpha*exp(-R*R)*(-R*beta);
							Wxx=alpha*exp(-R*R)*(R*R*beta*beta-delta);
							sq[i].W.push_back(W);
							sq[i].Wx.push_back(Wx);
							sq[i].Wxx.push_back(Wxx);
//							file4 << i << " " << j << " " << W << " " << Wx << " " << Wxx << '\n';
						}else{/*Aquí se toman los valores donde hay superposición*/
							W=alpha*exp(eta*eta);
							Wx=alpha*exp(-eta*eta)*(-eta*beta);
							Wxx=alpha*exp(-eta*eta)*(eta*eta*beta*beta-delta);
							sq[i].W.push_back(W);
							sq[i].Wx.push_back(Wx);
							sq[i].Wxx.push_back(Wxx);
//							file4 << i << " " << j << " " << W << " " << Wx << " " << Wxx << '\n';
						}
					}	
				}
			}
	}
	
		ofstream file5("rho.dat");
		for(int i=0; i<N; ++i){
			sq[i].Rho=0;
			for(int j=0; j<sq[i].in.size(); ++j){
				sq[i].Rho+=m*sq[i].W.at(j);
			}
			file5 << i << " " << sq[i].Rho << '\n'; 
		}
		file5.close();


}

















