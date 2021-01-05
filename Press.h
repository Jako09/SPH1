#include<bits/stdc++.h>
#include<vector>
#define NPV 100.0

using namespace std;

double Press(int j, vector<double> W, vector<double> Wk, vector<double> Wxk, vector<double> Wxxk){
	double m=1.0/NPV;
		double Dx, Dxx, Q;
		Dx=Denx(Wxk);
		Dxx=Denx(Wxxk);
		double Rhok=Den(Wk);
		Q=Dx*Dx/Rhok-Dxx;
		
		return	(m/(Rhok*4.0))*W.at(j)*(Q);
}
