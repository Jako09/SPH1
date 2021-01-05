#include<bits/stdc++.h>
#include<vector>


double Ace(int j, double Pi, double Pj, double rj, vector<double> Wi, vector<double> Wj, vector<double> Wx){
	double m=1.0/NPV;
	double Di, Dj;
	Di=Den(Wi);
	Dj=Den(Wj);
	double Vj=(1.0/2.0)*(rj*rj);
	
	return m*(Pi/(Di*Di)+Pj/(Dj*Dj)+Vj/Dj)*Wx.at(j);
}
