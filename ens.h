#include<bits/stdc++.h>
#define flag 2 /*flag1/flag2/flag3*/
#define h 0.4 /*0.5/0.4/0.1*/
#define N 100.0

using namespace std;
//-------------------------Busqueda de vecinos----------------------------------------------------------------------------
int cell(double ri){
	if(flag == 1){
		return (int)(ri/h);
	}
	if(flag == 2){
		return (int)((ri+4.0)/h);
	}
	if(flag == 3){
		return 0;
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
