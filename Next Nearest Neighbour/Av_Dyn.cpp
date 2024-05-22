#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <array>
#include <string> 

using namespace std;

    int NCol = 2; 			//Number of columns 
    int Runs = 10;			//Number of runs 
    
int main(){
	
    double t[2];     		//Auxiliary variable 
    double Val[NCol-1];	    //Readout array short
    double Vall[NCol-1];    //Readout array long
    int SPL;		    	//Sample counter 
    FILE *Data;				//Input file
	FILE *Output; 			//Output files
	  
    Data = fopen("e_trajectories/tau=4/Spin.dat","r");
    Output = fopen("e_trajectories/tau=4/Av.dat","w+");
    
    fscanf(Data,"%lf", &t[0]); 
	for(int r=0; r<NCol-1; r++){fscanf(Data,"%lf", &Val[r]);}
		
	t[1]=t[0]+1; 
	SPL=0; 
	
	while(t[1]>t[0]){
        fscanf(Data,"%lf", &t[1]); 
        for(int r=0; r<NCol-1; r++){fscanf(Data,"%lf", &Val[r]);}
        SPL++;
        }
    
    double AvrArray[SPL][NCol]; 
    for(int r=0; r<SPL; r++){for(int u=0; u<NCol; u++){AvrArray[r][u]=0;}}
		
	fclose(Data); 
	Data = fopen("e_trajectories/tau=4/Spin.dat","r");
	
	for(int k=0; k<Runs; k++){
	  for(int p=0; p<SPL; p++){
		for(int r=0; r<NCol; r++){fscanf(Data,"%lf", &Vall[r]);}
		for(int r=0; r<NCol; r++){AvrArray[p][r] = AvrArray[p][r] + Vall[r];}
	  }
	}
	
    std::cout << SPL;
    
    for(int r=0; r<SPL; r++){for(int u=0; u<NCol; u++){AvrArray[r][u] = AvrArray[r][u]/Runs;}}
    
    for(int r=0; r<SPL; r++){
		for(int u=0; u<NCol; u++){
			fprintf(Output,"%lf   ", AvrArray[r][u]);
		}
		fprintf(Output,"\n"); 
	}
	
	fclose(Data);
	fclose(Output);
} 
