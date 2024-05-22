#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

int main()
{   //1.Set variables
	  int Col=  2; 	//Column to be processed 
	  int NCol= 2;  //Number of columns 
	  
	  //Histogram parameters 
	  double BucketSize = 0.01;
	  double Min = -1;
	  double Max =  1;
	  
	  //Number of buckets
	  int NBucket = ceil( (Max-Min)/BucketSize );
	  
	  //Counting arrays
	  double Cnt[NBucket];
	  
	  //Counting variable
	  int Bucket;
	  
      //Readout variables for a list of type {Index,Value} 
	  double Index; 
	  double Value; 
	  
	  //Normalization;
	  double Norm;
	
      //Input file
	  FILE *Data;
	  //Output files
	  FILE *Hist;
	  
	//2.Open files 
	  Data = fopen("Spin.dat","r");
	  Hist = fopen("Hist.dat","w+");
	  
	//3.Initialize counting arrays
	  for( int j=0; j<NBucket; j++){
	    Cnt[j] = 0;
	    }
	
	//4.Read input and create counting array
	  while (!feof(Data)) {
		
		for(int p=0; p<NCol; p++){
	    fscanf(Data,"%lf", &Index);
	    if(p==Col-1){Value=Index;}
	    }
	    
	    Bucket = floor((Value-Min)/BucketSize);
	    if(Bucket>-1 && Bucket<NBucket){
	    Cnt[Bucket] =  Cnt[Bucket] + 1;
	    }
      }
     
	fclose(Data);
	
	//4.Normalize counting arrays, fill LDF arrays, write output
	
	  Norm = 0;
	  for( int j=0; j<NBucket; j++){Norm = Norm + Cnt[j]*BucketSize; }
	  for( int j=0; j<NBucket; j++){Cnt[j] = Cnt[j]/Norm;}
	  
	  for( int j=0; j<NBucket; j++){
	  fprintf(Hist,"%lf %lf \n",  Min + BucketSize*(j+0.5), Cnt[j]);
      }
      
	fclose(Hist);
}
