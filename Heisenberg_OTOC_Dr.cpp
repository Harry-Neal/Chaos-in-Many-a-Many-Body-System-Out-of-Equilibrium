#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <random>
#include <fstream>
#include <string>

using namespace std;

//0.0 Declare global variables 

  //0.01 Global system parameters 
    int    ssize=     100;
    double lambda=    1;
    double Jvar=      0.5;
    double HField[3]= {0,0,0};
    double Beta=	  2.888; 
    double dt=        0.02;
    int	   MCSmp=	  5000;
    double MCVar=	  0.25; 
    double T=         200;
    double trel=	  100; 
    double tau=       1.25;
    int    Runs=      100;
    double epsilon=	  0.01;
	double T_init= 	  100*tau;

  //0.02 Numerical constants 
    double Pi=3.141592653589793;


//0.1 Declare external driving function
  double Hext(double t, int k, double tau){
  
    switch(k){
	  case 0  : return  cos( (2*Pi*t)/tau );
	  break;
	  
	  case 1  : return  sin( (2*Pi*t)/tau );
	  break; 
	  
	  default : return 0; 
    }
  }
  
//0.2 Declare rotation function
  double MRot(double veco[], double axis[], double angle, int l){
	//Rotated vector
	double vecr[3];
	//Axis parallel component
	double axp;
	//Axis orthogonal component
	double axo[3];
	
	axp=veco[0]*axis[0] + veco[1]*axis[1] + veco[2]*axis[2];
	
	axo[0]=axis[1]*veco[2] - axis[2]*veco[1];
	axo[1]=axis[2]*veco[0] - axis[0]*veco[2];
	axo[2]=axis[0]*veco[1] - axis[1]*veco[0];
	
	for( int k=0; k<3; k++ ){
	vecr[k]=veco[k]*cos(angle) + axis[k]*axp*(1-cos(angle)) - axo[k]*sin(angle);
	}
    
    return vecr[l];
  }

//0.3 Declare local field function on A  
  double AField(double SpinB[][3], double HFieldA[][3], double JCplA[][3], double JCplB[][3], int j, int k){
    if(j==0){
		return JCplA[j][k]*SpinB[j][k] + JCplB[ssize-1][k]*SpinB[ssize-1][k] - HFieldA[j][k];
    }else{
		return JCplA[j][k]*SpinB[j][k] + JCplB[j-1][k]*SpinB[j-1][k] - HFieldA[j][k];  
    }
  }

//0.4 Declare local field function on B
  double BField(double SpinA[][3], double HFieldB[][3],  double JCplA[][3], double JCplB[][3], int j, int k){
    if(j==ssize-1){
		return JCplA[j][k]*SpinA[j][k] + JCplB[j][k]*SpinA[0][k] - HFieldB[j][k];
    }else{
		return JCplA[j][k]*SpinA[j][k] + JCplB[j][k]*SpinA[j+1][k] - HFieldB[j][k];
    }
  }

 
//0.5 Declare system energy function
   double Esys(double SpinA[][3], double SpinB[][3], double HFieldA[][3], double HFieldB[][3], double JCplA[][3], double JCplB[][3]){
	  
    double EInc = 0; //Cumulative total energy 
    
      for( int j=0; j<ssize-1; j++){
        for( int k=0; k<3; k++){
	      EInc = EInc - SpinA[j][k]*JCplA[j][k]*SpinB[j][k] - SpinB[j][k]*JCplB[j][k]*SpinA[j+1][k] + SpinA[j][k]*HFieldA[j][k] + SpinB[j][k]*HFieldB[j][k];         
	    }
      }
      
      for( int k=0; k<3; k++){
	      EInc = EInc - SpinA[ssize-1][k]*JCplA[ssize-1][k]*SpinB[ssize-1][k] - SpinB[ssize-1][k]*JCplB[ssize-1][k]*SpinA[0][k] + SpinA[ssize-1][k]*HFieldA[ssize-1][k] + SpinB[ssize-1][k]*HFieldB[ssize-1][k];         
	    }
	       
    return EInc;
  }
  
//0.6 Declare system magnetization function
  double Msys(double SpinA[][3], double SpinB[][3], int k){
	  
	double MInc = 0; //Cumulative total magnetization
	
	  for( int j=0; j<ssize; j++){
	    MInc = MInc + SpinA[j][k];
	    MInc = MInc + SpinB[j][k];
	  }          
	return MInc;
  }

//0.7 Declare cross product function
  void crossProduct(double v_A[], double v_B[], double c_P[]) {
   c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

int main(){
	std::clock_t c_start = std::clock();  	
  //1 Declare variables
  
    //1.1 Files
    FILE *Output;
	FILE *Energy_out;
    FILE *Parameters;
    
    //1.2 Declare random number generator and distributions 
    std::random_device rd;  
    std::mt19937 e2(rd());
    
    std::normal_distribution<double> Jdist(0.0,Jvar);			//Dist of coupling constants
    std::uniform_real_distribution<double> MCRxdist(0.0,1.0);	//MC Dist for random rotation axis 
	std::normal_distribution<double> MCAgdist(0.0,MCVar);		//MC Dist for random rotation angle 
	
	//1.3 Declare internal variables
	double SpinA[ssize][3]; //Spin degrees of freedom sublattice A
	double SpinB[ssize][3]; //Spin degrees of freedom sublattice B

	double SpinA1[ssize][3]; //Spin degrees of freedom sublattice A
	double SpinB1[ssize][3]; //Spin degrees of freedom sublattice B

	double SpinA2[ssize][3]; //Spin degrees of freedom sublattice A
	double SpinB2[ssize][3]; //Spin degrees of freedom sublattice B
	
	double HFieldA[ssize][3]; //External field on sublattice A
    double HFieldB[ssize][3]; //External field on sublattice B
    
    double JCplA[ssize][3]; //Spin-spin coupling constants for sublattice A
    double JCplB[ssize][3]; //Spin-spin coupling constants for sublattice B
	
	double FieldL[3]; //Local field vector
	double StrthL;    //Local field strength
	double SpinL[3];  //Local spin vector
	
	double phi; 	//Angle for random unit vector
	double theta; 	//Angle for random unit vector 
	
	double MCSpinA[ssize][3]; //MC copy of Spin config on sublattice A
	double MCSpinB[ssize][3]; //MC copy of Spin config on sublattice B
	double MCField[3]; 		  //MC random magnetic field 
	double z_ax[3] = {0,0,1}; //z-axis vector
	double norm[3];			  //normal vector

	int max_t = 0;
	for(double t=0; t<=T; t=t+dt){
		if(fmod(t+dt/20,tau)<= dt/2){
			max_t = max_t  + 1;
		}
	}
	double corrA[ssize];	  //correlator for copy A
	double corrB[ssize];	  //correlator for copy B
	double corr[2*ssize][max_t];//correlator for both coppies at given t
	double MCEi; 			  //MC initial energy 	
	double MCEf; 			  //MC final energy 
	double MCInd;			  //MC update indicator 	
	
	double Usys;			  //System energy 
	
  //2 Create parameter file	
    Parameters = fopen("Parameters.dat","w+");
    fprintf(Parameters,"ssize:   	%i  \n", 2*ssize);
    fprintf(Parameters,"lambda:    	%lf \n", lambda);
    fprintf(Parameters,"Jvar:  		%lf \n", Jvar);
    fprintf(Parameters,"HField:   	%lf %lf %lf \n", HField[0], HField[1], HField[2]);
    fprintf(Parameters,"Beta:  		%lf \n", Beta);
    fprintf(Parameters,"dt:       	%lf \n", dt);
    fprintf(Parameters,"MCSmp:    	%i  \n", MCSmp);
    fprintf(Parameters,"MCVar:    	%lf \n", MCVar);
    fprintf(Parameters,"T:     		%lf \n", T);
	fprintf(Parameters,"Runs:     	%i \n",  Runs);
	fprintf(Parameters,"tau:     	%lf \n", tau);
    fprintf(Parameters,"trel:     	%lf \n", trel);
	fprintf(Parameters,"epsilon:    %lf \n", epsilon);
    fclose(Parameters);
  //Open output file and 
	Energy_out = fopen("Spin.dat","w+");

  //3 Start iteration over trajectories (Runs) 
	for( int u=0; u<Runs; u++){
	std::cout << "Run " << u <<" out of " << Runs <<" \n";
  //4 Initialize system variables 
	//4.1 Initialize magnetic field with time independent external field  
	for( int j=0; j<ssize; j++){
		for( int k=0; k<3; k++){
		  HFieldA[j][k] = HField[k];
		  HFieldB[j][k] = HField[k];
		}
	  }
	  
	//4.2 Initialize fluctuating coupling constants  
	for( int j=0; j<ssize; j++){
		for( int k=0; k<3; k++){
		  if(k!=2){		  
			JCplA[j][k] = 1 + Jdist(e2);
			JCplB[j][k] = 1 + Jdist(e2);        
		  }else{
			JCplA[j][k] = lambda + Jdist(e2);
			JCplB[j][k] = lambda + Jdist(e2);
		  }
		}
	  }  
	  
	//4.3 Initialize spin vectors
	//4.3.1 Initialize random configuration (infinite temperature)
	for( int j=0; j<ssize; j++){
		theta 		= 2*Pi*MCRxdist(e2);
		phi   		= acos(1-2*MCRxdist(e2)); 
		SpinA[j][0] = sin(phi)*cos(theta);
		SpinA[j][1] = sin(phi)*sin(theta); 
		SpinA[j][2] = cos(phi);
		
		theta 		= 2*Pi*MCRxdist(e2);
		phi   		= acos(1-2*MCRxdist(e2)); 
		SpinB[j][0] = sin(phi)*cos(theta);
		SpinB[j][1] = sin(phi)*sin(theta); 
		SpinB[j][2] = cos(phi);
	}
	
	//4.3.2 MC sampling of thermal state at inverse temperature Beta 
	for( int n=0; n<MCSmp+1; n++){
	
	//4.3.2.1 Take initial energy, write output and copy spin config 
	MCEi = Esys(SpinA, SpinB, HFieldA, HFieldB, JCplA, JCplB);
	
	for( int j=0; j<ssize; j++){
		for( int k=0; k<3; k++){
			MCSpinA[j][k] = SpinA[j][k]; 
			MCSpinB[j][k] = SpinB[j][k]; 
		}
	}
	
	//4.3.2.2.Perform random change of spin config 
	for( int j=0; j<ssize; j++){
		theta 		= 2*Pi*MCRxdist(e2);
		phi   		= acos(1-2*MCRxdist(e2)); 
		MCField[0] 	= sin(phi)*cos(theta); 
		MCField[1] 	= sin(phi)*sin(theta); 
		MCField[2] 	= cos(phi);
		StrthL = MCAgdist(e2); 
		if(StrthL==0){}else{
			for( int k=0; k<3; k++){SpinL[k] = SpinA[j][k];}
			for( int k=0; k<3; k++){SpinA[j][k] = MRot(SpinL, MCField, StrthL, k);}
		}
		
		theta 		= 2*Pi*MCRxdist(e2);
		phi   		= acos(1-2*MCRxdist(e2)); 
		MCField[0] 	= sin(phi)*cos(theta); 
		MCField[1] 	= sin(phi)*sin(theta); 
		MCField[2] 	= cos(phi);
		StrthL = MCAgdist(e2);  
		if(StrthL==0){}else{
			for( int k=0; k<3; k++){SpinL[k] = SpinB[j][k];}
			for( int k=0; k<3; k++){SpinB[j][k] = MRot(SpinL, MCField, StrthL, k);}
		}
	}
	
	//4.3.2.3 Take final energy, check update condition, restore initial config if necessary 
	MCEf = Esys(SpinA, SpinB, HFieldA, HFieldB, JCplA, JCplB);
	if(MCEf > MCEi){
		MCInd = MCRxdist(e2);
		if(MCInd > exp(-Beta*(MCEf-MCEi))){
			for( int j=0; j<ssize; j++){
				for( int k=0; k<3; k++){
					SpinA[j][k] = MCSpinA[j][k]; 
					SpinB[j][k] = MCSpinB[j][k]; 
				}
			}
		}
	}		
	}

	// 5.1 Evolve spin configiration for T_init time steps
	for(double t=0; t<=T_init; t=t+dt){
	//5.1 Set external field and evaluate observables (stroboscopically)
	if(fmod(t+dt/20,tau)<= dt/2){
		if(t>=trel){
			for( int j=0; j<ssize; j++){
				for( int k=0; k<3; k++){
					HFieldA[j][k] = HField[k] + Hext(t-trel,k,tau);
					HFieldB[j][k] = HField[k] + Hext(t-trel,k,tau);
				}
			}   
		}
	Usys = Esys(SpinA, SpinB, HFieldA, HFieldB, JCplA, JCplB);
	fprintf(Energy_out,"%lf %lf \n", t/tau, Usys/(2*ssize)); 
	}
								   
	 //5.2 Update external field 
	 if(t>=trel){
		for( int j=0; j<ssize; j++){
		  for( int k=0; k<3; k++){
			HFieldA[j][k] = HField[k] + Hext(t-trel+dt/2,k,tau);
			HFieldB[j][k] = HField[k] + Hext(t-trel+dt/2,k,tau);
		  }
		}
	}
	  //5.3 Propagate spin configuration on A 
	  for( int j=0; j<ssize; j++){
  
		for( int k=0; k<3; k++){FieldL[k] = AField(SpinB, HFieldA, JCplA, JCplB, j, k);}
		StrthL = sqrt(pow(FieldL[0],2) + pow(FieldL[1],2) + pow(FieldL[2],2));
	  
		if(StrthL==0){}else{
		for( int k=0; k<3; k++){FieldL[k] = FieldL[k]/StrthL;
								SpinL[k] = SpinA[j][k];}
		for( int k=0; k<3; k++){SpinA[j][k] = MRot(SpinL, FieldL, StrthL*dt/2, k);}
		}
	  }     
	
	  //5.4 Propagate spin configuration on B
	  for( int j=0; j<ssize; j++){
  
		for( int k=0; k<3; k++){FieldL[k] = BField(SpinA, HFieldB, JCplA, JCplB, j, k);}
		StrthL = sqrt(pow(FieldL[0],2) + pow(FieldL[1],2) + pow(FieldL[2],2));
		
		if(StrthL==0){}else{
		for( int k=0; k<3; k++){FieldL[k] = FieldL[k]/StrthL;
								SpinL[k] = SpinB[j][k];}
		for( int k=0; k<3; k++){SpinB[j][k] = MRot(SpinL, FieldL, StrthL*dt, k);}
		}
	  }
	  
	  //5.5 Propagate spin configuration on A 
	  for( int j=0; j<ssize; j++){
  
		for( int k=0; k<3; k++){FieldL[k] = AField(SpinB, HFieldA, JCplA, JCplB, j, k);}
		StrthL = sqrt(pow(FieldL[0],2) + pow(FieldL[1],2) + pow(FieldL[2],2));
	  
		if(StrthL==0){}else{
		for( int k=0; k<3; k++){FieldL[k] = FieldL[k]/StrthL;
								SpinL[k] = SpinA[j][k];}
		for( int k=0; k<3; k++){SpinA[j][k] = MRot(SpinL, FieldL, StrthL*dt/2, k);}
		}
	  }
  }
	//6.1 Duplicate spins 
	for( int j=0; j<ssize; j++){
		for( int k=0; k<3; k++){
			SpinA1[j][k] = SpinA[j][k];
			SpinB1[j][k] = SpinB[j][k];

			SpinA2[j][k] = SpinA[j][k];
			SpinB2[j][k] = SpinB[j][k];
		}
	}
	// 6.2 make small change to one configuration by rotating about an axis norm = Z x S_0 by an angle theta = epsilon
	theta 		= epsilon;
	phi   		= 0; 
	crossProduct(z_ax,SpinA2[0],norm);
	for( int k=0; k<3; k++){norm[k] = norm[k]/sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);}
	for( int k=0; k<3; k++){SpinL[k] = SpinA2[0][k];}
	for( int k=0; k<3; k++){SpinA2[0][k] = MRot(SpinL, norm, StrthL, k);}
  	//7 Evolve spin configurations 
	int t_step = 0;
	for(double t=0; t<=T; t=t+dt){

	//7.1 Set external field and evaluate observables (stroboscopically)
	if(fmod(t+dt/20,tau)<= dt){
		std::cout << t/tau;
		std::cout << "\n";
	for( int j=0; j<ssize; j++){
		for( int k=0; k<3; k++){
			HFieldA[j][k] = HField[k] + Hext(t+T_init,k,tau);
			HFieldB[j][k] = HField[k] + Hext(t+T_init,k,tau);
		}
	}   
	
	//7.2 evaluate correlators for each position in chain
  	for( int j=0; j<ssize; j++){
		//7.2.1 evaluate correlators at each even position (sublattice A) at time t
		corrA[j] = (SpinA1[j][0]*SpinA2[j][0] + SpinA1[j][1]*SpinA2[j][1] + SpinA1[j][2]*SpinA2[j][2]);
		//7.2.2 append them to the total correlator array at even positions 2*j
		corr[2*j][t_step] = corr[2*j][t_step] + corrA[j];

		//7.2.3 evaluate correlators at each odd position (sublattice B)
		corrB[j] = (SpinB1[j][0]*SpinB2[j][0] + SpinB1[j][1]*SpinB2[j][1] + SpinB1[j][2]*SpinB2[j][2]);
		//7.2.4 append them to the total correlator array odd positions 2*j +1
		corr[2*j+1][t_step] = corr[2*j+1][t_step] + corrB[j];
		}
	//7.2.5 increment time step
	t_step = t_step + 1;
	
	Usys = Esys(SpinA1, SpinB1, HFieldA, HFieldB, JCplA, JCplB);
	fprintf(Energy_out,"%lf %lf \n", t/tau + T_init/tau, Usys/(2*ssize));

	Usys = Esys(SpinA2, SpinB2, HFieldA, HFieldB, JCplA, JCplB);
	fprintf(Energy_out,"%lf %lf \n", t/tau + T_init/tau, Usys/(2*ssize));
	}
								   
	//7.3 Update external field 
	for( int j=0; j<ssize; j++){
		for( int k=0; k<3; k++){
		HFieldA[j][k] = HField[k] + Hext(t+T_init+dt/2,k,tau);
		HFieldB[j][k] = HField[k] + Hext(t+T_init+dt/2,k,tau);
		}
	}
	//7.4 Propogate copy 1 using Suzuki-Troter decomposition
	  //7.4.1 Propagate spin configuration on A1
	  for( int j=0; j<ssize; j++){
  
		for( int k=0; k<3; k++){FieldL[k] = AField(SpinB1, HFieldA, JCplA, JCplB, j, k);}
		StrthL = sqrt(pow(FieldL[0],2) + pow(FieldL[1],2) + pow(FieldL[2],2));
	  
		if(StrthL==0){}else{
		for( int k=0; k<3; k++){FieldL[k] = FieldL[k]/StrthL;
								SpinL[k] = SpinA1[j][k];}
		for( int k=0; k<3; k++){SpinA1[j][k] = MRot(SpinL, FieldL, StrthL*dt/2, k);}
		}
	  }     
	
	  //7.4.2 Propagate spin configuration on B1
	  for( int j=0; j<ssize; j++){
  
		for( int k=0; k<3; k++){FieldL[k] = BField(SpinA1, HFieldB, JCplA, JCplB, j, k);}
		StrthL = sqrt(pow(FieldL[0],2) + pow(FieldL[1],2) + pow(FieldL[2],2));
		
		if(StrthL==0){}else{
		for( int k=0; k<3; k++){FieldL[k] = FieldL[k]/StrthL;
								SpinL[k] = SpinB1[j][k];}
		for( int k=0; k<3; k++){SpinB1[j][k] = MRot(SpinL, FieldL, StrthL*dt, k);}
		}
	  }
	  
	  //7.4.3 Propagate spin configuration on A1
	  for( int j=0; j<ssize; j++){
  
		for( int k=0; k<3; k++){FieldL[k] = AField(SpinB1, HFieldA, JCplA, JCplB, j, k);}
		StrthL = sqrt(pow(FieldL[0],2) + pow(FieldL[1],2) + pow(FieldL[2],2));
	  
		if(StrthL==0){}else{
		for( int k=0; k<3; k++){FieldL[k] = FieldL[k]/StrthL;
								SpinL[k] = SpinA1[j][k];}
		for( int k=0; k<3; k++){SpinA1[j][k] = MRot(SpinL, FieldL, StrthL*dt/2, k);}
		}
	  }

	//7.5 Propogate copy 2 using Suzuki-Troter decomposition
	  //7.5.1 Propagate spin configuration on A2
	  for( int j=0; j<ssize; j++){
  
		for( int k=0; k<3; k++){FieldL[k] = AField(SpinB2, HFieldA, JCplA, JCplB, j, k);}
		StrthL = sqrt(pow(FieldL[0],2) + pow(FieldL[1],2) + pow(FieldL[2],2));
	  
		if(StrthL==0){}else{
		for( int k=0; k<3; k++){FieldL[k] = FieldL[k]/StrthL;
								SpinL[k] = SpinA2[j][k];}
		for( int k=0; k<3; k++){SpinA2[j][k] = MRot(SpinL, FieldL, StrthL*dt/2, k);}
		}
	  }     
	
	  //7.5.2 Propagate spin configuration on B2
	  for( int j=0; j<ssize; j++){
  
		for( int k=0; k<3; k++){FieldL[k] = BField(SpinA2, HFieldB, JCplA, JCplB, j, k);}
		StrthL = sqrt(pow(FieldL[0],2) + pow(FieldL[1],2) + pow(FieldL[2],2));
		
		if(StrthL==0){}else{
		for( int k=0; k<3; k++){FieldL[k] = FieldL[k]/StrthL;
								SpinL[k] = SpinB2[j][k];}
		for( int k=0; k<3; k++){SpinB2[j][k] = MRot(SpinL, FieldL, StrthL*dt, k);}
		}
	  }
	  
	  //7.5.3 Propagate spin configuration on A2
	  for( int j=0; j<ssize; j++){
  
		for( int k=0; k<3; k++){FieldL[k] = AField(SpinB2, HFieldA, JCplA, JCplB, j, k);}
		StrthL = sqrt(pow(FieldL[0],2) + pow(FieldL[1],2) + pow(FieldL[2],2));
	  
		if(StrthL==0){}else{
		for( int k=0; k<3; k++){FieldL[k] = FieldL[k]/StrthL;
								SpinL[k] = SpinA2[j][k];}
		for( int k=0; k<3; k++){SpinA2[j][k] = MRot(SpinL, FieldL, StrthL*dt/2, k);}
		}
	  }

     }
	}
  fclose(Energy_out);
  //8 Output correlator to file
  Output = fopen("OTOC_Dr.dat","w+");
  //8.1 loop over times and positions
  for(int time = 0; time<max_t; time++){
    for(int pos = 0; pos<2*ssize; pos++){
	  //8.2 average correlator by dividing by number of runs 
	  corr[pos][time] = 1 - corr[pos][time]/Runs;
	  //8.3 Output to file
	  fprintf(Output, "%lf	", corr[pos][time]);
	  }
	//8.4 create newline in output for next time step
	fprintf(Output,  "\n");
	}
  fclose(Output);
  
  //9 Determine CPU time
  std::clock_t c_end = std::clock();
  double time_elapsed_ms = 1000.0*(c_end-c_start)/CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms/1000 << " s\n";
  
  Parameters = fopen("Parameters.dat","a");
  fprintf(Parameters,"CPU[s]:   	%lf \n", time_elapsed_ms/1000);
  fclose(Parameters);
}