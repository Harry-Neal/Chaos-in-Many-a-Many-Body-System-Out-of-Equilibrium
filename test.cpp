#include <iostream>
#include <math.h>
int main() {
    int T =200;
    double dt = 0.02;
    int t_step = 0;
    int w_indx = 0;
	for(double t=0; t<=T; t=t+dt){
	    if(t_step%100==0){
            std::cout << w_indx;
            w_indx = w_indx+1;
        }
    t_step=t_step+1;
    }
}