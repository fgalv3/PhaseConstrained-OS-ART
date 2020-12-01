/********************************************************************/
//  Take SAMPLE, generate signal from dense pixels, then reconstruct.
//  Reconst. algorithms:
//     * PWS (equivalent to DFT for dt=1/2BW, BW=gradient bandwith
//     * ARTabs (project RHO to real values)
//     * ARTlinear (do not project)
//     * ARTtv (use l1-approx TV penalty)
//
//  COMPILATION: icpc -march=native -mtune=native -Ofast -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG -I /usr/include/ -L /usr/lib -lblas -llapack -fopenmp -std=c++11 program.c kSpace.c quality.c reconst.c  -o thing.x 
//  LIBRARIES REQUIRED: armadillo-cpp
/********************************************************************/


#include "cDEF.h"

int main(int argc,char *argv[]){
    double tTOT=atof(argv[1]);
    double DT=atof(argv[2]);
    double lambda=atof(argv[3]);
    int nITERATIONS=atoi(argv[4]);
    double BETA=atof(argv[5]);
    int ACCEL=atoi(argv[6]);
    
    double dt=DT*pow(10,-6);
    
    int tSTEPS=(int)(tTOT/dt);
    
    int nx=256;
    int nn=nx*nx;
    
    vec SAMPLE(nn);
    FILE *archivo;
    archivo=fopen("SHEPP256x256.dat","r");
   
    for(int i=0; i<nn;i++){
        fscanf(archivo, "%lf",&SAMPLE(i));
    }
    
    //generate EPI single shot k-Space
    vec kS(tSTEPS*2,fill::zeros);
    kS=SPI(tTOT,dt,ACCEL);
    FILE *kspace;
    kspace=fopen("kSPACE.dat","w");
    for(int i=0;i<tSTEPS;i++){
        fprintf(kspace,"%2.20f %2.20f\n",kS[i*2],kS[i*2+1]);
    }
    fclose(kspace);
    
    //generate signal from SAMPLE
    cx_vec signal(tSTEPS,fill::zeros);
    signal=sigGen(kS,tSTEPS,SAMPLE,nx);
    FILE *sig;
    sig=fopen("SIGNAL.dat","w");
    for(int i=0;i<tSTEPS;i++){
        fprintf(sig,"%2.20f %2.20f\n",real(signal(i)),imag(signal(i)));
    }
    fclose(sig);
    
    //reconstruct using plane wave sum
    int NX=120;
    int NN=NX*NX;
    vec outRHO(NN,fill::zeros);
    
//     outRHO=PWS(kS,signal,tSTEPS, NX);//DFT
//     outRHO=ARTabs(kS,signal,tSTEPS,NX,lambda,nITERATIONS);
//     outRHO=ARTlinear(kS,signal,tSTEPS,NX,lambda,nITERATIONS);
    outRHO=ARTtv(kS,signal,tSTEPS,NX,lambda,nITERATIONS,BETA);
    
    //normalize reconstructed image to be at levels similar to phantom20red
    double AVERAGE=0.0;
    for(int i=0; i<NN; i++){
                AVERAGE+=outRHO[i]/NN;
        }
    for(int t=0; t < NN; t++){
        outRHO[t]=outRHO[t]/AVERAGE;
    }
    
    //print out reconstructed image
    FILE *out;
    out=fopen("outRHO.dat","w");
    for(int i=0;i<NN;i++){
        fprintf(out,"%2.20f\n",outRHO[i]);
    }
    fclose(out);
    
    return 0;
}
    
    
