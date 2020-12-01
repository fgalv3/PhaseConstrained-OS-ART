#include "cDEF.h"

cx_vec sigGen(vec kSPACE, int tSTEPS, vec SAMPLE, int nx){
    cx_vec signal(tSTEPS,fill::zeros);
    double SINC;
    double x,y,GX,GY;
    double dx=FoV/nx;
    
    #pragma omp parallel for shared(signal) private(GX,GY) num_threads(n_threads) 
    for(int t=0;t<tSTEPS;t++){
        GX=kSPACE[t*2];
        GY=kSPACE[t*2+1];
        SINC=2.0*sin(dx*GX/2.0)*sin(dx*GY/2.0)/(GX*GY);
        for(int i = 0; i <nx; i++){
            x=-FoV/2.0+i*FoV/(nx-1);
            for(int j =0; j < nx; j++){
                y=-FoV/2.0+j*FoV/(nx-1);
                signal(t)+=SAMPLE(i + j*nx)*exp((GX*x+GY*y)*II)*SINC;
            }
        }
    }
    return signal;
}

vec PWS(vec kSPACE, cx_vec signal, int tSTEPS, int NX){
    int nn=NX*NX;
    cx_vec RHO(nn,fill::zeros);
    double x,y,GX,GY;
    
    #pragma omp parallel for shared(RHO,signal) num_threads(n_threads)
    for(int i = 0; i <NX; i++){
        x=-FoV/2.0+i*FoV/(NX-1);
        for(int j =0; j < NX; j++){
            y=-FoV/2.0+j*FoV/(NX-1);
            for(int t=0;t<tSTEPS;t++){
                GX=kSPACE[t*2];
                GY=kSPACE[t*2+1];
                RHO(i + j*NX)+=exp(-(GX*x+GY*y)*II)*signal(t);
            }
        }
    }
    
    return abs(RHO);
}

vec ARTabs(vec kSPACE, cx_vec signal, int tSTEPS, int NX, double lambda, int nITERATIONS){
    int nn=NX*NX;
    cx_vec RHO(nn,fill::zeros);
    cx_vec Mt(nn);
    double x,y,GX,GY;
    std::complex<double> error=0.0;
    
    for(int indice=0; indice<nITERATIONS;indice++){        
        for(int t=0;t<tSTEPS;t++){
            GX=kSPACE[t*2];
            GY=kSPACE[t*2+1];
                                   
            #pragma omp parallel for shared(Mt) num_threads(n_threads)
            for(int i = 0; i <NX; i++){
                x=-FoV/2.0+i*FoV/(NX-1);
                for(int j =0; j < NX; j++){
                    y=-FoV/2.0+j*FoV/(NX-1);
                    Mt(i + j*NX )=exp((GX*x+GY*y)*II);
                }
            }
            error=(signal(t)-dot(Mt,RHO))/nn;
            //ART update step
            RHO+=(lambda*error)*conj(Mt);
            //projection to real values
            for(int n=0;n<nn;n++){RHO(n)=abs(RHO(n));}
        }
    }
    
    return abs(RHO);
}

vec ARTlinear(vec kSPACE, cx_vec signal, int tSTEPS, int NX, double lambda, int nITERATIONS){
    int nn=NX*NX;
    cx_vec RHO(nn,fill::zeros);
    cx_vec Mt(nn);
    double x,y,GX,GY;
    std::complex<double> error=0.0;
    
    for(int indice=0; indice<nITERATIONS;indice++){        
        for(int t=0;t<tSTEPS;t++){
            GX=kSPACE[t*2];
            GY=kSPACE[t*2+1];
                                   
            #pragma omp parallel for shared(Mt) num_threads(n_threads)
            for(int i = 0; i <NX; i++){
                x=-FoV/2.0+i*FoV/(NX-1);
                for(int j =0; j < NX; j++){
                    y=-FoV/2.0+j*FoV/(NX-1);
                    Mt(i + j*NX )=exp((GX*x+GY*y)*II);
                }
            }
            error=(signal(t)-dot(Mt,RHO))/nn;
            //ART update step
            RHO+=(lambda*error)*conj(Mt);
        }
    }
    
    return abs(RHO);
}

vec ARTtv(vec kSPACE, cx_vec signal, int tSTEPS, int NX, double lambda, int nITERATIONS, double beta){
    int nn=NX*NX;
    cx_vec RHO(nn,fill::zeros);
    cx_vec TV(nn,fill::zeros);
    cx_vec Mt(nn);
    double x,y,GX,GY;
    std::complex<double> error=0.0;
    
    double eps=pow(10,-20);
    
    for(int indice=0; indice<nITERATIONS;indice++){        
        for(int t=0;t<tSTEPS;t++){
            GX=kSPACE[t*2];
            GY=kSPACE[t*2+1];
                                   
            #pragma omp parallel for shared(Mt) num_threads(n_threads)
            for(int i = 0; i <NX; i++){
                x=-FoV/2.0+i*FoV/(NX-1);
                for(int j =0; j < NX; j++){
                    y=-FoV/2.0+j*FoV/(NX-1);
                    Mt(i + j*NX )=exp((GX*x+GY*y)*II);
                }
            }
            error=(signal(t)-dot(Mt,RHO))/nn;
            //ART update step
            RHO+=(lambda*error)*conj(Mt);
            
            //TV l1-norm: sqrt(x^2+eps)
            #pragma omp parallel for shared(TV) num_threads(n_threads)
            for(int i = 1; i <NX-1; i++){
                for(int j =1; j < NX-1; j++){
                    TV(i + j*NX )= (RHO(i+j*NX)-RHO((i-1)+j*NX))/sqrt(pow(RHO(i+j*NX)-RHO((i-1)+j*NX),2.0)+pow(RHO(i-1+(j+1)*NX)-RHO((i-1)+j*NX),2.0)+eps*eps);
                    TV(i + j*NX )+=(RHO(i+j*NX)-RHO(i+(j-1)*NX))/sqrt(pow(RHO(i+1+(j-1)*NX)-RHO(i+(j-1)*NX),2.0)+pow(RHO(i+j*NX)-RHO(i+(j-1)*NX),2.0)+eps*eps);
                    TV(i + j*NX )-=(RHO(i+1+j*NX)+RHO(i+(j+1)*NX)-2.0*RHO(i+j*NX))/sqrt(pow(RHO(i+1+j*NX)-RHO(i+j*NX),2.0)+pow(RHO(i+(j+1)*NX)-RHO(i+j*NX),2.0)+eps*eps);         
                }
            }
            RHO-=beta*max(abs(RHO))*TV;
            
            //projection to real values
            for(int n=0;n<nn;n++){RHO(n)=abs(RHO(n));}
        }
    }
    
    return abs(RHO);
}
