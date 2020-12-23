#include "cDEF.h"

vec EPI(double tTOT,double dt, int ACCEL){
        double dk=ACCEL*2.0*PI/FoV;
        int NkY=floor( pow(GAMMA*B1*tTOT/dk,1.0/2.0) );
        
        double tX=tTOT/NkY;
        double KMAX=GAMMA*B1*tX;
        
        int NkX=floor(tX/dt);
        
        double dkX=KMAX/NkX;
        double dkY=KMAX/NkY;
        
        int tSTEPS=NkX*NkY;
        vec kS(tSTEPS*2,fill::zeros);
        
        int t=0;
        for(int kx=0;kx<NkX;kx+=2){
            for(int ky=0;ky<NkY;ky++){
                kS[t*2]=(-(NkX-1)/2+kx)*dkX;
                kS[t*2+1]=(-(NkY-1)/2+ky)*dkY;
                t++;
            }
            for(int ky=0;ky<NkY;ky++){
                kS[t*2]=(-(NkX-1)/2+kx+1)*dkX;
                kS[t*2+1]=((NkY-1)/2-ky)*dkY;
                t++;
            }
        }
        
         return kS;
}

vec SPI(double tTOT,double dt, int ACCEL){
    
        double G0=B1/ACCEL;//this keeps gradients independent of acceleration factor
        double LAMBDA=1.0;
        double S=1000.0/ACCEL;//this keeps slew rate independent of acceleration factor
        double lam=1.0/(2.0*PI*FoV);
        double beta=GAMMA*S/(2.0*PI*lam);
        double a2=pow(9.0*beta/4.0,1.0/3.0);
        double Ts=pow(3.0*GAMMA*G0/(4.0*PI*lam*a2*a2),3.0);
        
        double theta=0.0;
        double origin=0.0;
        
        int tSTEPS=(int)(tTOT/dt);
        vec kS(tSTEPS*2);
        //Fill k-space
        for(int t=0; t < tSTEPS; t++){
            
            if(t*dt<=Ts){
                theta=(beta*t*dt*t*dt/2.0)/(LAMBDA+beta*pow(t*dt,4.0/3.0)/(2.0*a2));
            }
            else{
                origin=(beta*Ts*Ts/2.0)/(LAMBDA+beta*pow(Ts,4.0/3.0)/(2.0*a2));
                theta=sqrt(pow(origin,2.0)+GAMMA*G0*(t*dt-Ts)/(PI*lam) );
            }
            kS[t*2]  =ACCEL*theta*cos(theta)/FoV;
            kS[t*2+1]=ACCEL*theta*sin(theta)/FoV;
        }
        
        return kS;
}
