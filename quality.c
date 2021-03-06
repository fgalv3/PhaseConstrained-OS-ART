#include "cDEF.h"

double SSIM(vec RHO,vec SAMPLE, int nn){
        
        double AVphantom=0.0;
        double AVERAGE=0.0;
        for(int i=0; i<nn; i++){
            AVERAGE+=abs(RHO[i])/nn;
            AVphantom+=SAMPLE[i]/nn;
        }
        
        //normalize reconstructed image to be at levels similar to phantom20red
        for(int t=0; t < nn; t++){
            RHO[t]=abs(RHO[t])*(AVphantom/AVERAGE);///MAX;
        }
        
        //now recalculate average (should be equal to AVphantom, but just to check..)
        AVERAGE=0.0;
        for(int i=0; i<nn; i++){
                AVERAGE+=abs(RHO[i])/nn;
        }    
        double xi,xiPH;
        double PEARSON=0.0;
        double sigma=0.0;
        double sigmaPH=0.0;
        double sigmaXY=0.0;
        
        for(int t=0; t < nn; t++){
            xi=abs(RHO[t])-AVERAGE;
            xiPH=SAMPLE[t]-AVphantom;
            sigma+=pow(xi,2.0);
            sigmaPH+=pow(xiPH,2.0);
        }
        sigma=sqrt(sigma/nn);
        sigmaPH=sqrt(sigmaPH/nn);
        
        for(int t=0; t < nn; t++){
            xi=abs(RHO[t])-AVERAGE;
            xiPH=SAMPLE[t]-AVphantom;
            sigmaXY+=xi*xiPH/nn;
        }

        double C1=0.01;double C2=0.03;
        
        //SSIM
        double SSIM=(2.0*AVERAGE*AVphantom+C1)/(pow(AVERAGE,2.0)+pow(AVphantom,2.0)+C1);
        SSIM*=(2.0*sigmaXY+C2)/(pow(sigma,2.0)+pow(sigmaPH,2.0)+C2);
        
        return SSIM;
        }

void print2Dquality(double tTOT, vec RHO,vec SAMPLE, int nn){
        
        double AVphantom=0.0;
        double AVERAGE=0.0;
        for(int i=0; i<nn; i++){
            AVERAGE+=abs(RHO[i])/nn;
            AVphantom+=SAMPLE[i]/nn;
        }
        
        //normalize reconstructed image to be at levels similar to phantom20red
        for(int t=0; t < nn; t++){
            RHO[t]=abs(RHO[t])*(AVphantom/AVERAGE);///MAX;
        }
        
        //now recalculate average (should be equal to AVphantom, but just to check..)
        AVERAGE=0.0;
        for(int i=0; i<nn; i++){
                AVERAGE+=abs(RHO[i])/nn;
        }    
        double xi,xiPH;
        double PEARSON=0.0;
        double sigma=0.0;
        double sigmaPH=0.0;
        double sigmaXY=0.0;
        
        for(int t=0; t < nn; t++){
            xi=abs(RHO[t])-AVERAGE;
            xiPH=SAMPLE[t]-AVphantom;
            sigma+=pow(xi,2.0);
            sigmaPH+=pow(xiPH,2.0);
        }
        sigma=sqrt(sigma/nn);
        sigmaPH=sqrt(sigmaPH/nn);
        
        for(int t=0; t < nn; t++){
            xi=abs(RHO[t])-AVERAGE;
            xiPH=SAMPLE[t]-AVphantom;
            sigmaXY+=xi*xiPH/nn;
        }
        
        PEARSON=sigmaXY/(sigma*sigmaPH);
        
        //RMS
        double MSE=0.0;
        for(int t=0; t < nn; t++){
            MSE+=pow(abs(RHO[t])-SAMPLE[t],2.0)/nn;
        }
        double RMS=sqrt(MSE);
        
        //PSNR
        double PSNR=10*log10(1.0/MSE);
        
        
        //abs error
        double ABS=0.0;
        for(int t=0; t < nn; t++){
            ABS+=abs(abs(RHO[t])-SAMPLE[t]);
        }
        ABS=100.0*ABS/nn;
        
        //count voxels above 50%
        double c70=0.0;
        for(int i=0; i < nn; i++){
            if(abs(abs(RHO[i])-SAMPLE[i])>0.5){c70++;}
        }
        
        c70=100.0*c70/nn;
        
        double C1=0.01;double C2=0.03;
        double contr=(2.0*sigma*sigmaPH+C2)/(pow(sigma,2.0)+pow(sigmaPH,2.0)+C2);
        
        //SSIM
        double SSIM=(2.0*AVERAGE*AVphantom+C1)/(pow(AVERAGE,2.0)+pow(AVphantom,2.0)+C1);
        SSIM*=(2.0*sigmaXY+C2)/(pow(sigma,2.0)+pow(sigmaPH,2.0)+C2);
        
        printf("%2.10f  %2.5f  %2.5f  %2.5f  %2.5f  %2.5f  %2.5f\n",tTOT,PEARSON*contr,c70,ABS,RMS,PSNR,SSIM);

        
        } 
