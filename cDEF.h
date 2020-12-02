#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

#define n_threads 16

const double PI=3.1415926535897932385;
const double GAMMA=2.67538*pow(10.0,8.0);//gyromagnetic ratio in rad/s·T
const double B1=0.1;//gradient in T/m
const double FoV=0.02;//FoV in m
const std::complex<double> II(0,1);

using namespace arma;

double SSIM(vec RHO,vec SAMPLE2, int nn2);
void print2Dquality(double tTOT, vec RHO,vec SAMPLE2, int nn2);

vec EPI(double tTOT,double dt, int ACCEL);
vec SPI(double tTOT,double dt, int ACCEL);

cx_vec sigGen(vec kSPACE, int tSTEPS, vec SAMPLE, int nx);

vec PWS(vec kSPACE, cx_vec signal, int tSTEPS, int NX);
vec ARTabs(vec kSPACE, cx_vec signal, int tSTEPS, int NX, double lambda, int nITERATIONS);
vec ARTlinear(vec kSPACE, cx_vec signal, int tSTEPS, int NX, double lambda, int nITERATIONS);
vec ARTtv(vec kSPACE, cx_vec signal, int tSTEPS, int NX, double lambda, int nITERATIONS, double beta);

