#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libdcl.h"
#include <gsl/gsl_sf_bessel.h>

int n;

double f(double x,void *serviceData){
  return LegendreP(n,x);
}

double g(double x,void *serviceData){
  return x*x*x-3*x*x+2*x;
}

double f1(double x1[],void *serviceData){
  double x,y,z;
  x=x1[0];
  y=x1[1];
  z=x1[2];
  return x*x*x+x*sin(y*z);
}

double complex g1(double x1[],void *sd){
  return cexp(I*x1[0]);
}
void istart(){
  double sep,x0,x1;
  double *ksi;
  double *an;
  int i,nf;
  x0=0.0;
  x1=1.0;
  sep=1.0/(n*n);
  ksi=(double *)malloc(n*sizeof(double));
  an=(double *)malloc(n*sizeof(double));
  for(i=0;i<n;i++){
    ksi[i]=FindNZero(f,x0,x1,NULL,sep,i+1,&nf);
    an[i]=2.0/((1.0-ksi[i]*ksi[i])*(DLegendreP(n,ksi[i])*DLegendreP(n,ksi[i])));
//    printf("ksi[%d]=%.16lg\n",i,ksi[i]);
//    printf("an[%d]=%.16lg\n",i,an[i]);
  }
  printf("{");
  for(i=0;i<n;i++){
    printf("%.16lg,",ksi[i]);
  }
  printf("}\n");
  printf("{");
  for(i=0;i<n;i++){
    printf("%.16lg,",an[i]);
  }
  printf("}\n");
}

void F(int neq,double t,double * y,double rv[],void *serviceData){
  rv[0]=y[0]*y[0]+1;
}

int main(){
  int i,j;
  double x,y,t0,nsteps,dt;
  double complex z,p1,q1,res;
  double *a,*q,*r;
  int N;
  printf("Gamma()=%.16lg\n",log(Gamma(100.0)));
  printf("Ei()=%.16lg\n",ExpIntegralEi(0.01));
  printf("LogE=%.16lg\n",log(e_const));
  printf("C()=%.16lg\n",BinomCoeff(4,2));
  printf("Ylm()=(%.16lg,%.16lg\n",creal(SphericalHarmonicY(10,5,1,1)),cimag(SphericalHarmonicY(10,5,1,1)));
  printf("Hn=%.16lg\n",HermiteH(10,0.2));
  printf("Lnm=%.16lg\n",LaguerreL(10,5,0.2));
  printf("K()=%.16lg\n",EllipticE(0.1));
  printf("C=%.16lg\n",ClebschGordan(0.5,0.5,0.5,-0.5,1,0));
  print_complex("sf=",besselk(-2,2));
  return 0;
}