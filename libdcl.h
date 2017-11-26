#ifndef _LIB_DCL
#define _LIB_DCL
//Constants
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <openssl/rand.h>
#include <stdint.h>
#include <complex.h>
#include <omp.h>
#include "libamos.h"
#define pi          3.14159265358979323846264338327
#define e_const     2.71828182845904523536028747135
#define EulerGamma 0.577215664901532860606512090082
#define e0 4.803242e-10
#define m0 9.10938188e-28
#define hbar 1.05457159642e-27
#define kB 1.3806503e-16
#define evolt 1.602176462e-12
#define angstrem 1e-8
#define c_light 2.99792458e10

int increment_add(int *array,int * add,int N,int nmax);
double ClebschGordan(double j1,double m1,double j2,double m2,double j3,double m3);
int KronekerDelta(int l1,int l2);
double min(double a,double b);
double max(double a,double b);
double ChebyshevT(int n,double x);
double ChebyshevU(int n,double x);
double HermiteH(int n,double x);
double LaguerreL(int n,int m,double x);
double EllipticE(double m);
double EllipticK(double m);
double complex SphericalHarmonicY(int l,int m,double theta,double phi);
double Gamma(double x);
double ExpIntegralEi(double x);
double BinomCoeff(unsigned int n,unsigned int k);
double LegendrePlm(int l,int m,double x);
double GaussIntegrateElem(double (*f)(double[],void *),void * serviceData,int ndim,double a[],double b[]);
double complex ZGaussIntegrateElem(double complex(*f)(double[],void *),void * serviceData,int ndim,double a[],double b[]);
double GaussIntegrate(double (*f)(double[],void *),void * serviceData,int ndim,double a[],double b[],int m);
double complex ZGaussIntegrate(double complex(*f)(double[],void *),void * serviceData,int ndim,double a[],double b[],int m);
int IsNaN(double x);
int sign(double x);
double FindZero(double (*f)(double,void*),double a,double b,void *serviceData,int *errcode);
double FindNZero(double(*f)(double,void *),double x0,double x1,void *serviceData,double sep,int i,int * nf);
double LegendreP(int n,double x);
double DLegendreP(int n,double x);
void MatrixMatrixMultiply(double *r,double *a,double *b,int m,int n,int k);
double SQR(double x);
double pythag(double a,double b );
void QR_decompose_rotation(double *a,double *q,double *r,int n);
void QR_decompose_reflection(double *a,double *q,double *r,int n);
void printMatrix(double *a,char * title,int m,int n);
double random_double();
int rk4(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t0,double dt,int nsteps,void *serviceData);
int rk4_step(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t,double dt,void *serviceData);
int rk5(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t0,double dt,int nsteps,void *serviceData);
int rk5_step(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t,double dt,void *serviceData);
double ipow(double x,int k);
double complex ZGammaIncomplete(double s,double complex z);
double complex ZExpIntegralE(double s,double complex z);
void print_complex(char * msg,double complex z);
double complex besselj(double nu, double complex z);
double complex bessely(double nu, double complex z);
double complex besseli(double nu, double complex z);
double complex besselk(double nu, double complex z);
double dbesselj(double nu, double z);
double dbessely(double nu, double z);
double dbesseli(double nu, double z);
double dbesselk(double nu, double z);
double sj(int l,double x);
double sy(int l,double x);
double si(int l,double x);
double sk(int l,double x);

double complex BesselJ(double nu,double complex z);
double ifact(int k);
double complex rk(int k,double complex nu);
double complex B2k(unsigned int k2,double complex p);
double complex Bk(double k,double complex p);
double complex Ipq(double complex p,double complex q);


#endif