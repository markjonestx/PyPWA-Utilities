#ifndef plibH_INCLUDED
#define plibH_INCLUDED
#include <Vec.h>
#include <matrix.h>
matrix<double> lowq2DenMat(fourVec ein,fourVec eout);
double breit(double m,double mr,double gamma);
double gauss(double m,double mr,double sigma);
double breitW(double m,double mr,double gamma);
double betaCM(double pbeam,double mbeam,double mtarget);
double gammaCM(double pbeam,double mbeam,double mtarget);
double plab(double m,double pbeam,double mbeam,double mtarget,double thetarCM);
double ener(double p,double m);
double pprime(double w,double m,double m3);
double s(double pbeam,double mbeam,double mtarget);
double tmin(double plab,double ma,double mb,double mc,double md);
double randm (double low, double high);
double multscat(double p,double m,double x);
double Qsq(double E,double Ep,double theta);
double epsilon(double E,double Ep,double theta);
double epsilonL(double E,double Ep,double theta);
double  gammq(double a,double x);
void gcf(double *gammcf,double a,double x,double *gln);
void gser(double *gamser,double a,double x,double *gln);
double gammln(double xx);
int gpoisson(double nu);
void nrerror(char *);
double expDist(double slope);
double Fdist(double f,int nu1,int nu2);
double gammaFnc(double x);
double factrl(int n);
#endif
