#ifndef __PARAMETERS_HH__
#define __PARAMETERS_HH__
#include <cmath>

// number of points for grid, max value, and integrand calls
#ifdef NNET
#define NGRID  1000
#define LNEVOLMAX 4
#define LNEVOLMIN -8
#define SUDAMAX 7.0
#endif
#define NEVLGRID  500
#define NSUDGRID  200
#define NCALLSUDA 1e4

// medium stat parameters
// #define NEVLGRID  400
// #define NSUDGRID  600
// #define NCALLSUDA 1e5
// high stat parameters
// #define NEVLGRID  200
// #define NSUDGRID  500
// #define NCALLSUDA 1e6

// rapidity and min kt bounds
#define RAPMAX    5.0
#define KTSOFT    1.0e-8

// evolution cutoff
#define EVOLCUT   1.0

// QCD parameters
#define CA        3.0
#define TF        0.5
// Large Nc
//#define CF        1.5
// Full Nc
#define CF        1.33333333333333333333
#define NF        5.0

// constants
static const double b0 = (11*CA - 2*NF)/(12*M_PI);
static const double b1 = (17*CA*CA - 5*CA*NF - 3*CF*NF)/(24*M_PI*M_PI);
static const double as = 0.118;
static const double KCMW = CA*(67.0/18.0 - M_PI*M_PI/6.0) - TF*NF*10.0/9.0;
static const double H1 = CF*(-8.0 + 7.0/6.0*M_PI*M_PI);
  
// first order evolution of alphas 
double alphas1(double lnkt);

// second order evolution of alphas
double alphas2(double lnkt);

#endif // __PARAMETERS_HH__
