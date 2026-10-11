#ifndef __PARAMETERS_HH__
#define __PARAMETERS_HH__
#include <cmath>

// evolution cutoff
#define EVOLCUT   1.0

// QCD parameters
#define CA        3.0
#define TF        0.5
// Large Nc
#define CF        1.5
//#define NF        0.
// Full Nc
//#define CF        1.33333333333333333333
#define NF        5.0

// constants
static const double b0 = (11*CA - 2*NF)/(12*M_PI);
static const double b1 = (17*CA*CA - 5*CA*NF - 3*CF*NF)/(24*M_PI*M_PI);
static const double KCMW = CA*(67.0/18.0 - M_PI*M_PI/6.0) - TF*NF*10.0/9.0;
static const double H1 = CF*(-8.0 + 7.0/6.0*M_PI*M_PI);

// the following flag is for private use only.
// it decides whether for the expanded NLL corrections
// the algorithm uses two extra branches in addition
// to the default algorithm, or just one
static const bool insertion_is_part_of_NLL_ensemble = true;

// external variables
extern double lnktmax;
extern double RAPMAX;
extern double as;
// Stable arithmetic and precision-aware observable/histogram paths.
extern bool PRECISE_KERNELS;
void set_precise_kernels(bool enabled);
// flag to run only the hard contributions
extern bool HARD_ONLY;
// flag to run only the single insertion contributions
extern bool SINGLE_INSERTION_ONLY;
// flag to run only the double insertion contributions
extern bool DOUBLE_INSERTION_ONLY;
// flag to match to NLO
extern bool MATCHING;
// flag to compute observable in SL approximation
extern bool SL_OBSERVABLE;

// set infrared cutoff of the evolution
void set_lnktmax(double lnkt);

// set collinear cutoff for dipole emissions
void set_rapmax(double etamax);

// set alphas at Q (hard scale)  
void set_alphas_at_Q(double alphas);

// Preserve the original analytic-grid endpoint veto as a scalar cutoff.
// This is not an evolution kernel or an inverse interpolation table.
double legacy_evolution_limit(double xmur, double xQ);

// set HARD_ONLY flag 
void set_hard_only(bool hard_only);

// set SINGLE_INSERTION_ONLY flag
void set_single_insertion_only(bool single_insertion_only);

// set DOUBLE_INSERTION_ONLY flag
void set_double_insertion_only(bool double_insertion_only);

// set MATCHING flag 
void set_matching(bool matching);

// compute observable in SL approximation
void set_sl_observable(bool sl_observable);

// first order evolution of alphas 
double alphas1(double lnkt);

// second order evolution of alphas
double alphas2(double lnkt);

#endif // __PARAMETERS_HH__
