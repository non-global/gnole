#include "Parameters.hh"
#include <iostream>

bool PRECISE_KERNELS = true;
void set_precise_kernels(bool enabled) { PRECISE_KERNELS = enabled; }

// set infrared cutoff of the evolution
double lnktmax = 0;
void set_lnktmax(double lnkt) {lnktmax = lnkt;};

// set collinear cutoff for dipole emissions
double RAPMAX = 0;
void set_rapmax(double etamax) {RAPMAX = etamax;};

// set alphas at Q (hard scale)
double as = 0; // Default value
void set_alphas_at_Q(double alphas) {as = alphas;};

// set HARD_ONLY flag
bool HARD_ONLY = false;
void set_hard_only(bool hard_only) {HARD_ONLY = hard_only;};

// set SINGLE_INSERTION_ONLY flag
bool SINGLE_INSERTION_ONLY = false;
void set_single_insertion_only(bool single_insertion_only) {SINGLE_INSERTION_ONLY = single_insertion_only;};

// set DOUBLE_INSERTION_ONLY flag
bool DOUBLE_INSERTION_ONLY = false;
void set_double_insertion_only(bool double_insertion_only) {DOUBLE_INSERTION_ONLY = double_insertion_only;};

bool DOUBLE_INSERTION_GAUGE_ONLY = false;
void set_double_insertion_gauge_only(bool gauge_only) {DOUBLE_INSERTION_GAUGE_ONLY = gauge_only;}
bool DOUBLE_INSERTION_FERMION_ONLY = false;
void set_double_insertion_fermion_only(bool fermion_only) {DOUBLE_INSERTION_FERMION_ONLY = fermion_only;}

// set MATCHING flag
bool MATCHING = false;
void set_matching(bool matching) {MATCHING = matching;};

// set SL_OBSERVABLE
bool SL_OBSERVABLE = false;
void set_sl_observable(bool sl_observable) {SL_OBSERVABLE = sl_observable;};

// first order evolution of alphas 
double alphas1(double lnkt) {
  if (2.0 * as * b0 * lnkt >= 1.0) return 0.0;
  return as / (1.0 - 2.0 * as * b0 * lnkt);
}

// second order evolution of alphas
double alphas2(double xmur) {
  double lnoxmur = -log(xmur);
  return as/(1.-2.*as*b0*lnoxmur)*(1.-as*b1/b0*log(1.-2.*as*b0*lnoxmur)
				/(1.0 - 2.0*as*b0*lnoxmur));
}

// The old expanded implementation also vetoed emissions beyond the endpoint
// of its unused inverse grid. Preserve that cutoff, including at unusual
// scale choices where it is below EVOLCUT, without allocating a grid.
double legacy_evolution_limit(double xmur, double xQ) {
  const int npoints = 1000;
  double Lmin = log(xQ);
  double Lmax = 1.0/(2.0*alphas2(xmur)*b0);
  double last_time = 0.0;
  for (int i = 0; i < npoints; ++i) {
    double L = Lmin + i*(Lmax-Lmin)/npoints;
    double asmur = alphas2(xmur);
    double time = (2*asmur*asmur*b0*b1*L + (b0 + asmur*b1 - 2*asmur*b0*b0*L)*
                   log(1 - 2*asmur*b0*L))/(4.*b0*b0*(-1 + 2*asmur*b0*L)*M_PI)
      + (asmur*KCMW*asmur*b0*L)/(4*b0*M_PI*M_PI - 8*b0*asmur*b0*L*M_PI*M_PI)
      + 2*asmur*b0*b0*M_PI*(-2*asmur*b0*L*log(xmur) + log(xQ))/
        (4.*b0*b0*(-1 + 2*asmur*b0*L)*M_PI*M_PI);
    if (i > 0 && time < last_time) break;
    last_time = time;
  }
  return last_time;
}

// // second order evolution of alphas
// double alphas2(double lnkt) {
//   if (2.0 * as * b0 * lnkt >= 1.0) return 0.0;
//   return as/(1.-2.*as*b0*lnkt)*(1.-as*b1/b0*log(1.-2.*as*b0*lnkt)
// 				/(1.0 - 2.0*as*b0*lnkt));
// }
