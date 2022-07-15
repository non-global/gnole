#include "Parameters.hh"
#include <iostream>

// set infrared cutoff of the evolution
double lnktmax = 0;
void set_lnktmax(double lnkt) {lnktmax = lnkt;};

// set collinear cutoff for dipole emissions
double RAPMAX = 0;
void set_rapmax(double etamax) {RAPMAX = etamax;};

// set alphas at Q (hard scale)
double as = 0; // Default value
void set_alphas_at_Q(double alphas) {as = alphas;};


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

// // second order evolution of alphas
// double alphas2(double lnkt) {
//   if (2.0 * as * b0 * lnkt >= 1.0) return 0.0;
//   return as/(1.-2.*as*b0*lnkt)*(1.-as*b1/b0*log(1.-2.*as*b0*lnkt)
// 				/(1.0 - 2.0*as*b0*lnkt));
// }
