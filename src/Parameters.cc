#include "Parameters.hh"

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
