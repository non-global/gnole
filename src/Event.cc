#include "Event.hh"
#include "Parameters.hh"

//----------------------------------------------------------------------
/// reset shower to two jet configuration
void Event::reset_twojet() {
  dipoles_.clear();
  Momentum q (0.0, 0.0, -0.5, 0.5);
  Momentum qb(0.0, 0.0,  0.5, 0.5);
  q=(1/q.E())*q;
  qb=(1/qb.E())*qb;
  add(Dipole(DipoleEnd(q,  true),
	     DipoleEnd(qb, true),
	     -1, -1));
  eta_tot = dipoles_.back().delta_rap();
  thrust_axis_= Momentum(0,0,1,1);
  weight = 1.0;
}

//----------------------------------------------------------------------
// reset shower to three jet configuration and return lnkt
double Event::reset_threejet_LL(double xmur, double xQ, double r1, double r2, double r3, Momentum& gluon) {
  dipoles_.clear();
  add(Dipole(DipoleEnd(Momentum(0.0,0.0,-0.5,0.5), true),
	     DipoleEnd(Momentum(0.0,0.0, 0.5,0.5), true),
	     -1, -1));
  eta_tot = dipoles_.back().delta_rap();
  double asmur = alphas2(xmur);
  double lnkt = r1 * 1.0/(2.0*asmur*b0);
  gluon = dipoles_.back().radiate(lnkt-log(xQ), r2, r3);
  thrust_axis_= Momentum(0,0,1,1);
  weight = eta_tot/(2.0*asmur*b0);
  weight *= 4 * CF * asmur / (2.0 * M_PI);
  return lnkt;
}
  
//----------------------------------------------------------------------
/// reset shower to three jet configuration and return lnkt, taken from Dokshitzer paper
double Event::reset_threejet(double xmur, double xQ, double r1, double r2, double r3, Momentum& gluon) {
  dipoles_.clear();
  double jacobian = 1.0;
  double asmur = alphas2(xmur);
  
  // angle between emission and quark 1
  double t13 = M_PI*r1;
  jacobian *= M_PI*sin(t13);
  double a13 = 1.0 - cos(t13);
  //double a13;
  //if (abs(t13) > pow(10,-3)) {
  //  a13 = 1.0 - cos(t13);
  //} else {
  //  a13 = pow(t13,2)/2.;
  //}
  if (abs(1.0-a13)==1.0) jacobian=0.0;
  
  // generate x3 and phi13
  double x3 = r2;
  double phi13 = r3 * 2.0 * M_PI;
  jacobian *= 2.0 * M_PI;

  // p3 now fixed by the above three and requiring on-shell (needs phi13 rotation)
  // PM: pg = Q/two*x3*(/sqrt(one-cost3**2)*sin(phi3),sqrt(one-cost3**2)*cos(phi3),cost3,one/)
  //double p3xy = sqrt(1.0 - (1.0 - a13)*(1.0 - a13));
  double p3xy = sqrt((2.0 - a13)*a13);
  Momentum p3(p3xy*sin(phi13), p3xy*cos(phi13), 1.0 - a13, 1.0);
  p3 = 0.5 * x3 * p3;
    
  // fix x1 using B.28
  double x1 = 2.0 * (1.0 - x3)/(2.0 - x3*a13);
  // fix p1
  Momentum p1(0.0, 0.0, 0.5 * x1, 0.5 * x1);

  // fix x2 using delta function
  double x2 = 2.0 - x1 - x3;
  // fix p2
  Momentum p2(-p1.px() - p3.px(), -p1.py() - p3.py(),
	      -p1.pz() - p3.pz(), 1.0 - p1.E() - p3.E());
  // from B.29
  double J3 = 4.0 * (1.0 - x3)/((2.0 - x3*a13)*(2.0 - x3*a13));
  // so that B.26
  jacobian *= x3 /(2.0 * M_PI) * J3;

  gluon = p3;
    
  double dot13 = dot3(p1,p3);
  double cost3_sq = dot13*dot13/(p1.modpsq()*p3.modpsq());
  double lnkt = -log(0.5*x3*sqrt(1.0 - cost3_sq));
  // pink book eq 3.31
  double matelm = 2 * CF * asmur / (2.0 * M_PI);
  matelm *= 0.5*(x1*x1 + x2*x2)/((1.-x1)*(1.-x2));
  // factor two error in Dokshitzer et al
  matelm /= 2.0;

  // various safety cutoffs
  if (std::abs(gluon.rap()) > RAPMAX) {
    weight = 0.0;
    bad = true;
    return 0.0;
  }
  // add a collinear regulator to reduce instabilities (helpful at small as)
  if (((1.-x1) < pow(10,-6)) or ((1.-x2) < pow(10,-6))) {
    weight = 0.0;
    bad = true;
    return 0.0;
  }
  // avoid return nans below Landau pole
  if (2.0*asmur*b0*(lnkt+log(xQ))>=1.0) {
    weight = 0.0;
    bad = true;
  }
  else weight = matelm * jacobian;
  // remove rare events that have nans
  if (weight!=weight) {
    weight = 0.0;
    bad = true;
    return 0.0;
  }

  add(Dipole(DipoleEnd(p1, true),  DipoleEnd(p2, true), -1, -1));
  eta_tot = dipoles_.back().delta_rap();
  // finally compute the thrust axis
  if ((p3.E() > p1.E()) and (p3.E() > p2.E()))      thrust_axis_=(1.0/p3.E())*p3;
  else if ((p2.E() > p1.E()) and (p2.E() > p3.E())) thrust_axis_=(1.0/p2.E())*p2;
  else                                              thrust_axis_=(1.0/p1.E())*p1;
  return lnkt;
}

