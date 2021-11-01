/// Matrix elements for dipole emissions.
#ifndef __MATRIX_ELEMENTS_HH__
#define __MATRIX_ELEMENTS_HH__
#include <gsl/gsl_sf_dilog.h>

inline double integrated_counterterm(double delta_eta) {
  double c = (exp(delta_eta) - 1.0) / (exp(delta_eta) + 1.0);
  double li2OneMinusC2 = gsl_sf_dilog((1.0 - c)/2.0);
  double li2OnePlusC2  = gsl_sf_dilog((1.0 + c)/2.0);
  double counter_term = (4*(57 - 7*M_PI*M_PI) - (3*(64*c - 64* c*c*c - 192*c*log(2.) - 64*c*c*c*log(2.) -
	8*log(1 - c) + 96*c*log(1 - c) - 144*c*c*log(1 - c) +
        32*c*c*c*log(1 - c) + 24*c*c*c*c*log(1 - c) + 40*log(2.)*log(1 - c) -
        80*c*c*log(2.)*log(1 - c) + 40*c*c*c*c*log(2.)*log(1 - c) -
        3*log(256.)*log(1 - c) + 6*c*c*log(256.)*log(1 - c) -
        3*c*c*c*c*log(256.)*log(1 - c) - 8*log(1 - c)*log(1 - c) + 16*c*c*log(1 - c)*log(1 - c) -
        8*c*c*c*c*log(1 - c)*log(1 - c) + 8*log(1 + c) + 96*c*log(1 + c) +
        144*c*c*log(1 + c) + 32*c*c*c*log(1 + c) - 24*c*c*c*c*log(1 + c) -
        40*log(2.)*log(1 + c) + 80*c*c*log(2.)*log(1 + c) -
        40*c*c*c*c*log(2.)*log(1 + c) + 3*log(256.)*log(1 + c) -
        6*c*c*log(256.)*log(1 + c) + 3*c*c*c*c*log(256.)*log(1 + c) +
        8*log(1 + c)*log(1 + c) - 16*c*c*log(1 + c)*log(1 + c) + 8*c*c*c*c*log(1 + c)*log(1 + c) +
        16*(-1 + c*c)*(-1 + c*c)*li2OneMinusC2 -
	16*(-1 + c*c)*(-1 + c*c)*li2OnePlusC2))/((-1 + c*c)*(-1 + c*c)))/24.;
  return counter_term;
}

  // !---------------------------------------------------------------------------
  // ! this function returns the integrated counterterm for the slice problems
  // ! It is the integral of the subtraction term used in the 3-jet contribution
  // ! See constants_calculation.nb for the analytic calculation
  // ! delta_eta is the total width of the slice
  // subroutine integrated_counterterm(delta_eta, counter_term)
  //   real(dp), intent(in)  :: delta_eta
  //   real(dp), intent(out) :: counter_term
  //   !--------------------------------
  //   real(dp) :: c
  //   ! The following declarations are necessary for hplog
  //   integer, parameter :: n1=-1
  //   integer, parameter :: n2= 1
  //   integer, parameter :: nw= 2
  //   real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
  //        & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
  //   real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
  //        & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
  //   complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
  //        & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)
  //   complex(dp) :: Li2OneMinusC2, Li2OnePlusC2
    
  //   c = (exp(delta_eta) - one) / (exp(delta_eta) + one)

  //   call hplog((one-c)/two,nw,Hc1,Hc2,Hc3,Hc4, &
  //        &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
  //   Li2OneMinusC2 = Hc2(0,1)

  //   call hplog((one+c)/two,nw,Hc1,Hc2,Hc3,Hc4, &
  //        &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
  //   Li2OnePlusC2 = Hc2(0,1)

  //   counter_term = (4*(57 - 7*Pi**2) - (3*(64*c - 64*c**3 - 192*c*log(2._dp) - 64*c**3*log(2._dp) - &
  //        &         8*log(1 - c) + 96*c*log(1 - c) - 144*c**2*log(1 - c) + &
  //        &         32*c**3*log(1 - c) + 24*c**4*log(1 - c) + 40*log(2._dp)*log(1 - c) - &
  //        &         80*c**2*log(2._dp)*log(1 - c) + 40*c**4*log(2._dp)*log(1 - c) - &
  //        &         3*log(256._dp)*log(1 - c) + 6*c**2*log(256._dp)*log(1 - c) - &
  //        &         3*c**4*log(256._dp)*log(1 - c) - 8*log(1 - c)**2 + 16*c**2*log(1 - c)**2 - &
  //        &         8*c**4*log(1 - c)**2 + 8*log(1 + c) + 96*c*log(1 + c) + &
  //        &         144*c**2*log(1 + c) + 32*c**3*log(1 + c) - 24*c**4*log(1 + c) - &
  //        &         40*log(2._dp)*log(1 + c) + 80*c**2*log(2._dp)*log(1 + c) - &
  //        &         40*c**4*log(2._dp)*log(1 + c) + 3*log(256._dp)*log(1 + c) - &
  //        &         6*c**2*log(256._dp)*log(1 + c) + 3*c**4*log(256._dp)*log(1 + c) + &
  //        &         8*log(1 + c)**2 - 16*c**2*log(1 + c)**2 + 8*c**4*log(1 + c)**2 + &
  //        &         16*(-1 + c**2)**2*Li2OneMinusC2 - &
  //        &         16*(-1 + c**2)**2*Li2OnePlusC2))/(-1 + c**2)**2)/24._dp
  //   return    
  // end subroutine integrated_counterterm
  
inline double Mgg_strongly_ordered_large_Nc(const Momentum& p, const Momentum& pbar,
					    const Momentum& k1,const Momentum& k2) {
  double W12 = dot_product(p,pbar)/(dot_product(p,k1)*dot_product(k1,k2)*dot_product(k2,pbar));
  double W21 = dot_product(p,pbar)/(dot_product(p,k2)*dot_product(k2,k1)*dot_product(k1,pbar));
  // double W1 = dot_product(p,pbar)/(dot_product(p,k1)*dot_product(k1,pbar));
  // double W2 = dot_product(p,pbar)/(dot_product(p,k2)*dot_product(k2,pbar));
  double S = W12 + W21; // cancels against CF^2: (- W1*W2);
  double m = CA*(2*S);
  return m;
}

inline double Mgg_full_large_Nc(const Momentum& p, const Momentum& pbar,
				const Momentum& k1,const Momentum& k2) {
  double W12 = dot_product(p,pbar)/(dot_product(p,k1)*dot_product(k1,k2)*dot_product(k2,pbar));
  double W21 = dot_product(p,pbar)/(dot_product(p,k2)*dot_product(k2,k1)*dot_product(k1,pbar));
  // double W1 = dot_product(p,pbar)/(dot_product(p,k1)*dot_product(k1,pbar));
  // double W2 = dot_product(p,pbar)/(dot_product(p,k2)*dot_product(k2,pbar));
  double S = W12 + W21; // cancels against CF^2: (- W1*W2);
  double R = (dot_product(p,k1)*dot_product(k2,pbar)+dot_product(p,k2)*dot_product(k1,pbar))
    /(dot_product(p,k1+k2)*dot_product(k1+k2,pbar));
  double J2 = (dot_product(p,k1)*dot_product(k2,pbar)-dot_product(p,k2)*dot_product(k1,pbar))
    /(dot_product(p,k1+k2)*dot_product(k1+k2,pbar)*dot_product(k1,k2));
  J2 *= 2*J2;
  double Hg = -S*R+J2-4*dot_product(p,pbar)/(dot_product(p,k1+k2)*dot_product(k1+k2,pbar)*dot_product(k1,k2));
  
  double m = CA*(2*S + Hg);
  return m;
}


inline double single_emsn_antenna(const Momentum& a, const Momentum& b, const Momentum& c) {
  return 2*dot_product(a,c)/(dot_product(a,b)*dot_product(b,c));
}

inline double double_emsn_antenna_strongly_ordered(const Momentum& a, const Momentum& b,
						   const Momentum& c, const Momentum& d) {
  // A(a,b,c,d) (c<<b) = A(a,b,d)A(b,c,d) = [2(a d)/((a b)(b d))] * [2(b d)/((b c) (c d))
  //                   = 4 (a d) / ((a b) (b c) (c d))
  return 4*dot_product(a,d)/(dot_product(a,b)*dot_product(b,c)*dot_product(c,d));
}

inline double double_emsn_antenna(const Momentum& a, const Momentum& b,
				  const Momentum& c, const Momentum& d) {
  // A(a,b,c,d) = (D-2)/(b c)^2 * (1 - (a b)/(a b+c) - (c d)/(b+c d))^2
  //             + 2 (a d)^2/((a b)(c d)(a b+c)(b+c d))
  //             + 2 (a d)/(b c) * (1/((a b)(b+c d)) + 1/((a b)(c d))
  //                       + 1/((a b+c)(c d)) - 4/((a b+c)(b+c d)))
  double ab=dot_product(a,b);
  double ac=dot_product(a,c);
  double ad=dot_product(a,d);
  double bc=dot_product(b,c);
  double bd=dot_product(b,d);
  double cd=dot_product(c,d);
  double res = (1 - ab/(ab+ac) - cd/(bd+cd));
  res*=res*2/(bc*bc);
  res+=2*ad*ad/(ab*cd*(ab+ac)*(bd+cd));
  res+=(2*ad/bc)*(1/(ab*(bd+cd)) + 1/(ab*cd) + 1/((ab+ac)*cd) - 4/((ab+ac)*(bd+cd)));
  return res;
}

inline double double_soft(double zeta, double cdphi, double deta) {
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  return (-2*(3 + 2*cdphi*cdphi)*zeta*(1 + zeta2) + 4*cdphi*(1 + 4*zeta2 + zeta4) +
	  cosh(deta)*(-4*(-1 + 2*cdphi*cdphi)*zeta2 - (4 + zeta2)*(1 + 4*zeta2) + 
		      16*cdphi*(zeta + zeta2*zeta) + 4*zeta*(1 - 2*cdphi*zeta + zeta2)*sinh(deta)) + 
	  zeta*((-6 + 8*cdphi*zeta - 6*zeta2)*cosh(2*deta) - 4*cdphi*(1 + zeta2)*sinh(deta) + 
		    zeta*(-3*cosh(3*deta) + 5*sinh(deta) + sinh(3*deta))))/ 
    (4.*(cdphi - cosh(deta))*(1 + zeta2 + 2*zeta*cosh(deta))*(1 + zeta2 + 2*zeta*cosh(deta)));
}
#endif //__MATRIX_ELEMENTS_HH__
