#include "Momentum.hh"
#include <cassert>


//----------------------------------------------------------------------
/// rapidity of the jet
double Momentum::rap(const Momentum* thrust_axis) const  {
  double rap, dot;
  if (thrust_axis) dot = dot3(*this, *thrust_axis);
  else             dot = pz_;
  if (E_ == std::abs(dot)) {
    if (dot >= 0.0) {
      rap = maxrap_;
    } else {
      rap = -maxrap_;
    }
  } else {
    rap = 0.5*log((E_ + dot)/(E_ - dot));
  }
  return rap;
}

//----------------------------------------------------------------------
/// phi of the jet
double Momentum::phi() const {
  double phi;
  if (px_*px_ + py_*py_ == 0.0) {
    phi = 0.0;
  } else {
    phi = atan2(py_,px_);
  }
  if (phi < 0.0) {phi += 2.0*M_PI;}
  if (phi >= 2.0*M_PI) {phi -= 2.0*M_PI;}
  return phi;
}

//----------------------------------------------------------------------
/// rotate jet according to theta and phi
void Momentum::rotate(double thetaIn, double phiIn) {

  double cthe = cos(thetaIn);
  double sthe = sin(thetaIn);
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  /* rotation:
    cthe*cphi     cthe*sphi     -sthe
      -sphi         cphi          0
    sthe*cphi     sthe*sphi      cthe
  */
  double tmpx =  cthe * cphi * px_ + cthe * sphi * py_ - sthe * pz_;
  double tmpy =  -sphi *       px_ + cphi        * py_;
  double tmpz =  sthe * cphi * px_ + sthe * sphi * py_ + cthe * pz_;
  px_         = tmpx;
  py_         = tmpy;
  pz_         = tmpz;
}

//----------------------------------------------------------------------
/// perform inverse rotation of jet
void Momentum::unrotate(double thetaIn, double phiIn) {

  double cthe = cos(thetaIn);
  double sthe = sin(thetaIn);
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  /* inverse rotation:
    cthe*cphi     -sphi      sthe*cphi
    cthe*sphi      cphi      sthe*sphi
      -sthe          0         cthe
  */
  double tmpx =  cthe * cphi * px_ -    sphi * py_ + sthe * cphi * pz_;
  double tmpy =  cthe * sphi * px_ +    cphi * py_ + sthe * sphi * pz_;
  double tmpz = -sthe *        px_ +                cthe *        pz_;
  px_         = tmpx;
  py_         = tmpy;
  pz_         = tmpz;
}

//----------------------------------------------------------------------
/// perform rotation of jet according to cos and sin of angles
void Momentum::rotate2(double cthe, double sthe, double cphi, double sphi) {
  /* rotation:
    cthe*cphi     -sphi      sthe*cphi
    cthe*sphi      cphi      sthe*cphi
      -sthe          0         cthe
  */
  // double tmpx =  cthe * cphi * px_ -    sphi * py_ + sthe * cphi * pz_;
  // double tmpy =  cthe * sphi * px_ +    cphi * py_ + sthe * sphi * pz_;
  // double tmpz = -sthe *        px_ +                cthe *        pz_;
  /* inverse rotation:
    cthe*cphi     cthe*sphi     -sthe
      -sphi         cphi          0
    sthe*cphi     sthe*sphi      cthe
  */
  double tmpx =  cthe * px_ - sthe * sphi * py_ - sthe * cphi * pz_;
  double tmpy =  -sphi *       pz_ + cphi        * py_;
  double tmpz =  sthe * px_ + cthe * sphi * py_ + cthe * cphi * pz_;
  px_         = tmpx;
  py_         = tmpy;
  pz_         = tmpz;
}
//----------------------------------------------------------------------
/// multiply the jet's momentum by the coefficient
Momentum & Momentum::operator*=(double coeff) {
  px_ *= coeff;
  py_ *= coeff;
  pz_ *= coeff;
  E_  *= coeff;
  return *this;
}

//----------------------------------------------------------------------
/// transform this jet (given in the rest frame of prest) into a jet
/// in the lab frame 
//
// NB: code adapted from that in herwig f77 (checked how it worked
// long ago)
Momentum & Momentum::boost(const Momentum & prest) {
  
  if (prest.px() == 0.0 && prest.py() == 0.0 && prest.pz() == 0.0) 
    return *this;

  double m_local = prest.m();
  assert(m_local != 0);

  double pf4  = (  px()*prest.px() + py()*prest.py()
                 + pz()*prest.pz() + E()*prest.E() )/m_local;
  double fn   = (pf4 + E()) / (prest.E() + m_local);
  px_ +=  fn*prest.px();
  py_ +=  fn*prest.py();
  pz_ +=  fn*prest.pz();
  E_ = pf4;

  return *this;
}


//----------------------------------------------------------------------
/// transform this jet (given in lab) into a jet in the rest
/// frame of prest  
//
// NB: code adapted from that in herwig f77 (checked how it worked
// long ago)
Momentum & Momentum::unboost(const Momentum & prest) {
  
  if (prest.px() == 0.0 && prest.py() == 0.0 && prest.pz() == 0.0) 
    return *this;

  double m_local = prest.m();
  //assert(m_local != 0);

  double pf4  = ( -px()*prest.px() - py()*prest.py()
                 - pz()*prest.pz() + E()*prest.E() )/m_local;
  double fn   = (pf4 + E()) / (prest.E() + m_local);
  px_ -=  fn*prest.px();
  py_ -=  fn*prest.py();
  pz_ -=  fn*prest.pz();
  E_ = pf4;

  return *this;
}

//----------------------------------------------------------------------
// return "sum" of two pseudojets
Momentum operator+ (const Momentum & jet1, const Momentum & jet2) {
  //return Momentum(jet1.four_mom()+jet2.four_mom());
  return Momentum(jet1.px()+jet2.px(),
		  jet1.py()+jet2.py(),
		  jet1.pz()+jet2.pz(),
		  jet1.E() +jet2.E()  );
} 

//----------------------------------------------------------------------
// return difference of two pseudojets
Momentum operator- (const Momentum & jet1, const Momentum & jet2) {
  //return Momentum(jet1.four_mom()-jet2.four_mom());
  return Momentum(jet1.px()-jet2.px(),
		  jet1.py()-jet2.py(),
		  jet1.pz()-jet2.pz(),
		  jet1.E() -jet2.E()  );
}

//----------------------------------------------------------------------
// return the product, coeff * jet
Momentum operator* (double coeff, const Momentum & jet) {
  Momentum coeff_times_jet(jet);
  coeff_times_jet *= coeff;
  return coeff_times_jet;
}

//----------------------------------------------------------------------
/// Returns the 3-vector cross-product of p1 and p2. If lightlike is false
/// then the energy component is zero; if it's true the the energy component
/// is arranged so that the vector is lighlike
Momentum cross(const Momentum & p1, const Momentum & p2, bool lightlike) {
  double px = p1.py() * p2.pz() - p2.py() * p1.pz();
  double py = p1.pz() * p2.px() - p2.pz() * p1.px();
  double pz = p1.px() * p2.py() - p2.px() * p1.py();

  double E;
  if (lightlike) {
    E = sqrt(px*px + py*py + pz*pz);
  } else {
    E = 0.0;
  }
  return Momentum(px, py, pz, E);
}

//----------------------------------------------------------------------
/// Returns two 4-vectors, each with square = -1, that have a zero
/// (4-vector) dot-product with dir1 and dir2.
///
/// When dir1 and dir2 form a plane, then perp1 will be out of that
/// plane and perp2 will be in the plane.
void twoPerp(const Momentum & dir1, const Momentum & dir2,
            Momentum & perp1, Momentum & perp2) {

  // First get a 3-vector that is perpendicular to both dir1 and dir2.
  perp1 = cross(dir1,dir2);

  // for now test for exact zero -- later we will have to be more sophisticated...
  if (perp1.px()*perp1.px()+perp1.py()*perp1.py() + perp1.pz()*perp1.pz() == 0.0) {
    if (abs(dir1.px()) < abs(dir1.py()) && abs(dir1.px()) < abs(dir1.pz())) {
      // x is smallest direction, so use that as a starting point
      // for a cross product
      perp2 = cross(dir1, Momentum(1,0,0,0));
    } else if (abs(dir1.py()) < abs(dir1.pz())) {
      perp2 = cross(dir1, Momentum(0,1,0,0));
    } else {
      perp2 = cross(dir1, Momentum(0,0,1,0));
    }
    perp1 = cross(dir1,perp2);
  } else {
    // requirements for perp2:
    // - 3-vector should be perpendicular to perp1
    // - 3-vector dot-product with dir1 and dir2 should be zero
    // - norm = -1
    //-----
    // go to centre-of-mass frame of the dipole
    Momentum dir12 = dir1 + dir2;
    Momentum dir1rest = dir1; dir1rest.unboost(dir12);
    // get something perpendicular to perp1 and either of the dipole directions
    // (since perp1 is perpendicular to both, it doesn't change with the boost)
    perp2 = cross(dir1rest,perp1);
    
    // Momentum dir2rest = dir2; dir2rest.unboost(dir12);
    // cout << dir1rest << "  xxxd1" << endl;
    // cout << dir2rest << "  xxxd2" << endl;
    // cout << perp2 << " xxxp2" << endl;
    // boost back
    perp2.boost(dir12);
  }

  // Momentum dir12 = dir1 + dir2;
  // Momentum dir1rest = dir1; dir1rest.unboost(dir12);
  // if (dir12.modp2() < 1e-16 * dir12.m2) {
  //   // we are already in the centre of mass
  // }
  
  // arrange norm -1 for perp1 and perp2
  perp1 *= 1.0/sqrt((-perp1.m2()));
  perp2 *= 1.0/sqrt((-perp2.m2()));

}
