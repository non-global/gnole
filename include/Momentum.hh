#ifndef __MOMENTUM_HH__
#define __MOMENTUM_HH__

#include <iostream>
#include <iomanip>
#include <cmath>

// #include "fastjet/PseudoJet.hh"
// typedef fastjet::PseudoJet Momentum;

//----------------------------------------------------------------------
/// \class Momentum
/// contains a four vector, adapted from FastJet's PseudoJet
class Momentum {
public:
  Momentum() : px_(0),py_(0),pz_(0),E_(0), stored_E_(0.0){}
  /// constructor
  Momentum(double px, double py, double pz, double E)
    : px_(px),py_(py),pz_(pz),E_(E), stored_E_(E) {}

  /// mass squared
  double m2() const {return E_*E_ - px_*px_ - py_*py_ - pz_*pz_;}

  /// mass
  double m() const {
    double mm = m2();
    return mm < 0.0 ? -std::sqrt(-mm) : std::sqrt(mm);
  }

  // squared norm of three-vector
  double modpsq() const {
    return px_*px_ + py_*py_ + pz_*pz_;
  }

  // norm of three vector
  double modp() const {
    return sqrt(modpsq());
  }

  /// rapidity
  double rap(const Momentum* thrust_axis = 0) const;

  /// phi
  double phi() const;

  /// px
  double px() const {return px_;}

  /// py
  double py() const {return py_;}

  /// pz
  double pz() const {return pz_;}

  /// energy
  double E () const {return E_;}

  /// full energy without normalisation
  double stored_E() const {return stored_E_;}

  /// modify full energy without normalisation
  void stored_E(double E) {stored_E_=E;}

  /// theta
  double theta() const {return atan2(sqrt(px_*px_ + py_*py_), pz_);}

  /// rotate
  void rotate(double theta, double phi);
  /// rotate
  void unrotate(double theta, double phi);
  /// rotate
  void rotate2(double ctht, double stht, double cphi, double sphi);
  
  /// multiplication by double
  Momentum & operator*=(double);
  
  /// transform this jet (given in the rest frame of prest) into a jet
  /// in the lab frame
  Momentum & boost(const Momentum & prest);

  /// transform this jet (given in lab) into a jet in the rest
  /// frame of prest
  Momentum & unboost(const Momentum & prest);

private:
  /// four momentum components
  double px_, py_ ,pz_, E_, stored_E_;

  /// maximum rapidity value
  static constexpr double maxrap_ = 1e5;
};

Momentum operator*(double, const Momentum &);
Momentum operator+(const Momentum &, const Momentum &);
Momentum operator-(const Momentum &, const Momentum &);

/// returns the 4-vector dot product of a and b
inline double dot_product(const Momentum & a, const Momentum & b) {
  return a.E()*b.E() - a.px()*b.px() - a.py()*b.py() - a.pz()*b.pz();
}

/// returns the dot product of a and b vector
inline double dot3(const Momentum & a, const Momentum & b) {
  return a.px()*b.px() + a.py()*b.py() + a.pz()*b.pz();
}

//----------------------------------------------------------------------
/// output a momentum
inline std::ostream & operator<<(std::ostream & ostr, const Momentum & p) {
  ostr << std::setw(13) << p.px() << " "
       << std::setw(13) << p.py() << " "
       << std::setw(13) << p.pz() << " "
       << std::setw(13) << p.E();
  return ostr;
}

//----------------------------------------------------------------------
/// Returns the 3-vector cross-product of p1 and p2. If lightlike is false
/// then the energy component is zero; if it's true the the energy component
/// is arranged so that the vector is lighlike
Momentum cross(const Momentum & p1, const Momentum & p2, bool lightlike = false);

//----------------------------------------------------------------------
/// Returns two 4-vectors, each with square = -1, that have a zero
/// (4-vector) dot-product with dir1 and dir2.
///
/// When dir1 and dir2 form a plane, then perp1 will be out of that
/// plane and perp2 will be in the plane.
void twoPerp(const Momentum & dir1, const Momentum & dir2,
             Momentum & perp1, Momentum & perp2);

#endif //  __MOMENTUM_HH__
