#ifndef __DIPOLE_HH__
#define __DIPOLE_HH__

#include "Momentum.hh"
#include "Parameters.hh"

//----------------------------------------------------------------------
/// \class DipoleEnd
/// contains information about a dipole end, such as its momentum and
/// whether it is Born like.
class DipoleEnd{
public:
  /// constructor taking a four-momentum and boolean (true for born
  /// dipole ends)
  DipoleEnd(Momentum mom, bool isBorn = false)
    : momentum_(mom), isBorn_(isBorn) {}

  /// copy constructor
  DipoleEnd(const DipoleEnd& de) : momentum_(de.momentum()), isBorn_(de.isBorn()) {}

  /// return true if dipole end is a born particle
  bool isBorn() const {return isBorn_;}

  /// return the momentum of the dipole end
  Momentum & momentum() {return momentum_;}

  /// return the (unmodifiable) momentum of the dipole end
  const Momentum & momentum() const {return momentum_;}
  
private:
  /// Four momentum
  Momentum momentum_;
  /// bool keeping track of born particles
  bool isBorn_;
};


//----------------------------------------------------------------------
/// \class Dipole
/// contains information about a dipole
class Dipole {
public:
  /// constructor taking left and right dipole ends as well as indices
  /// of the left and right neighbours
  Dipole(DipoleEnd left, DipoleEnd right, int left_neighbour,
	 int right_neighbour)
    : left_(left), right_(right),
      left_neighbour_(left_neighbour),
      right_neighbour_(right_neighbour) {
    reset();
  }

  /// copy constructor
  Dipole(const Dipole& dip) : left_(dip.left()), right_(dip.right()),
			      left_neighbour_(dip.left_neighbour()),
			      right_neighbour_(dip.right_neighbour()),
			      m2_(dip.m2()), rap_left_(dip.rap_left()),
			      rap_right_(dip.rap_right()), delta_rap_(dip.delta_rap()) {}

  //----------------------------------------------------------------------
  /// radiate a gluon off the dipole using two uniform random
  /// numbers ran1,ran2 \in [0,1].
  /// NOTE: This does NOT modify the dipole itself, and only returns
  /// the direction of the emitted parton.
  Momentum radiate(double lnkt, double ran1, double ran2) {
    double phi = ran1*2.0*M_PI;
    double eta = ran2*delta_rap_ + rap_right_;
    Momentum cm2orig = left_.momentum() + right_.momentum();
    Momentum lboost = left_.momentum();
    lboost.unboost(cm2orig);

    // produce a gluon relative to a fictional emitter along z axis
    double exp_eta = std::exp(eta);
    // std::cout<<" "  << exp_eta << std::endl;
    // need to add kt here for NLL
    double alpha = std::exp(-lnkt) * exp_eta;
    double beta  = std::exp(-lnkt) / exp_eta;
    double gluon_px = std::exp(-lnkt) * sin(phi);
    double gluon_py = std::exp(-lnkt) * cos(phi);
    double gluon_pz = 0.5 * (alpha - beta);
    double gluon_E  = 0.5 * (alpha + beta);

    // rotate gluon so that it is as if emitted from particle axis
    double mdls2 = sqrt(lboost.py()*lboost.py() + lboost.pz()*lboost.pz());
    double mdls3 = sqrt(mdls2*mdls2 + lboost.px()*lboost.px());
    double ctht = mdls2/mdls3;
    double stht = lboost.px()/mdls3;
    double temp1 =  gluon_px * ctht + gluon_pz * stht;
    double temp2 = -gluon_px * stht + gluon_pz * ctht;
    gluon_px = temp1;
    gluon_pz = temp2;
    double cphi = lboost.pz()/mdls2;
    double sphi = lboost.py()/mdls2;
    double temp3 =  gluon_py * cphi + gluon_pz * sphi;
    double temp4 = -gluon_py * sphi + gluon_pz * cphi;
    gluon_py = temp3;
    gluon_pz = temp4;
    Momentum gluon(gluon_px, gluon_py, gluon_pz, gluon_E);
    // boost back to original frame
    gluon.boost(cm2orig);
    gluon.stored_E(gluon.E());
    gluon = (1.0/gluon.E())*gluon;
    return gluon;
  }

  /// return left dipole end
  const DipoleEnd & left()  const {return left_;}

  /// change left dipole end
  void left(const DipoleEnd & left) {
    left_ = left;
    reset();
  }

  /// return right dipole end
  const DipoleEnd & right() const {return right_;}

  /// change right dipole end
  void right(const DipoleEnd & right) {
    right_ = right;
    reset();
  }

  /// return index of left neighbour
  int left_neighbour() const {return left_neighbour_;}
  
  /// change index of left neighbour
  void left_neighbour(int i) {left_neighbour_ = i;}
  
  /// return index of right neighbour
  int right_neighbour() const {return right_neighbour_;}

  /// change index of right neighbour
  void right_neighbour(int i) {right_neighbour_ = i;}

  /// return squared invariant mass
  double m2() const {return m2_;}

  /// return the rapidity separation
  double delta_rap() const {return delta_rap_;}
  /// return the left rapidity
  double rap_left() const {return rap_left_;}
  /// return the right rapidity
  double rap_right() const {return rap_right_;}

private:

  //----------------------------------------------------------------------
  /// set or reset the rapidity separation of the dipole and its
  /// invariant mass
  void reset() {
    // We set the relative boundaries wrt to dipole ends for new emissions in the lab frame 
    rap_left_  =  std::min(rapmax_, log(left_. momentum().E()/ktsoft_));
    rap_right_ = -std::min(rapmax_, log(right_.momentum().E()/ktsoft_));
    Momentum sum = left_.momentum() + right_.momentum();
    // std::cout << std::setw(8) << rap_left_ << " <> "
    // 	      << std::setw(8) << rap_right_
    // 	     << " => " << delta_rap_ << " (dipole frame)" <<std::endl; 

    // then we boost to the dipole frame and compute the rapidity
    // range for new emissions in that frame
    rap_left_  = rap_left_  + log(0.5*sum.m()/left_.momentum().E());
    rap_right_ = rap_right_ - log(0.5*sum.m()/right_.momentum().E());
    delta_rap_ = std::max(rap_left_ - rap_right_, 0.0);
    m2_        = sum.m2();
    // std::cout << std::setw(8) << rap_left_ << " <> "
    // 	      << std::setw(8) << rap_right_
    // 	      << " => " << delta_rap_ <<std::endl; 
  }

  /// left and right dipole ends
  DipoleEnd left_, right_;
  /// index of left and right neighbours
  int left_neighbour_, right_neighbour_;
  /// values of invariant mass and rapidities
  double m2_, rap_left_, rap_right_, delta_rap_;
  /// maximum rapidity value
  static constexpr double rapmax_ = RAPMAX;
  /// maximum kt value
  static constexpr double ktsoft_ = KTSOFT;
};

#endif // __DIPOLE_HH__
