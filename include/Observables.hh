#ifndef __OBSERVABLES_HH__
#define __OBSERVABLES_HH__
#include <cmath>
#include <string>
#include "Parameters.hh"
#include "Momentum.hh"
#include "SimpleHist.hh"
#include <math.h>
                                                                   
//----------------------------------------------------------------------
/// \class Observable
/// abstract observable base class
class Observable {
public:
  /// constructor setting up the histograms
  Observable(double p, int nbin = 100, double maxlnkt = 10.0, double maxt = EVOLCUT)
    : p_(p),  dSdlambda_(0.0, 0.5, 10), dSdlnET_(0.0, maxlnkt, nbin), dSdt_(0.0, maxt, nbin) {}
  
  /// description
  virtual std::string description() const = 0;

  /// add entries to the histograms if emission is valid, otherwise return false
  bool add_entries_in_region(const Momentum & emsn, double t, double lnkt,
			     double weight = 1.0,
			     const Momentum* thrust_axis = 0) {
    if (this->in_region(emsn, thrust_axis)) {
      dSdt_.add_entry(t, weight);
      // dSdlnkt_.add_entry(lnkt, weight);
      // dSdlnE_. add_entry(-log(emsn.E()), weight);

      double ET2;
      if (SL_OBSERVABLE) {
        ET2 = exp(-2.*lnkt); // LL approximation for the observable
      } else {
        if (thrust_axis) {
	  //double costh = dot3(emsn, *thrust_axis)/emsn.E();
	  //ET2 = emsn.E() * emsn.E() * (1.0 - costh*costh);
	  ET2 = cross(emsn, *thrust_axis, true).E();
	  ET2 *= ET2;
        } else ET2 = emsn.px()*emsn.px() + emsn.py()*emsn.py();
      }
      dSdlnET_. add_entry(-Ltilde(sqrt(ET2)), weight);
      dSdlambda_. add_entry(-as*Ltilde(sqrt(ET2)), weight);
      // consistency test: calculate the observable (soft) correction as 
      // as*dS/dlambda*ln(1/sqrt(f2)) 
      // Comment call to dSdlambda_.add_entry above and increase the binning of dSdlambda_
      // This gives the pure NLL correction
      //if (SL_OBSERVABLE) {
      //  Momentum p1(0., 0., 1., 1.);
      //  Momentum p2(0., 0., -1., 1.);
      //  double f2= 2.*dot_product(emsn, p1)*dot_product(emsn, p2)/dot_product(p1, p2)
      //          /exp(-2.*lnkt);
      //  dSdlambda_. add_entry(-as*Ltilde(sqrt(ET2)), - weight * log(f2) * 0.5);
      //  // compare to expression for one-loop soft function by Becher et al.
      //  //dSdlambda_. add_entry(-as*Ltilde(sqrt(ET2)), - weight * log(f2/4./pow(sin(emsn.stored_phi_dip()),2)) * 0.5);
      //}
      return true;
    }
    return false;
  }

  // /// add entries to the histograms 
  // bool add_entries(const Momentum & emsn, double t, double lnkt,
  // 		   double weight = 1.0, Momentum* thrust_axis = 0) {
  //   dSdt_.   add_entry(t, weight);
  //   dSdlnkt_.add_entry(lnkt, weight);
  //   dSdlnE_. add_entry(log(emsn.E()), weight);
  //   return true;
  // }

  /// write the histograms to file or cout with proper normalisation
  void write(int nev, std::ostream * ostr = (&std::cout)) {
    *ostr << "# dSdt (differential)" << std::endl;
    output_differential(dSdt_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdt (cumulative)" << std::endl;
    output_inverse_cumulative(dSdt_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdlnET (differential)" << std::endl;
    output_differential(dSdlnET_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdlnET (cumulative)" << std::endl;
    output_inverse_cumulative(dSdlnET_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdlambda (differential)" << std::endl;
    output_differential(dSdlambda_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdlambda (cumulative)" << std::endl;
    output_inverse_cumulative(dSdlambda_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
  }

  /// returns true if emission is within region of interest
  virtual bool in_region(const Momentum & emsn, const Momentum* thrust_axis) const = 0;

  /// parameter of the observable
  virtual double parameter() const = 0;
protected:
  double p_;
  // SimpleHist dSdt_, dSdlnkt_, dSdlnE_;
  SimpleHist dSdlambda_, dSdlnET_, dSdt_;
  double Ltilde(double v) {
    if (p_==-1) return log(v);
    return -1/p_*log(1/pow(v, p_) + 1);
  }
};

//----------------------------------------------------------------------
/// \class Slice
/// contains the slice observable
class Slice : public Observable {
public:
  /// constructor
  Slice(double delta_rap, double p, int nbin = 100, double maxlnkt = lnktmax, double maxt = EVOLCUT)
    : Observable(p, nbin, maxlnkt, maxt), delta_rap_(delta_rap) {}

  /// description
  virtual std::string description() const {
    return "Slice with delta rap = " + std::to_string(delta_rap_);
  }

  /// return true if emission is in the slice
  virtual bool in_region(const Momentum & emsn, const Momentum* thrust_axis) const {
    if (std::abs(emsn.rap(thrust_axis)) < delta_rap_/2.0) return true;
    return false;
  }

  /// return delta rap
  virtual double parameter() const {return delta_rap_;}
  
private:
  double delta_rap_;
};
                                         
#endif // __OBSERVABLES_HH__
