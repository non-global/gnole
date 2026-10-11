#ifndef __OBSERVABLES_HH__
#define __OBSERVABLES_HH__
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <memory>
#include <stdexcept>
#include "Parameters.hh"
#include "Momentum.hh"
#include "SimpleHist.hh"
#include "EventStatistics.hh"
#include <math.h>
                                                                   
//----------------------------------------------------------------------
/// \class Observable
/// abstract observable base class
class Observable {
public:
  /// constructor setting up the histograms
  Observable(double p, int nbin = 100, double maxlnkt = 10.0, double maxt = EVOLCUT)
    : p_(p),  dSdlambda_(0.0, 0.5, 10), dSdlnET_(0.0, maxlnkt, nbin), dSdt_(0.0, maxt, nbin),
      mc_errors_(false), event_score_(false) {}

  /// Enable event-wise errors. Set before running; toggling resets statistics.
  void enable_mc_errors(bool enabled) {
    mc_errors_ = enabled;
    configure_statistics();
  }

  /// Track net event-bin norm for adaptive sampling without printing errors.
  void enable_event_score(bool enabled) {
    if (enabled && correction_)
      throw std::invalid_argument("Observable correction does not support adaptive event scoring");
    event_score_ = enabled;
    configure_statistics();
  }

  /// Extract the linearised measurement correction on unchanged LL histories.
  /// Set before running. A mode change discards previous fills and statistics.
  void enable_obs_correction_only(bool enabled) {
    if (enabled) {
      validate_correction_mode();
      if (event_score_)
        throw std::invalid_argument("Observable correction does not support adaptive event scoring");
      correction_.reset(new CorrectionHistograms(dSdlnET_, dSdlambda_));
    } else {
      correction_.reset();
    }
    reset_histograms();
  }

  bool obs_correction_only() const { return static_cast<bool>(correction_); }

  void begin_event() {
    if (correction_) {
      if (mc_errors_) {
        correction_->l_statistics.begin_event();
        correction_->lambda_statistics.begin_event();
      }
      return;
    }
    if (mc_errors_) {
      lambda_statistics_.begin_event();
      time_statistics_.begin_event();
    }
    if (mc_errors_ || event_score_) et_statistics_.begin_event();
  }

  void end_event() {
    if (correction_) {
      if (mc_errors_) {
        correction_->l_statistics.end_event();
        correction_->lambda_statistics.end_event();
      }
      return;
    }
    if (mc_errors_) {
      lambda_statistics_.end_event();
      time_statistics_.end_event();
    }
    if (mc_errors_ || event_score_) et_statistics_.end_event();
  }

  double current_event_squared_norm() const {
    if (correction_)
      return mc_errors_ ? correction_->l_statistics.current_squared_norm() : 0.;
    return (mc_errors_ || event_score_) ? et_statistics_.current_squared_norm() : 0.;
  }
  double event_squared_norm() const { return current_event_squared_norm(); }
  double last_event_squared_norm() const {
    if (correction_)
      return mc_errors_ ? correction_->l_statistics.last_squared_norm() : 0.;
    return (mc_errors_ || event_score_) ? et_statistics_.last_squared_norm() : 0.;
  }

  /// Discard a pilot phase, including its histograms and error accumulators.
  void reset_histograms() {
    dSdlambda_.reset();
    dSdlnET_.reset();
    dSdt_.reset();
    if (correction_) {
      correction_->l.reset();
      correction_->lambda.reset();
    }
    configure_statistics();
  }
  
  /// description
  virtual std::string description() const = 0;

  /// add entries to the histograms if emission is valid, otherwise return false
  bool add_entries_in_region(const Momentum & emsn, double t, double lnkt,
			     double weight = 1.0,
			     const Momentum* thrust_axis = 0) {
    if (correction_) return add_correction_in_region(emsn, t, lnkt, weight, thrust_axis);
    // Keeping the SL observable in logarithmic form avoids losing a valid
    // flow entry when exp(-lnkt) underflows.
    if (PRECISE_KERNELS && SL_OBSERVABLE) {
      if (!std::isfinite(lnkt) || !std::isfinite(t) || !std::isfinite(weight))
        throw std::invalid_argument("Nonfinite observable coordinate or event weight");
      if (!in_region(emsn, thrust_axis)) return false;
      return add_transformed_entries(-Ltilde_from_log(-lnkt), t, weight);
    }
    double ET = transverse_energy(emsn, lnkt, thrust_axis);
    return add_entries(ET, t, weight);    
  }

  /// add entries to the histograms 
  bool add_entries(const double obs, double t, double weight = 1.0) {
    if (correction_)
      throw std::logic_error("Observable correction requires an emission momentum and its SL scale");
    if (!std::isfinite(obs) || obs < 0. || !std::isfinite(t) || !std::isfinite(weight))
      throw std::invalid_argument("Invalid observable value or event weight");
    if (obs == 0.) return false;
    return add_transformed_entries(-Ltilde(obs), t, weight);
  }

  /// compute observable
  double transverse_energy(const Momentum & emsn, double lnkt,
			     const Momentum* thrust_axis = 0) {
    if (PRECISE_KERNELS) {
      if (!in_region(emsn, thrust_axis)) return 0.;
      if (SL_OBSERVABLE) return std::exp(-lnkt);
      if (thrust_axis) return cross(emsn, *thrust_axis, true).E();
      return std::hypot(emsn.px(), emsn.py());
    }
    double ET2 = 0.;
    if (this->in_region(emsn, thrust_axis)) {
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
    }
    return sqrt(ET2);
  }
  
  /// write the histograms to file or cout with proper normalisation
  void write(int nev, std::ostream * ostr = (&std::cout)) {
    if (correction_) {
      if (nev <= 0) throw std::invalid_argument("Observable correction requires a positive event count");
      validate_correction_mode();
      if (mc_errors_)
        *ostr << "# Monte Carlo errors: event-wise standard errors of the mean; "
              << correction_->l_statistics.events() << " complete events" << std::endl;
      *ostr << "# Linearised LL measurement correction: c = log(v_SL/v_full); "
            << "tables are DeltaSigma, not densities to cumulatively sum" << std::endl;
      *ostr << "# observable correction (L)" << std::endl;
      write_correction(correction_->l, correction_->l_statistics, ostr, nev, 1.);
      *ostr << std::endl << std::endl;
      *ostr << "# observable correction (lambda)" << std::endl;
      write_correction(correction_->lambda, correction_->lambda_statistics, ostr, nev, as);
      *ostr << std::endl << std::endl;
      return;
    }
    if (mc_errors_)
      *ostr << "# Monte Carlo errors: event-wise standard errors of the mean; "
            << et_statistics_.events() << " complete events" << std::endl;
    *ostr << "# dSdt (differential)" << std::endl;
    write_differential(dSdt_, time_statistics_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdt (cumulative)" << std::endl;
    write_inverse_cumulative(dSdt_, time_statistics_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdlnET (differential)" << std::endl;
    write_differential(dSdlnET_, et_statistics_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdlnET (cumulative)" << std::endl;
    write_inverse_cumulative(dSdlnET_, et_statistics_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdlambda (differential)" << std::endl;
    write_differential(dSdlambda_, lambda_statistics_, ostr, 1.0/nev);
    *ostr << std::endl << std::endl;
    *ostr << "# dSdlambda (cumulative)" << std::endl;
    write_inverse_cumulative(dSdlambda_, lambda_statistics_, ostr, 1.0/nev);
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
  bool mc_errors_, event_score_;
  EventStatistics lambda_statistics_, et_statistics_, time_statistics_;

  struct CorrectionHistograms {
    SimpleHist l, lambda;
    EventStatistics l_statistics, lambda_statistics;
    CorrectionHistograms(const SimpleHist &ordinary_l, const SimpleHist &ordinary_lambda)
      : l(ordinary_l.min() - ordinary_l.binsize()/2.,
          ordinary_l.max() + ordinary_l.binsize()/2., ordinary_l.size() + 1),
        lambda(ordinary_lambda.min() - ordinary_lambda.binsize()/2.,
               ordinary_lambda.max() + ordinary_lambda.binsize()/2., ordinary_lambda.size() + 1) {}
  };
  // No correction grids or their event buffers are allocated in ordinary mode.
  std::unique_ptr<CorrectionHistograms> correction_;

  void validate_correction_mode() const {
    if (!SL_OBSERVABLE || p_ != -1.)
      throw std::invalid_argument("Observable correction requires the SL observable and p = -1");
    if (!std::isfinite(as) || as <= 0.)
      throw std::invalid_argument("Observable correction requires a finite positive coupling");
  }

  bool add_correction_in_region(const Momentum &emsn, double t, double lnkt,
                               double weight, const Momentum *thrust_axis) {
    validate_correction_mode();
    if (!std::isfinite(lnkt) || !std::isfinite(t) || !std::isfinite(weight) ||
        !std::isfinite(emsn.px()) || !std::isfinite(emsn.py()) ||
        !std::isfinite(emsn.pz()) || !std::isfinite(emsn.E()))
      throw std::invalid_argument("Nonfinite observable correction momentum, scale or weight");
    // This first implementation is an LL Born-dipole calculation. A rotated
    // H3 thrust axis would require a separate hard measurement prescription.
    if (thrust_axis && (thrust_axis->px() != 0. || thrust_axis->py() != 0. ||
                        std::abs(thrust_axis->pz()) != 1.))
      throw std::invalid_argument("Observable correction requires the LL Born z axis");
    if (!in_region(emsn, thrust_axis)) return false;
    const double transverse_scale = std::max(std::abs(emsn.px()), std::abs(emsn.py()));
    if (transverse_scale == 0.)
      throw std::invalid_argument("Observable correction requires nonzero physical transverse momentum");
    // log(hypot(px,py)) evaluated with a scaled norm remains finite even if
    // the norm itself overflows, and never constructs exp(-lnkt) or f^2.
    const double log_kt = std::log(transverse_scale) +
      std::log(std::hypot(emsn.px()/transverse_scale, emsn.py()/transverse_scale));
    const double c = -lnkt - log_kt;
    const double correction_weight = weight * c;
    const double lambda = as * lnkt;
    if (!std::isfinite(correction_weight) || !std::isfinite(lambda))
      throw std::overflow_error("Nonfinite observable correction coordinate or weight");
    const unsigned l_bin = correction_->l.bin(lnkt);
    const unsigned lambda_bin = correction_->lambda.bin(lambda);
    correction_->l.add_entry_at_bin(l_bin, lnkt, correction_weight);
    correction_->lambda.add_entry_at_bin(lambda_bin, lambda, correction_weight);
    if (mc_errors_) {
      correction_->l_statistics.add(l_bin, correction_weight);
      correction_->lambda_statistics.add(lambda_bin, correction_weight);
    }
    // A zero correction is still the first accepted in-slice emission and
    // must stop the LL shower; do not test correction_weight here.
    return true;
  }

  bool add_transformed_entries(double log_coordinate, double t, double weight) {
    if (!std::isfinite(log_coordinate) || !std::isfinite(t) || !std::isfinite(weight))
      throw std::invalid_argument("Nonfinite observable coordinate or event weight");
    const double lambda = as * log_coordinate;
    const unsigned time_bin = dSdt_.bin(t);
    const unsigned et_bin = dSdlnET_.bin(log_coordinate);
    const unsigned lambda_bin = dSdlambda_.bin(lambda);
    dSdt_.add_entry_at_bin(time_bin, t, weight);
    dSdlnET_.add_entry_at_bin(et_bin, log_coordinate, weight);
    dSdlambda_.add_entry_at_bin(lambda_bin, lambda, weight);
    if (mc_errors_) {
      time_statistics_.add(time_bin, weight);
      lambda_statistics_.add(lambda_bin, weight);
    }
    if (mc_errors_ || event_score_) et_statistics_.add(et_bin, weight);
    return true;
  }

  void configure_statistics() {
    if (correction_) {
      correction_->l_statistics.configure(mc_errors_ ? correction_->l.size() : 0, mc_errors_);
      correction_->lambda_statistics.configure(mc_errors_ ? correction_->lambda.size() : 0, mc_errors_);
      return;
    }
    lambda_statistics_.configure(mc_errors_ ? dSdlambda_.size() : 0, mc_errors_);
    time_statistics_.configure(mc_errors_ ? dSdt_.size() : 0, mc_errors_);
    et_statistics_.configure((mc_errors_ || event_score_) ? dSdlnET_.size() : 0,
                             mc_errors_);
  }

  void write_correction(const SimpleHist &hist, const EventStatistics &statistics,
                        std::ostream *ostr, int nev, double coefficient) const {
    *ostr << "# columns: lower center upper delta_cumulant";
    if (mc_errors_) *ostr << " standard_error";
    *ostr << std::endl;
    *ostr << "# centers align with ordinary cumulant edges; shifted bins are finite-bin averages"
          << std::endl;
    *ostr << "# underflow raw_weight_per_event " << hist.underflow()/nev;
    if (mc_errors_) *ostr << " standard_error_of_mean " << statistics.standard_error(hist.underflow_bin());
    *ostr << " (not DeltaSigma)" << std::endl;
    for (unsigned i = 0; i < hist.size(); ++i) {
      const double width = hist.binhi(i) - hist.binlo(i);
      *ostr << hist.binlo(i) << " " << hist.binmid(i) << " " << hist.binhi(i)
            << " " << coefficient * (hist[i]/nev) / width;
      if (mc_errors_) *ostr << " " << coefficient * statistics.standard_error(i) / width;
      *ostr << std::endl;
    }
    *ostr << "# overflow raw_weight_per_event " << hist.overflow()/nev;
    if (mc_errors_) *ostr << " standard_error_of_mean " << statistics.standard_error(hist.overflow_bin());
    *ostr << " (not DeltaSigma)" << std::endl;
  }

  void write_differential(const SimpleHist &hist, const EventStatistics &statistics,
                          std::ostream *ostr, double norm) const {
    if (!mc_errors_) { output_differential(hist, ostr, norm); return; }
    *ostr << "# columns: low midpoint high density standard_error" << std::endl;
    *ostr << "# under flow bin " << hist.underflow()
          << " standard_error_of_mean " << statistics.standard_error(hist.underflow_bin())
          << std::endl;
    for (unsigned i = 0; i < hist.size(); ++i) {
      const double width = hist.binhi(i) - hist.binlo(i);
      *ostr << hist.binlo(i) << " " << hist.binmid(i) << " " << hist.binhi(i)
            << " " << hist[i] * norm / width
            << " " << statistics.standard_error(i) / width << std::endl;
    }
    *ostr << "over flow bin " << hist.overflow()
          << " standard_error_of_mean " << statistics.standard_error(hist.overflow_bin())
          << std::endl;
  }

  void write_inverse_cumulative(const SimpleHist &hist, const EventStatistics &statistics,
                                std::ostream *ostr, double norm) const {
    if (!mc_errors_) { output_inverse_cumulative(hist, ostr, norm); return; }
    *ostr << "# columns: edge cumulative standard_error" << std::endl;
    double cumulative = 0.;
    double correction = 0.;
    hist.accumulate_bin_weight(hist.overflow_bin(), norm, cumulative, correction);
    *ostr << hist.binhi(hist.size() - 1) << " " << cumulative + correction << " "
          << statistics.inverse_cumulative_standard_error(0) << std::endl;
    unsigned step = 1;
    for (int i = hist.size() - 1; i >= 0; --i, ++step) {
      hist.accumulate_bin_weight(i, norm, cumulative, correction);
      *ostr << hist.binlo(i) << " " << cumulative + correction << " "
            << statistics.inverse_cumulative_standard_error(step) << std::endl;
    }
    hist.accumulate_bin_weight(hist.underflow_bin(), norm, cumulative, correction);
    *ostr << "# with_underflow " << cumulative + correction << " "
          << statistics.inverse_cumulative_standard_error(step) << std::endl;
  }
  double Ltilde(double v) {
    if (PRECISE_KERNELS) return Ltilde_from_log(std::log(v));
    if (p_==-1) return log(v);
    return -1/p_*log(1/pow(v, p_) + 1);
  }
  double Ltilde_from_log(double log_v) const {
    if (p_ == -1.) return log_v;
    const double argument = -p_ * log_v;
    // Factoring the leading logarithm also avoids overflowing p*log(v)
    // when p is large, while keeping a small softplus correction.
    if (argument > 0.) return log_v - std::log1p(std::exp(-argument)) / p_;
    return -std::log1p(std::exp(argument)) / p_;
  }
};

//----------------------------------------------------------------------
/// \class Slice
/// contains the slice observable
class Slice : public Observable {
public:
  /// constructor
  Slice(double delta_rap, double p, int nbin = 100, double maxlnkt = lnktmax, double maxt = EVOLCUT)
    : Observable(p, nbin, maxlnkt, maxt), delta_rap_(delta_rap),
      rapidity_threshold_(std::tanh(delta_rap / 2.)) {}

  /// description
  virtual std::string description() const {
    return "Slice with delta rap = " + std::to_string(delta_rap_);
  }

  /// return true if emission is in the slice
  virtual bool in_region(const Momentum & emsn, const Momentum* thrust_axis) const {
    if (PRECISE_KERNELS) {
      const double energy = emsn.E();
      const double projection = thrust_axis ? dot3(emsn, *thrust_axis) : emsn.pz();
      if (!std::isfinite(energy) || !std::isfinite(projection))
        throw std::invalid_argument("Nonfinite slice momentum or axis");
      const double absolute_projection = std::abs(projection);
      if (energy > 0. && absolute_projection < energy) {
        const double boundary = energy * rapidity_threshold_;
        const double tolerance = 8. * std::numeric_limits<double>::epsilon() * energy;
        // Preserve the original strict boundary convention near a rounding
        // ambiguity, including non-unit axes and exactly collinear momenta.
        if (std::abs(absolute_projection - boundary) > tolerance)
          return absolute_projection < boundary;
      }
    }
    if (std::abs(emsn.rap(thrust_axis)) < delta_rap_/2.0) return true;
    return false;
  }

  /// return delta rap
  virtual double parameter() const {return delta_rap_;}
  
private:
  double delta_rap_;
  double rapidity_threshold_;
};
                                         
#endif // __OBSERVABLES_HH__
