#ifndef __EVENTSTATISTICS_HH__
#define __EVENTSTATISTICS_HH__

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
#include "Parameters.hh"

/// Statistics for complete Monte Carlo events, including signed partners.
/// Histogram bins are [0,n), followed by underflow and overflow. Each event
/// contributes the SUM of all its fills to each bin before any square is taken.
class EventStatistics {
public:
  EventStatistics() : nbins_(0), events_(0), moments_enabled_(false),
                      active_(false), precise_(PRECISE_KERNELS),
                      last_squared_norm_(0.) {}

  void configure(std::size_t nbins, bool moments_enabled) {
    nbins_ = nbins;
    moments_enabled_ = moments_enabled;
    precise_ = PRECISE_KERNELS;
    weights_.assign(nbins + 2, 0.);
    weights_correction_.assign(nbins + 2, 0.);
    touched_.assign(nbins + 2, false);
    touched_bins_.clear();
    mean_.assign(moments_enabled ? nbins + 2 : 0, 0.);
    m2_.assign(moments_enabled ? nbins + 2 : 0, 0.);
    mean_correction_.assign(moments_enabled ? nbins + 2 : 0, 0.);
    m2_correction_.assign(moments_enabled ? nbins + 2 : 0, 0.);
    cumulative_mean_.assign(moments_enabled ? nbins + 2 : 0, 0.);
    cumulative_m2_.assign(moments_enabled ? nbins + 2 : 0, 0.);
    cumulative_mean_correction_.assign(moments_enabled ? nbins + 2 : 0, 0.);
    cumulative_m2_correction_.assign(moments_enabled ? nbins + 2 : 0, 0.);
    events_ = 0;
    active_ = false;
    last_squared_norm_ = 0.;
  }

  void reset() { configure(nbins_, moments_enabled_); }

  void begin_event() {
    assert(!active_);
    for (std::size_t bin : touched_bins_) {
      weights_[bin] = 0.;
      weights_correction_[bin] = 0.;
      touched_[bin] = false;
    }
    touched_bins_.clear();
    active_ = true;
  }

  void add(std::size_t bin, double weight) {
    assert(active_ && bin < weights_.size());
    if (!std::isfinite(weight))
      throw std::invalid_argument("Nonfinite complete-event statistical weight");
    if (!touched_[bin]) {
      touched_[bin] = true;
      touched_bins_.push_back(bin);
    }
    accumulate(weight, weights_[bin], weights_correction_[bin]);
  }

  double current_squared_norm() const {
    long double result = 0.;
    long double correction = 0.;
    for (std::size_t bin : touched_bins_) {
      const long double weight = value(bin);
      accumulate(weight * weight, result, correction);
    }
    result += correction;
    if (!std::isfinite(result) || result > std::numeric_limits<double>::max())
      throw std::overflow_error("Nonfinite complete-event squared norm");
    return static_cast<double>(result);
  }

  void end_event() {
    assert(active_);
    active_ = false;
    ++events_;
    last_squared_norm_ = current_squared_norm();
    if (!moments_enabled_) return;
    for (std::size_t bin = 0; bin < weights_.size(); ++bin)
      update(value(bin), mean_[bin], mean_correction_[bin],
             m2_[bin], m2_correction_[bin]);

    // Match inverse-cumulative output: first overflow alone at the upper
    // edge, then descending physical bins, finally the total with underflow.
    long double cumulative = weights_[nbins_ + 1], correction = weights_correction_[nbins_ + 1];
    update(cumulative + correction, cumulative_mean_[0], cumulative_mean_correction_[0],
           cumulative_m2_[0], cumulative_m2_correction_[0]);
    for (std::size_t step = 1; step <= nbins_; ++step) {
      accumulate(weights_[nbins_ - step], cumulative, correction);
      if (precise_) accumulate(weights_correction_[nbins_ - step], cumulative, correction);
      update(cumulative + correction, cumulative_mean_[step],
             cumulative_mean_correction_[step], cumulative_m2_[step],
             cumulative_m2_correction_[step]);
    }
    accumulate(weights_[nbins_], cumulative, correction);
    if (precise_) accumulate(weights_correction_[nbins_], cumulative, correction);
    update(cumulative + correction, cumulative_mean_[nbins_ + 1],
           cumulative_mean_correction_[nbins_ + 1],
           cumulative_m2_[nbins_ + 1], cumulative_m2_correction_[nbins_ + 1]);
  }

  std::size_t events() const { return events_; }
  double last_squared_norm() const { return last_squared_norm_; }
  double standard_error(std::size_t bin) const {
    assert(moments_enabled_ && bin < m2_.size());
    return error(m2_[bin] + m2_correction_[bin]);
  }
  double inverse_cumulative_standard_error(std::size_t step) const {
    assert(moments_enabled_ && step < cumulative_m2_.size());
    return error(cumulative_m2_[step] + cumulative_m2_correction_[step]);
  }

private:
  long double value(std::size_t bin) const {
    return weights_[bin] + weights_correction_[bin];
  }
  void accumulate(long double increment, long double &sum,
                  long double &correction) const {
    if (!precise_) { sum += increment; return; }
    const long double next = sum + increment;
    correction += std::abs(sum) >= std::abs(increment)
                ? (sum - next) + increment : (increment - next) + sum;
    sum = next;
  }
  void update(long double value, long double &mean, long double &mean_correction,
              long double &m2, long double &m2_correction) {
    const long double delta = (value - mean) - mean_correction;
    accumulate(delta / static_cast<long double>(events_), mean, mean_correction);
    accumulate(delta * ((value - mean) - mean_correction), m2, m2_correction);
    if (!std::isfinite(mean) || !std::isfinite(mean_correction) ||
        !std::isfinite(m2) || !std::isfinite(m2_correction))
      throw std::overflow_error("Nonfinite complete-event variance accumulator");
  }
  double error(long double m2) const {
    if (!std::isfinite(m2))
      throw std::overflow_error("Nonfinite complete-event variance accumulator");
    if (events_ < 2) return std::numeric_limits<double>::quiet_NaN();
    // Welford accumulation avoids subtracting two large nearly equal sums.
    const long double denominator = static_cast<long double>(events_) *
                                    static_cast<long double>(events_ - 1);
    return static_cast<double>(std::sqrt(std::max(0.L, m2) / denominator));
  }

  std::size_t nbins_, events_;
  bool moments_enabled_, active_, precise_;
  double last_squared_norm_;
  std::vector<long double> weights_, mean_, m2_;
  std::vector<long double> weights_correction_, mean_correction_, m2_correction_;
  std::vector<long double> cumulative_mean_, cumulative_m2_;
  std::vector<long double> cumulative_mean_correction_, cumulative_m2_correction_;
  std::vector<bool> touched_;
  std::vector<std::size_t> touched_bins_;
};

#endif // __EVENTSTATISTICS_HH__
