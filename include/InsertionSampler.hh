#ifndef GNOLE_INSERTION_SAMPLER_HH
#define GNOLE_INSERTION_SAMPLER_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

/// A draw of the base variable r in ln(k_t) = ln(x_Q) + r^power * range.
/// The existing insertion Jacobian must be multiplied by inverse_density.
struct InsertionSample {
  double r;
  double inverse_density;
  std::size_t bin;
  double density;
};

/// Piecewise constant importance sampling of the original uniform r variable.
/// Updates use completed events only; the density saved with each draw remains
/// its exact density even when record() updates the proposal for the next draw.
/// A positive uniform mixture guarantees support and bounds 1/q(r).
class InsertionSampler {
public:
  explicit InsertionSampler(bool adaptive = false, std::size_t bins = 20,
                            std::size_t update_interval = 5000,
                            double uniform_fraction = 0.2)
      : adaptive_(adaptive), bins_(bins), update_interval_(update_interval),
        uniform_fraction_(uniform_fraction), records_(0), updates_(0) {
    if (bins_ == 0 || update_interval_ == 0 ||
        !std::isfinite(uniform_fraction_) || uniform_fraction_ <= 0.0 ||
        uniform_fraction_ > 1.0)
      throw std::invalid_argument("Invalid insertion sampler configuration");
    reset(adaptive);
  }

  void reset(bool adaptive) {
    adaptive_ = adaptive;
    records_ = updates_ = 0;
    probability_.assign(bins_, 1.0 / static_cast<double>(bins_));
    cdf_.resize(bins_);
    square_sum_.assign(bins_, 0.0L);
    count_.assign(bins_, 0);
    rebuild_cdf();
  }

  InsertionSample sample(double uniform) const {
    if (!std::isfinite(uniform) || uniform < 0.0 || uniform >= 1.0)
      throw std::invalid_argument("Insertion uniform variate must be in [0,1)");
    if (!adaptive_) {
      // Preserve the original draw exactly, including its floating-point bits.
      const std::size_t bin = std::min(
          static_cast<std::size_t>(uniform * bins_), bins_ - 1);
      return {uniform, 1.0, bin, 1.0};
    }
    const std::size_t bin = static_cast<std::size_t>(
        std::upper_bound(cdf_.begin(), cdf_.end(), uniform) - cdf_.begin());
    const double left = bin == 0 ? 0.0 : cdf_[bin - 1];
    // Use the actual CDF width, so floating-point normalization is included.
    const double width = cdf_[bin] - left;
    const double density = static_cast<double>(bins_) * width;
    const double fraction = (uniform - left) / width;
    double r = (static_cast<double>(bin) + fraction) /
               static_cast<double>(bins_);
    // Guard only endpoint rounding; never clip or otherwise alter weights.
    if (r >= 1.0) r = std::nextafter(1.0, 0.0);
    return {r, 1.0 / density, bin, density};
  }

  /// squared_norm is the sum of squared *net* observable-bin contributions
  /// from this insertion after summing its signed S2 branches. For pure
  /// importance-weighted contributions, q^2 recovers the unsampled squared
  /// integrand. The caller weights its entire S2 estimator, including the
  /// LL/H1 baseline. Proposal learning never alters current event weights.
  void record(const InsertionSample &draw, double squared_norm) {
    if (!adaptive_) return;
    if (draw.bin >= bins_ || !std::isfinite(draw.density) || draw.density <= 0.0)
      throw std::invalid_argument("Invalid insertion draw passed to record");
    if (squared_norm < 0.0)
      throw std::invalid_argument("Insertion squared norm must be nonnegative");
    ++records_;
    // A nonfinite learning score carries no useful proposal information. This
    // does not truncate/reweight the observable; the caller owns that estimate.
    if (std::isfinite(squared_norm)) {
      const long double score = static_cast<long double>(squared_norm) *
                                draw.density * draw.density;
      if (std::isfinite(score)) {
        square_sum_[draw.bin] += score;
        ++count_[draw.bin];
      }
    }
    if (records_ % update_interval_ == 0) update();
  }

  const std::vector<double> &probabilities() const { return probability_; }
  std::size_t records() const { return records_; }
  std::size_t updates() const { return updates_; }
  bool adaptive() const { return adaptive_; }

private:
  void rebuild_cdf() {
    double cumulative = 0.0;
    for (std::size_t i = 0; i < bins_; ++i) {
      cumulative += probability_[i];
      cdf_[i] = cumulative;
    }
    cdf_.back() = 1.0;
    // Expose the probabilities actually used by sample(), including rounding.
    for (std::size_t i = 0; i < bins_; ++i)
      probability_[i] = cdf_[i] - (i == 0 ? 0.0 : cdf_[i - 1]);
  }

  void update() {
    long double total_square = 0.0L;
    std::size_t total_count = 0;
    for (std::size_t i = 0; i < bins_; ++i) {
      total_square += square_sum_[i];
      total_count += count_[i];
    }
    if (total_count == 0 || total_square <= 0.0L ||
        !std::isfinite(total_square))
      return;
    const long double global_mean = total_square / total_count;
    // A modest common prior protects sparsely sampled bins during warmup.
    // Cumulative conditional means are valid because q is constant within
    // each r bin, even when its selection probability changes between events.
    constexpr long double prior_count = 32.0L;
    std::vector<long double> target(bins_);
    long double total_target = 0.0L;
    for (std::size_t i = 0; i < bins_; ++i) {
      target[i] = std::sqrt((square_sum_[i] + prior_count * global_mean) /
                            (count_[i] + prior_count));
      total_target += target[i];
    }
    if (total_target <= 0.0L || !std::isfinite(total_target)) return;
    const double floor = uniform_fraction_ / static_cast<double>(bins_);
    for (std::size_t i = 0; i < bins_; ++i) {
      const double learned = floor + (1.0 - uniform_fraction_) *
          static_cast<double>(target[i] / total_target);
      // Smooth updates to avoid a noisy block replacing the entire proposal.
      probability_[i] = 0.5 * probability_[i] + 0.5 * learned;
    }
    rebuild_cdf();
    ++updates_;
  }

  bool adaptive_;
  std::size_t bins_;
  std::size_t update_interval_;
  double uniform_fraction_;
  std::size_t records_;
  std::size_t updates_;
  std::vector<double> probability_;
  std::vector<double> cdf_;
  std::vector<long double> square_sum_;
  std::vector<std::size_t> count_;
};

#endif
