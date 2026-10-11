#ifndef __GSLRANDOM_HH__
#define __GSLRANDOM_HH__

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

//---- get a shared pointer
#if __APPLE_CC__ >= 6000
#include <memory>
#define SHARED_PTR std::shared_ptr
#else
#include <tr1/memory>
#define SHARED_PTR std::tr1::shared_ptr
#endif // end of getting shared pointer

#include <string>
#include <vector>
#include <cstring>
#include <cstdint>
#include <stdexcept>

/// A (partial) C++ interface to the GSL random number generators;
///
/// - to link remember to use -lgsl -lgslcblas
///
/// - names of generators are to be found at
///     http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
///   & include gsl_rng_mt19937 (default) 
///             gsl_rng_ranlxs0 / s1 / s2 (second generation ranlux)
///             gsl_rng_ranlxd1 / d2      (48 bits, double precision output)
///             gsl_rng_ranlux / ranlux289 (orig. ranlux; 389: 24 decorr bits)
///             gsl_rng_cmrg (combined multiple recursive gen [l'Ecuyer])
///      etc.... (see docs for more info)
///
/// - some enquiry functions haven't been implemented (e.g. for state, range).
///
class GSLRandom {
public:
  struct Snapshot {
    std::vector<unsigned char> bytes;
    std::uint64_t uniform_draws = 0;
    const gsl_rng_type* type = nullptr;
  };
  inline GSLRandom() {
    reset(gsl_rng_alloc(gsl_rng_mt19937));
  }
  inline GSLRandom(unsigned long int s) {
    reset(gsl_rng_alloc(gsl_rng_mt19937));
    set(s);
  }
  inline GSLRandom(const gsl_rng_type * T) {reset(gsl_rng_alloc(T));}
  inline GSLRandom(const gsl_rng_type * T, unsigned long int s) {
    reset(gsl_rng_alloc(T));
    set(s);
  }

  inline void reset(gsl_rng * r) {
    require_no_replay("reset");
    _r.reset(r, gsl_rng_free);
    reset_counters();
  }

  /// set the seed
  inline void set(unsigned long int s) {
    require_no_replay("set");
    gsl_rng_set(_r.get(), s);
    reset_counters();
  }

  /// Deep checkpoints used to couple shower continuations. GSLRandom copies
  /// themselves share ownership and must not be used as independent snapshots.
  void save(Snapshot& snapshot) const {
    require_no_replay("save");
    snapshot.bytes.resize(gsl_rng_size(_r.get()));
    std::memcpy(snapshot.bytes.data(), gsl_rng_state(_r.get()), snapshot.bytes.size());
    snapshot_bytes_copied_ += snapshot.bytes.size();
    snapshot.uniform_draws = uniform_draws_;
    snapshot.type = _r.get()->type;
  }
  void restore(const Snapshot& snapshot) {
    require_no_replay("restore");
    if (snapshot.type != _r.get()->type || snapshot.bytes.size() != gsl_rng_size(_r.get()))
      throw std::invalid_argument("Incompatible random-generator checkpoint");
    std::memcpy(gsl_rng_state(_r.get()), snapshot.bytes.data(), snapshot.bytes.size());
    snapshot_bytes_copied_ += snapshot.bytes.size();
    uniform_draws_ = snapshot.uniform_draws;
  }
  std::uint64_t uniform_draws() const {return uniform_draws_;}

  /// Uniform-only continuation replay. The generator is never rewound:
  /// partners consume the same recorded prefix, appending additional draws
  /// only if they pass the longest earlier continuation. Once committed,
  /// the next independent draw therefore follows the furthest endpoint.
  /// Groups cannot nest. The vector retains its capacity between groups.
  void begin_replay_group() {
    require_no_replay("begin_replay_group");
    replay_start_ = uniform_draws_;
    replay_cursor_ = 0;
    replay_tape_.clear();
    replay_active_ = true;
  }
  void rewind_replay_group() {
    require_replay("rewind_replay_group");
    replay_cursor_ = 0;
    uniform_draws_ = replay_start_;
  }
  void end_replay_group() {
    require_replay("end_replay_group");
    uniform_draws_ = replay_start_ + replay_tape_.size();
    replay_cursor_ = replay_tape_.size();
    replay_active_ = false;
  }
  bool replay_group_active() const {return replay_active_;}
  std::size_t replay_tape_capacity() const {return replay_tape_.capacity();}
  std::uint64_t generator_uniform_draws() const {return generator_uniform_draws_;}
  std::uint64_t snapshot_bytes_copied() const {return snapshot_bytes_copied_;}

  /// returns in range [0,1) (includes 0, excludes 1)
  inline double uniform() const {
    if (replay_active_) {
      if (replay_cursor_ == replay_tape_.size()) {
        // Allocate before advancing GSL so an allocation failure cannot lose
        // a generator value. Normal groups reuse this storage after warmup.
        if (replay_tape_.size() == replay_tape_.capacity())
          replay_tape_.reserve(replay_tape_.capacity() ? 2*replay_tape_.capacity() : 256);
        replay_tape_.push_back(gsl_rng_uniform(_r.get()));
        ++generator_uniform_draws_;
      }
      const double value = replay_tape_[replay_cursor_++];
      uniform_draws_ = replay_start_ + replay_cursor_;
      return value;
    }
    ++uniform_draws_;
    ++generator_uniform_draws_;
    return gsl_rng_uniform(_r.get());
  }

  /// returns in range (0,1) (excludes 0, excludes 1)
  inline double uniform_pos() const {
    double value;
    do {value = uniform();} while (value == 0.0);
    return value;
  }

  /// returns in range [xmin,xmax) (includes xmin, excludes xmax -- modulo rounding)
  inline double uniform(double xmin, double xmax) const {return xmin + (xmax-xmin)*uniform();}

  /// returns +-1 with equal probability
  inline double sign() const {return uniform() < 0.5 ? -1 : 1;}

  /// returns a Gaussian distributed random number with standard
  /// deviation sigma and mean zero
  inline double gaussian(double sigma) const {
    require_no_replay("gaussian");
    return gsl_ran_gaussian(_r.get(),sigma);
  }
  /// Gaussian with sigma=1
  inline double gaussian() const {
    require_no_replay("gaussian");
    return gsl_ran_ugaussian(_r.get());
  }

  /// return random x according to p(x) dx = {1 \over \mu} \exp(-x/\mu) dx
  inline double exponential(double mu) const {
    require_no_replay("exponential");
    return gsl_ran_exponential(_r.get(),mu);
  }

  // returns random k according to p(k) = {\mu^k \over k!} \exp(-\mu)
  inline unsigned int poisson(double mu) const {
    require_no_replay("poisson");
    return gsl_ran_poisson(_r.get(),mu);
  }

  // returns a number generated according to the gamma distribution
  inline double gamma(double a, double b) const {
    require_no_replay("gamma");
    return gsl_ran_gamma(_r.get(),a,b);
  }

  /// returns integer in range [0,n-1]
  inline unsigned long int uniform_int(unsigned long int n) const {
    if (n <= 1) return 0;
    require_no_replay("uniform_int");
    return gsl_rng_uniform_int(_r.get(), n);
  }

  inline std::string name() const {return std::string(gsl_rng_name (_r.get()));}

  inline gsl_rng * gsl_generator() {
    require_no_replay("gsl_generator");
    return _r.get();
  }
  
  /// no specific destructor
  //inline ~GSLRandom() {}
protected:
  SHARED_PTR<gsl_rng> _r;
  // Shower continuations use only uniform()/uniform_pos(), allowing exact
  // comparison of consumed stream positions, including rejected zero draws.
  mutable std::uint64_t uniform_draws_ = 0;
  mutable std::uint64_t generator_uniform_draws_ = 0;
  mutable std::uint64_t snapshot_bytes_copied_ = 0;
  mutable std::vector<double> replay_tape_;
  mutable std::size_t replay_cursor_ = 0;
  std::uint64_t replay_start_ = 0;
  bool replay_active_ = false;

  void require_no_replay(const char* operation) const {
    if (replay_active_)
      throw std::logic_error(std::string(operation) + " during an active uniform replay group");
  }
  void require_replay(const char* operation) const {
    if (!replay_active_)
      throw std::logic_error(std::string(operation) + " without an active uniform replay group");
  }
  void reset_counters() {
    uniform_draws_ = 0;
    generator_uniform_draws_ = 0;
    snapshot_bytes_copied_ = 0;
    replay_cursor_ = 0;
    replay_start_ = 0;
    replay_tape_.clear();
  }
};

#endif // __GSLRANDOM_HH__
