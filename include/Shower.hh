/// Code for the core elements of a shower.
#ifndef __SHOWER_HH__
#define __SHOWER_HH__
#include "Event.hh"
#include "GSLRandom.hh"
#include "Observables.hh"
#include "matrix_elements.hh"
#include "InsertionSampler.hh"
#include <string>
#include <cassert>
#include <iostream>
#include <limits>

struct VarianceOptions {
  enum Mode {Legacy, Merge, Coupled};
  Mode mode = Coupled;
  bool adaptive_insertion = false;
  double insertion_power = 5.0;
  double second_power = 1.0;
  bool optimize = true;
  bool replay_tape = true;
  bool diagnostics = false;
  int output_precision = 17;
};

struct CoupledRandomState {
  GSLRandom::Snapshot start, furthest;
  bool valid = false;
};

//----------------------------------------------------------------------
/// \class Shower
/// contains the core shower evolution
class Shower {
public:
  /// constructor
  Shower(Observable& obs,
	 double xmur,
	 double xQ,
	 int order_evl,
	 std::string header = "",
	 int seed           = 0,
         VarianceOptions variance = VarianceOptions())
    : xmur_(xmur), xQ_(xQ), order_evl_(order_evl), obs_(&obs),
      evolution_limit_(std::numeric_limits<double>::infinity()),
      include_hard_factor_(false), header_(header),
      rng(seed), gluon_(0,0,0,0), event_cache_(new Event()),
      variance_(variance), insertion_sampler_(variance.adaptive_insertion) {
    assert((order_evl_==0) || (order_evl_==1));
    asmur_=alphas2(xmur_);
    log_xQ_ = log(xQ_);
    event_.cache_coupling(asmur_, variance_.optimize);
    // Preserve the former grid's endpoint veto for expanded NLL.
    if (order_evl_==1) {
      evolution_limit_ = legacy_evolution_limit(xmur_, xQ_);
      // set the integrated coefficient
      integrated_counterterm_ = integrated_counterterm(obs_->parameter());
      hard_factor_ = 1.0 + asmur_/(2.0*M_PI)*(CF*integrated_counterterm_+H1);
    }
  }

  Shower(const Shower&) = delete;
  Shower& operator=(const Shower&) = delete;

  /// destructor
  virtual ~Shower() {
    delete event_cache_;
  }

  /// run the shower with nev events
  virtual void run(int nev, const std::string& fn = "");
  
protected:
  /// run the threejet piece
  virtual void run_threejet(double tstart, bool soft = false, bool use_cached_variables = false);
  void run_threejet_combined(double tstart);
  
  /// reset the dipoles to an initial qqbar pair along the z axis
  void reset(bool threejet = false, bool soft = false, bool use_cached_variables = false);

  /// choose an emitting dipole
  int choose_emitter() const;

  /// decide whether to split the dipole using new emission
  bool do_split(int idip, Momentum&  emsn);

  /// evolve the shower between two scales
  void evolve_scale(double t, double tend = evol_cutoff_, bool include_as_constant = true);
  /// Rebuild a damaged incremental rapidity sum before sampling radiation.
  bool ensure_rapidity_sum();
  
  /// evolve the shower including insertion (NLL fully expanded out)
  void evolve_insertion_expanded(double tstart);

  /// singles out single insertion corrections (NLL fully expanded out)
  void evolve_single_insertion_corrections(double tstart);

  /// singles out double insertion corrections (NLL fully expanded out)
  void evolve_double_insertion_corrections(double tstart);

  /// perform the Z^{(0)} evolution
  void perform_branch0(double ta, const Momentum& ka);

  /// perform the perturbative NLL insertion of the Z^{(0)} evolution
  void perform_branch_single_insertion(double t_insertion, int ibranch, const Momentum& ka);

  /// perform the Z^{(1)} evolution
  void perform_branch_double_insertion(double t_insertion, int idipa, int ibranch, const Momentum& ka);

  /// perform the Z^{(1)} evolution (fermionic case)
  void perform_branch_double_insertion_fermion(double t_insertion, int idipa, int ibranch, const Momentum& ka);

  /// generate and insert the first insertion
  Momentum generate_first_insertion(double& t_insertion, int& idip_insertion);
  
  /// generate the second insertion
  Momentum generate_second_insertion(double t_insertion, int idip, int& idip_insertion, bool dipole_kt_ordering);
  
  /// return the t scale for a given ln kt value
  double t_scale(double lnkt) const;
  
  /// return the ln(kt) for a given evolution scale
  double ln_kt(double t) const;

  /// write current observable to ostream
  void write(int nev, const std::string& fn) const;

  /// reconstruct momentum of parent
  void reconstruct_parent(const Momentum& spec_left, const Momentum& spec_right, Momentum& kab) const;
  double log_scale_interval(double later, double earlier) const;
  double gauge_antenna_weight(const Momentum& a, const Momentum& b,
                             const Momentum& c, const Momentum& d, int branch);
  double fermion_antenna_weight(const Momentum& a, const Momentum& b,
                               const Momentum& c, const Momentum& d);
  double log_xQ() const { return variance_.optimize ? log_xQ_ : log(xQ_); }
  void collect_antenna_weight(double value) {
    antenna_weight_cache_ = value;
    antenna_weight_valid_ = true;
  }
  
  /// event containing dipole configuration
  Event event_;
  /// Evolution time of the last emission (for truncated shower)
  double tlast_ = 0.;
  /// renormalisation scale
  double xmur_;
  /// resummation scale
  double xQ_;
  /// alphas at renormalisation scale
  double asmur_;
  /// order of the shower (LL=0, expanded NLL=1)
  int order_evl_;
  /// pointer to observable
  Observable * obs_;
  /// legacy endpoint cutoff; inactive (infinite) for LL
  double evolution_limit_;
  /// include the H2 hard factor during the truncated S2 evolution
  bool include_hard_factor_;
  /// integrated counterterm
  double integrated_counterterm_;
  double log_xQ_ = 0.0, hard_factor_ = 1.0;
  /// header with description of run
  std::string header_;
  /// random number generator
  GSLRandom rng;
  /// momentum of hard gluon
  Momentum gluon_;
  /// flag for kt ordering in the parent dipole frame
  bool lab_ordering_;
  /// flag for kt ordering in the (respective) emitting dipole frame
  bool dip_ordering_;
  
  /// cache second insertion info
  double t_second_insertion_; 
  Momentum momentum_cache_;
  int  idip_cache_;
  bool event_bad_cache_;
  void cache_second_insertion(Momentum kb, int idipb, bool event_bad) {
    momentum_cache_  = kb;
    idip_cache_      = idipb;
    event_bad_cache_ = event_bad;
    antenna_weight_valid_ = false;
  }
  void retrieve_second_insertion(Momentum& kb, int& idipb, bool& event_bad) {
    kb        = momentum_cache_;
    idipb     = idip_cache_;
    event_bad = event_bad_cache_;
  }
  bool   weighted_second_insertion = true; // if true the second insertion is generated at FO with a weight (instead of with a Sudakov)
  double second_insertion_weight_;

  VarianceOptions variance_;
  InsertionSampler insertion_sampler_;
  InsertionSample insertion_sample_{};
  bool sampled_insertion_ = false;
  double insertion_base_weight_ = 1.0;
  CoupledRandomState single_random_, double_random_, hard_random_;
  Event soft_event_scratch_, hard_event_scratch_;
  double antenna_weight_cache_ = 0.0;
  bool antenna_weight_valid_ = false;
  struct Metrics {
    std::uint64_t emitter_calls = 0, emitter_dipoles = 0;
    unsigned maximum_dipoles = 0;
    std::uint64_t zero_weight_skips = 0, antenna_evaluations = 0, antenna_cache_hits = 0;
    std::uint64_t nonfinite_branches = 0;
    std::uint64_t rapidity_sum_rebuilds = 0;
  };
  mutable Metrics metrics_;
  double sample_insertion_log(double& jacobian);
  
  /// cached information for insertion
  Event * event_cache_;
  /// cutoffs for the evolution
  static constexpr double evol_cutoff_ = EVOLCUT;
  static constexpr double landau_pole_tolerance_ = 1.;
};

#endif //__SHOWER_HH__
