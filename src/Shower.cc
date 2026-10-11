#include "Shower.hh"
#include "Parameters.hh"
#include "Dipole.hh"
#include "helpers.hh"
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <iomanip>

namespace {
// Replay a pair/group from one uniform stream, then commit the furthest
// consumed position. Future events never reuse draws inspected by a partner.
class ContinuationGuard {
public:
  ContinuationGuard(bool enabled, bool first, bool final,
                    GSLRandom& rng, CoupledRandomState& state, bool tape = false)
    : enabled_(enabled), final_(final), tape_(tape), rng_(rng), state_(state) {
    if (!enabled_) return;
    if (tape_) {
      if (first) rng_.begin_replay_group();
      else rng_.rewind_replay_group();
      return;
    }
    if (first) {
      rng_.save(state_.start);
      rng_.save(state_.furthest);
      state_.valid = true;
    } else {
      if (!state_.valid) throw std::logic_error("Missing coupled continuation");
      rng_.restore(state_.start);
    }
  }
  ~ContinuationGuard() {
    if (!enabled_) return;
    if (tape_) {
      if (final_) rng_.end_replay_group();
      return;
    }
    if (rng_.uniform_draws() > state_.furthest.uniform_draws)
      rng_.save(state_.furthest);
    if (final_) {
      rng_.restore(state_.furthest);
      state_.valid = false;
    }
  }
private:
  bool enabled_, final_, tape_;
  GSLRandom& rng_;
  CoupledRandomState& state_;
};
}

//----------------------------------------------------------------------
/// run the shower with nev events
void Shower::run(int nev, const std::string& fn) {
  if (nev <= 0) throw std::invalid_argument("Number of events must be positive");
  for(int i=0; i<nev; ++i) {
    obs_->begin_event();
    sampled_insertion_ = false;
    reset();
    double tstart = 0.0;
    if (order_evl_ == 0) 
      evolve_scale(tstart);
    else if (order_evl_ == 1) {
      if (HARD_ONLY) {
        include_hard_factor_ = true;
        evolve_scale(tstart);
      } else if (SINGLE_INSERTION_ONLY) {
        // run only corrections encoded in the single insertion
        evolve_single_insertion_corrections(tstart);
      } else if (DOUBLE_INSERTION_ONLY) {
        // run only corrections encoded in the double insertion
        evolve_double_insertion_corrections(tstart);
      } else {
        // run the soft factor S2 with one-loop evolution convoluted with H2 at NLO
        // followed by the two loop corrections to S2 convoluted with H2 at LO
        evolve_insertion_expanded(tstart);
      }

      if (variance_.adaptive_insertion && sampled_insertion_)
        insertion_sampler_.record(insertion_sample_, obs_->current_event_squared_norm());

      if ((!SINGLE_INSERTION_ONLY) && (!DOUBLE_INSERTION_ONLY)) {
	// run the soft factor S3 with one-loop evolution convoluted with H3 at LO
        if (variance_.mode == VarianceOptions::Legacy) {
          run_threejet(tstart, true);
          run_threejet(tstart, false, true);
        } else {
          run_threejet_combined(tstart);
        }
      }
    }
    obs_->end_event();
    if (periodic_output()) {
      std::cerr << "# " << i+1 << " out of " << nev << " events." << std::endl;
      if (!fn.empty()) write(i+1, fn);
    }
  }
  write(nev, fn);
}

//----------------------------------------------------------------------
/// run the threejet piece
void Shower::run_threejet(double tstart, bool soft, bool use_cached_variables) {
  include_hard_factor_ = false; // S3 evolves at LL, also after the hard-only S2 branch.
  reset(true, soft, use_cached_variables);
  // flip sign for LL term
  if (soft) event_.weight = -event_.weight;
  evolve_scale(tstart);
  double w = -event_.weight;
  reset(false);
  event_.weight = w;
  // now veto configurations where gluon is in the slice as computed in 
  // the two jet limit within the integrated counterterm
  if (!MATCHING) {
    if (obs_->in_region(gluon_, &event_.axis())) event_.bad=true;
  }
  evolve_scale(tstart);
}

//----------------------------------------------------------------------
// Combine the two identical S2 counterterms before their common evolution.
// S3 supports and twojet gluon vetoes are deliberately evaluated separately.
void Shower::run_threejet_combined(double tstart) {
  include_hard_factor_ = false;
  reset(true, true);
  event_.weight = -event_.weight;
  soft_event_scratch_ = event_;
  Momentum soft_gluon = gluon_;
  const double soft_counter_weight = -event_.weight;

  // The cached hard-phase coordinates survive reset and are shared with soft.
  // Matching deposits made by reset() remain part of this same outer event.
  reset(true, false, true);
  hard_event_scratch_ = event_;
  Momentum hard_gluon = gluon_;
  const double hard_counter_weight = -event_.weight;
  const bool coupled = variance_.mode == VarianceOptions::Coupled;

  event_ = soft_event_scratch_;
  {
    ContinuationGuard continuation(coupled, true, false, rng, hard_random_, variance_.replay_tape);
    evolve_scale(tstart);
  }
  event_ = hard_event_scratch_;
  {
    ContinuationGuard continuation(coupled, false, false, rng, hard_random_, variance_.replay_tape);
    evolve_scale(tstart);
  }

  reset(false);
  double weight = 0.0;
  if (MATCHING || !obs_->in_region(soft_gluon, &event_.axis()))
    weight += soft_counter_weight;
  if (MATCHING || !obs_->in_region(hard_gluon, &event_.axis()))
    weight += hard_counter_weight;
  event_.weight = weight;
  {
    ContinuationGuard continuation(coupled, false, true, rng, hard_random_, variance_.replay_tape);
    if (weight != 0.0) evolve_scale(tstart);
  }
  gluon_ = hard_gluon;
}

//----------------------------------------------------------------------
double Shower::sample_insertion_log(double& jacobian) {
  insertion_sample_ = insertion_sampler_.sample(variance_.insertion_power < 5.0
                                               ? rng.uniform_pos() : rng.uniform());
  sampled_insertion_ = true;
  const double r = insertion_sample_.r;
  const double n = variance_.insertion_power;
  const double range = landau_pole_tolerance_*(1.0/(2.0*asmur_*b0) - log_xQ());
  // Preserve the entire original estimator, including its boundary terms:
  // its canonical variable was r0=u^(n/5), uniform for n=5. All S2 terms,
  // including the LL prefix and real '1', get dr0/du divided by q(u).
  // This remains valid when xQ<1 gives negative sampled insertion times.
  insertion_base_weight_ = insertion_sample_.inverse_density *
    (n == 5.0 ? 1.0 : n/5.0*pow(r,n/5.0-1.0));
  event_.weight = insertion_base_weight_;
  jacobian = n*pow(r,n-1)*range*insertion_sample_.inverse_density;
  return log_xQ() + pow(r,n)*range;
}

//----------------------------------------------------------------------
/// reset the dipoles to an initial qqbar pair along the z axis
void Shower::reset(bool threejet, bool soft, bool use_cached_variables) {
  event_.bad=false;
  if (!threejet) {
    event_.reset_twojet();
  } else if (soft){
    bool in_region = false;
    //event_.reset_threejet_LL(xmur_, xQ_, rng.uniform(), rng.uniform(), rng.uniform(), gluon_);
    event_.reset_threejet_soft(xmur_, xQ_, rng.uniform(), rng.uniform(), rng.uniform(), gluon_);
    if (event_.bad) return;
    Shower::do_split(0,gluon_);
    if (std::abs(gluon_.rap()) > RAPMAX) event_.bad = true;
    if (!MATCHING) {
      in_region = obs_->in_region(gluon_, &event_.axis());
      if (in_region or event_.weight==0.0) event_.bad = true;
    } else if (!event_.bad) {
      // fill histogram if a particle is in the slice
      if (obs_->add_entries_in_region(gluon_.stored_E()*gluon_, 0., 0., -event_.weight, &event_.axis())) {
        event_.bad = true; // setting this to avoid starting the evolution later
      }
      if (event_.weight==0.0) event_.bad = true;
    }
  } else {
    bool in_region = false;
    event_.reset_threejet(xmur_, xQ_, rng.uniform(), rng.uniform(), rng.uniform(), gluon_, use_cached_variables);
    if (event_.bad) return;
    Shower::do_split(0,gluon_);

    // make sure we don't have any hard partons in the observed region, otherwise set weight to 0
    if (!MATCHING) {
      in_region = obs_->in_region(event_[0].left().momentum(), &event_.axis());
      in_region = in_region || (obs_->in_region(event_[0].right().momentum(), &event_.axis()));
      in_region = in_region || (obs_->in_region(event_[1].right().momentum(), &event_.axis()));

      if (in_region or event_.weight==0.0) event_.bad = true;
    } else if (!event_.bad) {
      // fill histogram if a particle is in the slice (additive observables)
      double observable;
      observable  = obs_->transverse_energy(event_[0].left() .momentum().stored_E()*event_[0].left() .momentum(), 0., &event_.axis());
      observable += obs_->transverse_energy(event_[0].right().momentum().stored_E()*event_[0].right().momentum(), 0., &event_.axis());
      observable += obs_->transverse_energy(event_[1].right().momentum().stored_E()*event_[1].right().momentum(), 0., &event_.axis());

      if (obs_->add_entries(observable, 0., event_.weight)) {event_.bad = true;}
      if (event_.weight==0.0) event_.bad = true;
    }
  }
  return;
}

//----------------------------------------------------------------------
/// choose an emitting dipole
int Shower::choose_emitter() const {
  if (variance_.diagnostics) {
    ++metrics_.emitter_calls;
    metrics_.emitter_dipoles += event_.size();
    metrics_.maximum_dipoles = std::max(metrics_.maximum_dipoles, event_.size());
  }
  double rand_rap = rng.uniform() * event_.eta_tot;
  double curr_rap = 0.0;
  for (unsigned int idip=0; idip < event_.size(); ++idip) {
    curr_rap += event_[idip].delta_rap();
    if (rand_rap < curr_rap) return idip;
  }
  return -1;
}

//----------------------------------------------------------------------
/// split a dipole and insert the emission into the event
/// The algorithm in this function was originally written by G. Salam
bool Shower::do_split(int idip, Momentum& emsn) {
  // now insert the emission
  Dipole & emitter = event_[idip];
  
  // remove current value of delta_rap of dipole
  // this is to be added back at the end with updated value
  event_.eta_tot -= event_[idip].delta_rap();
  
  Dipole new_dipole(DipoleEnd(emsn), emitter.right(),
		    idip, emitter.right_neighbour());
  if (emitter.right_neighbour() >= 0)
    event_[emitter.right_neighbour()].left_neighbour(event_.size());
  emitter.right_neighbour(event_.size());
  emitter.right(DipoleEnd(emsn));
  // adding newly created dipole to end of dipoles vector
  event_.add(new_dipole);    
  // finally update the eta_tot value
  event_.eta_tot += event_.dipoles().back().delta_rap();
  event_.eta_tot += event_[idip].delta_rap();
  return true;
}

//----------------------------------------------------------------------
/// evolve the shower down to specified scale or down to cutoff scale
bool Shower::ensure_rapidity_sum() {
  // This is a trigger for exact recomputation, never an emission cutoff:
  // any genuinely positive rebuilt range, however small, remains active.
  const double rounding_scale = 64. * std::numeric_limits<double>::epsilon()
    * (2. * RAPMAX + 20.) * event_.size();
  if (event_.eta_tot > rounding_scale && std::isfinite(event_.eta_tot)) return true;
  double sum = 0., correction = 0.;
  for (unsigned int i = 0; i < event_.size(); ++i) {
    const double range = event_[i].delta_rap();
    if (!std::isfinite(range) || range < 0.)
      throw std::runtime_error("Invalid dipole rapidity range while rebuilding emission support");
    histogram_accumulate(range, sum, correction, true);
  }
  event_.eta_tot = sum + correction;
  if (!std::isfinite(event_.eta_tot))
    throw std::runtime_error("Nonfinite rebuilt dipole rapidity sum");
  if (variance_.diagnostics) ++metrics_.rapidity_sum_rebuilds;
  return event_.eta_tot > 0.;
}

void Shower::evolve_scale(double t, double tend, bool include_as_constant) {
  //if (event_.weight == 0.0 or event_.bad) return;
  if (event_.bad) return;
  tlast_ = t;
  if (variance_.optimize && variance_.mode != VarianceOptions::Legacy && event_.weight == 0.0) {
    if (variance_.diagnostics) ++metrics_.zero_weight_skips;
    return;
  }

  while (true) {
    // Subtracting obsolete dipole ranges can leave a negative rounding
    // residual when the exact emission support is zero. Check after every
    // split, before drawing a waiting time, while keeping reference replay.
    if (PRECISE_KERNELS && !ensure_rapidity_sum()) return;
    t += -log(rng.uniform_pos()) / (2.0 * CA * event_.eta_tot);
    if (!(t < tend)) break;
    // increase by amount chosen with distribution e^(-2*CA*sum_rap)    
    int idip = choose_emitter();
    // now set all emission properties
    double lnkt = ln_kt(t);
    if ((2.0*asmur_*b0*lnkt >= 1.0) or (t >= evolution_limit_) or (lnkt > lnktmax)) {
      event_.bad = true; // setting this to avoid starting another evolution later
      break;
    }

    //fixed coupling: if (event_.eta_tot < 0.0 or event_.eta_tot!=event_.eta_tot) break;
    assert(idip>=0);
    // the momentum of the gluon is xQ*exp(-lnkt)=exp[-(lnkt+ln(xQ))]

    Momentum emsn = event_[idip].radiate(lnkt - log_xQ(), rng.uniform_pos(),
					 rng.uniform_pos());
    if (!do_split(idip, emsn)) continue;

    tlast_ = t;
    // if emission is in observed region, add to histogram and stop the evolution
    double C1 = ((include_as_constant and include_hard_factor_) ?
                 (variance_.optimize ? hard_factor_ : 1.0 + asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1)) : 1.0);
    if (obs_->add_entries_in_region(emsn.stored_E()*emsn, t, lnkt-log_xQ(),
	  		    C1 * event_.weight, &event_.axis())) {
      event_.bad = true; // setting this to avoid starting another evolution later
      break;
    }
  }
  return;
}

//----------------------------------------------------------------------
/// evolve the shower down to cutoff scale while fully expanding 
/// out NLL corrections
void Shower::evolve_insertion_expanded(double t) {  
  //if (event_.weight == 0.0 or event_.bad) return;
  if (event_.bad) return;
  include_hard_factor_ = true;

  // generate log(Q/kt) in [log(xQ), 1/(2 as b0)]
  // generate scale t from log(Q/kt)
  double insertion_jacobian;
  // introduce partition of unity to identify the scale of the insertion
  double lnkt_insertion = sample_insertion_log(insertion_jacobian);
  double t_insertion = t_scale(lnkt_insertion);
  // Endpoint rounding in an extreme power map can reach the Landau pole.
  // Interpret this as an insertion beyond the cutoff, retaining its LL prefix.
  if (!std::isfinite(t_insertion)) t_insertion = std::numeric_limits<double>::infinity();

  // evolve S2 at one loop kernel and include H2 at one loop
  evolve_scale(t, std::min(t_insertion, evol_cutoff_));
  int idipa = -1;
  if (t_insertion >= evol_cutoff_) {
    include_hard_factor_ = false;
    return;
  }
  // generate the insertion
  Momentum ka = generate_first_insertion(t_insertion, idipa);
  if (event_.bad) {
    include_hard_factor_ = false;
    return;
  }

  // branch0 = Z^{(0)} evolution at LL
  if (!insertion_is_part_of_NLL_ensemble) {
    event_cache_->copy(event_);
    perform_branch0(t_insertion, ka);
    event_.retrieve(event_cache_);
  }

  // now handle extra NLL branches
  // add jacobian to account for the sampling of one of the LL emissions
  if (insertion_is_part_of_NLL_ensemble) {
    event_.weight  = insertion_jacobian / log_scale_interval(t_insertion, tlast_);
  } else {
    event_.weight  = insertion_jacobian * event_.eta_tot;
  }

  // branch -2 = NLL virtual corrections to Z^{(0)}
  event_cache_->copy(event_);
  {
    ContinuationGuard continuation(variance_.mode == VarianceOptions::Coupled,
                                   true, false, rng, single_random_, variance_.replay_tape);
    perform_branch_single_insertion(t_insertion, -2, ka);
  }
  event_.retrieve(event_cache_);

  // split the dipole that emits ka
  // branch -1 = NLL real corrections to Z^{(0)}
  do_split(idipa, ka);

  // cache the event with the first insertion in;
  // to be used later for the second insertion.  
  // Then evolve and retrieve the event before carrying on
  event_cache_->copy(event_);
  {
    ContinuationGuard continuation(variance_.mode == VarianceOptions::Coupled,
                                   false, true, rng, single_random_, variance_.replay_tape);
    perform_branch_single_insertion(t_insertion, -1, ka);
  }
  event_.retrieve(event_cache_);

  // branch1 = Z^{(1)} evolution double real with kta > ktb'
  perform_branch_double_insertion(t_insertion, idipa, 1, ka);
  event_.retrieve(event_cache_);

  // branch2 = Z^{(1)} evolution collinear counterterm with kta > ktb'
  perform_branch_double_insertion(t_insertion, idipa, 2, ka);
  event_.retrieve(event_cache_);
  
  // branch3 = Z^{(1)} evolution double real with kta > ktb
  perform_branch_double_insertion(t_insertion, idipa, 3, ka);
  event_.retrieve(event_cache_);
  
  // branch4 = Z^{(1)} evolution collinear counterterm with kta > ktb
  perform_branch_double_insertion(t_insertion, idipa, 4, ka);

  if (NF > 0.) {
    event_.retrieve(event_cache_);
    // branch1 (fermion) = Z^{(1)} evolution double real with kta > ktb'
    perform_branch_double_insertion_fermion(t_insertion, idipa, 1, ka);
    event_.retrieve(event_cache_);

    // branch2 (fermion) = Z^{(1)} evolution collinear counterterm with kta > ktb'
    perform_branch_double_insertion_fermion(t_insertion, idipa, 2, ka);
  }
  include_hard_factor_ = false;
}


//----------------------------------------------------------------------
/// singles out single insertion corrections as defined in perform_branch_single_insertion
void Shower::evolve_single_insertion_corrections(double t) {
  //if (event_.weight == 0.0 or event_.bad) return;
  if (event_.bad) return;
  include_hard_factor_ = false;

  assert(insertion_is_part_of_NLL_ensemble
         && "Single-insertion selector requires insertion_is_part_of_NLL_ensemble=true");

  // generate log(Q/kt) in [log(xQ), 1/(2 as b0)]
  // generate scale t from log(Q/kt)
  double insertion_jacobian;
  // introduce partition of unity to identify the scale of the insertion
  double lnkt_insertion = sample_insertion_log(insertion_jacobian);
  double t_insertion = t_scale(lnkt_insertion);
  if (!std::isfinite(t_insertion)) t_insertion = std::numeric_limits<double>::infinity();

  // evolve S2 at one loop kernel (without O(as) constants)
  evolve_scale(t, std::min(t_insertion, evol_cutoff_), include_hard_factor_);
  int idipa = -1;
  if (t_insertion >= evol_cutoff_) {return;}

  // generate the insertion
  Momentum ka = generate_first_insertion(t_insertion, idipa);
  if (event_.bad) {return;}

  // now handle extra NLL branches
  event_.weight  = insertion_jacobian / log_scale_interval(t_insertion, tlast_);

  // branch -2 = NLL virtual corrections to Z^{(0)}
  event_cache_->copy(event_);
  {
    ContinuationGuard continuation(variance_.mode == VarianceOptions::Coupled,
                                   true, false, rng, single_random_, variance_.replay_tape);
    perform_branch_single_insertion(t_insertion, -2, ka);
  }
  event_.retrieve(event_cache_);

  // split the dipole that emits ka
  // branch -1 = NLL real corrections to Z^{(0)}
  do_split(idipa, ka);

  // cache the event with the first insertion in;
  // to be used later for the second insertion.
  // Then evolve and retrieve the event before carrying on
  event_cache_->copy(event_);
  {
    ContinuationGuard continuation(variance_.mode == VarianceOptions::Coupled,
                                   false, true, rng, single_random_, variance_.replay_tape);
    perform_branch_single_insertion(t_insertion, -1, ka);
  }
}


//----------------------------------------------------------------------
/// singles out double insertion corrections as defined in perform_branch_double_insertion
void Shower::evolve_double_insertion_corrections(double t) {
  //if (event_.weight == 0.0 or event_.bad) return;
  if (event_.bad) return;
  include_hard_factor_ = false;

  assert(insertion_is_part_of_NLL_ensemble
         && "Double-insertion selector requires insertion_is_part_of_NLL_ensemble=true");

  // generate log(Q/kt) in [log(xQ), 1/(2 as b0)]
  // generate scale t from log(Q/kt)
  double insertion_jacobian;
  // introduce partition of unity to identify the scale of the insertion
  double lnkt_insertion = sample_insertion_log(insertion_jacobian);
  double t_insertion = t_scale(lnkt_insertion);
  if (!std::isfinite(t_insertion)) t_insertion = std::numeric_limits<double>::infinity();

  // evolve S2 at one loop kernel (without O(as) constants)
  evolve_scale(t, std::min(t_insertion, evol_cutoff_));
  int idipa = -1;
  if (t_insertion >= evol_cutoff_) {return;}

  // generate the insertion
  Momentum ka = generate_first_insertion(t_insertion, idipa);
  if (event_.bad) {return;}

  // weight associated with the first insertion
  event_.weight  = insertion_jacobian / log_scale_interval(t_insertion, tlast_);

  // split the dipole that emits ka
  do_split(idipa, ka);
  event_cache_->copy(event_);

  // branch0 = Z^{(0)} evolution at LL
  event_.weight = insertion_base_weight_;
  perform_branch0(t_insertion, ka);

  // Current selection leaves the shared LL prefix and branch zero intact.
  // Keep each complete real/counterterm pair so its replay group always ends.
  if (!DOUBLE_INSERTION_FERMION_ONLY) {
    // gauge current contribution
    event_.retrieve(event_cache_);

    // branch1 = Z^{(1)} evolution double real with kta > ktb'
    perform_branch_double_insertion(t_insertion, idipa, 1, ka);
    event_.retrieve(event_cache_);

    // branch2 = Z^{(1)} evolution collinear counterterm with kta > ktb'
    perform_branch_double_insertion(t_insertion, idipa, 2, ka);
    event_.retrieve(event_cache_);

    // branch3 = Z^{(1)} evolution double real with kta > ktb
    perform_branch_double_insertion(t_insertion, idipa, 3, ka);
    event_.retrieve(event_cache_);

    // branch4 = Z^{(1)} evolution collinear counterterm with kta > ktb
    perform_branch_double_insertion(t_insertion, idipa, 4, ka);
  }

  // fermion current contribution
  if (NF > 0. && !DOUBLE_INSERTION_GAUGE_ONLY) {
    event_.retrieve(event_cache_);
    // branch1 (fermion) = Z^{(1)} evolution double real with kta > ktb'
    perform_branch_double_insertion_fermion(t_insertion, idipa, 1, ka);
    event_.retrieve(event_cache_);

    // branch2 (fermion) = Z^{(1)} evolution collinear counterterm with kta > ktb'
    perform_branch_double_insertion_fermion(t_insertion, idipa, 2, ka);
  }
}


//----------------------------------------------------------------------
/// perform the Z0 evolution
void Shower::perform_branch0(double ta, const Momentum& ka) {
  if (event_.bad) return;
  // if emission is in observed region, add to histogram and stop the evolution
  if (insertion_is_part_of_NLL_ensemble) {
    double C1 = (include_hard_factor_ ?
                 (variance_.optimize ? hard_factor_ : 1.0 + asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1)) : 1.0);
    if (obs_->add_entries_in_region(ka.stored_E()*ka, ta, ln_kt(ta)-log_xQ(),
		  		  C1*event_.weight, &event_.axis())) return;
  }
  evolve_scale(ta, evol_cutoff_, include_hard_factor_);
  return;
}

//----------------------------------------------------------------------
/// perform the perturbative NLL insertion of the Z0 evolution
void Shower::perform_branch_single_insertion(double t_insertion,  int ibranch, const Momentum& ka) {
  // calculate weight of the insertion
  // ibranch = -1 => reals
  // ibranch = -2 => virtuals
  double L     = ln_kt(t_insertion);
  double rho   = 2.0*asmur_*b0*L;
  const double log_one_minus_rho = PRECISE_KERNELS ? std::log1p(-rho) : log(1.0-rho);
  double scale_correction = 2.0*asmur_*b0*log(xmur_/xQ_)/(1.0-rho);

  // update event weight and return
  if (!insertion_is_part_of_NLL_ensemble) {
    // Version 1: add two extra branches describing real and
    // virtual expansions of the NLL correction to the Sudakov
    // Derivation:
    // the ll evolution takes place separately in branch0,
    // so we only need to correct for the pure NLL corrections.
    // accept the emission with weight w = weight_nll,
    // and reject it with weight -w
    // The virtual corrections coming from the Sudakov are now
    // multiplying a purely NLL correction and therefore cancel
    // in the difference real - virtual up to NNLL corrections
    double w = 2. * CA * asmur_ / (2.0*M_PI) / (1.0 - rho);
    w *= - asmur_*(-(b0*KCMW) + 2.0*M_PI*b1*log_one_minus_rho)/(2.0*b0*M_PI*(1 - rho)) + scale_correction;
    // update event weight
    event_.weight *= (ibranch == -2 ? -w : w);
  } else {
    // Version 2 (slightly more efficient): 
    // correct directly the LL evolution with a probability weight,
    // and incorporate directly branch0 in the reals (branch -1)
    // this effectively adds a single branch to the existing algorithm
    // Derivation:
    // accept emission with a probability w = (weight_ll+weight_nll)/weight_ll
    // reject it (virtuals) with probability 1-w
    // In this way, the total virtual probability will be 
    // w_virt = 1 - weight_ll (from the LL Sudakov) + (1 - w) * weight_ll
    //        = 1 - (weight_ll + weight_nll)
    // Finally, include the H2 hard factor in the weight of the real to multiply
    // the sole LL contribution
    double w = (- asmur_*(-(b0*KCMW) + 2.0*M_PI*b1*log_one_minus_rho)/(2.0*b0*M_PI*(1 - rho)) + scale_correction) * event_.weight;
    if (!SINGLE_INSERTION_ONLY) {
      if (ibranch == -1) w += insertion_base_weight_ * asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1);
    }
    // reset event weight
    event_.weight = (ibranch == -2 ? -w : insertion_base_weight_ + w);
  }

  // for real corrections, check if the insertion is in the slice
  // and if so fill the histogram
  if (ibranch == -1) {
    if (obs_->add_entries_in_region(ka.stored_E()*ka, 
      t_insertion, L-log_xQ(), event_.weight, &event_.axis())) return;
  }
  // now complete the evolution until the cutoff scale
  evolve_scale(t_insertion, evol_cutoff_, false);
  return;
}

//----------------------------------------------------------------------
/// perform the Z1 evolution pieces (gauge currents)
double Shower::gauge_antenna_weight(const Momentum& a, const Momentum& b,
                                  const Momentum& c, const Momentum& d, int branch) {
  if (variance_.optimize && antenna_weight_valid_) {
    if (variance_.diagnostics) ++metrics_.antenna_cache_hits;
    return antenna_weight_cache_;
  }
  if (variance_.diagnostics) ++metrics_.antenna_evaluations;
  double weight = 0.0;
  if (branch == 1 || branch == 2) {
    weight = PRECISE_KERNELS ? double_emsn_antenna_residual(a,b,c,d)
      : double_emsn_antenna(a,b,c,d)/double_emsn_antenna_strongly_ordered(a,b,c,d)-1.0;
  }
  if (!(lab_ordering_ && dip_ordering_)) {
    // The original numerator and denominator are the same ordered antenna.
    weight += PRECISE_KERNELS ? 1.0
      : double_emsn_antenna_strongly_ordered(a,b,c,d)/double_emsn_antenna_strongly_ordered(a,b,c,d);
  }
  if (variance_.optimize) collect_antenna_weight(weight);
  return weight;
}

double Shower::fermion_antenna_weight(const Momentum& a, const Momentum& b,
                                    const Momentum& c, const Momentum& d) {
  if (variance_.optimize && antenna_weight_valid_) {
    if (variance_.diagnostics) ++metrics_.antenna_cache_hits;
    return antenna_weight_cache_;
  }
  if (variance_.diagnostics) ++metrics_.antenna_evaluations;
  double weight = PRECISE_KERNELS ? double_emsn_antenna_fermion_ratio(a,b,c,d)
    : double_emsn_antenna_fermion(a,b,c,d)/double_emsn_antenna_strongly_ordered(a,c,d,b);
  // Retain the old fermion-kernel rejection prescription, now observable in
  // diagnostics rather than silently hiding invalid evaluations.
  if (!std::isfinite(weight)) {
    if (variance_.diagnostics) ++metrics_.nonfinite_branches;
    weight = 0.0;
  }
  weight *= NF/CA;
  if (variance_.optimize) collect_antenna_weight(weight);
  return weight;
}

void Shower::perform_branch_double_insertion(double t_insertion, int idipa, int ibranch, const Momentum& ka) {
  int idipb = -1;
  Momentum kb;
  if (ibranch==1 or ibranch==3) {
    kb = generate_second_insertion(t_insertion, idipa, idipb, (ibranch==3 or ibranch==4));
    cache_second_insertion(kb, idipb, event_.bad);
  } else {
    retrieve_second_insertion(kb, idipb, event_.bad);
  }

  ContinuationGuard continuation(variance_.mode == VarianceOptions::Coupled,
                                 ibranch==1 || ibranch==3,
                                 ibranch==2 || ibranch==4, rng, double_random_, variance_.replay_tape);

  if (event_.bad) return;
  // add in the weight of the second insertion
  event_.weight*=2.0;
  const Momentum* emitter = &ka;
  const Momentum* spec_left;
  const Momentum* spec_right;
  const Momentum* emission = &kb;
  Momentum parent_storage; // lifetime covers every return path
  double w = 0.;
  event_.weight *= second_insertion_weight_;
  if (variance_.optimize && variance_.mode != VarianceOptions::Legacy && event_.weight == 0.0) {
    if (variance_.diagnostics) ++metrics_.zero_weight_skips;
    return;
  }

  bool thetaIn_ka = obs_->in_region(emitter->stored_E()*(*emitter), &event_.axis());
  bool thetaIn_kb = obs_->in_region(emission->stored_E()*(*emission), &event_.axis());

  // avoid calculating zeroes
  if (ibranch==1 or ibranch==2) {
    if (thetaIn_ka and thetaIn_kb) return;
  } else if (ibranch==3 or ibranch==4) {
    if (thetaIn_ka) return;
  }
  //if (thetaIn_ka and thetaIn_kb) return;
  // select specific configurations
  //if (thetaIn_ka and (!thetaIn_kb)) return;
  //if ((!thetaIn_ka) and thetaIn_kb) return;
  //if ((!thetaIn_ka) and (!thetaIn_kb)) return;

  if (idipb == idipa) { // emitter is to the right
    // here (1a) emitted b and split into (1b) and (ba)
    // with the index of (1a), idipa, now corresponding to
    // the new index of (1b), idipb
    // (1(b)a)(a2) : index of (1a) is idipb
    spec_left = &event_[idipb].left().momentum();
    spec_right = &event_[event_[idipb].right_neighbour()].right().momentum();
    w = gauge_antenna_weight(*spec_left, emission->stored_E()*(*emission),
                             emitter->stored_E()*(*emitter), *spec_right, ibranch);
    if (variance_.optimize && variance_.mode != VarianceOptions::Legacy && w == 0.0) {
      if (variance_.diagnostics) ++metrics_.zero_weight_skips;
      return;
    }

    // replace emitter with massless version of parent for ibranch 2
    if (ibranch==2) {
      parent_storage = emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission);
      Momentum* kab = &parent_storage;
      *kab = (1.0/kab->E())*(*kab);
      reconstruct_parent(*spec_left, *spec_right, *kab);
      emitter = kab;
      // now updated the dipoles containing ka
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
      event_[idipb].right(DipoleEnd(*kab));
      event_[event_[idipb].right_neighbour()].left(DipoleEnd(*kab));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
    }
    // split the dipole into two for branches 1 and 3
    if (ibranch==1 or ibranch==3) 
      do_split(idipb, kb);

  } else { // emitter is to the left
    // here (a2), with index event[idipa].right_neighbour,
    // emitted b, with the index of (ab) being idipb
    // (1a)(a(b)2) : index of (a2) is idipb
    spec_right = &event_[idipb].right().momentum();
    spec_left = &event_[event_[idipb].left_neighbour()].left().momentum();
    w = gauge_antenna_weight(*spec_left, emitter->stored_E()*(*emitter),
                             emission->stored_E()*(*emission), *spec_right, ibranch);
    if (variance_.optimize && variance_.mode != VarianceOptions::Legacy && w == 0.0) {
      if (variance_.diagnostics) ++metrics_.zero_weight_skips;
      return;
    }

    // replace emitter with massless version of parent for ibranch 2
    if (ibranch==2) {
      parent_storage = emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission);
      Momentum* kab = &parent_storage;
      *kab = (1.0/kab->E())*(*kab);
      reconstruct_parent(*spec_left, *spec_right, *kab);
      emitter = kab;
      // now updated the dipoles containing ka
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
      event_[idipb].left(DipoleEnd(*kab));
      event_[event_[idipb].left_neighbour()].right(DipoleEnd(*kab));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
    }
    // split the dipole into two for branches 1 and 3
    if (ibranch==1 or ibranch==3) {
      do_split(idipb, kb);
      // update spec_right with the new dipole indices. Only if
      // emitter is to the left (new dipole is added at the end of the chain)
      spec_right = &event_[event_[idipb].right_neighbour()].right().momentum();
    }
  }

  // flip the sign of the weight for branches 2 and 3
  event_.weight *= ((ibranch==1 or ibranch==4) ? w : -w);
  if (isnan(event_.weight) or isinf(event_.weight)) {
    if (variance_.diagnostics) ++metrics_.nonfinite_branches;
    event_.bad = true;
    return;
  }

  // now check if any emission is in the slice and fill the histogram accordingly
  if (ibranch==1 or ibranch==3) {
    // both emissions ka and kb are inside the slice 
    if (thetaIn_ka and thetaIn_kb) {
      if (ibranch==1) {
        // ==> bin the parent (defined in the massless scheme @ NLL)
        Momentum kab = emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission);
        kab = (1.0/kab.E())*kab;
        reconstruct_parent(*spec_left, *spec_right, kab);
        // Dipole-frame massless reconstruction can slightly modify lab rapidity.
        if (obs_->add_entries_in_region(kab.stored_E()*kab, t_insertion,
                                       ln_kt(t_insertion) - log_xQ(), event_.weight, &event_.axis())) return;
      } else {
        // ==> bin the emitter (much harder than the emission in LL kinematics)
	obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
					   ln_kt(t_insertion) - log_xQ(), event_.weight, &event_.axis());
	return;
      }
    // only one of the two emissions ka and kb is inside the slice 
    } else if (thetaIn_ka or thetaIn_kb) {
      obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
                                 ln_kt(t_insertion) - log_xQ(), event_.weight, &event_.axis());
      obs_->add_entries_in_region(emission->stored_E()*(*emission), t_insertion,
                                 ln_kt(t_insertion) - log_xQ(), event_.weight, &event_.axis());
      return;
    }

  // in branches 2 and 4 bin the emitter  
  } else if (ibranch==2 or ibranch==4) {
    if (obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
                                   ln_kt(t_insertion) - log_xQ(), event_.weight, &event_.axis())) return;
  }

  // now complete the evolution until the cutoff scale
  evolve_scale(t_insertion, evol_cutoff_, false);
  // equivalently (NNLL difference) one can start from the scale of the second insertion
  //evolve_scale(t_second_insertion_, evol_cutoff_, false);
}

//----------------------------------------------------------------------
/// perform the Z1 evolution pieces (fermion currents)
void Shower::perform_branch_double_insertion_fermion(double t_insertion, int idipa, int ibranch, const Momentum& ka) {
  int idipb = -1;
  Momentum kb;

  if (!(ibranch == 1 || ibranch == 2)) {
    event_.bad = true;
    return;
  }
  //>> identify the emission either with the quark (idipb == idipa) or
  //>> with the antiquark (idipb != idipa). Since the squared
  //>> amplitude is fully symmetric this is purely conventional
  if (ibranch==1) {
    kb = generate_second_insertion(t_insertion, idipa, idipb, (ibranch==3 or ibranch==4));
    cache_second_insertion(kb, idipb, event_.bad);
  } else {
    retrieve_second_insertion(kb, idipb, event_.bad);
  }

  ContinuationGuard continuation(variance_.mode == VarianceOptions::Coupled,
                                 ibranch==1, ibranch==2, rng, double_random_, variance_.replay_tape);

  if (event_.bad) return;
  // add in the weight of the second insertion
  event_.weight*=2.0;
  const Momentum* emitter = &ka;
  const Momentum* spec_left;
  const Momentum* spec_right;
  const Momentum* emission = &kb;
  Momentum parent_storage; // lifetime covers every return path
  double w = 0.;
  event_.weight *= second_insertion_weight_;
  if (variance_.optimize && variance_.mode != VarianceOptions::Legacy && event_.weight == 0.0) {
    if (variance_.diagnostics) ++metrics_.zero_weight_skips;
    return;
  }

  bool thetaIn_ka = obs_->in_region(emitter->stored_E()*(*emitter), &event_.axis());
  bool thetaIn_kb = obs_->in_region(emission->stored_E()*(*emission), &event_.axis());

  // avoid calculating zeroes
  if (ibranch==1 or ibranch==2) {
    if (thetaIn_ka and thetaIn_kb) return;
  } else if (ibranch==3 or ibranch==4) {
    if (thetaIn_ka) return;
  }
  //if (thetaIn_ka and thetaIn_kb) return;
  // select specific configurations
  //if (thetaIn_ka and (!thetaIn_kb)) return;
  //if ((!thetaIn_ka) and thetaIn_kb) return;
  //if ((!thetaIn_ka) and (!thetaIn_kb)) return;

  if (idipb == idipa) { // emitter is to the right
    // here (1a) emitted b and split into (1b) and (ba)
    // with the index of (1a), idipa, now corresponding to
    // the new index of (1b), idipb
    // (1(b)a)(a2) : index of (1a) is idipb
    spec_left = &event_[idipb].left().momentum();
    spec_right = &event_[event_[idipb].right_neighbour()].right().momentum();
    w = fermion_antenna_weight(*spec_left, *spec_right,
                               emission->stored_E()*(*emission), emitter->stored_E()*(*emitter));
    // to bring the normalisation of the NF term consistent with the SO CA term
    // one needs to multiply the former by 2 (see eq. 2.4 of 9707532). 
    // Since we are generating two identical NF contributions 
    // (i.e. sum over two colour dipoles), we further need to multiply by 1/2.
    if (variance_.optimize && variance_.mode != VarianceOptions::Legacy && w == 0.0) {
      if (variance_.diagnostics) ++metrics_.zero_weight_skips;
      return;
    }
    
    // replace emitter with massless version of parent for ibranch 2
    if (ibranch==2) {
      parent_storage = emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission);
      Momentum* kab = &parent_storage;
      *kab = (1.0/kab->E())*(*kab);
      reconstruct_parent(*spec_left, *spec_right, *kab);
      emitter = kab;
      // now updated the dipoles containing ka
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
      event_[idipb].right(DipoleEnd(*kab));
      event_[event_[idipb].right_neighbour()].left(DipoleEnd(*kab));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
    }
    // split the dipole into two for branches 1 and 3
    if (ibranch==1) {
      // split dipole without adding a new one
      // >> emitter is to the right => ff dipole is at the end of the chain
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
      event_[idipb].right(DipoleEnd(*emission));
      event_[event_[idipb].right_neighbour()].left(DipoleEnd(*emitter));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
    }

  } else { // emitter is to the left
    // here (a2), with index event[idipa].right_neighbour,
    // emitted b, with the index of (ab) being idipb
    // (1a)(a(b)2) : index of (a2) is idipb
    spec_right = &event_[idipb].right().momentum();
    spec_left = &event_[event_[idipb].left_neighbour()].left().momentum();
    w = fermion_antenna_weight(*spec_left, *spec_right,
                               emitter->stored_E()*(*emitter), emission->stored_E()*(*emission));
    // to bring the normalisation of the NF term consistent with the SO CA term
    // one needs to multiply the former by 2 (see eq. 2.4 of 9707532). 
    // Since we are generating two identical NF contributions 
    // (i.e. sum over two colour dipoles), we further need to multiply by 1/2.
    if (variance_.optimize && variance_.mode != VarianceOptions::Legacy && w == 0.0) {
      if (variance_.diagnostics) ++metrics_.zero_weight_skips;
      return;
    }

    // replace emitter with massless version of parent for ibranch 2
    if (ibranch==2) {
      parent_storage = emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission);
      Momentum* kab = &parent_storage;
      *kab = (1.0/kab->E())*(*kab);
      reconstruct_parent(*spec_left, *spec_right, *kab);
      emitter = kab;
      // now updated the dipoles containing ka
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
      event_[idipb].left(DipoleEnd(*kab));
      event_[event_[idipb].left_neighbour()].right(DipoleEnd(*kab));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
    }
    // split the dipole into two for branches 1 and 3
    if (ibranch==1) {
      // split dipole without adding a new one
      // >> emitter is to the left => ff dipole is at idipb
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
      event_[idipb].left(DipoleEnd(*emission));
      event_[event_[idipb].left_neighbour()].right(DipoleEnd(*emitter));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
      // no need to update spec_right with the new dipole indices 
      // since no new dipole was created
    }
  }

  // flip the sign of the weight for branches 2 and 3
  event_.weight *= ((ibranch==1) ? w : -w);
  if (isnan(event_.weight) or isinf(event_.weight)) {
    if (variance_.diagnostics) ++metrics_.nonfinite_branches;
    event_.bad = true;
    return;
  }

  // now check if any emission is in the slice and fill the histogram accordingly
  if (ibranch==1) {
    // both emissions ka and kb are inside the slice 
    if (thetaIn_ka and thetaIn_kb) {
      // ==> bin the parent (defined in the massless scheme @ NLL)
      Momentum kab = emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission);
      kab = (1.0/kab.E())*kab;
      reconstruct_parent(*spec_left, *spec_right, kab);
      // Dipole-frame massless reconstruction can slightly modify lab rapidity.
      if (obs_->add_entries_in_region(kab.stored_E()*kab, t_insertion,
                                     ln_kt(t_insertion) - log_xQ(), event_.weight, &event_.axis())) return;
    // only one of the two emissions ka and kb is inside the slice 
    } else if (thetaIn_ka or thetaIn_kb) {
      obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
                                 ln_kt(t_insertion) - log_xQ(), event_.weight, &event_.axis());
      obs_->add_entries_in_region(emission->stored_E()*(*emission), t_insertion,
                                 ln_kt(t_insertion) - log_xQ(), event_.weight, &event_.axis());
      return;
    }

  // in branches 2 and 4 bin the emitter  
  } else if (ibranch==2) {
    if (obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
                                   ln_kt(t_insertion) - log_xQ(), event_.weight, &event_.axis())) return;
  }

  // now complete the evolution until the cutoff scale
  evolve_scale(t_insertion, evol_cutoff_, false);
  // equivalently (NNLL difference) one can start from the scale of the second insertion
  //evolve_scale(t_second_insertion_, evol_cutoff_, false);
}

//----------------------------------------------------------------------
/// generate first insertion
Momentum Shower::generate_first_insertion(double& t_insertion, int& idip_insertion) {
  if (insertion_is_part_of_NLL_ensemble) t_insertion += - log(rng.uniform_pos()) / (2.0 * CA * event_.eta_tot);
  double lnkt = ln_kt(t_insertion);

  if (event_.bad or (2.0*asmur_*b0*lnkt >= 1.0) or (t_insertion >= evolution_limit_) or (lnkt > lnktmax)) {
    event_.bad=true;
    return Momentum();
  }
  idip_insertion = choose_emitter();
  assert(idip_insertion>=0);
  
  Momentum insertion = event_[idip_insertion].radiate(lnkt - log_xQ(), rng.uniform_pos(),
					    rng.uniform_pos());
  return insertion;
}

//----------------------------------------------------------------------
/// generate second insertion
Momentum Shower::generate_second_insertion(double t_insertion, int idip, int& idip_insertion, bool dipole_kt_ordering) {
  if (rng.uniform() < 0.5) {
    //idip_insertion = event_[idip].delta_rap()!=0.0 ? idip : event_[idip].right_neighbour();
    // Alternative: set weight to zero if selected dipole has no phase space for radiating
    if (event_[idip].delta_rap() != 0.0) {
      idip_insertion = idip;
      // emitter in to the right - update its rapidity bounds
      event_.eta_tot -= event_[idip_insertion].delta_rap();
      event_[idip_insertion].reset(RAPMAX,-10.);
      event_.eta_tot += event_[idip_insertion].delta_rap();
    } else {
      event_.bad = true;
      return Momentum();
    }
  } else {
    assert(event_[idip].right_neighbour()>=0);
    //idip_insertion = (event_[event_[idip].right_neighbour()].delta_rap() != 0.0) ?
    //event_[idip].right_neighbour() : idip;
    // Alternative: set weight to zero if selected dipole has no phase space for radiating
    if (event_[event_[idip].right_neighbour()].delta_rap() != 0.0) {
      idip_insertion = event_[idip].right_neighbour();
      // emitter in to the left - update its rapidity bounds
      event_.eta_tot -= event_[idip_insertion].delta_rap();
      event_[idip_insertion].reset(10.,-RAPMAX);
      event_.eta_tot += event_[idip_insertion].delta_rap();
    } else {
      event_.bad = true;
      return Momentum();
    }  
  }

  // generate insertion with unit dipole transverse momentum
  Momentum insertion = event_[idip_insertion].radiate(0.0, rng.uniform_pos(), rng.uniform_pos());
  double f2 = 1.0;
  int rn = event_[idip].right_neighbour();
  Momentum emitter = event_[idip].right().momentum();
  if (!dipole_kt_ordering) 
    f2 = 0.5*dot_product(event_[idip].left().momentum(),event_[rn].right().momentum())
      / dot_product(event_[idip].left().momentum(),insertion.stored_E()*insertion)
      / dot_product(insertion.stored_E()*insertion,event_[rn].right().momentum());
  double ln_buffer = 0.5*log(f2);
  double tb = t_scale(ln_kt(t_insertion) - ln_buffer);

  // Cut events with a t above the starting scale t=0.
  // Moreover, the starting tb will be nan if 
  // ln_kt(t_insertion)-ln_buffer is above the landau pole
  // (this effectively acts as a collinear cutoff)
  double lnkt;
  if (!weighted_second_insertion) {
    // generate second insertion with a Sudakov (unweighted)
    while (true) {
      tb += - log(rng.uniform_pos()) / (2.0 * CA * event_[idip_insertion].delta_rap());
      if ((tb > 0.) or (tb != tb)) break;
    }
    lnkt = ln_kt(tb);
    second_insertion_weight_ = 1.;
    // reweigh scale of coupling to that of first emission
    double rho_second = 2.*asmur_*b0*lnkt;
    double rho_first  = 2.*asmur_*b0*ln_kt(t_insertion);
    second_insertion_weight_ *= (1. - rho_second)/(1. - rho_first);
  } else {
    // generate second insertion at fixed order (weighted)
    const double lower = ln_kt(t_insertion) - ln_buffer;
    const double range = 1.0/(2.0*asmur_*b0) - lower;
    const double u = rng.uniform_pos();
    const double power = variance_.second_power;
    const double mapped = power == 1.0 ? u : pow(u, power);
    lnkt = mapped*range + lower;
    const double second_jacobian = power == 1.0 ? range : power*pow(u, power-1)*range;

    // calculate weight
    double rho = 2.*asmur_*b0*ln_kt(t_insertion); // scale of first insertion (NNLL difference)
    //double rho = 0; // fixed coupling (for fixed order tests)

    second_insertion_weight_  = second_jacobian;
    second_insertion_weight_ *= 2. * CA * asmur_ / (2.0*M_PI) / (1.0 - rho) * event_[idip_insertion].delta_rap();
    tb = t_scale(lnkt);
  }
  // update emission's energy
  insertion.stored_E(insertion.stored_E()*exp(-lnkt+log_xQ()));

  if (isnan(lnkt) or isinf(lnkt) or (2.0*asmur_*b0*lnkt >= 1.0)
      or (tb >= evolution_limit_) or (lnkt > lnktmax)) {
    event_.bad = true;
    return Momentum();
  }

  // check kt ordering in parent dipole frame
  double kta2 = 2./dot_product(event_[idip].left().momentum(),event_[rn].right().momentum())
      * dot_product(event_[idip].left().momentum(),emitter.stored_E() * emitter)
      * dot_product(emitter.stored_E() * emitter,event_[rn].right().momentum());
  double ktb2 = 2./dot_product(event_[idip].left().momentum(),event_[rn].right().momentum())
      * dot_product(event_[idip].left().momentum(),insertion.stored_E()*insertion)
      * dot_product(insertion.stored_E()*insertion,event_[rn].right().momentum());

  if (kta2 > ktb2) {
    lab_ordering_ = true;
  } else {
    lab_ordering_ = false;
  }

  t_second_insertion_ = tb; // define emission's time in the emitting dipole frame
  if (isnan(t_second_insertion_) or isinf(t_second_insertion_)) {
    t_second_insertion_ = EVOLCUT;
    event_.bad = true;
    return Momentum();
  }

  if (ln_kt(t_insertion) < ln_kt(t_second_insertion_)) {
    dip_ordering_ = true;
  } else {
    dip_ordering_ = false;
  }

  return insertion;
}

//----------------------------------------------------------------------
/// return the t scale for a given ln kt value
double Shower::t_scale(double lnkt) const {  
  if (PRECISE_KERNELS) return -std::log1p(-2*asmur_*b0*lnkt)/(4.*b0*M_PI);
  // fixed coupling:
  // at LL:
  //   t = Log[1/(2 - 2 b0 alphas L)]/(4 pi b0)
  return -log(1 - 2*asmur_*b0*lnkt)/(4.*b0*M_PI);
  // fixed coupling: return as*lnkt/(2.0*M_PI);
  //return asmur_/(2.*M_PI)*lnkt;
}

//----------------------------------------------------------------------
/// return the ln(kt) for a given evolution scale
double Shower::ln_kt(double t) const {
  if (PRECISE_KERNELS) return -std::expm1(-4*M_PI*t*b0)/(2.0*asmur_*b0);
  // fixed coupling
  /// at LL:
  ///    t = Log[1/(1 - 2 b0 alphas L)]/(4 pi b0)
  /// => ln(kt) = (1 - exp(-4 pi t b0))/(2 alphas b0)
  return (1 - std::exp(-4*M_PI*t*b0))/(2.0*asmur_*b0);
  // fixed coupling: return 2.0*M_PI*t/as;
  //return (2.*M_PI)/asmur_*t;
}

double Shower::log_scale_interval(double later, double earlier) const {
  if (!PRECISE_KERNELS) return ln_kt(later)-ln_kt(earlier);
  return -std::exp(-4*M_PI*b0*earlier)*std::expm1(-4*M_PI*b0*(later-earlier))
    /(2.0*asmur_*b0);
}

//----------------------------------------------------------------------
/// write current observable to file (or to cout)
void Shower::write(int nev, const std::string& fn) const {
  std::ofstream file;
  std::ostream * ostr = &std::cout;
  if (!fn.empty()) {
    file.open(fn);
    if (!file) throw std::runtime_error("Cannot open output file: " + fn);
    ostr = &file;
  }
  *ostr << header_ << std::endl;
  *ostr << std::setprecision(variance_.output_precision);
  if (variance_.diagnostics) {
    *ostr << "# emitter_calls = " << metrics_.emitter_calls << std::endl;
    *ostr << "# mean_dipoles = " << (metrics_.emitter_calls ?
      double(metrics_.emitter_dipoles)/metrics_.emitter_calls : 0.0) << std::endl;
    *ostr << "# maximum_dipoles = " << metrics_.maximum_dipoles << std::endl;
    *ostr << "# zero_weight_skips = " << metrics_.zero_weight_skips << std::endl;
    *ostr << "# antenna_evaluations = " << metrics_.antenna_evaluations << std::endl;
    *ostr << "# antenna_cache_hits = " << metrics_.antenna_cache_hits << std::endl;
    *ostr << "# nonfinite_weight_rejections = " << metrics_.nonfinite_branches << std::endl;
    *ostr << "# rapidity_sum_rebuilds = " << metrics_.rapidity_sum_rebuilds << std::endl;
    *ostr << "# rng_generator_draws = " << rng.generator_uniform_draws() << std::endl;
    *ostr << "# rng_snapshot_bytes = " << rng.snapshot_bytes_copied() << std::endl;
    *ostr << "# rng_tape_capacity = " << rng.replay_tape_capacity() << std::endl;
  }
  if (variance_.adaptive_insertion) {
    *ostr << "# adaptive_updates = " << insertion_sampler_.updates() << std::endl;
    *ostr << "# adaptive_probabilities =";
    for (double probability : insertion_sampler_.probabilities()) *ostr << " " << probability;
    *ostr << std::endl;
  }
  *ostr << "# output with nev = " << nev << std::endl;
  obs_->write(nev, ostr);
  if (!*ostr) throw std::runtime_error("Failed writing shower output: " + fn);
}

//----------------------------------------------------------------------
/// reconstruct the parent kab such that it is massless
void Shower::reconstruct_parent(const Momentum& spec_left, const Momentum& spec_right, Momentum& kab) const {
  // Create a massless parent in the dipole frame, aligned along the z axis.
  Momentum k12 = spec_left+spec_right;
  Momentum k1 = spec_left;
  kab *= kab.stored_E();
  kab.unboost(k12);
  k1.unboost(k12);
  double theta=k1.theta();
  double phi=k1.phi();
  kab.rotate(theta, phi);
  double nab = kab.rap();
  double ktab = sqrt(kab.px()*kab.px()+kab.py()*kab.py());
  kab = Momentum(kab.px(), kab.py(), ktab*sinh(nab), ktab*cosh(nab));
  kab.unrotate(theta,phi);
  kab.boost(k12);
  kab.stored_E(kab.E());
  kab = (1.0/kab.E())*kab;
  // Expanded measurements use the first-insertion time, not a parent time.
  return;
}
