#include "Shower.hh"
#include "Parameters.hh"
#include "Dipole.hh"
#include "helpers.hh"
#include <cmath>

//----------------------------------------------------------------------
/// run the shower with nev events
void Shower::run(int nev, const std::string& fn) {
  for(int i=0; i<nev; ++i) {
    reset();
    double tstart = 0.0;
    if (periodic_output()) {
      std::cerr << "# " << i+1 << " out of " << nev << " events." << std::endl;
      if (!fn.empty()) write(i+1, fn);
    }
    if (order_evl_ == 0) 
      evolve_scale(tstart);
    else if (order_evl_ == 1) {
      if (!NLL_EXPANDED){
        // run the soft factor S2 with two-loop evolution convoluted with H2 at NLO
        evolve_insertion(tstart);
      } else {
        // run the soft factor S2 with one-loop evolution convoluted with H2 at NLO
        // followed by the two loop corrections to S2 convoluted with H2 at LO
        evolve_insertion_expanded(tstart);
      }
      // run the soft factor S3 with one-loop evolution convoluted with H3 at LO
      run_threejet(tstart, true);
      run_threejet(tstart, false, true);
    }
  }
  write(nev, fn);
}

//----------------------------------------------------------------------
/// run the threejet piece
void Shower::run_threejet(double tstart, bool soft, bool use_cached_variables) {
  reset(true, soft, use_cached_variables);
  // flip sign for LL term
  if (soft) event_.weight = -event_.weight;
  evolve_scale(tstart);
  double w = -event_.weight;
  reset(false);
  event_.weight = w;
  // now veto configurations where gluon is in the slice as computed in 
  // the two jet limit within the integrated counterterm
  if (obs_->in_region(gluon_, &event_.axis())) event_.bad=true;
  evolve_scale(tstart);
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
    Shower::do_split(0,gluon_);
    if (std::abs(gluon_.rap()) > RAPMAX) event_.bad = true;
    in_region = obs_->in_region(gluon_, &event_.axis());
    if (in_region or event_.weight==0.0) event_.bad = true;
  } else {
    bool in_region = false;
    event_.reset_threejet(xmur_, xQ_, rng.uniform(), rng.uniform(), rng.uniform(), gluon_, use_cached_variables);
    Shower::do_split(0,gluon_);

    // make sure we don't have any hard partons in the observed region, otherwise set weight to 0
    in_region = obs_->in_region(event_[0].left().momentum(), &event_.axis());
    in_region = in_region || (obs_->in_region(event_[0].right().momentum(), &event_.axis()));
    in_region = in_region || (obs_->in_region(event_[1].right().momentum(), &event_.axis()));

    if (in_region or event_.weight==0.0) event_.bad = true;
  }
  return;
}

//----------------------------------------------------------------------
/// choose an emitting dipole
int Shower::choose_emitter() const {
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
void Shower::evolve_scale(double t, double tend, bool include_as_constant) {
  //if (event_.weight == 0.0 or event_.bad) return;
  if (event_.bad) return;
  tlast_ = t;

  while ((t += - log(rng.uniform_pos()) / (2.0 * CA * event_.eta_tot)) < tend) {
    // increase by amount chosen with distribution e^(-2*CA*sum_rap)    
    int idip = choose_emitter();
    // now set all emission properties
    double lnkt = ln_kt(t);
    if ((2.0*asmur_*b0*lnkt >= 1.0) or (evl_grid_ and t>=evl_grid_->xlim()) or (lnkt > lnktmax)) {
      event_.bad = true; // setting this to avoid starting another evolution later
      break;
    }

    //fixed coupling: if (event_.eta_tot < 0.0 or event_.eta_tot!=event_.eta_tot) break;
    assert(idip>=0);
    // the momentum of the gluon is xQ*exp(-lnkt)=exp[-(lnkt+ln(xQ))]

    Momentum emsn = event_[idip].radiate(lnkt - log(xQ_), rng.uniform_pos(),
					 rng.uniform_pos());
    if (!do_split(idip, emsn)) continue;

    tlast_ = t;
    // if emission is in observed region, add to histogram and stop the evolution
    double C1 = ((include_as_constant and NLL_evolution_) ?  1.0 + asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1) : 1.0);
    if (obs_->add_entries_in_region(emsn.stored_E()*emsn, t, lnkt-log(xQ_),
	  		    C1 * event_.weight, &event_.axis())) {
      event_.bad = true; // setting this to avoid starting another evolution later
      break;
    }
  }
  return;
}

//----------------------------------------------------------------------
/// evolve the shower down to cutoff scale
void Shower::evolve_insertion(double t) {
  //if (event_.weight == 0.0 or event_.bad) return;
  if (event_.bad) return;
  NLL_evolution_ = true;

  // generate log(Q/kt) in [log(xQ), 1/(2 as b0)]
  // generate scale t from log(Q/kt)
  double r = rng.uniform();
  double n = 1.;
  double lnkt_insertion = log(xQ_) + pow(r,n)*landau_pole_tolerance_*(1.0/(2.0*asmur_*b0) - log(xQ_));
  double t_insertion = t_scale(lnkt_insertion);
  
  evolve_scale(t, std::min(t_insertion, evol_cutoff_));
  int idipa = -1;
  if (t_insertion >= evol_cutoff_) {
    NLL_evolution_ = false;
    return;
  }
  Momentum ka = generate_first_insertion(t_insertion, idipa);
  if (event_.bad) {
    NLL_evolution_ = false;
    return;
  }
  assert(do_split(idipa, ka));

  // jacobian associated to the insertion
  double w = n*pow(r,n-1) * landau_pole_tolerance_*(1.0/(2.0*asmur_*b0) - log(xQ_)) / (ln_kt(t_insertion) - ln_kt(tlast_));
  event_cache_->copy(event_);
  
  // branch0 = Z^{(0)} evolution
  perform_branch0(t_insertion, ka);
  event_.retrieve(event_cache_);

  // branch1 = Z^{(1)} evolution double real with kta > ktb'
  event_.weight = w;
  perform_branch_double_insertion(t_insertion, idipa, 1, ka);
  event_.retrieve(event_cache_);

  // branch2 = Z^{(1)} evolution collinear counterterm with kta > ktb'
  event_.weight = w;
  perform_branch_double_insertion(t_insertion, idipa, 2, ka);
  event_.retrieve(event_cache_);
  
  // branch3 = Z^{(1)} evolution double real with kta > ktb
  event_.weight = w;
  perform_branch_double_insertion(t_insertion, idipa, 3, ka);
  event_.retrieve(event_cache_);
  
  // branch4 = Z^{(1)} evolution collinear counterterm with kta > ktb
  event_.weight = w;
  perform_branch_double_insertion(t_insertion, idipa, 4, ka);
  NLL_evolution_ = false;
}


//----------------------------------------------------------------------
/// evolve the shower down to cutoff scale while fully expanding 
/// out NLL corrections
void Shower::evolve_insertion_expanded(double t) {  
  //if (event_.weight == 0.0 or event_.bad) return;
  if (event_.bad) return;
  NLL_evolution_ = true;

  // generate log(Q/kt) in [log(xQ), 1/(2 as b0)]
  // generate scale t from log(Q/kt)
  double r = rng.uniform();
  double n = 1.;
  // introduce partition of unity to identify the scale of the insertion
  double lnkt_insertion = log(xQ_) + pow(r,n)*landau_pole_tolerance_*(1.0/(2.0*asmur_*b0) - log(xQ_));
  double t_insertion = t_scale(lnkt_insertion);

  // evolve S2 at one loop kernel and include H2 at one loop
  evolve_scale(t, std::min(t_insertion, evol_cutoff_));
  int idipa = -1;
  if (t_insertion >= evol_cutoff_) {
    NLL_evolution_ = false;
    return;
  }
  // generate the insertion
  Momentum ka = generate_first_insertion(t_insertion, idipa);
  if (event_.bad) {
    NLL_evolution_ = false;
    return;
  }

  // branch0 = Z^{(0)} evolution at LL
  if (!insertion_is_part_of_NLL_ensemble) {
    event_cache_->copy(event_);
    perform_branch0(t_insertion, ka);
    event_.retrieve(event_cache_);
  }

  // now handle extra NLL branches
  if (insertion_is_part_of_NLL_ensemble) {
    event_.weight  = n*pow(r,n-1) * landau_pole_tolerance_*(1.0/(2.0*asmur_*b0) - log(xQ_)) / (ln_kt(t_insertion) - ln_kt(tlast_));
  } else {
    event_.weight  = n*pow(r,n-1) * landau_pole_tolerance_*(1.0/(2.0*asmur_*b0) - log(xQ_)) * event_.eta_tot;
  }

  // branch -2 = NLL virtual corrections to Z^{(0)}
  event_cache_->copy(event_);
  perform_branch_single_insertion(t_insertion, -2, ka);
  event_.retrieve(event_cache_);

  // split the dipole that emits ka
  // branch -1 = NLL real corrections to Z^{(0)}
  assert(do_split(idipa, ka));

  // cache the event with the first insertion in;
  // to be used later for the second insertion.  
  // Then evolve and retrieve the event before carrying on
  event_cache_->copy(event_);
  perform_branch_single_insertion(t_insertion, -1, ka);
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
  NLL_evolution_ = false;
}


//----------------------------------------------------------------------
/// perform the Z0 evolution
void Shower::perform_branch0(double ta, const Momentum& ka) {
  if (event_.bad) return;
  // if emission is in observed region, add to histogram and stop the evolution
  if (!(NLL_EXPANDED and (!insertion_is_part_of_NLL_ensemble))) {
    double C1 = (NLL_evolution_ ?  1.0 + asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1) : 1.0);
    if (obs_->add_entries_in_region(ka.stored_E()*ka, ta, ln_kt(ta)-log(xQ_),
		  		  C1*event_.weight, &event_.axis())) return;
  }
  evolve_scale(ta, evol_cutoff_);
}

//----------------------------------------------------------------------
/// perform the perturbative NLL insertion of the Z0 evolution (when NLL_EXPANDED=true)
void Shower::perform_branch_single_insertion(double t_insertion,  int ibranch, const Momentum& ka) {
  // calculate weight of the insertion
  // ibranch = -1 => reals
  // ibranch = -2 => virtuals
  double L     = ln_kt(t_insertion);
  double rho   = 2.0*asmur_*b0*L;

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
    w *= - asmur_*(-(b0*KCMW) + 2.0*M_PI*b1*log(1.0 - rho))/(2.0*b0*M_PI*(1 - rho));
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
    double w = - asmur_*(-(b0*KCMW) + 2.0*M_PI*b1*log(1.0 - rho))/(2.0*b0*M_PI*(1 - rho)) * event_.weight;
    if (ibranch == -1) w += asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1);
    // reset event weight
    event_.weight = (ibranch == -2 ? -w : 1. + w);
  }

  // for real corrections, check if the insertion is in the slice
  // and if so fill the histogram
  if (ibranch == -1) {
    if (obs_->add_entries_in_region(ka.stored_E()*ka, 
      t_insertion, L-log(xQ_), event_.weight, &event_.axis())) return;
  }
  // now complete the evolution until the cutoff scale
  evolve_scale(t_insertion, evol_cutoff_, !NLL_EXPANDED);
  return;
}

//----------------------------------------------------------------------
/// perform the Z1 evolution pieces
void Shower::perform_branch_double_insertion(double t_insertion, int idipa, int ibranch, const Momentum& ka) {
  int idipb = -1;
  Momentum kb;
  // branches 3 & 4 cancel with the SO part of branches 1 & 2 up to subleading corrections
  if (ibranch==3 or ibranch==4) return;
  if (ibranch==1 or ibranch==3) {
    kb = generate_second_insertion(t_insertion, idipa, idipb, (ibranch==3 or ibranch==4));
    cache_second_insertion(kb, idipb, event_.bad);
  } else {
    retrieve_second_insertion(kb, idipb, event_.bad);
  }

  if (event_.bad) return;
  // add in the weight of the second insertion
  event_.weight*=2.0;
  const Momentum* emitter = &ka;
  const Momentum* spec_left;
  const Momentum* spec_right;
  const Momentum* emission = &kb;
  double tab;
  double w = 1.0;

  if (idipb == idipa) { // emitter is to the right
    // here (1a) emitted b and split into (1b) and (ba)
    // with the index of (1a), idipa, now corresponding to
    // the new index of (1b), idipb
    // (1(b)a)(a2) : index of (1a) is idipb
    spec_left = &event_[idipb].left().momentum();
    spec_right = &event_[event_[idipb].right_neighbour()].right().momentum();
    if (ibranch==1 or ibranch==2) {
      w = double_emsn_antenna(*spec_left, (emission->stored_E())*(*emission),
			      (emitter->stored_E())*(*emitter), *spec_right)
	          /double_emsn_antenna_strongly_ordered(*spec_left, emission->stored_E()*(*emission),
					  emitter->stored_E()*(*emitter), *spec_right);
      w -= 1.; // remove SO limit
    }
    // replace emitter with massless version of parent for ibranch 2
    if (ibranch==2) {
      Momentum* kab = new Momentum(emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission));
      *kab = (1.0/kab->E())*(*kab);
      reconstruct_parent(*spec_left, *spec_right, *kab, tab);
      emitter = kab;
      // now updated the dipoles containing ka
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
      event_[idipb].right(DipoleEnd(*kab));
      event_[event_[idipb].right_neighbour()].left(DipoleEnd(*kab));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
    }
    // split the dipole into two for branches 1 and 3
    if (ibranch==1 or ibranch==3) 
      assert(do_split(idipb, kb));

  } else { // emitter is to the left
    // here (a2), with index event[idipa].right_neighbour,
    // emitted b, with the index of (ab) being idipb
    // (1a)(a(b)2) : index of (a2) is idipb
    spec_right = &event_[idipb].right().momentum();
    spec_left = &event_[event_[idipb].left_neighbour()].left().momentum();
    if (ibranch==1 or ibranch==2) {
      w = double_emsn_antenna(*spec_left, (emitter->stored_E())*(*emitter),
			      (emission->stored_E())*(*emission), *spec_right)
	          /double_emsn_antenna_strongly_ordered(*spec_left, emitter->stored_E()*(*emitter),
					  emission->stored_E()*(*emission), *spec_right);
      w -= 1.; // remove SO limit
    }
    // replace emitter with massless version of parent for ibranch 2
    if (ibranch==2) {
      Momentum* kab = new Momentum(emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission));
      *kab = (1.0/kab->E())*(*kab);
      reconstruct_parent(*spec_left, *spec_right, *kab, tab);
      emitter = kab;
      // now updated the dipoles containing ka
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
      event_[idipb].left(DipoleEnd(*kab));
      event_[event_[idipb].left_neighbour()].right(DipoleEnd(*kab));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
    }
    // split the dipole into two for branches 1 and 3
    if (ibranch==1 or ibranch==3) {
      assert(do_split(idipb, kb));
      // update spec_right with the new dipole indices. Only if
      // emitter is to the left (new dipole is added at the end of the chain)
      spec_right = &event_[event_.size()-1].right().momentum();
    }
  }

  // flip the sign of the weight for branches 2 and 3
  event_.weight *= ((ibranch==1 or ibranch==4) ? w : -w);
  if (event_.weight != event_.weight) {
    event_.bad = true;
    return;
  }

  // now check if any emission is in the slice and fill the histogram accordingly
  double C1 = (((!NLL_EXPANDED) and NLL_evolution_) ?  1.0 + asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1) : 1.0);
  if (ibranch==1 or ibranch==3) {
    bool thetaIn_ka = obs_->in_region(emitter->stored_E()*(*emitter), &event_.axis());
    bool thetaIn_kb = obs_->in_region(emission->stored_E()*(*emission), &event_.axis());

    // both emissions ka and kb are inside the slice 
    if (thetaIn_ka and thetaIn_kb) {
      if (ibranch==1) {
        // ==> bin the parent (defined in the massless scheme @ NLL)
        Momentum kab = emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission);
        kab = (1.0/kab.E())*kab;
        reconstruct_parent(*spec_left, *spec_right, kab, tab);
        // sanity check: the following assert is only satisfied with Option 1 & 3 in reconstruct_parent
        // since otherwise the rapidity of the parent in the lab frame is slightly modified
        //assert(obs_->add_entries_in_region(kab.stored_E()*kab, tab,
				//	ln_kt(tab) - log(xQ_), C1*event_.weight, &event_.axis()));
        return;
      } else {
        // ==> bin the emitter (much harder than the emission in LL kinematics)
	      assert(obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
				  ln_kt(t_insertion) - log(xQ_), C1*event_.weight, &event_.axis()));
	      return;
      }
    // only one of the two emissions ka and kb is inside the slice 
    } else if (thetaIn_ka or thetaIn_kb) {
      obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
				  ln_kt(t_insertion) - log(xQ_), C1*event_.weight, &event_.axis());
      //obs_->add_entries_in_region(emission->stored_E()*(*emission), t_insertion,
			//	  ln_kt(t_insertion) - log(xQ_), C1*event_.weight, &event_.axis());
      obs_->add_entries_in_region(emission->stored_E()*(*emission), t_second_insertion_,
				  ln_kt(t_second_insertion_) - log(xQ_), C1*event_.weight, &event_.axis());
      return;
    }

  // in branches 2 and 4 bin the emitter  
  } else if (ibranch==2 or ibranch==4) {
    if (ibranch==2) {
      if (obs_->add_entries_in_region(emitter->stored_E()*(*emitter), tab,
				    ln_kt(tab) - log(xQ_), C1*event_.weight, &event_.axis())) {
        delete emitter; //clean up from the allocation above
        return;
      }
      delete emitter;
    } else {
      if (obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
              ln_kt(t_insertion) - log(xQ_), C1*event_.weight, &event_.axis())) return;
    }
  }

  // now complete the evolution until the cutoff scale
  evolve_scale(t_insertion, evol_cutoff_, !NLL_EXPANDED);
}

//----------------------------------------------------------------------
/// generate first insertion
Momentum Shower::generate_first_insertion(double& t_insertion, int& idip_insertion) {
  if (!(NLL_EXPANDED and (!insertion_is_part_of_NLL_ensemble))) t_insertion += - log(rng.uniform_pos()) / (2.0 * CA * event_.eta_tot);
  double lnkt = ln_kt(t_insertion);

  if (event_.bad or (2.0*asmur_*b0*lnkt >= 1.0) or (evl_grid_ and t_insertion >= evl_grid_->xlim()) or (lnkt > lnktmax)) {
    event_.bad=true;
    return Momentum();
  }
  idip_insertion = choose_emitter();
  assert(idip_insertion>=0);
  
  Momentum insertion = event_[idip_insertion].radiate(lnkt - log(xQ_), rng.uniform_pos(),
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
    } else {
      event_.bad = true;
      return Momentum();
    }
  } else {
    assert(event_[idip].right_neighbour()>=0);
    //idip_insertion = (event_[event_[idip].right_neighbour()].delta_rap() != 0.0) ?
    //  event_[idip].right_neighbour() : idip;
    // Alternative: set weight to zero if selected dipole has no phase space for radiating
    if (event_[event_[idip].right_neighbour()].delta_rap() != 0.0) {
      idip_insertion = event_[idip].right_neighbour();
    } else {
      event_.bad = true;
      return Momentum();
    }  
  }

  // generate insertion with unit dipole transverse momentum
  Momentum insertion = event_[idip_insertion].radiate(0.0, rng.uniform_pos(),
						      rng.uniform_pos());
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
  while (true) {
    tb += - log(rng.uniform_pos()) / (2.0 * CA * event_[idip_insertion].delta_rap());
    if ((tb > 0.) or (tb != tb)) break;
  }
  double lnkt = ln_kt(tb);

  if (lnkt!=lnkt or (2.0*asmur_*b0*lnkt >= 1.0)
      or (evl_grid_ and tb >= evl_grid_->xlim()) or (lnkt > lnktmax)) {
    event_.bad = true;
    return Momentum();
  }
  // update emission's energy
  insertion.stored_E(insertion.stored_E()*exp(-lnkt+log(xQ_)));
  // sanity check: kt ordering in parent dipole frame
  double kta2 = 2./dot_product(event_[idip].left().momentum(),event_[rn].right().momentum())
      * dot_product(event_[idip].left().momentum(),emitter.stored_E() * emitter)
      * dot_product(emitter.stored_E() * emitter,event_[rn].right().momentum());
  double ktb2 = 2./dot_product(event_[idip].left().momentum(),event_[rn].right().momentum())
      * dot_product(event_[idip].left().momentum(),insertion.stored_E()*insertion)
      * dot_product(insertion.stored_E()*insertion,event_[rn].right().momentum());
  if ((!dipole_kt_ordering) and (kta2 < ktb2)) {
    event_.bad = true;
    return Momentum();
  }

  //t_second_insertion_ = tb; // define emission's time in the emitting dipole frame
  t_second_insertion_ = t_scale(-0.5*log(ktb2)); // define emission's time in the parent dipole frame
  if (t_second_insertion_ != t_second_insertion_) t_second_insertion_ = EVOLCUT;
  return insertion;
}

//----------------------------------------------------------------------
/// return the t scale for a given ln kt value
double Shower::t_scale(double lnkt) const {  
  if ((!NLL_EXPANDED) and evl_grid_ and NLL_evolution_) return evl_grid_->t(lnkt, xmur_, xQ_);
  // fixed coupling:
  //if ((!NLL_EXPANDED) and evl_grid_ and NLL_evolution_) return asmur_/(2.*M_PI)*lnkt * (1.0 + asmur_ / (2.*M_PI) * KCMW);
  // at LL:
  //   t = Log[1/(2 - 2 b0 alphas L)]/(4 pi b0)
  return -log(1 - 2*asmur_*b0*lnkt)/(4.*b0*M_PI);
  // fixed coupling: return as*lnkt/(2.0*M_PI);
  //return asmur_/(2.*M_PI)*lnkt;
}

//----------------------------------------------------------------------
/// return the ln(kt) for a given evolution scale
double Shower::ln_kt(double t) const {
  if ((!NLL_EXPANDED) and evl_grid_ and NLL_evolution_) return evl_grid_->ln_kt(t);
  // fixed coupling
  //if ((!NLL_EXPANDED) and evl_grid_ and NLL_evolution_) return (2.*M_PI)/(asmur_*(1.0 + asmur_ / (2.*M_PI) * KCMW))*t;
  /// at LL:
  ///    t = Log[1/(1 - 2 b0 alphas L)]/(4 pi b0)
  /// => ln(kt) = (1 - exp(-4 pi t b0))/(2 alphas b0)
  return (1 - std::exp(-4*M_PI*t*b0))/(2.0*asmur_*b0);
  // fixed coupling: return 2.0*M_PI*t/as;
  //return (2.*M_PI)/asmur_*t;
}

//----------------------------------------------------------------------
/// write current observable to file (or to cout)
void Shower::write(int nev, const std::string& fn) const {
  ostream * ostr;
  if (!fn.empty()) ostr = new ofstream(fn);
  else             ostr = &std::cout;
  *ostr << header_ << std::endl;
  *ostr << "# output with nev = " << nev << std::endl;
  obs_->write(nev, ostr);
}

//----------------------------------------------------------------------
/// reconstruct the parent kab such that it is massless
void Shower::reconstruct_parent(const Momentum& spec_left, const Momentum& spec_right, Momentum& kab, double& tab) const {

  // Option 1: keep parent massive.
  // Recombination is possible in any frame due to linearity of Lorentz transformations
  // compute lnkt=ln(Q/kt) in the parent dipole frame to calculate evolution time
  double m2 = dot_product(kab.stored_E()*kab, kab.stored_E()*kab);
  double lnktab = -log(sqrt(2.*dot_product(spec_left, kab.stored_E()*kab)
                *dot_product(spec_right, kab.stored_E()*kab)/dot_product(spec_left, spec_right) - m2));

  // Option 2: create massless parent in dipole frame (with parent dipole aligned along z axis)
  //Momentum k12 = spec_left+spec_right;
  //Momentum k1 = spec_left;
  //kab *= kab.stored_E();
  //kab.unboost(k12);
  //k1.unboost(k12);
  //double theta=k1.theta();
  //double phi=k1.phi();
  //kab.rotate(theta, phi);
  //double nab = kab.rap();
  //double ktab = sqrt(kab.px()*kab.px()+kab.py()*kab.py());
  //kab = Momentum(kab.px(), kab.py(), ktab*sinh(nab), ktab*cosh(nab));
  //kab.unrotate(theta,phi);
  //kab.boost(k12);
  //kab.stored_E(kab.E());
  //kab = (1.0/kab.E())*kab;
  //// compute lnkt=ln(Q/kt) in the parent dipole frame to calculate evolution time
  //double lnktab = -log(ktab);

  // Option 3: create massless parent directly in lab frame
  //kab = kab.stored_E()*kab;
  //double nab = kab.rap();
  //double ktab = sqrt(kab.px()*kab.px()+kab.py()*kab.py());
  //kab = Momentum(kab.px(), kab.py(), ktab*sinh(nab), ktab*cosh(nab));
  //kab.stored_E(kab.E());
  //kab = (1.0/kab.stored_E())*kab;
  //// compute lnkt=ln(Q/kt) in the parent dipole frame to calculate evolution time
  //double lnktab = -log(sqrt(2.*dot_product(spec_left, kab.stored_E()*kab)
  //              *dot_product(spec_right, kab.stored_E()*kab)/dot_product(spec_left, spec_right)));

  // finally compute evolution time
  tab = t_scale(lnktab);
  if (tab != tab) tab = EVOLCUT;

  return;
}
