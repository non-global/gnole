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
      evolve_insertion(tstart);
      run_threejet(tstart, false);
      run_threejet(tstart, true);
    }
  }
  write(nev, fn);
}

//----------------------------------------------------------------------
/// run the threejet piece
void Shower::run_threejet(double tstart, bool soft) {
  reset(true, soft);
  // flip sign for LL term
  if (soft) event_.weight = -event_.weight;
  evolve_scale(tstart);
  double w = -event_.weight;
  reset(false);
  event_.weight = w;
  // now veto configurations where gluon is in the side as computed in the two jet limit
  if (obs_->in_region(gluon_, &event_.axis())) event_.bad=true;
  evolve_scale(tstart);
}

//----------------------------------------------------------------------
/// reset the dipoles to an initial qqbar pair along the z axis
void Shower::reset(bool threejet, bool soft) {
  event_.bad=false;
  if (!threejet) {
    event_.reset_twojet();
  } else if (soft){
    bool in_region = false;
    event_.reset_threejet_LL(xmur_, xQ_, rng.uniform(), rng.uniform(), rng.uniform(), gluon_);
    Shower::do_split(0,gluon_);
    in_region = obs_->in_region(gluon_, &event_.axis());
    if (in_region or event_.weight==0.0) event_.bad = true;
  } else {
    bool in_region = false;
    event_.reset_threejet(xmur_, xQ_, rng.uniform(), rng.uniform(), rng.uniform(), gluon_);
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
void Shower::evolve_scale(double t, double tend) {
  if (event_.weight == 0.0 or event_.bad) return;
  while ((t += - log(rng.uniform_pos()) / (2.0 * CA * event_.eta_tot)) < tend) {
    // increase by amount chosen with distribution e^(-2*CA*sum_rap)    
    int idip = choose_emitter();
    // now set all emission properties
    double lnkt = ln_kt(t);
    if (2.0*asmur_*b0*lnkt>=1.0 or (evl_grid_ and t>=evl_grid_->xlim())) {
      event_.bad = true; // setting this to avoid starting another evolution later
      break;
    }
    //fixed coupling: if (event_.eta_tot < 0.0 or event_.eta_tot!=event_.eta_tot) break;
    assert(idip>=0);
    // the momentum of the gluon is xQ*exp(-lnkt)=exp[-(lnkt+ln(xQ))]
      
    Momentum emsn = event_[idip].radiate(lnkt - log(xQ_), rng.uniform_pos(),
					 rng.uniform_pos());

    if (!do_split(idip, emsn)) continue;
    // if emission is in observed region, add to histogram and stop the evolution
    double C1 = (NLL_counterterm_ ?  1.0 + asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1) : 1.0);
    if (obs_->add_entries_in_region(emsn.stored_E()*emsn, t, lnkt-log(xQ_),
				    C1*event_.weight, &event_.axis())) {
      event_.bad = true; // setting this to avoid starting another evolution later
      break;
    }
  }
  return;
}

//----------------------------------------------------------------------
/// evolve the shower down to cutoff scale
void Shower::evolve_insertion(double t) {
  if (event_.weight == 0.0 or event_.bad) return;
  NLL_counterterm_ = true;

  // generate log(Q/kt) in [log(xQ), 1/(2 as b0)]
  // generate scale t from log(Q/kt)
  double lnkt_insertion = log(xQ_) + rng.uniform()*(1.0/(2.0*asmur_*b0) - log(xQ_));
  double t_insertion = t_scale(lnkt_insertion);
  
  evolve_scale(t, min(t_insertion, evol_cutoff_));
  int idipa = -1;
  if (t_insertion >= evol_cutoff_) {
    NLL_counterterm_ = false;
    return;
  }
  Momentum ka = generate_first_insertion(t_insertion, idipa);
  if (event_.bad) {
    NLL_counterterm_ = false;
    return;
  }

  event_cache_->copy(event_);
  
  // branch0 = Z^{(0)} evolution 
  perform_branch0(t_insertion, ka);
  event_.retrieve(event_cache_);
  
  // branch1 = Z^{(1)} evolution double real with kta > ktb'
  perform_branch(t_insertion, idipa, 1, ka);
  event_.retrieve(event_cache_);

  // branch2 = Z^{(1)} evolution collinear counterterm with kta > ktb'
  perform_branch(t_insertion, idipa, 2, ka);
  event_.retrieve(event_cache_);
  
  // branch3 = Z^{(1)} evolution double real with kta > ktb
  perform_branch(t_insertion, idipa, 3, ka);
  event_.retrieve(event_cache_);
  
  // branch4 = Z^{(1)} evolution collinear counterterm with kta > ktb
  perform_branch(t_insertion, idipa, 4, ka);
  NLL_counterterm_ = false;
}

void Shower::perform_branch0(double ta, const Momentum& ka) {
  if (event_.bad) return;
  // if emission is in observed region, add to histogram and stop the evolution
  double C1 = (NLL_counterterm_ ?  1.0 + asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1) : 1.0);
  if (obs_->add_entries_in_region(ka.stored_E()*ka, ta, ln_kt(ta)-log(xQ_),
				  C1*event_.weight, &event_.axis())) return;
  evolve_scale(ta, evol_cutoff_);
}

void Shower::perform_branch(double t_insertion, int idipa, int ibranch, const Momentum& ka) {
  int idipb = -1;
  Momentum kb = generate_second_insertion(t_insertion, idipa, idipb, (ibranch==3 or ibranch==4));
  if (event_.bad) return;
  // add in the weight of the second insertion
  event_.weight*=2.0;
  const Momentum* emitter = &ka;
  const Momentum* spec_left;
  const Momentum* spec_right;
  const Momentum* emission = &kb;
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
    }
    // replace emitter with massless version of parent for ibranch 2
    if (ibranch==2) {
      Momentum* kab = new Momentum(emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission));
      *kab = (1.0/kab->E())*(*kab);
      reconstruct_parent(*spec_left, *spec_right, *kab);
      emitter = kab;
      // now updated the dipoles containing ka
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
      event_[idipb].right(DipoleEnd(*kab));
      event_[event_[idipb].right_neighbour()].left(DipoleEnd(*kab));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].right_neighbour()].delta_rap();
    }
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
    }
    // replace emitter with massless version of parent for ibranch 2
    if (ibranch==2) {
      Momentum* kab = new Momentum(emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission));
      *kab = (1.0/kab->E())*(*kab);
      reconstruct_parent(*spec_left, *spec_right, *kab);
      emitter = kab;
      // now updated the dipoles containing ka
      event_.eta_tot -= event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
      event_[idipb].left(DipoleEnd(*kab));
      event_[event_[idipb].left_neighbour()].right(DipoleEnd(*kab));
      event_.eta_tot += event_[idipb].delta_rap() + event_[event_[idipb].left_neighbour()].delta_rap();
    }
  }
  
  // update the sign, and flip sign for branches 2 and 3
  event_.weight*= ((ibranch==1 or ibranch==4) ? w : -w);
  if (ibranch==1 or ibranch == 3)
    assert(do_split(idipb, kb));

  double C1 = (NLL_counterterm_ ?  1.0 + asmur_/(2.0*M_PI) * (CF*integrated_counterterm_ + H1) : 1.0);
  if (ibranch==1 or ibranch==3) {
    bool thetaIn_ka = obs_->in_region(emitter->stored_E()*(*emitter), &event_.axis());
    bool thetaIn_kb = obs_->in_region(emission->stored_E()*(*emission), &event_.axis());
    if (thetaIn_ka and thetaIn_kb) {
      if (ibranch==1) {
	Momentum kab = emitter->stored_E()*(*emitter) + emission->stored_E()*(*emission);
	kab = (1.0/kab.E())*kab;
	reconstruct_parent(*spec_left, *spec_right, kab);
	if (obs_->add_entries_in_region(kab.stored_E()*kab, t_insertion,
					ln_kt(t_insertion) - log(xQ_),
					C1*event_.weight, &event_.axis()))
	  return;
      } else {
	assert(obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
					   ln_kt(t_insertion) - log(xQ_),
					   C1*event_.weight, &event_.axis()));
	return;
      }
    } else if (thetaIn_ka or thetaIn_kb) {
      obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
				  ln_kt(t_insertion) - log(xQ_),
				  C1*event_.weight, &event_.axis());
      obs_->add_entries_in_region(emission->stored_E()*(*emission), t_insertion,
				  ln_kt(t_insertion) - log(xQ_),
				  C1*event_.weight, &event_.axis());
      return;
    }
  } else if (ibranch==2 or ibranch==4) {
    if (obs_->add_entries_in_region(emitter->stored_E()*(*emitter), t_insertion,
				    ln_kt(t_insertion) - log(xQ_),
				    C1*event_.weight, &event_.axis())) {
      if (ibranch==2) delete emitter; //clean up from the allocation above
      return;
    } else if (ibranch==2) delete emitter;
  }
  evolve_scale(t_insertion, evol_cutoff_);
}

Momentum Shower::generate_first_insertion(double& t_insertion, int& idip_insertion) {
  t_insertion += - log(rng.uniform_pos()) / (2.0 * CA * event_.eta_tot);
  double lnkt = ln_kt(t_insertion);
  if (event_.bad or 2.0*asmur_*b0*lnkt>=1.0 or (evl_grid_ and t_insertion>=evl_grid_->xlim())) {
    event_.bad=true;
    return Momentum();
  }
  idip_insertion = choose_emitter();
  assert(idip_insertion>=0);
  Momentum insertion = event_[idip_insertion].radiate(lnkt - log(xQ_), rng.uniform_pos(),
					    rng.uniform_pos());
  assert(do_split(idip_insertion, insertion));
  return insertion;
}
Momentum Shower::generate_second_insertion(double t_insertion, int idip, int& idip_insertion, bool dipole_kt_ordering) {
  if (rng.uniform() < 0.5) {
    idip_insertion = event_[idip].delta_rap()!=0.0 ? idip : event_[idip].right_neighbour();
  } else {
    assert(event_[idip].right_neighbour()>=0);
    idip_insertion = (event_[event_[idip].right_neighbour()].delta_rap() != 0.0) ?
      event_[idip].right_neighbour() : idip;
  }
  Momentum insertion = event_[idip_insertion].radiate(0.0, rng.uniform_pos(),
						      rng.uniform_pos());
  int f2 = 1.0;
  int rn = event_[idip].right_neighbour();
  if (!dipole_kt_ordering) 
    f2 = 0.5*dot_product(event_[idip].left().momentum(),event_[rn].right().momentum())
      / dot_product(event_[idip].left().momentum(),insertion)
      / dot_product(insertion,event_[rn].right().momentum());
  double tb = t_scale(ln_kt(t_insertion)-0.5*log(f2));
  tb += - log(rng.uniform_pos()) / (2.0 * CA * event_[idip_insertion].delta_rap());
  double lnkt = ln_kt(tb);
  
  if (lnkt!=lnkt or 2.0*asmur_*b0*lnkt>=1.0
      or (evl_grid_ and t_insertion>=evl_grid_->xlim())) {
    event_.bad=true;
    return Momentum();
  }
  insertion.stored_E(exp(-lnkt+log(xQ_)));
  return insertion;
}

//----------------------------------------------------------------------
/// return the t scale for a given ln kt value
double Shower::t_scale(double lnkt) const {  
  if (evl_grid_) return evl_grid_->t(lnkt, xmur_, xQ_);
  // at LL:
  //   t = Log[1/(2 - 2 b0 alphas L)]/(4 pi b0)
  return -log(1 - 2*asmur_*b0*lnkt)/(4.*b0*M_PI);
  // fixed coupling: return as*lnkt/(2.0*M_PI);
}

//----------------------------------------------------------------------
/// return the ln(kt) for a given evolution scale
double Shower::ln_kt(double t) const {
  if (evl_grid_) return evl_grid_->ln_kt(t);
  /// at LL:
  ///    t = Log[1/(1 - 2 b0 alphas L)]/(4 pi b0)
  /// => ln(kt) = (1 - exp(-4 pi t b0))/(2 alphas b0)
  return (1 - std::exp(-4*M_PI*t*b0))/(2.0*asmur_*b0);
  // fixed coupling: return 2.0*M_PI*t/as;
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

void Shower::reconstruct_parent(const Momentum& spec_left, const Momentum& spec_right, Momentum& kab) const {
  // Create parent in dipole frame
  Momentum k12 = spec_left+spec_right;
  Momentum k1 = spec_left;
  double storedE = kab.stored_E();
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
  kab.stored_E(kab.E()*storedE);
  kab = (1.0/kab.E())*kab;
  return;
}
