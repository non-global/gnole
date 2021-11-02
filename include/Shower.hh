/// Code for the core elements of a shower.
#ifndef __SHOWER_HH__
#define __SHOWER_HH__
#include "Event.hh"
#include "GSLRandom.hh"
#include "EvolGrid.hh"
#include "Observables.hh"
#include "matrix_elements.hh"
#include <string>
#include <cassert>
#include <iostream>

//----------------------------------------------------------------------
/// \class Shower
/// contains the core shower evolution
class Shower {
public:
  /// constructor
  Shower(Observable& obs,
	 double xmur,
	 double xQ,
#ifdef NNET
	 std::string fn_evl_nn,
#endif
	 int order_evl,
	 std::string header = "",
	 int seed           = 0)
    : xmur_(xmur), xQ_(xQ), order_evl_(order_evl), obs_(&obs), NLL_counterterm_(false),
      evl_grid_(0), header_(header), rng(seed), gluon_(0,0,0,0), event_cache_(new Event()) {
    assert((order_evl_==0) || (order_evl_==1));
    asmur_=alphas2(xmur_);
    // if order = 1, set up the grid for ln kt
    if (order_evl_==1) {
      std::cerr << "# Initialization of grids" << std::flush;
#ifdef NNET
      if (fn_evl_nn!="") evl_grid_ = new EvolGrid(fn_evl_nn);
      else
#endif
	evl_grid_ = new EvolGrid(NEVLGRID, xmur_, xQ_);
      std::cerr << ": done." << std::endl;
      // set the integrated coefficient
      integrated_counterterm_ = integrated_counterterm(obs_->parameter());
    }
  }

  /// destructor
  virtual ~Shower() {
    if (evl_grid_) delete evl_grid_;
    delete event_cache_;
  }

  /// run the shower with nev events
  virtual void run(int nev, const std::string& fn = "");
  
protected:
  /// run the threejet piece
  virtual void run_threejet(double tstart, bool soft = false);
  
  /// reset the dipoles to an initial qqbar pair along the z axis
  void reset(bool threejet = false, bool soft = false);

  /// choose an emitting dipole
  int choose_emitter() const;

  /// decide whether to split the dipole using new emission
  bool do_split(int idip, Momentum&  emsn);
  
  /// evolve the shower between two scales
  void evolve_scale(double t, double tend = evol_cutoff_);
  
  /// evolve the shower including insertion
  void evolve_insertion(double tstart);

  /// perform the Z^{(0)} evolution
  void perform_branch0(double ta, const Momentum& ka);

  /// perform the Z^{(1)} evolution
  void perform_branch(double t_insertion, int idipa, int ibranch, const Momentum& ka);
  
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
  
  /// event containing dipole configuration
  Event event_;
  /// renormalisation scale
  double xmur_;
  /// resummation scale
  double xQ_;
  /// alphas at renormalisation scale
  double asmur_;
  /// order of the shower (LL=0, NLL=1)
  int order_evl_;
  /// use virtual corrections counteterm
  bool NLL_counterterm_;
  /// integrated counterterm
  double integrated_counterterm_;
  /// pointer to observable
  Observable * obs_;
  /// pointer to evolution grid
  EvolGrid * evl_grid_;
  /// header with description of run
  std::string header_;
  /// random number generator
  GSLRandom rng;
  /// momentum of hard gluon
  Momentum gluon_;
  /// cached information for insertion
  Event * event_cache_;
  /// cutoff for the evolution
  static constexpr double evol_cutoff_ = EVOLCUT;
};

#endif //__SHOWER_HH__
