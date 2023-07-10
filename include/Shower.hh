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
    : xmur_(xmur), xQ_(xQ), order_evl_(order_evl), obs_(&obs),
      evl_grid_(0), NLL_evolution_(false),  header_(header),
      rng(seed), gluon_(0,0,0,0), event_cache_(new Event()) {
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
  virtual void run_threejet(double tstart, bool soft = false, bool use_cached_variables = false);
  
  /// reset the dipoles to an initial qqbar pair along the z axis
  void reset(bool threejet = false, bool soft = false, bool use_cached_variables = false);

  /// choose an emitting dipole
  int choose_emitter() const;

  /// decide whether to split the dipole using new emission
  bool do_split(int idip, Momentum&  emsn);
  
  /// evolve the shower between two scales
  void evolve_scale(double t, double tend = evol_cutoff_, bool include_as_constant = true);
  
  /// evolve the shower including insertion
  void evolve_insertion(double tstart);

  /// evolve the shower including insertion (NLL fully expanded out)
  void evolve_insertion_expanded(double tstart);

  /// perform the Z^{(0)} evolution
  void perform_branch0(double ta, const Momentum& ka);

  /// perform the perturbative NLL insertion of the Z^{(0)} evolution (used when NLL_EXPANDED=true)
  void perform_branch_single_insertion(double t_insertion, int ibranch, const Momentum& ka);

  /// perform the Z^{(1)} evolution
  void perform_branch_double_insertion(double t_insertion, int idipa, int ibranch, const Momentum& ka);
  
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
  void reconstruct_parent(const Momentum& spec_left, const Momentum& spec_right, Momentum& kab, double& tab) const;
  
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
  /// order of the shower (LL=0, NLL=1)
  int order_evl_;
  /// pointer to observable
  Observable * obs_;
  /// pointer to evolution grid
  EvolGrid * evl_grid_;
  /// use virtual corrections counterterm and NLL evolution
  bool NLL_evolution_;
  /// integrated counterterm
  double integrated_counterterm_;
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
  }
  void retrieve_second_insertion(Momentum& kb, int& idipb, bool& event_bad) {
    kb        = momentum_cache_;
    idipb     = idip_cache_;
    event_bad = event_bad_cache_;
  }
  bool   weighted_second_insertion = false; // if true the second insertion is generated at FO with a weight (instead of with a Sudakov)
  double second_insertion_weight_;
  
  /// cached information for insertion
  Event * event_cache_;
  /// cutoffs for the evolution
  static constexpr double evol_cutoff_ = EVOLCUT;
  static constexpr double landau_pole_tolerance_ = 1.;
};

#endif //__SHOWER_HH__
