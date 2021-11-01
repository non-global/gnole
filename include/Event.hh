#ifndef __EVENT_HH__
#define __EVENT_HH__

#include "Momentum.hh"
#include "Dipole.hh"
#include <vector>

/// class containing a current event with a configuration of dipoles
class Event {
public:
  Event() : weight(1.0), bad(false), eta_tot(0.0), thrust_axis_(0,0,1,1) {}
  // Event(const Event& event) : weight(event.weight), bad(event.bad),
  // 			      eta_tot(event.eta_tot),
  // 			      dipoles_(event.dipoles()),
  // 			      thrust_axis(event.thrust_axis) {}
  /// number of dipoles
  unsigned int size() const {
    return dipoles_.size();
  }
  
  /// access a specific dipole in the event
  Dipole & operator[](const int i) {
    return dipoles_[i];
  }

  /// access a specific dipole in the event
  const Dipole & operator[](const int i) const {
    return dipoles_[i];
  }

  /// return vector of dipoles
  std::vector<Dipole>& dipoles() {return dipoles_;}

  /// add dipole to vector of dipoles
  void add(const Dipole& dipole) {
    dipoles_.push_back(dipole);
  }
  
  /// return the thrust axis
  const Momentum& axis() const {return thrust_axis_;}
  
  /// copy an event from a cached pointer
  void copy(const Event& event) {
    weight = event.weight;
    bad = event.bad;
    eta_tot = event.eta_tot;
    dipoles_.clear();
    for(unsigned int i = 0; i<event.size(); ++i) {
      dipoles_.push_back(event[i]);
    }
  }
  /// retrieve an event from a cached pointer
  void retrieve(Event* cache) {
    weight = cache->weight;
    bad = cache->bad;
    eta_tot = cache->eta_tot;
    dipoles_.clear();
    for(Dipole dip : cache->dipoles()) {
      dipoles_.push_back(dip);
    }
  }
  
  /// event weight
  double weight;

  /// event bool
  bool bad;
  
  /// total rapidity of dipoles
  double eta_tot;

  /// reset shower to two jet configuration
  void reset_twojet();
  
  /// reset shower to three jet LL configuration and return lnkt
  double reset_threejet_LL(double xmur, double xq, double r1, double r2, double r3, Momentum& gluon);

  /// reset shower to three jet configuration and return ln kt, taken from Dokshitzer paper
  double reset_threejet(double xmur, double xq, double r1, double r2, double r3, Momentum& gluon);
  
private:
  /// vector of dipoles
  std::vector<Dipole> dipoles_;
  
  /// thrust axis
  Momentum thrust_axis_;
};
#endif //__EVENT_HH__
