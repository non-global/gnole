#ifndef __EVOLGRID_HH__
#define __EVOLGRID_HH__

#include "Grid.hh"
#include "Parameters.hh"
#include <iostream>
//----------------------------------------------------------------------
/// \class EvolGrid
/// a grid containing the relation between ln kt and
/// evolution scale for fast interpolation
class EvolGrid : public Grid {
public:
  /// constructor
  EvolGrid(int npoints, double xmur, double xQ) : Grid(npoints) {
    double Lmin = log(xQ);
    double Lmax = 1.0/(2.0*alphas2(xmur)*b0); // Landau pole
    for(int i=0; i<npoints_; ++i) {
      double L=Lmin + i*(Lmax-Lmin)/(npoints_);
      double tval=t(L, xmur, xQ);
      if (!grid.empty() and tval<grid.back().first) break;
      grid.push_back(std::make_pair(tval,L));
    }
    npoints_=grid.size();
    xlim_=grid.back().first;
    ylim_=grid.back().second;
  }
#ifdef NNET
  /// constructor for neural net
  EvolGrid(std::string fn) : Grid(fn) {
    npoints_=NGRID;
    grid.push_back(std::make_pair(0.0,0.0));
    for(int i=1; i<npoints_; ++i) {
      double lnt=LNEVOLMIN+double(i)*(LNEVOLMAX-LNEVOLMIN)/(npoints_-1);
      // because the neural network is fitted as log(x),log(y)
      // we need to take the appropriate exp and log accordingly
      grid.push_back(std::make_pair(exp(lnt), exp(ln_kt(lnt))));
    }
    xlim_=grid.back().first;
    ylim_=grid.back().second;
    nnet_=false;
  }
#endif
  
  virtual ~EvolGrid() {}
  
  //----------------------------------------------------------------------
  /// get the t value corresponding to a given ln kt
  double t(double L, double xmur, double xQ) const {
    double asmur=alphas2(xmur);
    return (2*asmur*asmur*b0*b1*L + (b0 + asmur*b1 - 2*asmur*b0*b0*L)*
	    log(1 - 2*asmur*b0*L))/(4.*b0*b0*(-1 + 2*asmur*b0*L)*M_PI)
      + (asmur*KCMW*asmur*b0*L)/(4*b0*M_PI*M_PI - 8*b0*asmur*b0*L*M_PI*M_PI)
      + 2*asmur*b0*b0*M_PI*(-2*asmur*b0*L*log(xmur) + log(xQ))/
      (4.*b0*b0*(-1 + 2*asmur*b0*L)*M_PI*M_PI);
    // at LL:
    //   return -log(1 - 2*asmur*b0*L)/(4.*b0*M_PI);
  }

  //----------------------------------------------------------------------
  /// get the ln kt value corresponding to a given t scale
  double ln_kt(double tval) {
    return Grid::eval(tval);
  }
  
  //----------------------------------------------------------------------
  /// find the bin corresponding to a t value
  virtual int bin(double tval) const {
    int start=0, end=grid.size()-1;
    while(start<=end) {
      int mid=start+(end-start)/2;
      if (grid[mid].first < tval) {
	     start = mid+1;
      } else if (grid[mid].first >= tval) {
	     end = mid-1;
      } else {
	break;
      }
    }
    return start;
  }

};

#endif //  __EVOLGRID_HH__
