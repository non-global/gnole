#ifndef __GRID_HH__
#define __GRID_HH__

#include <vector>
#include <string>

#ifdef NNET
#include "torch/script.h"
#endif

//======================================================================
/// class with a grid containing precomputed values
class Grid {
public:
  /// constructor
  Grid(int npoints, double xlim=1e100, double ylim=0.0)
    : npoints_(npoints), xlim_(xlim), ylim_(ylim), clip_(false)
#ifdef NNET
    , nnet_(false)
#endif
  {}
  
#ifdef NNET
  /// constructor with neural net
  Grid(std::string fn, double xlim=1e100, double ylim=0.0)
    : npoints_(0), xlim_(xlim), ylim_(ylim), clip_(false), nnet_(true) {
    module_ = torch::jit::load(fn.c_str());
  }
#endif

  virtual ~Grid() {}
  
  double eval(double x) {
    if (x >= xlim_) return ylim_;
#ifdef NNET
    if (nnet_) {
      std::vector<torch::jit::IValue> inputs;
      inputs.push_back(torch::ones({1,1})*x);
      auto output = module_.forward(inputs).toTensor();
      double res = output[0][0].item<double>();
      if (clip_) res = res>clipmax_ ? clipmax_ : res;
      if (clip_) res = res<clipmin_ ? clipmin_ : res;
      return res;
    }
#endif
    unsigned int i = bin(x);
    double f = (x-grid[i].first)/(grid[i+1].first-grid[i].first);
    double y = grid[i].second + f*(grid[i+1].second - grid[i].second);
    if (clip_) y = y>clipmax_ ? clipmax_ : y;
    if (clip_) y = y<clipmin_ ? clipmin_ : y;
    return y;
  }

  void set_clip(double min, double max) {
    clip_ = true;
    clipmin_ = min;
    clipmax_ = max;
  }
  
  virtual int bin(double tval) const = 0;

  double xlim() const {return xlim_;}
  double ylim() const {return ylim_;}
protected:
  int npoints_;
  double xlim_, ylim_;
  bool clip_;
  double clipmin_, clipmax_;
#ifdef NNET
  bool nnet_;
  torch::jit::script::Module module_;
#endif
  /// grid containing the pairs of values
  std::vector< std::pair<double,double> > grid;
};
#endif //  __GRID_HH__
