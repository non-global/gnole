#ifndef __SIMPLEHIST_HH__
#define __SIMPLEHIST_HH__

//STARTHEADER
// $Id: SimpleHist.hh 3392 2012-08-01 10:42:51Z salam $
//
// Copyright (c) 2007-2011, Matteo Cacciari, Gavin Salam and Gregory Soyez
//
//----------------------------------------------------------------------
//
//  SimpleHist is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  SimpleHist is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//----------------------------------------------------------------------
//ENDHEADER


#include<valarray>
#include<string>
#include<cmath>
#include<iostream>
#include<cassert>
#include<vector>
#include<stdexcept>
#include "Parameters.hh"

// Neumaier summation preserves a small signed residual even when a later
// opposite contribution cancels the leading sum (e.g. 1e16 + 1 - 1e16).
inline void histogram_accumulate(double value, double &sum, double &correction,
                                 bool precise) {
  if (!precise) { sum += value; return; }
  const double next = sum + value;
  correction += std::abs(sum) >= std::abs(value)
              ? (sum - next) + value : (value - next) + sum;
  sum = next;
}

class SimpleHist {
public:
  SimpleHist() : _minv(0.), _maxv(0.), _dv(0.),
                 _precise(PRECISE_KERNELS) { reset(); }
  SimpleHist(double minv, double maxv, int n): 
    _minv(minv), _maxv(maxv), _dv((maxv-minv)/n), _weights(n+2),
    _corrections(n+2), _view(n+2), _precise(PRECISE_KERNELS) {
    reset();
  }

  SimpleHist(double minv, double maxv, unsigned int n) {
    declare(minv, maxv, n);
  }

  SimpleHist(double minv, double maxv, double bin_size) {
    declare(minv, maxv, int(0.5+(maxv-minv)/bin_size));
  }

  // declare (or redeclare) the histogram
  void declare(double minv, double maxv, double bin_size) {
    declare(minv, maxv, int(0.5+(maxv-minv)/bin_size));
  }

  void declare(double minv, double maxv, int n) {
    declare(minv, maxv, unsigned(n));
  }

  // declare (or redeclare) the histogram
  void declare(double minv, double maxv, unsigned int n) {
    _minv = minv; _maxv = maxv; _dv = (maxv-minv)/n; 
    _weights.resize(n+2); // two extra bins for the underflow and overflow
    _corrections.resize(n+2);
    _view.resize(n+2);
    _precise = PRECISE_KERNELS;
    reset();
  }

  /// reset the contents of the histogram to zero (does not
  /// modify the histogram bounds)
  void reset() {
    _weights = 0.0;
    _corrections = 0.0;
    _view = 0.0;
    _weight_v = 0.0;
    _weight_vsq = 0.0;
    _weight_v_correction = 0.0;
    _weight_vsq_correction = 0.0;
    _total_weight = 0.0;
    _have_total = false;
    _n_entries = 0.0;
  }

  double min() const {return _minv;};
  double max() const {return _maxv;};
  /// returns the size of the histogram proper
  unsigned int size() const {
    unsigned outflow_sz = _weights.size();
    assert(outflow_sz > 0); // help capture uninitialised histogram bugs
    return outflow_sz-2;
  }
  /// returns the size of the histogram plus outflow bin
  unsigned int outflow_size() const {return _weights.size();};

  double & operator[](int i) {
    _have_total = false;
    _weights[i] += _corrections[i];
    _corrections[i] = 0.;
    return _weights[i];
  }
  const double & operator[](int i) const {
    if (!_precise) return _weights[i];
    _view[i] = _weights[i] + _corrections[i];
    return _view[i];
  }
  
  /// handle underflow and overflow bins
  unsigned int underflow_bin() const {return size();}
  unsigned int overflow_bin()  const {return size()+1;}
  double & underflow() {return (*this)[underflow_bin()];}
  double & overflow()  {return (*this)[overflow_bin()];}
  const double & underflow() const {return (*this)[underflow_bin()];}
  const double & overflow()  const {return (*this)[overflow_bin()];}

  /// returns the outflow bin
  double outflow() {
    const SimpleHist &hist = *this;
    return hist.underflow() + hist.overflow();
  }
  //double & outflow() {return _weights[size()];};
  //const double & outflow() const {return _weights[size()];};

  double binlo (int i) const {return i*_dv + _minv;};
  double binhi (int i) const {return (i+1)*_dv + _minv;};
  double binmid(int i) const {return (i+0.5)*_dv + _minv;};
  double binsize()     const {return _dv;};

  unsigned int bin(double v) const {
    if (!std::isfinite(v))
      throw std::invalid_argument("Nonfinite histogram coordinate");
    if (v >= _minv && v < _maxv) {
      int i = int((v-_minv)/_dv); 
      if (i >= 0) {
        // do unsigned comparisons to avoid issues with huge
        // histograms
        unsigned unsigned_i = unsigned(i);
        if (unsigned_i < size()) {return unsigned_i;} 
      }
    }
    // otherwise handle outflow bins
    if (v <  _minv) {
      return underflow_bin();
    } else {
      assert(v >= _maxv);
      return overflow_bin();
    }
    //return size();
  }

  /// return the mean value of all events given to histogram
  /// including those that were outside the histogram edges
  double mean() const {return (_weight_v + _weight_v_correction) / total_weight();}

  /// return the total weight in the histogram (inefficient)...
  double total_weight() const {
    if (!_have_total) {
      _total_weight = 0.0;
      double correction = 0.;
      for (unsigned i = 0; i < _weights.size(); i++) {
        histogram_accumulate(_weights[i], _total_weight, correction, _precise);
        if (_precise)
          histogram_accumulate(_corrections[i], _total_weight, correction, true);
      }
      _total_weight += correction;
      _have_total = true;
    }
    return _total_weight;
  }

  double n_entries() const {return _n_entries;}

  /// Include the residual separately when summing bins into a cumulative
  /// quantity, so an unresolved bin-level residual can survive cancellation.
  void accumulate_bin_weight(unsigned int index, double norm, double &sum,
                             double &correction) const {
    // Keep this multiply/add in one full expression for the original
    // compiler contraction policy in numerical-reference mode.
    if (!_precise) { sum += _weights[index] * norm; return; }
    histogram_accumulate(_weights[index] * norm, sum, correction, _precise);
    if (_precise)
      histogram_accumulate(_corrections[index] * norm, sum, correction, true);
  }

  void add_entry(double v, double weight = 1.0) {
    add_entry_at_bin(bin(v), v, weight);
  }

  /// Fill an already computed bin, avoiding duplicate observable/bin work.
  void add_entry_at_bin(unsigned int index, double v, double weight = 1.0) {
    if (!std::isfinite(v) || !std::isfinite(weight))
      throw std::invalid_argument("Nonfinite histogram value or event weight");
    if (index >= _weights.size()) throw std::out_of_range("Histogram bin");
    _have_total = false;
    histogram_accumulate(weight, _weights[index], _corrections[index], _precise);
    histogram_accumulate(weight * v, _weight_v, _weight_v_correction, _precise);
    histogram_accumulate(weight * v * v, _weight_vsq, _weight_vsq_correction, _precise);
    if (!std::isfinite(_weights[index]) || !std::isfinite(_corrections[index]) ||
        !std::isfinite(_weight_v) || !std::isfinite(_weight_v_correction) ||
        !std::isfinite(_weight_vsq) || !std::isfinite(_weight_vsq_correction))
      throw std::overflow_error("Nonfinite histogram weight or moment accumulator");
    _n_entries += 1.0;
  };

  // Operations with constants ---------------------------------------
  SimpleHist & operator*=(double fact) {
    _have_total = false;
    for (unsigned i = 0; i < outflow_size(); i++) {
      _weights[i] *= fact;
      _corrections[i] *= fact;
    }
    _weight_v *= fact;
    _weight_vsq *= fact;
    _weight_v_correction *= fact;
    _weight_vsq_correction *= fact;
    _total_weight *= fact;
    return *this;
  };
  SimpleHist & operator/=(double fact) {
    *this *= 1.0/fact;
    return *this;
  };

  // Operations with another histogram -------------------------------
  SimpleHist & operator*=(const SimpleHist & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] *= other[i];
    return *this;
  };

  SimpleHist & operator/=(const SimpleHist & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] /= other[i];
    return *this;
  };

  SimpleHist & operator+=(const SimpleHist & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) {
      histogram_accumulate(other._weights[i], _weights[i], _corrections[i], _precise);
      histogram_accumulate(other._corrections[i], _weights[i], _corrections[i], _precise);
    }
    histogram_accumulate(other._weight_v, _weight_v, _weight_v_correction, _precise);
    histogram_accumulate(other._weight_v_correction, _weight_v, _weight_v_correction, _precise);
    histogram_accumulate(other._weight_vsq, _weight_vsq, _weight_vsq_correction, _precise);
    histogram_accumulate(other._weight_vsq_correction, _weight_vsq, _weight_vsq_correction, _precise);
    _have_total = false;
    return *this;
  };

  SimpleHist & operator-=(const SimpleHist & other) {
    assert(other.outflow_size() == outflow_size());
    _have_total = false;
    for (unsigned i = 0; i < outflow_size(); i++) {
      histogram_accumulate(-other._weights[i], _weights[i], _corrections[i], _precise);
      histogram_accumulate(-other._corrections[i], _weights[i], _corrections[i], _precise);
    }
    return *this;
  };


  friend SimpleHist operator*(const SimpleHist & hist, double fact);
  friend SimpleHist operator/(const SimpleHist & hist, double fact);

private:
  double _minv, _maxv, _dv;
  std::valarray<double> _weights;
  std::valarray<double> _corrections;
  mutable std::valarray<double> _view;
  bool _precise;
  std::string _name;
  double _weight_v, _weight_vsq;
  double _weight_v_correction, _weight_vsq_correction;
  mutable double _total_weight;
  mutable bool   _have_total;
  double _n_entries;
};



// Binary operations with constants -----------------------------
inline SimpleHist operator*(const SimpleHist & hist, double fact) {
  SimpleHist result(hist.min(), hist.max(), hist.size());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = hist[i] * fact;
  result._weight_v = (hist._weight_v + hist._weight_v_correction) * fact;
  result._weight_vsq = (hist._weight_vsq + hist._weight_vsq_correction) * fact;
  return result;
}

inline SimpleHist operator/(const SimpleHist & hist, double fact) {
  SimpleHist result(hist.min(), hist.max(), hist.size());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = hist[i] / fact;
  result._weight_v = (hist._weight_v + hist._weight_v_correction) / fact;
  result._weight_vsq = (hist._weight_vsq + hist._weight_vsq_correction) / fact;
  return result;
}

inline SimpleHist operator*(double fact, const SimpleHist & hist) {
  return hist*fact;
}

inline SimpleHist operator/(double fact, const SimpleHist & hist) {
  return hist/fact;
}


// Binary operations with other histograms ------------------------
inline SimpleHist operator*(const SimpleHist & hista, const SimpleHist & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist result(hista.min(), hista.max(), hista.size());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] * histb[i];
  return result;
}
inline SimpleHist operator/(const SimpleHist & hista, const SimpleHist & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist result(hista.min(), hista.max(), hista.size());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] / histb[i];
  return result;
}
inline SimpleHist operator+(const SimpleHist & hista, const SimpleHist & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist result(hista.min(), hista.max(), hista.size());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] + histb[i];
  return result;
}
inline SimpleHist operator-(const SimpleHist & hista, const SimpleHist & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist result(hista.min(), hista.max(), hista.size());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] - histb[i];
  return result;
}


// Unary mathematical functions
inline SimpleHist sqrt(const SimpleHist & hist) {
  SimpleHist result(hist.min(), hist.max(), hist.size());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = sqrt(hist[i]);
  return result;
}

// Unary mathematical functions
inline SimpleHist pow2(const SimpleHist & hist) {
  SimpleHist result(hist.min(), hist.max(), hist.size());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = hist[i]*hist[i];
  return result;
}

// cumulative output of histograms ----------------------------------------------
/// output the differential histogram
inline void output_differential(const SimpleHist & hist0, 
                        std::ostream * ostr = (&std::cout),
                        double norm = 1.0) {
  *ostr << "# under flow bin " << hist0.underflow() << std::endl;
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i)  << " "
          << hist0[i] * norm / (hist0.binhi(i) - hist0.binlo(i)) << std::endl;
  }
  *ostr << "over flow bin " << hist0.overflow() << std::endl;
}

/// output the cumulative histogram.
inline void output_cumulative(const SimpleHist & hist0,
                        std::ostream * ostr = (&std::cout),
                        double norm = 1.0) {
  double cumul = 0., correction = 0.;
  hist0.accumulate_bin_weight(hist0.underflow_bin(), norm, cumul, correction);
  *ostr << hist0.binlo(0) << " " << cumul + correction << std::endl;
  for (unsigned i = 0; i < hist0.size(); i++) {
    hist0.accumulate_bin_weight(i, norm, cumul, correction);
    *ostr << hist0.binhi(i) << " " << cumul + correction << std::endl;
  }
  hist0.accumulate_bin_weight(hist0.overflow_bin(), norm, cumul, correction);
  *ostr << "# with_overflow " << cumul + correction << std::endl;
}

/// output the inverse cumulative histogram.
inline void output_inverse_cumulative(const SimpleHist & hist0, 
                        std::ostream * ostr = (&std::cout),
                        double norm = 1.0) {
  double cumul = 0., correction = 0.;
  hist0.accumulate_bin_weight(hist0.overflow_bin(), norm, cumul, correction);
  *ostr << hist0.binhi(hist0.size() - 1) << " " << cumul + correction << std::endl;
  for (int i = hist0.size() - 1; i >= 0; i--) {
    hist0.accumulate_bin_weight(i, norm, cumul, correction);
    *ostr << hist0.binlo(i) << " " << cumul + correction << std::endl;
  }
  hist0.accumulate_bin_weight(hist0.underflow_bin(), norm, cumul, correction);
  *ostr << "# with_underflow " << cumul + correction << std::endl;
}
// -------------------------------------------------------------------

/// output the histogram to standard output -- an operator<< might
/// have seemed nice, but less easy to generalize to multiple
/// histograms; the output is multipled by the factor norm.
inline void output(const SimpleHist & hist0, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i]*norm << std::endl;
  }
}


// std::ostream & operator<<(std::ostream & ostr, const SimpleHist & hist) {
//   output(hist, &ostr);
//   return ostr;
// }

inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << std::endl;
  }
}

inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   const SimpleHist & hist2, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  assert(hist0.size() == hist2.size() && 
         hist0.min()  == hist2.min() &&
         hist0.max()  == hist2.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << hist2[i] * norm << " "
          << std::endl;
  }
}


inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   const SimpleHist & hist2, 
                   const SimpleHist & hist3, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  assert(hist0.size() == hist2.size() && 
         hist0.min()  == hist2.min() &&
         hist0.max()  == hist2.max());
  assert(hist0.size() == hist3.size() && 
         hist0.min()  == hist3.min() &&
         hist0.max()  == hist3.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << hist2[i] * norm << " "
          << hist3[i] * norm << " "
          << std::endl;
  }
}

inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   const SimpleHist & hist2, 
                   const SimpleHist & hist3, 
                   const SimpleHist & hist4, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  assert(hist0.size() == hist2.size() && 
         hist0.min()  == hist2.min() &&
         hist0.max()  == hist2.max());
  assert(hist0.size() == hist3.size() && 
         hist0.min()  == hist3.min() &&
         hist0.max()  == hist3.max());
  assert(hist0.size() == hist4.size() && 
         hist0.min()  == hist4.min() &&
         hist0.max()  == hist4.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << hist2[i] * norm << " "
          << hist3[i] * norm << " "
          << hist4[i] * norm << " "
          << std::endl;
  }
}

inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   const SimpleHist & hist2, 
                   const SimpleHist & hist3, 
                   const SimpleHist & hist4, 
                   const SimpleHist & hist5, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  assert(hist0.size() == hist2.size() && 
         hist0.min()  == hist2.min() &&
         hist0.max()  == hist2.max());
  assert(hist0.size() == hist3.size() && 
         hist0.min()  == hist3.min() &&
         hist0.max()  == hist3.max());
  assert(hist0.size() == hist4.size() && 
         hist0.min()  == hist4.min() &&
         hist0.max()  == hist4.max());
  assert(hist0.size() == hist5.size() && 
         hist0.min()  == hist5.min() &&
         hist0.max()  == hist5.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i)  << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << hist2[i] * norm << " "
          << hist3[i] * norm << " "
          << hist4[i] * norm << " "
          << hist5[i] * norm << " "
          << std::endl;
  }
}

inline void output(const SimpleHist *hists, 
                   const unsigned int nb_hist, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  unsigned int ih;
  for (ih=1;ih<nb_hist;ih++)
    assert(hists[0].size() == hists[ih].size() && 
	   hists[0].min()  == hists[ih].min() &&
	   hists[0].max()  == hists[ih].max());

  for (unsigned i = 0; i < hists[0].size(); i++) {
    *ostr << hists[0].binlo(i)  << " " 
          << hists[0].binmid(i) << " "
          << hists[0].binhi(i)  ;
    for (ih=0;ih<nb_hist;ih++)
      *ostr << " " << hists[ih][i] * norm;
    *ostr << std::endl;
  }
}

inline void output(const std::vector<SimpleHist*> &hists, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  unsigned int ih;
  unsigned int nb_hist = hists.size();
  for (ih=1;ih<nb_hist;ih++)
    assert(hists[0]->size() == hists[ih]->size() && 
	   hists[0]->min()  == hists[ih]->min() &&
	   hists[0]->max()  == hists[ih]->max());

  for (unsigned i = 0; i < hists[0]->size(); i++) {
    *ostr << hists[0]->binlo(i)  << " " 
          << hists[0]->binmid(i) << " "
          << hists[0]->binhi(i)  ;
    for (ih=0;ih<nb_hist;ih++)
      *ostr << " " << (*hists[ih])[i] * norm;
    *ostr << std::endl;
  }
}

inline void output(const std::vector<SimpleHist> &hists, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  unsigned int ih;
  unsigned int nb_hist = hists.size();
  for (ih=1;ih<nb_hist;ih++)
    assert(hists[0].size() == hists[ih].size() && 
	   hists[0].min()  == hists[ih].min() &&
	   hists[0].max()  == hists[ih].max());

  for (unsigned i = 0; i < hists[0].size(); i++) {
    *ostr << hists[0].binlo(i)  << " " 
          << hists[0].binmid(i) << " "
          << hists[0].binhi(i)  ;
    for (ih=0;ih<nb_hist;ih++)
      *ostr << " " << hists[ih][i] * norm;
    *ostr << std::endl;
  }
}

#endif // __SIMPLEHIST_HH__
