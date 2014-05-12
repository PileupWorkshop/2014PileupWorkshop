#ifndef __PROFILEHIST_HH__
#define __PROFILEHIST_HH__

#include "SimpleHist.hh"

/// class which will provide average of entries going into each bin
class ProfileHist {
public:

  ProfileHist() {}

  ProfileHist(double minv, double maxv, int n) {
    declare(minv, maxv, n);
  }

  ProfileHist(double minv, double maxv, unsigned int n) {
    declare(minv, maxv, n);
  }

  ProfileHist(double minv, double maxv, double bin_size) {
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
    _sum.declare(minv, maxv, n);
    _sum2.declare(minv, maxv, n);
    _weights.declare(minv, maxv, n);
    _nentries.declare(minv, maxv, n);
  }

  void add_entry(double v, double contents, double weight = 1.0) {
    _sum.add_entry (v, contents * weight);
    _sum2.add_entry(v, contents*contents * weight);
    _weights.add_entry(v, weight);
    _nentries.add_entry(v, 1.0);
  };

  double binlo (int i) const {return _sum.binlo(i);}
  double binhi (int i) const {return _sum.binhi(i);}
  double binmid(int i) const {return _sum.binmid(i);}
  double binsize()     const {return _sum.binsize();}
  unsigned int bin(double v) const {return _sum.bin(v);}

  double min() const {return _sum.min();};
  double max() const {return _sum.max();};
  /// returns the size of the histogram proper
  unsigned int size() const {return _sum.size();};
  /// returns the size of the histogram plus outflow bin
  unsigned int outflow_size() const {return _sum.outflow_size();};
  
  /// return the average of the current bin
  double average(int i) const {
    return _sum[i] / _weights[i];
  }

  /// return the average of the current bin
  double average_of_squares(int i) const {
    return _sum2[i] / _weights[i];
  }


  /// return the standard deviation of the current bin
  double stddev(int i) const {
   if (_nentries[i] > 1) {
     double avg = average(i);
     double result = 
          sqrt(abs(_sum2[i]/_weights[i] - avg*avg));
     return result;
   } else {
      return 0.0;
   }   
  }

  /// return the error on the average of the current bin
  double error(int i) const {
      return stddev(i) / sqrt(_nentries[i]);
  }

  /// return the number of entries in the current bin
  int nentries(int i) const {
    return int(_nentries[i]);
  }


private:
  SimpleHist _sum, _sum2, _weights, _nentries;
};

//----------------------------------------------------------------------
inline void output_noNaN(const ProfileHist hist, 
                   std::ostream * ostr = (&std::cout)) {
  
  for (unsigned i = 0; i < hist.size(); i++) {
    *ostr << hist.binlo(i)  << " " 
          << hist.binmid(i) << " "
          << hist.binhi(i) << " ";
    if (hist.nentries(i) > 0) {
      *ostr << hist.average(i) << " " << hist.stddev(i) << " " << hist.error(i) << " " <<  hist.average_of_squares(i) << std::endl;
    } else {
      *ostr << "0.0 0.0 0.0 0.0" << std::endl;
    }

  }
}


//----------------------------------------------------------------------
inline void output_avg_and_sqr_noNaN(const ProfileHist hist, 
                   std::ostream * ostr = (&std::cout)) {
  
  for (unsigned i = 0; i < hist.size(); i++) {
    *ostr << hist.binlo(i)  << " " 
          << hist.binmid(i) << " "
          << hist.binhi(i) << " ";
    if (hist.nentries(i) > 0) {
      *ostr << hist.average(i) << " " << hist.average_of_squares(i);
    } else {
      *ostr << "0.0 0.0";
    }

    *ostr << std::endl;
  }
}


#endif //  __PROFILEHIST_HH__
