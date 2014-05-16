#ifndef __CORRELATIONCOEFFICIENT_HH__
#define __CORRELATIONCOEFFICIENT_HH__

#include<cmath>

/// micro class to calculate the Pearson correlation coefficient between x and y
class CorrelationCoefficient {
public:

  /// default constructor
  CorrelationCoefficient() { _sumx = 0.0; _sumy = 0.0; _sumxx = 0.0; _sumyy = 0.0;
                             _sumxy = 0.0, _n=0;}
  
  /// add one entry
  inline void add(double x, double y) { _sumx += x;
                                        _sumy += y; 
                                        _sumxx += x*x;
                                        _sumyy += y*y;
                                        _sumxy += x*y;
                                        _n += 1;
                                      }
  
  /// alternative way to add one or more entries, for consistency with 
  /// SimpleHist, AveragingHist, etc
  inline void add_entry(double x, double y) { add(x,y); }

  /// return sum of x
  inline double sumx() const { return _sumx; }

  /// return sum of y
  inline double sumy() const { return _sumy; }

  /// return sumx2, second way for consistency with AveragingHist
  inline double sumxx() const { return _sumxx; }
  inline double sum_of_xsquares() const { return _sumxx; }

  /// return sumy2, second way for consistency with AveragingHist
  inline double sumyy() const { return _sumyy; }
  inline double sum_of_ysquares() const { return _sumyy; }

  /// return sumxy, second way for consistency with AveragingHist
  inline double sumxy() const { return _sumxy; }
  
  /// return number of events
  inline int n() const { return _n; }
  /// alternative way to return number of events, for consistency with
  /// SimpleHist, AveragingHist, etc
  inline int n_entries() const { return n(); }
  
  /// allow the user to reset the effective value of n
  inline void set_n(int n_in) { _n = n_in;}

  /// calculate and return averages (notation consistent with CorrelationHist)
  inline double avgx () const {return _sumx/_n;}
  inline double avgy () const {return _sumy/_n;}
  inline double avgxx () const {return _sumxx/_n;}
  inline double avgyy () const {return _sumyy/_n;}
  inline double avgxy() const {return _sumxy/_n;}
  inline double varx () const {return _sumxx/_n - pow(avgx(),2);}
  inline double vary () const {return _sumyy/_n - pow(avgy(),2);}
  
  /// calculate and return Pearson correlation coefficient
  // the line below should be multiplied by (n-1)/(n-2) to get an
  // unbiased estimate of the correlation when using a finite "sample"
  inline double r() const { if (_n <= 1) return 0.0;
                            double _avgx = avgx();
                            double _avgy = avgy();                            
			    return (avgxy()-_avgx*_avgy)/sqrt(std::abs(avgxx()-_avgx*_avgx)*std::abs(avgyy()-_avgy*_avgy));
			  }
  
double _sumx, _sumy, _sumxx, _sumyy, _sumxy;
int _n;

};

#endif // __CORRELATIONCOEFFICIENT_HH__
