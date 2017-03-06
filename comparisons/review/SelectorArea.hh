#ifndef __SELECTOR_AREA_HH__
#define __SELECTOR_AREA_HH__

#include <fastjet/Selector.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <sstream>
#include <iostream>

FASTJET_BEGIN_NAMESPACE

// implements a selector that keeps the jets that are above
//  cr * rho * A  +  cs * sigma * sqrt(A)
// where cr and cs are user-defined constants
// 
// rho and sigma will be estimated using a BackgroundEstiumator passed
// as an argument
class SW_AreaRhoThreshold : public SelectorWorker {
public:
  /// ctor with specification of
  ///   \param bge     a pointer to the background estimator to use (has
  ///                  to be available!)
  ///   \param crho    the factor of rho*A in the cut
  ///   \param csigms  the factor of sigma*sqrt(A) in the cut
  SW_AreaRhoThreshold(BackgroundEstimatorBase *bge, 
                      double crho, double csigma)
    : _bge(bge), _crho(crho), _csigma(csigma){}

  /// return true if the jet carries a large enough fraction of the reference.
  /// Throw an error if the reference is not initialised.
  virtual bool pass(const PseudoJet & jet) const {
    // otherwise, just call that method on the jet
    double a = jet.area();
    return (jet.perp() >= _crho  *_bge->rho(jet)  *a 
                        + _csigma*_bge->sigma(jet)*sqrt(a));
  }
  
  /// returns a description of the worker
  virtual std::string description() const {
    std::ostringstream ostr;
    ostr << "pt >= " << _crho << "* rho(jet)*area + " << _csigma << "* sigma(jet)*sqrt(area)";
    return ostr.str();
  }

protected:
  mutable BackgroundEstimatorBase * _bge;
  double _crho, _csigma;
};

// select objects that carry at least a fraction "fraction" of the reference jet
// (Note that this selectir takes a reference)
Selector SelectorAreaRhoThreshold(BackgroundEstimatorBase *bge, double crho, double csigma){
  return Selector(new SW_AreaRhoThreshold(bge, crho, csigma));
}

FASTJET_END_NAMESPACE

#endif // __SELECTOR_AREA_HH__
