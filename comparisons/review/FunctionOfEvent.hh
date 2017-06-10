#ifndef __FUNCTION_OF_EVENT_HH__
#define __FUNCTION_OF_EVENT_HH__

#include <fastjet/PseudoJet.hh>
#include <string>

//----------------------------------------------------------------------
/// abstract base class 
template<class T> class FunctionOfEvent {
public:
  virtual ~FunctionOfEvent(){}
  virtual T result(const std::vector<fastjet::PseudoJet> & ev) const = 0;
  T operator()(const std::vector<fastjet::PseudoJet> & ev) const {return result(ev);}
  virtual std::string description() const{ return "Unknown";}
};


//----------------------------------------------------------------------
/// wrapper that applies a * (b * event)
template<class T> class PairFoE : public FunctionOfEvent<T> {
public:
  PairFoE(FunctionOfEvent<T> * a, FunctionOfEvent<T> * b) : _a(a), _b(b) {}
  virtual T result(const std::vector<fastjet::PseudoJet> & ev) const {
    return (*_a)( (*_b)(ev) );
  }
  
private:
  FunctionOfEvent<T> * _a, * _b;
};

#endif  // __FUNCTION_OF_EVENT_HH__
