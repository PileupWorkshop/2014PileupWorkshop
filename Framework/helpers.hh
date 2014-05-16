/// This file is provides a collection of helpers:
///
/// MasslessTransformer   a FastJet Transformer that sets particle masses to 0
///
///

#ifndef __PU14_HELPERS_HH__
#define __PU14_HELPERS_HH__


/// \class MasslessTransformer
/// a FastJet transformer that sets the mass of the particles/jets to
/// 0 (preserving all the other info, e.g. rapidity)
class MasslessTransformer : public fastjet::Transformer{
public:
  MasslessTransformer(){};
  virtual std::string description() const{ return "makes a particle massless, co
nserving rapidity";}
  virtual fastjet::PseudoJet result(const fastjet::PseudoJet &jet) const{
    fastjet::PseudoJet j = jet;
    j.reset_momentum(fastjet::PtYPhiM(jet.pt(), jet.rap(), jet.phi()));
    return j;
  }
};


#endif  // __PU14_HELPERS_HH__
