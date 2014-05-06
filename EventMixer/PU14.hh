#ifndef __PU14_HH__
#define __PU14_HH__

#include "fastjet/PseudoJet.hh"
#include "HepPID/ParticleIDMethods.hh"

//----------------------------------------------------------------------
/// \class PU14
/// 
/// class that stores basic user information for the 2014 pileup
/// workshop.
///
class PU14 : public fastjet::PseudoJet::UserInfoBase {
public:
  /// construct 
  PU14(int pdg_id, int barcode = 0, int vertex = 0) {
    _pdg_id = pdg_id;
    _three_charge = HepPID::threeCharge(pdg_id);
    _barcode = barcode;
    _vertex = vertex;
  }
  
  /// return the pdg ID of the particle
  int pdg_id() const {return _pdg_id;}
  
  /// return its charge
  double charge() const {return _three_charge / 3.0;}

  /// return three times its charge (an integer)
  int three_charge() const {return _three_charge;}

private:
  int _pdg_id, _three_charge;
  int _barcode;
  int _vertex;
};



/// returns a Selector that is true for charged particles
fastjet::Selector SelectorIsCharged();
/// alternative name for charged-particle selector
fastjet::Selector SelectorCharged() {return SelectorIsCharged();}

/// returns a Selector that is true for particles that come from the
/// specified vertex
fastjet::Selector SelectorVertexNumber(int i);

/// returns a Selector that is true for particles from the hard vertex
fastjet::Selector SelectorIsHard()   {return SelectorVertexNumber(0);}
fastjet::Selector SelectorHard()   {return SelectorIsHard();}

/// returns a Selector that is true for particles from pileup vertices
fastjet::Selector SelectorIsPileup() {return !SelectorVertexNumber(0);}
fastjet::Selector SelectorPileup() {return !SelectorIsPileup();}

/// returns a Selector that is true for particles with a specific PDGId
fastjet::Selector SelectorPDGId(int i);

/// returns a Selector that is true for particles with a specific absolute PDGId
fastjet::Selector SelectorAbsPDGId(int i);


#endif
