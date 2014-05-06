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

#endif
