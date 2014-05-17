#ifndef __MATCHING_HH__
#define __MATCHING_HH__

#include <vector>
#include <limits>
#include <string>
#include <sstream>
#include "fastjet/PseudoJet.hh"


/// class to help match hard jets with full jets (e.g. pileup subtracted jets). 
///
/// Usage: 
///
/// - create an instance of the class with your chosen delta_R_max
///
/// - call set_full_jets(full_jets) to register a vector of full jets
///
/// - for each hard jet, call match(hard_jet); the class will return a
///   pointer to corresponding matched full jet, or 0 if there was no
///   match.
///
/// Matching is simple geometrical: for each hard jet it looks for the
/// nearest full jet.
///
class Matching {
public:
  /// constructor that sets the maximum allowed delta_R for the matching
  Matching(double delta_R_max = 0.3) : _delta_R_max(delta_R_max) {}

  /// set the list of full jets
  void set_full_jets(const std::vector<fastjet::PseudoJet> & full_jets) {_full_jets = full_jets;}

  /// Returns a pointer to the full jet that is closest in delta R to
  /// the supplied hard jet. If there is no full jet within
  /// delta_R_max, return a null pointer.
  ///
  /// The pointer is to a PseudoJet within a full-jet vector stored
  /// inside the matching class instance, i.e. it doesn't point inside
  /// the original full_jets vector that was supplied to
  /// set_ful_jets(..).
  const fastjet::PseudoJet * match(const fastjet::PseudoJet & hard_jet) const {
    double smallest_delta_R2 = std::numeric_limits<double>::max();
    const fastjet::PseudoJet * m = 0;
    for (unsigned i = 0; i < _full_jets.size(); i++) {
      double delta_R2 = _full_jets[i].squared_distance(hard_jet);
      if (delta_R2 < smallest_delta_R2) {
        m = & _full_jets[i];
        smallest_delta_R2 = delta_R2;
      }
    }
    if (smallest_delta_R2 < pow(_delta_R_max,2)) {return m;}
    else                                         {return 0;} 
  }

  /// returns the value of delta_R_max
  double delta_R_max() const {return _delta_R_max;}

  /// returns a textual description
  std::string description() const {
    std::ostringstream ostr;
    ostr << "Geometrical matching class with delta_R_max = " << delta_R_max();
    return ostr.str();
  }

private:
  double _delta_R_max;
  std::vector<fastjet::PseudoJet> _full_jets;
};

#endif // __MATCHING_HH__
