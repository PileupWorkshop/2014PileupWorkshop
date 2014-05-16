#ifndef __EVENTMIXER_HH__
#define __EVENTMIXER_HH__

#include "CmdLine.hh"
#include "EventSource.hh"

//----------------------------------------------------------------------
/// \class EventMixer
/// 
/// Class for that produces merged events from a mix of individual
/// hard and pileup interactions. It deduces the hard and pileup event
/// filenames from the command line, with options:
///
///   -hard HardFileName  -pileup PileupFileName  -npu  NPU
///
/// Currently the number of pileup events is fixed (not Poisson
/// distributed).
class EventMixer {
public:
  EventMixer(CmdLine * cmdline);
  ~EventMixer() {}

  /// causes the next event to be read in and mixed (hard + multiple pileup).
  /// Returns true if it successfully produced the event, false otherwise.
  bool next_event(); 

  /// returns a reference to vector of particles in the last event
  /// that was read in
  const std::vector<fastjet::PseudoJet> & particles() const {return _particles;}

  /// returns the number of pileup events generated in the last mixed event 
  int npu() const {return _npu;}

  std::string description() const;

private:
  CmdLine * _cmdline;
  std::string _hard_name, _pileup_name;
  std::auto_ptr<EventSource> _hard, _pileup;
  int _npu;
  std::vector<fastjet::PseudoJet> _particles;

};

#endif  // __EVENTMIXER_HH__

