#include "EventMixer.hh"

using namespace std;

//----------------------------------------------------------------------
EventMixer::EventMixer(CmdLine * cmdline) : _cmdline(cmdline) {

  _hard_name = _cmdline->value<string>("-hard");
  _pileup_name = _cmdline->value<string>("-pileup");
  _npu = _cmdline->value("-npu", 20);

  _hard  .reset(new EventSource(_hard_name  ));
  _pileup.reset(new EventSource(_pileup_name));
}

//----------------------------------------------------------------------
bool EventMixer::next_event() {
  _particles.resize(0);
  

  // first get the hard event
  if (! _hard->append_next_event(_particles,0)) return false;

  for (int i = 1; i <= _npu; i++) {
    if (! _pileup->append_next_event(_particles,i)) return false;
  }

  return true;
}


//----------------------------------------------------------------------
string EventMixer::description() const {
  ostringstream ostr;
  ostr << "Event mixer using hard events from " << _hard_name 
       << " and " << _npu << " pileup events from " << _pileup_name;
  return ostr.str();
}
