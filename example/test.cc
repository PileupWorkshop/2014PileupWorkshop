#include "EventSource.hh"
#include "CmdLine.hh"
#include "PU14.hh"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  EventSource events(cmdline.value<string>("-hard"));
  
  vector<PseudoJet> event;
  event.resize(0);
  bool out = events.append_next_event(event);
  cout << out << " " << event.size() << " " << endl;
  
  PseudoJet sum;
  double net_charge = 0;
  for (unsigned i = 0; i < event.size(); i++) {
    const PseudoJet & p  = event[i];
    cout << i << " "
         << p.px() << " " << p.py() << " " << p.pz() << " " << p.E() << " " << p.m() << " " 
         << p.user_info<PU14>().pdg_id() << " " << p.user_info<PU14>().charge() << endl;
    sum += p;
    net_charge += p.user_info<PU14>().charge();
  }
  cout << "total momentum: " << sum << endl;
  cout << "total charge: " << net_charge << endl;
}
