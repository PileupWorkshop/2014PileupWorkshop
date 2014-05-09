#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
//#include "SimpleHist.hh"

using namespace std;
using namespace fastjet;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Run this example using, for instance,
//
// ./example01 -hard ../sample-events/lhc14-pythia8-4C-dijet50-nev20.pu14.gz \
//             -pileup ../sample-events/lhc14-pythia8-4C-minbias-nev100.pu14.gz \
//             -npu 20 -nev 2 
//
// Other possible command line options are -R, -verbose (see code for role)
//
///////////////////////////////////////////////////////////////////////////////////////////////////


int main (int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nev = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  bool verbose = cmdline.present("-verbose");
  double R = cmdline.value<double>("-R",0.4);
  double rapmax = cmdline.value<double>("-rapmax",2.);
  
  // some definitions
  JetDefinition jet_def(antikt_algorithm,R);         // the jet definition....
  AreaDefinition area_def(active_area);              // the area definition
  cout << jet_def.description() << endl;            
  cout << area_def.description() << endl;           
  Selector sel_twohardest = SelectorNHardest(2);   // definition of a selector for the two hardest jets
  Selector sel_acceptance = SelectorAbsRapMax(rapmax);
  // define background estimator (grid type, does not need to recluster)
  GridMedianBackgroundEstimator * gmbge = new GridMedianBackgroundEstimator(rapmax+R,0.55);
  // define subtractor
  Subtractor sub(gmbge);
  
  // create mixer that will construct events by mixing hard and pileup
  // events read from files given from command line using 
  // -hard hard_events_file(.gz) -pilup pileup_events_file(.gz)
  EventMixer mixer(&cmdline);  
  
  // loop over events
  int iev = 0;
  while ( mixer.next_event() && iev < nev ) {
     // increment event number    
     iev++;
     
     // extract particles from event
     vector<PseudoJet> full_event = mixer.particles() ;
     cout << "\nEvent " << iev << endl;
     cout << "nPU = " << mixer.npu() << endl;
     
     // extract hard event
     vector<PseudoJet> hard_event, pileup_event;
     SelectorIsHard().sift(full_event, hard_event, pileup_event); 

     // cluster hard event only                                                     
     ClusterSequenceArea cs_hard(hard_event,jet_def,area_def);
     // cluster full event (hard + pileup)
     ClusterSequenceArea cs_full(full_event,jet_def,area_def);
          
     // write out two hardest jets in full event
     vector<PseudoJet> twohardest = sel_twohardest(sorted_by_pt(cs_full.inclusive_jets()));
     cout << "Full event" << endl;
     cout << "  jet 0: " << twohardest[0] << endl;
     cout << "  jet 1: " << twohardest[1] << endl;	 

     // subtract the two hardest jets
     // 1. feed particles to background estimator
     gmbge->set_particles(full_event);
     // 2. subtract
     vector<PseudoJet> subtracted_jets = sub(twohardest); 
     // write out jets in subtracted event
     cout << "Subtracted event" << endl;
     cout << "  jet 0: " << subtracted_jets[0] << endl;
     cout << "  jet 1: " << subtracted_jets[1] << endl;	 
         
     // write out two hardest jets in hard event
     twohardest = sel_twohardest(sorted_by_pt(cs_hard.inclusive_jets()));
     cout << "Hard event" << endl;
     cout << "  jet 0: " << twohardest[0] << endl;
     cout << "  jet 1: " << twohardest[1] << endl;	 
     
  }  
}
