///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This example02 is meant to illustrate the implementation of a simple pileup subtraction
// (using FastJet's GridMedianBackGroundEstimator and Subtractor classes, and the calculation
// and output of some quality measures, namely the average offset between subtracted and hard,
// its dispersion, and the correlation coefficient.
//
// Run this example using, for instance,
//
// ./example02 -hard ../sample-events/lhc14-pythia8-4C-dijet50-nev20.pu14.gz \
//             -pileup ../sample-events/lhc14-pythia8-4C-minbias-nev100.pu14.gz \
//             -npu 20 -nev 2 
//
// Other possible command line options are -R, -rapmax
//
// TODO: add rapidity rescaling, add matching between jets, add a -out option
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "AverageAndError.hh"
#include "CorrelationCoefficient.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include <iomanip>      // std::setprecision

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nev = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  int maxprintout = cmdline.value<int>("-maxprintout",1);  
  //bool verbose = cmdline.present("-verbose");
  double R = cmdline.value<double>("-R",0.4);
  double rapmax = cmdline.value<double>("-rapmax",3.);
  cout << "# " << cmdline.command_line() << "\n#" << endl;
  
  // some definitions
  JetDefinition jet_def(antikt_algorithm,R);         // the jet definition
  AreaDefinition area_def(active_area);              // the area definition
  cout << "# "  << jet_def.description() << endl;            
  cout << "# "  << area_def.description() << endl;           
  // selects two hardest jets in event, and THEN select only those within rapmax-R
  Selector sel_jets = SelectorAbsRapMax(rapmax-R)*SelectorNHardest(2);  
  AverageAndError npu, offset;
  CorrelationCoefficient subhardcorr;
  
  // define background estimator (grid type, does not need to recluster)
  GridMedianBackgroundEstimator * gmbge = new GridMedianBackgroundEstimator(rapmax,0.55);
  // define subtractor
  Subtractor sub(gmbge);
  cout << "# "  << sub.description() << endl;
  
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
     if ( iev <= maxprintout ) { cerr << "\nEvent " << iev << endl; }
     if ( iev <= maxprintout ) { cerr << "nPU = " << mixer.npu() << endl; }
     npu.add_entry(mixer.npu());
     
     // extract hard event
     vector<PseudoJet> hard_event, pileup_event;
     SelectorIsHard().sift(full_event, hard_event, pileup_event); 

     // cluster hard event only                                                     
     ClusterSequenceArea cs_hard(hard_event,jet_def,area_def);
     // cluster full event (hard + pileup)
     ClusterSequenceArea cs_full(full_event,jet_def,area_def);
          
     // write out the selected jets in full event
     vector<PseudoJet> full_jets = sel_jets(sorted_by_pt(cs_full.inclusive_jets()));
     if ( iev <= maxprintout ) { cerr << "Full event" << endl; }
     for (unsigned int i=0; i < full_jets.size(); i++) {
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": "  << full_jets[i] << endl; }
     }
     
     // subtract the selected full jets
     // 1. feed particles to background estimator
     gmbge->set_particles(full_event);
     // 2. subtract
     vector<PseudoJet> subtracted_jets = sub(full_jets); 
     // write out jets in subtracted event
     if ( iev <= maxprintout ) { cerr << "Subtracted event" << endl; }
     for (unsigned int i=0; i < subtracted_jets.size(); i++) {
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": "  << subtracted_jets[i] << endl; }
     }
        
     // select jets in hard event
     vector<PseudoJet> hard_jets = sel_jets(sorted_by_pt(cs_hard.inclusive_jets()));
     if ( iev <= maxprintout ) { cerr << "Hard event" << endl; }
     for (unsigned int i=0; i < hard_jets.size(); i++) {
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": " << hard_jets[i] << endl; }
     }
     
     // calculate quality measures, offset and dispersion and correlation coefficient,
     // for the jet transverse momentum
     for (unsigned int i=0; i < subtracted_jets.size(); i++) {
           offset.add_entry(subtracted_jets[i].pt() - hard_jets[i].pt());
	   subhardcorr.add_entry(subtracted_jets[i].pt(),hard_jets[i].pt());
     }
     
  }  // end loop over events

  // output quality measures as a function of <npu>
  cout << "# npu  <DeltaO>    sigma_DeltaO    corr.coeff." << endl;     
  cout << setprecision(9) << npu.average()    << "    " 
                          << offset.average() << "    " 
			  << offset.sd()      << "    " 
			  << subhardcorr.r() << endl;
}
