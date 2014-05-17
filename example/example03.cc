///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This example02 is meant to illustrate the implementation of a simple pileup subtraction
// (using FastJet's GridMedianBackGroundEstimator and Subtractor classes, and the calculation
// and output of some quality measures, namely the average offset between subtracted and hard,
// its dispersion, and the correlation coefficient.
//
// Rho rescaling using a pre-determined function of rapidity is also performed, and
// can be switched off using the command line option "-norescale"
//
// A simple matching between subtracted and hard, based on a delta_R distance <= 0.2, is
// performed before calculation of quality measures.
//
// Run this example using, for instance,
//
// ./example03 -hard ../sample-events/lhc14-pythia8-4C-dijet50-nev20.pu14.gz \
//             -pileup ../sample-events/lhc14-pythia8-4C-minbias-nev100.pu14.gz \
//             -massless -npu 5 -nev 20
//
// Other possible command line options are documented in the code below
//
// TODO: - add a -out option? 
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "AverageAndError.hh"
#include "CorrelationCoefficient.hh"
#include "ProfileHist.hh"
#include "Matching.hh"

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

#include <iomanip>      // std::setprecision
#include <cassert>
#include <cmath>

using namespace std;
using namespace fastjet;

////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  cout << "# " << cmdline.command_line() << "\n#" << endl;

  // inputs read from command line
  int nev = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  int maxprintout = cmdline.value<int>("-maxprintout",1);  

  // particle and jet selection & R
  double particle_rapmax = cmdline.value<double>("-particle-rapmax",4.);
  double jet_rapmax = cmdline.value<double>("-jet-rapmax",2.5);
  double jet_ptmin  = cmdline.value<double>("-jet-ptmin",20);
  double R = cmdline.value<double>("-R",0.4);

  // matching condition
  double maxdeltaR = cmdline.value<double>("-maxdeltaR",0.3);
  Matching matching(maxdeltaR);
  cout << "# matching: " << matching.description() << endl;

  // rapidity rescaling of rho
  bool rescale = ! cmdline.present("-norescale");
  cout << "# rapidity rescaling for rho = " << rescale << endl;

  
  // some definitions
  JetDefinition jet_def(antikt_algorithm,R);             // the jet definition
  AreaDefinition area_def(active_area_explicit_ghosts);  // the area definition
  cout << "# jet_def: "  << jet_def.description() << endl;            
  cout << "# area_def: "  << area_def.description() << endl;           

  // select 2 hardest > 20, 100, 500 in hard event, and then |y| < 2.5
  Selector sel_jets = SelectorAbsRapMax(jet_rapmax)*SelectorPtMin(jet_ptmin)*SelectorNHardest(2);  
  cout << "# sel_jets: " << sel_jets.description() << endl;

  // take particles only within |y|<4
  Selector sel_particles = SelectorAbsRapMax(particle_rapmax);
  cout << "# sel_particles: " << sel_particles.description() << endl;

  AverageAndError npu, njets, offset, matching_efficiency, pass_fraction;
  CorrelationCoefficient subhardcorr;
  ProfileHist offset_v_rapidity(-jet_rapmax,jet_rapmax,0.50);

  // define background estimator (grid type, does not need to recluster)
  GridMedianBackgroundEstimator * gmbge = new GridMedianBackgroundEstimator(particle_rapmax,0.55);
  // define (and then apply) function for rapidity rescaling of rho.
  // NB These parameters have actually been determined at 13 TeV, but they vary slowly
  // and should therefore also be aproproate for 14 TeV
  FunctionOfPseudoJet<double> * rescaling = new BackgroundRescalingYPolynomial(1.1685397, 0, -0.0246807, 0, 5.94119e-05);
  if ( rescale) { gmbge -> set_rescaling_class(rescaling); }
  // define subtractor
  Subtractor sub(gmbge);
  cout << "# subtractor: "  << sub.description() << endl;
  
  // create mixer that will construct events by mixing hard and pileup
  // events read from files given from command line using 
  // -hard hard_events_file(.gz) -pileup pileup_events_file(.gz)
  EventMixer mixer(&cmdline);  

  // make sure there are no unused command-line arguments
  cmdline.assert_all_options_used();
  
  // loop over events ----------------------------------------------------------------------
  int iev = 0;
  while ( mixer.next_event() && iev < nev ) {
     // increment event number    
     iev++;
     
     // extract particles from event
     vector<PseudoJet> full_event = sel_particles(mixer.particles());
     if ( iev <= maxprintout ) { cerr << "\nEvent " << iev << endl; }
     if ( iev <= maxprintout ) { cerr << "nPU = " << mixer.npu() << endl; }
     npu.add_entry(mixer.npu());
     
     // extract hard event
     vector<PseudoJet> hard_event, pileup_event;
     SelectorIsHard().sift(full_event, hard_event, pileup_event); 

     // cluster hard event only (Note: area may not be needed here)                                      
     ClusterSequenceArea cs_hard(hard_event,jet_def,area_def);
     // select two hardest jets in hard event
     vector<PseudoJet> hard_jets = sel_jets(sorted_by_pt(cs_hard.inclusive_jets()));
     if ( iev <= maxprintout ) { cerr << "Hard event" << endl; }
     for (unsigned int i=0; i < hard_jets.size(); i++) {
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": " << hard_jets[i] << endl; }
     }
     
     // if no hard jets pass the selection, continue on to the next event
     if (hard_jets.size() == 0) {
       pass_fraction += 0;
       continue;
     }
     pass_fraction += 1;

     // cluster full event (hard + pileup)
     ClusterSequenceArea cs_full(full_event,jet_def,area_def);
          
     // get all jets in full event
     vector<PseudoJet> full_jets = sorted_by_pt(cs_full.inclusive_jets());
     if ( iev <= maxprintout ) { cerr << "Full event" << endl; }
     for (unsigned int i=0; i < 4U && i < full_jets.size(); i++) {
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": "  << full_jets[i] << endl; }
     }
     
     // subtract the full jets
     // 1. feed particles to background estimator
     gmbge->set_particles(full_event);
     // 2. subtract
     vector<PseudoJet> subtracted_jets = sub(full_jets); 
     // write out jets in subtracted event
     if ( iev <= maxprintout ) { cerr << "Subtracted event" << endl; }
     for (unsigned int i=0; i < 4U && i < subtracted_jets.size(); i++) {
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": "  << subtracted_jets[i] << endl; }
     }
        
     // calculate quality measures, offset and dispersion and correlation coefficient,
     // for the jet transverse momentum.
     // Also fill an histogram of offset v. rapidity, to appreciate the effect
     // of rho rescaling

     // set up the set of full/subtracted jets from which to match
     matching.set_full_jets(subtracted_jets);
     
     // run over the hard jets
     for (unsigned int i=0; i < hard_jets.size(); i++) {
       // for each hard jet, find the corresponding full/subtracted jet that matches
       // (if any)
       const PseudoJet * match = matching.match(hard_jets[i]);
       if (match) {
         matching_efficiency += 1.0;
         double deltapt = match->pt() - hard_jets[i].pt();
         offset.add_entry(deltapt);
         subhardcorr.add_entry(match->pt(),hard_jets[i].pt());
         offset_v_rapidity.add_entry(hard_jets[i].rap(),deltapt);
       } else {
         matching_efficiency += 0.0;
       }
     }
     
     // keep track of number of jets above 20 GeV within the jet rapidity window
     njets += ( SelectorPtMin(20.0)
               && SelectorAbsRapMax(jet_rapmax) ).count(subtracted_jets);
     
     
  }  // end loop over events

  // output quality measures as a function of <npu>
  cout << "# fraction of events passing the basic jet cuts = " << pass_fraction.average() << " +- " << pass_fraction.error() << endl;
  cout << "# npu    jet_ptmin     <DeltaO>           sigma_DeltaO     corr.coeff.     njets>20GeV        match_eff " << endl;     
  cout << setprecision(4) 
       << setw(4) << npu.average()    << "    "
       << setw(6) << jet_ptmin << "    "
       << setw(6) << offset.average() << " +- " << setw(6) << offset.error() << "    "
       << setw(6) << offset.sd()      << " +- " << setw(6) << offset.error_on_sd() << "   "
       << setw(6) << subhardcorr.r()  << "    "
       << setw(6) << njets.average() << " +- " << setw(6) << njets.error() << "    "
       << setw(6) << matching_efficiency.average() << " +- " << setw(6) << matching_efficiency.error()
       << endl;

  // output histograms
  cout << "\n\n# offset_v_rapidity" << endl;
  cout << "# binlo binmid binhi avg stddev err avgsquares" << endl;
  output_noNaN(offset_v_rapidity);

  // free allocated memory
  delete rescaling;
  delete gmbge;
}

