///////////////////////////////////////////////////////////////////////////////////////////////////
//
// test of subtraction performance done by the FJ authors
//
// We include the following methods:
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "AverageAndError.hh"
#include "CorrelationCoefficient.hh"
#include "SimpleHist.hh"
#include "ProfileHist.hh"
#include "Matching.hh"

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

#include <iomanip>      // std::setprecision
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;
using namespace fastjet;

class OutputInfo{
public:
  OutputInfo(string outfile) : _outfile(outfile){}

  // what will get 
  string _outfile;
  ostringstream header;

  // the various infos
  map<string,SimpleHist> hists; //< normalised to the bin size
  map<string,ProfileHist> profiles;
  map<string,AverageAndError> averages;

  // output once in a while
  void write(const unsigned int iev){
    ofstream ostr(_outfile.c_str());

    cout << "Output of results for nev = " << iev << endl;

    ostr << header.str();
    ostr << "# nev = " << iev << endl;
  
    // write out averages
    for(map<string,AverageAndError>::iterator result = averages.begin(); 
        result != averages.end(); result++) {
      ostr << "# " << result->first << ": average = " << result->second.average()
           << " stddev = " << result->second.sd() << endl;
    }

    // histograms (normalised to the number of events) 
    ostr << "#------------------------------" << endl;
    for(map<string,SimpleHist>::iterator hist = hists.begin(); 
	hist != hists.end(); hist++) {
      ostr << "# " << hist->first << " [1/x-axis units] " << endl;
      output(hist->second, &ostr, 1.0/iev/hist->second.binsize());
      ostr << endl << endl;
    }

    // output profiles
    for(map<string,ProfileHist>::iterator hist = profiles.begin(); 
        hist != profiles.end(); hist++) {
      ostr << "# " << hist->first << " [average within each bin; stddev; error on avg; average of squares] " << endl;
      output_noNaN(hist->second, &ostr);
      ostr << endl << endl;
    }
  }
};

////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nev = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  int maxprintout = cmdline.value<int>("-maxprintout",1);  
  double R = cmdline.value<double>("-R",0.4);
  double rapmax = cmdline.value<double>("-rapmax",3.);
  double maxdeltaR = cmdline.value<double>("-maxdeltaR",0.3);
  bool rescale = ! cmdline.present("-norescale");

  // some containier for what we want to output
  OutputInfo output(cmdline.value<string>("-out"));
  int iev_periodic=10;

  cout << "# " << cmdline.command_line() << "\n#" << endl;
  output.header << "# " << cmdline.command_line() << "\n#" << endl;
  output.header << "# time now = " << cmdline.time_stamp() << endl;
  
  // some definitions
  JetDefinition jet_def(antikt_algorithm,R);             // the jet definition
  AreaDefinition area_def(active_area_explicit_ghosts);  // the area definition
  cout << "# jet_def: "  << jet_def.description() << endl;            
  cout << "# area_def: "  << area_def.description() << endl;           
  // selects two hardest jets in event, and THEN select only those within rapmax-R
  Selector sel_jets = SelectorAbsRapMax(rapmax-R)*SelectorNHardest(2);  
  cout << "# sel_jets: " << sel_jets.description() << endl;
  AverageAndError npu, offset, matching_efficiency;
  CorrelationCoefficient subhardcorr;
  ProfileHist offset_v_rapidity(-rapmax,rapmax,0.50);
  
  Matching matching(maxdeltaR);

  // define background estimator (grid type, does not need to recluster)
  GridMedianBackgroundEstimator * gmbge = new GridMedianBackgroundEstimator(rapmax,0.55);
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
  // -hard hard_events_file(.gz) -pilup pileup_events_file(.gz)
  EventMixer mixer(&cmdline);  

  // make sure there are no unused command-line arguments
  cmdline.assert_all_options_used();
  
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

     // cluster hard event only (Note: area may not be needed here)                                      
     ClusterSequenceArea cs_hard(hard_event,jet_def,area_def);
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
        
     // select two hardest jets in hard event
     vector<PseudoJet> hard_jets = sel_jets(sorted_by_pt(cs_hard.inclusive_jets()));
     if ( iev <= maxprintout ) { cerr << "Hard event" << endl; }
     for (unsigned int i=0; i < hard_jets.size(); i++) {
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": " << hard_jets[i] << endl; }
     }
     
     // calculate quality measures, offset and dispersion and correlation coefficient,
     // for the jet transverse momentum.
     // Also fill an histogram of offset v. rapidity, to appreciate the effect
     // of rho rescaling
     matching.set_full_jets(subtracted_jets);
     for (unsigned int i=0; i < hard_jets.size(); i++) {
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

     // output from time to time
     if (iev % iev_periodic == 0){
       if (iev == 15*iev_periodic) iev_periodic*=10;
       output.write(iev);
     }
     
  }  // end loop over events

  // output quality measures as a function of <npu>
  cout << "# npu  <DeltaO>    sigma_DeltaO    corr.coeff.     match_eff    match_eff_error" << endl;     
  cout << setprecision(9) 
       << npu.average()    << "    " 
       << offset.average() << "    " 
       << offset.sd()      << "    " 
       << subhardcorr.r()  << "    "
       << setw(9) << matching_efficiency.average()  << "    "
       << setw(9) << matching_efficiency.error()    << "    "
       << endl;

  // output histograms
  cout << "\n\n# offset_v_rapidity" << endl;
  cout << "# binlo binmid binhi avg stddev err avgsquares" << endl;
  output_noNaN(offset_v_rapidity);

  output.write(iev);

  // free allocated memory
  delete rescaling;
  delete gmbge;
}

