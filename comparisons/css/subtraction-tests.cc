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

/// a specific class that handles Matteo's format
///
/// Usage: 
///   declare it with a name
///   for a given entry add
///     data.add_entry();   // for un-matched 
///     data.add_entry(hard_value, full_value); // for matched
class MatteoAverageAndError{
public:
  MatteoAverageAndError(){}
  
  AverageAndError match_eff, delta, njets;
  CorrelationCoefficient correl;
  
  void add_entry(){
    match_eff += 0.0;
  }

  void add_entry(double hard, double full){
    match_eff += 1.0;
    delta += full-hard;
    correl.add_entry(hard, full);
  }

};

/// a class for output management
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

  map<string,MatteoAverageAndError> matteos;

  /// output once in a while
  /// Arguments are:
  ///   iev    the number of events so far
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

    // Mattteo's output of the form 
    //   name  <delta>  sigma_delta   corr.coef  match.eff  match.eff_error
    ostr << "#----------------------------------------------------------------------" << endl;
    ostr << "# name                         <delta>             sigma_delta       corr.coef     Njets(>20GeV)        match.eff   " << endl;
    for(map<string,MatteoAverageAndError>::iterator result = matteos.begin(); 
        result != matteos.end(); result++) {
      ostr << "# " << setw(20) << result->first << "    " 
           << setw(6) << result->second.delta.average()     << " +- " << setw(6) << result->second.delta.error()       << "    "
           << setw(6) << result->second.delta.sd()          << " +- " << setw(6) << result->second.delta.error_on_sd() << "    "
           << setw(6) << result->second.correl.r()          << "    "
           << setw(6) << result->second.njets.average()     << " +- " << setw(6) << result->second.njets.error()       << "    "
           << setw(6) << result->second.match_eff.average() << " +- " << setw(6) << result->second.match_eff.error()   << endl;
    }
    ostr << "#----------------------------------------------------------------------" << endl;

    // histograms (normalised to the number of hard jets) 
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
  cout << "# " << cmdline.command_line() << "\n#" << endl;

  // inputs read from command line
  int nev = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  int maxprintout = cmdline.value<int>("-maxprintout",1);  

  // particle and jet selection & R
  double particle_rapmax = cmdline.value<double>("-particle-rapmax",4.);
  double jet_rapmax = cmdline.value<double>("-jet-rapmax",2.5);
  double jet_ptmin  = cmdline.value<double>("-jet-ptmin",20);
  double R = cmdline.value<double>("-R",0.4);

  // some containier for what we want to output
  OutputInfo output(cmdline.value<string>("-out"));
  int iev_periodic=10;

  output.header << "# " << cmdline.command_line() << "\n#" << endl;
  output.header << "# time now = " << cmdline.time_stamp() << endl;

  // matching condition
  double maxdeltaR = cmdline.value<double>("-maxdeltaR",0.3);
  Matching matching(maxdeltaR);
  cout << "# matching: " << matching.description() << endl;
  output.header << "# matching: " << matching.description() << endl;

  // rapidity rescaling of rho
  bool rescale = ! cmdline.present("-norescale");
  cout << "# rapidity rescaling for rho = " << rescale << endl;
  output.header << "# rapidity rescaling for rho = " << rescale << endl;
  
  // some definitions
  JetDefinition jet_def(antikt_algorithm,R);             // the jet definition
  AreaDefinition area_def(active_area_explicit_ghosts);  // the area definition
  cout << "# jet_def: "  << jet_def.description() << endl;            
  cout << "# area_def: "  << area_def.description() << endl;           
  output.header << "# jet_def: "  << jet_def.description() << endl;            
  output.header << "# area_def: "  << area_def.description() << endl;           

  // select 2 hardest > 20, 100, 500 in hard event, and then |y| < 2.5
  Selector sel_jets = SelectorAbsRapMax(jet_rapmax)*SelectorPtMin(jet_ptmin)*SelectorNHardest(2);  
  cout << "# sel_jets: " << sel_jets.description() << endl;
  output.header << "# sel_jets: " << sel_jets.description() << endl;

  // take particles only within |y|<4
  Selector sel_particles = SelectorAbsRapMax(particle_rapmax);
  cout << "# sel_particles: " << sel_particles.description() << endl;
  output.header << "# sel_particles: " << sel_particles.description() << endl;

  AverageAndError npu;
  //ProfileHist offset_v_rapidity(-jet_rapmax,jet_rapmax,0.50);

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
  output.header << "# subtractor: "  << sub.description() << endl;
  
  // create mixer that will construct events by mixing hard and pileup
  // events read from files given from command line using 
  // -hard hard_events_file(.gz) -pileup pileup_events_file(.gz)
  EventMixer mixer(&cmdline);  
  output.header << "# Mixer: "  << mixer.description() << endl;

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

    // set up the set of full/subtracted jets from which to match
    matching.set_full_jets(subtracted_jets);
     
    // run over the hard jets
    for (unsigned int i=0; i < hard_jets.size(); i++) {
      // for each hard jet, find the corresponding full/subtracted jet that matches
      // (if any)
      const PseudoJet * match = matching.match(hard_jets[i]);
      if (match){
        output.matteos["pt_areasub"].add_entry(hard_jets[i].pt(), match->pt());
      } else {
        output.matteos["pt_areasub"].add_entry();
      }
    }
     
    // keep track of number of jets above 20 GeV within the jet rapidity window
    output.matteos["pt_areasub"].njets += ( SelectorPtMin(20.0)
                                            && SelectorAbsRapMax(jet_rapmax) ).count(subtracted_jets);
     

    // output from time to time
    if (iev % iev_periodic == 0){
       if (iev == 15*iev_periodic) iev_periodic*=10;
       output.write(iev);
    }
  }  // end loop over events

  //// output histograms
  //cout << "\n\n# offset_v_rapidity" << endl;
  //cout << "# binlo binmid binhi avg stddev err avgsquares" << endl;
  //output_noNaN(offset_v_rapidity);

  output.write(iev);

  // free allocated memory
  delete rescaling;
  delete gmbge;
}

