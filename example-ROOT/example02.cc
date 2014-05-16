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
// This example store some output into a TTree, 
// the output file name is example02.root
//
// A simple matching between subtracted and hard, based on a delta_R distance <= 0.2, is
// performed before calculation of quality measures.
//
// Run this example using, for instance,
//
// ./example02 -hard ../sample-events/lhc14-pythia8-4C-dijet50-nev20.pu14.gz \
//             -pileup ../sample-events/lhc14-pythia8-4C-minbias-nev100.pu14.gz \
//             -npu 20 -nev 2 
//
// Other possible command line options are -R, -rapmax, -maxdeltaR, -maxprintout
//
// TODO: - add a -out option? 
//       - think through jet selection
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "AverageAndError.hh"
#include "CorrelationCoefficient.hh"
#include "ProfileHist.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include <iomanip>      // std::setprecision
#include <cassert>

#include "TTree.h"
#include "TFile.h"


using namespace std;
using namespace fastjet;

/// very dumb matching routine, expected to work with at most two jets.
/// It returns true if it manages to return at least one valid
/// matching between subtracted and hard, using a deltaR <= maxdeltaR rule.
/// If needed, it performs approprate swaps.
/// If only one matching is found the sub vector is resized to 1.
/// If no matching is found it returns false.
///
bool match(vector<PseudoJet> & sub, vector<PseudoJet> & hard, double maxdeltaR);


////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nev = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  int maxprintout = cmdline.value<int>("-maxprintout",1);  
  double R = cmdline.value<double>("-R",0.4);
  double rapmax = cmdline.value<double>("-rapmax",3.);
  double maxdeltaR = cmdline.value<double>("-maxdeltaR",0.2);
  bool rescale = ! cmdline.present("-norescale");
  cout << "# " << cmdline.command_line() << "\n#" << endl;
  
  // some definitions
  JetDefinition jet_def(antikt_algorithm,R);         // the jet definition
  AreaDefinition area_def(active_area);              // the area definition
  cout << "# jet_def: "  << jet_def.description() << endl;            
  cout << "# area_def: "  << area_def.description() << endl;           
  // selects two hardest jets in event, and THEN select only those within rapmax-R
  Selector sel_jets = SelectorAbsRapMax(rapmax-R)*SelectorNHardest(2);  
  cout << "# sel_jets: " << sel_jets.description() << endl;
  AverageAndError npu, offset;
  CorrelationCoefficient subhardcorr;
  ProfileHist offset_v_rapidity(-rapmax,rapmax,0.25);
  
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
  


  // declare root Ttree and add branches ----------------------------------
  TFile * tF = new TFile("example02.root", "RECREATE");
  TTree * tT = new TTree("Events", "");

  const int MaxNJets = 20; // max number of jets stored in Ttree

  int   fTEventNumber   = 0;
  int   fTNPV           = 0;
  int   fTNJets         = 0;
  int   fTNJetsTruth    = 0;
  float fTJetPt      [MaxNJets];
  float fTJetEta     [MaxNJets];
  float fTJetPhi     [MaxNJets];
  float fTJetM       [MaxNJets];
  float fTJetPtTruth [MaxNJets];

  // initialization 
  for(int ij=0; ij<MaxNJets; ++ij){
      fTJetPt       [ij] = -999;
      fTJetEta      [ij] = -999;
      fTJetPhi      [ij] = -999;
      fTJetM        [ij] = -999;

      fTJetPtTruth  [ij] = -999;

  }
  tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
  tT->Branch("NPV",                       &fTNPV,                    "NPV/I");
  tT->Branch("NJets",                     &fTNJets,                  "NJets/I");
  tT->Branch("JetPt",                     &fTJetPt,                  "JetPt[NJets]/F");
  tT->Branch("JetEta",                    &fTJetEta,                 "JetEta[NJets]/F");
  tT->Branch("JetPhi",                    &fTJetPhi,                 "JetPhi[NJets]/F");
  tT->Branch("JetM",                      &fTJetM,                   "JetM[NJets]/F");
  tT->Branch("NJetsTruth",                &fTNJetsTruth,             "NJetsTruth/I");
  tT->Branch("JetPtTruth",                &fTJetPtTruth,             "JetPtTruth[NJets]/F");
  // done with tTree 
  
  // loop over events
  int iev = 0;
  while ( mixer.next_event() && iev < nev ) {
     // increment event number    
     iev++;

     cout << iev << " " << nev << endl;

     // 
     fTEventNumber = iev;
     fTNPV         = mixer.npu();
     
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
          
     // write out the selected jets in full event
     vector<PseudoJet> full_jets = sel_jets(sorted_by_pt(cs_full.inclusive_jets()));
     if ( iev <= maxprintout ) { cerr << "Full event" << endl; }
     for (unsigned int i=0; i < full_jets.size(); i++) {
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": "  << full_jets[i] << endl; }
        // fill ttree 
        if (fTNJets == MaxNJets) continue;
        fTJetPt [fTNJets] =  full_jets[i].pt();
        fTJetEta[fTNJets] =  full_jets[i].eta();
        fTJetPhi[fTNJets] =  full_jets[i].phi();
        fTJetM  [fTNJets] =  full_jets[i].m();
        fTNJets++;
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
        if (fTNJetsTruth == MaxNJets) continue;
        if( fTNJetsTruth == MaxNJets) continue;
        fTJetPtTruth [fTNJetsTruth] = hard_jets[i].pt();
        fTNJetsTruth++;
     }
     
     // calculate quality measures, offset and dispersion and correlation coefficient,
     // for the jet transverse momentum.
     // Also fill an histogram of offset v. rapidity, to appreciate the effect
     // of rho rescaling
     if ( match(subtracted_jets, hard_jets, maxdeltaR) ) {
        for (unsigned int i=0; i < subtracted_jets.size(); i++) {
           double deltapt = subtracted_jets[i].pt() - hard_jets[i].pt();
           offset.add_entry(deltapt);
	   subhardcorr.add_entry(subtracted_jets[i].pt(),hard_jets[i].pt());
           offset_v_rapidity.add_entry(hard_jets[i].rap(),deltapt);
        }
     }	
     
     // fill tTree
     tT->Fill();

  }  // end loop over events

  // output quality measures as a function of <npu>
  cout << "# npu  <DeltaO>    sigma_DeltaO    corr.coeff." << endl;     
  cout << setprecision(9) << npu.average()    << "    " 
                          << offset.average() << "    " 
			  << offset.sd()      << "    " 
			  << subhardcorr.r() << endl;

  // output histograms
  cout << "\n\n# offset_v_rapidity" << endl;
  cout << "# binlo binmid binhi avg std err avgsquares" << endl;
  output_noNaN(offset_v_rapidity);

  // write ttree
  tT->Write();
  tF->Close();

  // free allocated memory
  delete rescaling;
  delete gmbge;
}





/// very dumb matching routine, expected to work with at most two jets.
/// It returns true if it manages to return at least one valid
/// matching between subtracted and hard, using a deltaR <= maxdeltaR rule.
/// If needed, it performs approprate swaps.
/// If only one matching is found the sub vector is resized to 1.
/// If no matching is found it returns false.
///
bool match(vector<PseudoJet> & sub, vector<PseudoJet> & hard, double maxdeltaR) {
   // check size of vectors
   if(sub.size() > 2 || hard.size() > 2 ) {
      cerr << "\nError: match() only works with vectors of size <= 2 for the time being.\nAborting." << endl;
      assert(false);  
   }   
   if ( sub.size() == 2  && hard.size() == 2 ) {
       if ( sub[0].delta_R(hard[0]) <= maxdeltaR ) {
          if ( sub[1].delta_R(hard[1]) > maxdeltaR ) { sub.resize(1); }
	  return true;
       }
       if ( sub[0].delta_R(hard[1]) <= maxdeltaR ) {
          swap(hard[0],hard[1]);
          if ( sub[1].delta_R(hard[1]) > maxdeltaR ) { sub.resize(1); }
	  return true;
       }
       if ( sub[1].delta_R(hard[0]) <= maxdeltaR ) {
          swap(sub[0],sub[1]);
	  sub.resize(1);
	  return true;
       }      
       if ( sub[1].delta_R(hard[1]) <= maxdeltaR ) { 
          swap(hard[0],hard[1]);
          swap(sub[0],sub[1]);
	  sub.resize(1);       
          return true;
       }
   }
   else if (sub.size() == 1 && hard.size() == 2 ) {
      if ( sub[0].delta_R(hard[0]) > maxdeltaR  &&  sub[0].delta_R(hard[1]) <= maxdeltaR ) {
	   swap(hard[0],hard[1]);
         return true;
      }
   }      
   else if (sub.size() == 2 && hard.size() == 1 ) {
      if ( sub[0].delta_R(hard[0]) > maxdeltaR  &&  sub[1].delta_R(hard[0]) <= maxdeltaR ) {
	   sub[0] = sub[1];
      }
      sub.resize(1);
      return true;
   }
   else if (sub.size() == 1 && hard.size() == 1 ) {
      if ( sub[0].delta_R(hard[0]) <= maxdeltaR ) { return true; }
   }
   //cout << "matching false" << endl;
   return false; 
}
