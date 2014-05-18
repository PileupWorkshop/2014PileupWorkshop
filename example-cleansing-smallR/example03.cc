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
#include "fastjet/tools/Filter.hh"

#include <iomanip>      // std::setprecision
#include <cassert>
#include <cmath>

#include "fastjet/contrib/JetCleanser.hh"
#include "TTree.h"
#include "TFile.h"

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
  string cleansing_config = cmdline.value<string>("-cleansing","jvf");

  // cleansing mode
  cout <<" # cleansing config: " << cleansing_config << endl;
  contrib::JetCleanser::cleansing_mode mode;
  if      (cleansing_config == "jvf"   ) { mode = contrib::JetCleanser::jvf_cleansing;}
  else if (cleansing_config == "linear") { mode = contrib::JetCleanser::linear_cleansing;}
  else if (cleansing_config == "gaus"  ) { mode = contrib::JetCleanser::gaussian_cleansing;}

  // matching condition
  double maxdeltaR = cmdline.value<double>("-maxdeltaR",0.3);
  Matching matching(maxdeltaR);
  cout << "# matching: " << matching.description() << endl;

  // matching reco to truth jets
  Matching matching_to_truth(maxdeltaR);


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

  AverageAndError npu, njets, offset, matching_efficiency, pass_fraction, hard_pt;
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
  
  // selectors for cleansing
  Selector SelectorHStracks =    SelectorVertexNumber(0) * SelectorIsCharged()  ;
  Selector SelectorPUtracks = (!SelectorVertexNumber(0)) * SelectorIsCharged()  ;


  contrib::JetCleanser cleanser(JetDefinition(kt_algorithm, 0.2), mode, contrib::JetCleanser::input_nc_together);
  if       ( mode == contrib::JetCleanser::jvf_cleansing){
    cleanser.SetTrimming(0.);
  }else if ( mode == contrib::JetCleanser::linear_cleansing){
    cleanser.SetLinearParameters(0.65);
  }else if ( mode == contrib::JetCleanser::gaussian_cleansing){
    cleanser.SetGaussianParameters(0.67,0.62,0.20,0.25);
  }


  // trimming
  fastjet::Filter trimmer (fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.0));
  trimmer.set_subtractor(&sub);

  // declare root Ttree and add branches ----------------------------------

//  string treename = cmdline.value<double>("-",0.3);

  ostringstream filename;
  filename << cleansing_config << "_npv" << mixer.npu() << ".root";

  TFile * tF = new TFile(filename.str().c_str(), "RECREATE");
  TTree * tT = new TTree("Events", "");

  const int MaxNJets = 20; // max number of jets stored in Ttree

  int   fTEventNumber   = 0;
  int   fTNPV           = 0;
  int   fTNJets         = 0;
  int   fTNJetsJVFCl    = 0;
  int   fTNJetsTruth    = 0;
  int   fTNJetsTrimmed    = 0;
  float fTJetPt      [MaxNJets];
  float fTJetPtOffset[MaxNJets];
  float fTJetEta     [MaxNJets];
  float fTJetPhi     [MaxNJets];
  float fTJetM       [MaxNJets];
  float fTJetMOffset [MaxNJets];
  // cleansing
  float fTJetJVFClPt [MaxNJets];
  float fTJetJVFClPtOffset [MaxNJets];
  float fTJetJVFClEta[MaxNJets];
  float fTJetJVFClPhi[MaxNJets];
  float fTJetJVFClM  [MaxNJets];
  float fTJetJVFClMOffset  [MaxNJets];
  // truth 
  float  fTJetPtTruth [MaxNJets];
  float  fTJetPhiTruth[MaxNJets];
  float  fTJetEtaTruth[MaxNJets];
  float  fTJetMTruth  [MaxNJets];
  // trimmed
  float  fTJetTrimmedPt[MaxNJets];
  float  fTJetTrimmedPhi[MaxNJets];
  float  fTJetTrimmedEta[MaxNJets];
  float  fTJetTrimmedM[MaxNJets];
  float  fTJetTrimmedMOffset [MaxNJets];
  float fTJetTrimmedPtOffset[MaxNJets];

  tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
  tT->Branch("NPV",                       &fTNPV,                    "NPV/I");
  tT->Branch("NJets",                     &fTNJets,                  "NJets/I");
  tT->Branch("JetPt",                     &fTJetPt,                  "JetPt[NJets]/F");
  tT->Branch("JetPtOffset",               &fTJetPtOffset,            "JetPtOffset[NJets]/F");
  tT->Branch("JetEta",                    &fTJetEta,                 "JetEta[NJets]/F");
  tT->Branch("JetPhi",                    &fTJetPhi,                 "JetPhi[NJets]/F");
  tT->Branch("JetM",                      &fTJetM,                   "JetM[NJets]/F");
  tT->Branch("JetMOffset",                &fTJetMOffset,             "JetMOffset[NJets]/F");
  //
  tT->Branch("NJetsTruth",                &fTNJetsTruth,             "NJetsTruth/I");
  tT->Branch("JetPtTruth",                &fTJetPtTruth,             "JetPtTruth[NJetsTruth]/F");
  tT->Branch("JetEtaTruth",               &fTJetEtaTruth,            "JetEtaTruth[NJetsTruth]/F");
  tT->Branch("JetPhiTruth",               &fTJetPhiTruth,            "JetPhiTruth[NJetsTruth]/F");
  tT->Branch("JetMTruth",                 &fTJetMTruth,              "JetMTruth[NJetsTruth]/F");
  //
  tT->Branch("NJetsJVFCl",                &fTNJetsJVFCl,             "NJetsJVFCl/I");
  tT->Branch("JetJVFClPt",                &fTJetJVFClPt,             "JetJVFClPt[NJetsJVFCl]/F");
  tT->Branch("JetJVFClPtOffset",          &fTJetJVFClPtOffset,       "JetJVFClPtOffset[NJetsJVFCl]/F");
  tT->Branch("JetJVFClEta",               &fTJetJVFClEta,            "JetJVFClEta[NJetsJVFCl]/F");
  tT->Branch("JetJVFClPhi",               &fTJetJVFClPhi,            "JetJVFClPhi[NJetsJVFCl]/F");
  tT->Branch("JetJVFClM",                 &fTJetJVFClM,              "JetJVFClM[NJetsJVFCl]/F");
  tT->Branch("JetJVFClMOffset",           &fTJetJVFClMOffset,        "JetJVFClMOffset[NJetsJVFCl]/F");
  //
  tT->Branch("NJetsTrimmed",         &fTNJetsTrimmed,      "NJetsTrimmed/I");
  tT->Branch("JetTrimmedPt",         &fTJetTrimmedPt,      "JetTrimmedPt[NJetsTrimmed]/F");
  tT->Branch("JetTrimmedEta",        &fTJetTrimmedEta,     "JetTrimmedEta[NJetsTrimmed]/F");
  tT->Branch("JetTrimmedPhi",        &fTJetTrimmedPhi,     "JetTrimmedPhi[NJetsTrimmed]/F");
  tT->Branch("JetTrimmedM",          &fTJetTrimmedM,       "JetTrimmedM[NJetsTrimmed]/F");
  tT->Branch("JetTrimmedMOffset",          &fTJetTrimmedMOffset,       "JetTrimmedMOffset[NJetsTrimmed]/F");
  tT->Branch("JetTrimmedPtOffset",         &fTJetTrimmedPtOffset,      "JetTrimmedPtOffset[NJetsTrimmed]/F");
  //
  // done with tTree 
  
  // loop over events ----------------------------------------------------------------------
  int iev = 0;
  while ( mixer.next_event() && iev < nev ) {
     // increment event number    
     iev++;
     
     // Reset counters
     fTNJets         = 0;
     fTNJetsTruth    = 0;
     fTNJetsJVFCl    = 0;
     fTNJetsTrimmed  = 0;

     // 
     fTEventNumber = iev;
     fTNPV         = mixer.npu();



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
        if (fTNJetsTruth == MaxNJets) continue;
        fTJetPtTruth  [fTNJetsTruth]= hard_jets[i].pt();
        fTJetEtaTruth [fTNJetsTruth]= hard_jets[i].eta();
        fTJetPhiTruth [fTNJetsTruth]= hard_jets[i].phi();
        fTJetMTruth   [fTNJetsTruth]= hard_jets[i].m();
        fTNJetsTruth++;
        // record info about the average hard pt
        hard_pt += hard_jets[i].pt();
     }
     
     // if no hard jets pass the selection, continue on to the next event
     if (hard_jets.size() == 0) {
       pass_fraction += 0;
       continue;
     }
     pass_fraction += 1;
    


     // cluster full event (hard + pileup)
     ClusterSequenceArea cs_full(full_event,jet_def,area_def);


     // setup matching to truth
     matching_to_truth.set_full_jets(hard_jets);
          


     // get all jets in full event
     vector<PseudoJet> full_jets       = sorted_by_pt(cs_full.inclusive_jets());
     vector<PseudoJet> jvfcleansed_jets;
     if ( iev <= maxprintout ) { cerr << "Full event" << endl; }
     for (unsigned int i=0; i < 4U && i < full_jets.size(); i++) {
        PseudoJet JVFcleansedJet = cleanser( full_jets[i], SelectorHStracks(full_jets[i].constituents()),  SelectorPUtracks(full_jets[i].constituents()));
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": "  << full_jets[i]  << " cleansed " << JVFcleansedJet << endl; }
         // fill ttree 
        if (fTNJetsJVFCl == MaxNJets) continue;
        fTJetJVFClPt       [fTNJetsJVFCl] =  JVFcleansedJet.pt();
        fTJetJVFClEta      [fTNJetsJVFCl] =  JVFcleansedJet.eta();
        fTJetJVFClPhi      [fTNJetsJVFCl] =  JVFcleansedJet.phi();
        fTJetJVFClM        [fTNJetsJVFCl] =  JVFcleansedJet.m();

        const PseudoJet * truthmatch = matching_to_truth.match(JVFcleansedJet);
        if(truthmatch ){
            fTJetJVFClPtOffset[fTNJetsJVFCl]  = JVFcleansedJet.pt() - truthmatch->pt();
            fTJetJVFClMOffset [fTNJetsJVFCl]  = JVFcleansedJet.m()  - truthmatch->m();
        }else{
            fTJetJVFClPtOffset[fTNJetsJVFCl]  = -999;
            fTJetJVFClMOffset [fTNJetsJVFCl]  = -999;
        }
        fTNJetsJVFCl++;
        jvfcleansed_jets.push_back(JVFcleansedJet);
     }



     
     // subtract the full jets
     // 1. feed particles to background estimator
     gmbge->set_particles(full_event);
     // 2. subtract
     vector<PseudoJet> subtracted_jets = sub(full_jets); 
     // write out jets in subtracted event
     if ( iev <= maxprintout ) { cerr << "Subtracted event" << endl; }
     for (unsigned int i=0; i < 4U && i < subtracted_jets.size(); i++) {
         // fill ttree 
        if (fTNJets == MaxNJets) continue;
        fTJetPt       [fTNJets] =  subtracted_jets[i].pt();
        fTJetEta      [fTNJets] =  subtracted_jets[i].eta();
        fTJetPhi      [fTNJets] =  subtracted_jets[i].phi();
        fTJetM        [fTNJets] =  subtracted_jets[i].m();

        const PseudoJet * truthmatch = matching_to_truth.match( subtracted_jets[i]);
        if(truthmatch ){
            fTJetPtOffset[fTNJets]  =  subtracted_jets[i].pt() - truthmatch->pt();
            fTJetMOffset [fTNJets]  =  subtracted_jets[i].m()  - truthmatch->m();
        }else{
            fTJetPtOffset[fTNJets]  = -999;
            fTJetMOffset [fTNJets]  = -999;
        }
        fTNJets++;

        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": "  << subtracted_jets[i] << endl; }
     }
        
     // trimmed -----------------------
     if ( iev <= maxprintout ) { cerr << "Trimmed event" << endl; }
     for (unsigned int i=0; i < 4U && i < full_jets.size(); i++) {
         PseudoJet trimmed = trimmer(full_jets[i]);
        if ( iev <= maxprintout ) { cerr << "  jet " << i << ": "  << full_jets[i]  << "  trimmed " << trimmed << endl; }
         // fill ttree 
        if (fTNJetsTrimmed == MaxNJets) continue;
        fTJetTrimmedPt       [fTNJetsTrimmed] =  trimmed.pt();
        fTJetTrimmedEta      [fTNJetsTrimmed] =  trimmed.eta();
        fTJetTrimmedPhi      [fTNJetsTrimmed] =  trimmed.phi();
        fTJetTrimmedM        [fTNJetsTrimmed] =  trimmed.m();

        const PseudoJet * truthmatch = matching_to_truth.match( trimmed);
        if(truthmatch ){
            fTJetTrimmedPtOffset[fTNJetsTrimmed]  =  trimmed.pt() - truthmatch->pt();
            fTJetTrimmedMOffset [fTNJetsTrimmed]  =  trimmed.m()  - truthmatch->m();
        }else{
            fTJetTrimmedPtOffset[fTNJetsTrimmed]  = -999;
            fTJetTrimmedMOffset [fTNJetsTrimmed]  = -999;
        }
        fTNJetsTrimmed++;
     }
     
     // calculate quality measures, offset and dispersion and correlation coefficient,
     // for the jet transverse momentum.
     // Also fill an histogram of offset v. rapidity, to appreciate the effect
     // of rho rescaling

     // set up the set of full/subtracted jets from which to match
     matching.set_full_jets(jvfcleansed_jets);
     
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
     
     // fill tTree
     tT->Fill();
     
  }  // end loop over events


    // output quality measures as a function of <npu>
  cout << "# fraction of events passing the basic jet cuts = " << pass_fraction.average() << " +- " << pass_fraction.error() << endl;
  cout << "# npu    jet_ptmin     <DeltaO>         sigma_DeltaO   corr.coeff.     njets>20GeV           match_eff       <O_hard>           name " << endl;     
  cout << setprecision(4) 
       << setw(4) << npu.average()    << "    "
       << setw(6) << jet_ptmin << "    "
       << setw(6) << offset.average() << " +- " << setw(6) << offset.error() << "    "
       << setw(6) << offset.sd()      << " +- " << setw(6) << offset.error_on_sd() << "   "
       << setw(6) << subhardcorr.r()  << "    "
       << setw(6) << njets.average() << " +- " << setw(6) << njets.error() << "    "
       << setw(6) << matching_efficiency.average() << " +- " << setw(6) << matching_efficiency.error() << "    "
       << setw(6) << hard_pt.average() << " +- " << setw(6) << hard_pt.error()
       << "    # jet cleansing " << cleansing_config  // a label to say what the observable and tool were
       << endl;


  // output histograms
  cout << "\n\n# offset_v_rapidity" << endl;
  cout << "# binlo binmid binhi avg stddev err avgsquares" << endl;
  output_noNaN(offset_v_rapidity);
  
  // write ttree
  tT->Write();
  tF->Close();

  // free allocated memory
  delete rescaling;
  delete gmbge;
}

