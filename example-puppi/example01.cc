///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This example is meant to illustrate the following features:
//   - Use of the class EventMixer to read in signal and pileup events from
//     their respective files, and merge them
//   - Access to particle properties through PseudoJet's user_info method
//     and the PU14 class
//   - Use of SelectorIsHard to select, in the mixed event, particles
//     originating from the hard event only
//   - Use of the CmdLine class to pass optsions from the command line
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

#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "fastjet/ClusterSequence.hh"
#include "puppiContainer.hh"
#include "TTree.h"
#include "TFile.h"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {
    CmdLine cmdline(argc,argv);
    // inputs read from command line
    int nev = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
    bool verbose = cmdline.present("-verbose");
    double R = cmdline.value<double>("-R",0.4);
    
    // some definitions
    JetDefinition jet_def(antikt_algorithm,R);         // the jet definition....
    cout << jet_def.description() << endl;             // ....and the output of its description
    Selector sel_threehardest = SelectorNHardest(3);   // definition of a selector for the three hardest jets
    
    // create mixer that will construct events by mixing hard and pileup
    // events read from files given from command line using
    // -hard hard_events_file(.gz) -pileup pileup_events_file(.gz)
    EventMixer mixer(&cmdline);
    
    TFile * tF = new TFile("example01.root", "RECREATE");
    TTree * tT = new TTree("Events", "");
    const int NJets = 3;
    int   fTEventNumber   = 0;
    float fJet_PF_pt[NJets];
    float fJet_PF_ma[NJets];
    float fJet_PFCHS_pt[NJets];
    float fJet_PFCHS_ma[NJets];
    float fJet_PUPPI_pt[NJets];
    float fJet_PUPPI_ma[NJets];
    float fJet_GEN_pt[NJets];
    float fJet_GEN_ma[NJets];
    
    tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tT->Branch("Jet_PF_pt",                 &fJet_PF_pt,               "Jet_PF_pt[3]/F");
    tT->Branch("Jet_PF_ma",                 &fJet_PF_ma,               "Jet_PF_ma[3]/F");
    tT->Branch("Jet_PFCHS_pt",              &fJet_PFCHS_pt,            "Jet_PFCHS_pt[3]/F");
    tT->Branch("Jet_PFCHS_ma",              &fJet_PFCHS_ma,            "Jet_PFCHS_ma[3]/F");
    tT->Branch("Jet_PUPPI_pt",              &fJet_PUPPI_pt,            "Jet_PUPPI_pt[3]/F");
    tT->Branch("Jet_PUPPI_ma",              &fJet_PUPPI_ma,            "Jet_PUPPI_ma[3]/F");
    tT->Branch("Jet_GEN_pt",                &fJet_GEN_pt,              "Jet_GEN_pt[3]/F");
    tT->Branch("Jet_GEN_ma",                &fJet_GEN_ma,              "Jet_GEN_ma[3]/F");
    
    // loop over events
    int iev = 0;
    while ( mixer.next_event() && iev < nev ) {
        // increment event number
        iev++;
        
        // extract particles from event
        vector<PseudoJet> full_event = mixer.particles() ;
        cout << "\nEvent " << iev << endl;
        cout << "nPU = " << mixer.npu() << endl;
        
        // write out some details about the event and its particles, exploiting
        // PseudoJet::user_info() to extract information stored via the PU14 class
        PseudoJet sum;
        double net_charge = 0;
        for (unsigned i = 0; i < full_event.size(); i++) {
            const PseudoJet & p  = full_event[i];
            if ( verbose ) {
                cout << i << " "
                << p.px() << " " << p.py() << " " << p.pz() << " " << p.E() << " " << p.m() << " "
                << p.user_info<PU14>().pdg_id() << " " << p.user_info<PU14>().charge() << endl;
            }
            sum += p;
            net_charge += p.user_info<PU14>().charge();
        }
        cout << "total momentum: " << sum << endl;
        cout << "total charge: " << net_charge << endl;
        
        // cluster full event (hard + pileup)
        ClusterSequence cs_full(full_event,jet_def);
        
        // write out three hardest jets in full event
        vector<PseudoJet> threehardest = sel_threehardest(sorted_by_pt(cs_full.inclusive_jets()));
        cout << "Full event" << endl;
        cout << "  jet 0: " << threehardest[0] << endl;
        cout << "  jet 1: " << threehardest[1] << endl;
        cout << "  jet 2: " << threehardest[2] << endl;
        
        // cluster hard event only
        vector<PseudoJet> hard_event, pileup_event;
        SelectorIsHard().sift(full_event, hard_event, pileup_event); // this sifts the full event into two vectors
        // of PseudoJet, one for the hard event, one for the pileup
        ClusterSequence cs_hard(hard_event,jet_def);
        
        // write out three hardest jets in hard event
        threehardest = sel_threehardest(sorted_by_pt(cs_hard.inclusive_jets()));
        cout << "Hard event" << endl;
        cout << "  jet 0: " << threehardest[0] << endl;
        cout << "  jet 1: " << threehardest[1] << endl;
        cout << "  jet 2: " << threehardest[2] << endl;
        
        
        //////////////////////////////////////////////////////
        
        puppiContainer curEvent(hard_event, pileup_event);
        vector<PseudoJet> puppiParticles = curEvent.puppiFetch(20);
        
        ClusterSequence cs_pf(curEvent.pfFetch(),jet_def);
        ClusterSequence cs_pfchs(curEvent.pfchsFetch(-1),jet_def);
        ClusterSequence cs_puppi(puppiParticles,jet_def);
        ClusterSequence cs_gen(curEvent.genFetch(),jet_def);

        vector<PseudoJet> hard3_pf = sel_threehardest(sorted_by_pt(cs_pf.inclusive_jets()));
        vector<PseudoJet> hard3_pfchs = sel_threehardest(sorted_by_pt(cs_pfchs.inclusive_jets()));
        vector<PseudoJet> hard3_puppi = sel_threehardest(sorted_by_pt(cs_puppi.inclusive_jets()));
        vector<PseudoJet> hard3_gen = sel_threehardest(sorted_by_pt(cs_gen.inclusive_jets()));
        
        for (unsigned int i = 0; i < 3; i++){
            
            cout << hard3_pf[i].pt() << ", " << hard3_pfchs[i].pt() << ", " << hard3_puppi[i].pt() << ", " << hard3_gen[i].pt() << endl;
            
            // PF event
            fJet_PF_pt[i] = hard3_pf[i].pt();
            fJet_PF_ma[i] = hard3_pf[i].m();

            // PF CHS event
            fJet_PFCHS_pt[i] = hard3_pfchs[i].pt();
            fJet_PFCHS_ma[i] = hard3_pfchs[i].m();

            // PUPPI event
            fJet_PUPPI_pt[i] = hard3_puppi[i].pt();
            fJet_PUPPI_ma[i] = hard3_puppi[i].m();

            // GEN event
            fJet_GEN_pt[i] = hard3_gen[i].pt();
            fJet_GEN_ma[i] = hard3_gen[i].m();
            
        }
        
        tT->Fill();
        
        //////////////////////////////////////////////////////
        
    }
    
    tT->Write();
    tF->Close();

}
