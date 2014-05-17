#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "puppiContainer.hh"

using namespace std;
using namespace fastjet;
using namespace contrib;

struct JetInfo { 
  float pt;
  float ptraw;
  float ptclean;
  float pttrim;
  float eta;
  float phi;
  float m;
  float mraw;
  float mclean;
  float mtrim;
};
void setJet(PseudoJet &iJet,JetInfo &iJetI,std::vector<PseudoJet> &iParticles, GridMedianBackgroundEstimator* iGMBE,Subtractor &iSub,
	    std::vector<PseudoJet> &iCharged,std::vector<PseudoJet> &iNeutral,std::vector<PseudoJet> &iPileup
	    ) { 
  JetDefinition subjet_def(kt_algorithm,0.2);
  JetCleanser gsn_cleanser(subjet_def,JetCleanser::gaussian_cleansing,JetCleanser::input_nc_separate);
  gsn_cleanser.SetGaussianParameters(0.617,0.62,0.15,0.22);
  PseudoJet     lClean = gsn_cleanser(iNeutral,iCharged,iPileup);

  iGMBE->set_particles(iParticles);
  vector<PseudoJet> lCorrs = iSub(iParticles); 
  PseudoJet lCorr = lCorrs[0];
  
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), fastjet::SelectorPtFractionMin(0.05)));
  PseudoJet lTrim = (trimmer)(iJet);
  iJetI.pt      = lCorr .pt();  
  iJetI.ptraw   = iJet  .pt();
  iJetI.ptclean = lClean.pt();
  iJetI.pttrim  = lTrim .pt();
  iJetI.eta     = iJet.eta();
  iJetI.phi     = iJet.phi();
  iJetI.mraw    = iJet.m();  
  iJetI.m       = lCorr .m();  
  iJetI.mclean  = lClean.m();  
  iJetI.mtrim   = lTrim .m();
}
void setupTree(TTree *iTree,JetInfo &iJet,std::string iName) { 
  iTree->Branch((iName+"pt"     ).c_str(),&iJet.pt     ,(iName+"pt/F"   ).c_str());
  iTree->Branch((iName+"ptraw"  ).c_str(),&iJet.ptraw  ,(iName+"ptraw/F").c_str());
  iTree->Branch((iName+"ptclean").c_str(),&iJet.ptclean,(iName+"ptclean/F").c_str());
  iTree->Branch((iName+"pttrim" ).c_str(),&iJet.pttrim ,(iName+"pttrim/F").c_str());
  iTree->Branch((iName+"eta"    ).c_str(),&iJet.eta    ,(iName+"eta/F"   ).c_str());
  iTree->Branch((iName+"phi"    ).c_str(),&iJet.phi    ,(iName+"phi/F"   ).c_str());
  iTree->Branch((iName+"m"      ).c_str(),&iJet.m      ,(iName+"m/F"     ).c_str());
  iTree->Branch((iName+"mraw"   ).c_str(),&iJet.mraw   ,(iName+"mraw/F"  ).c_str());
  iTree->Branch((iName+"mtrim"  ).c_str(),&iJet.mtrim  ,(iName+"mtrim/F" ).c_str());
  iTree->Branch((iName+"mclean" ).c_str(),&iJet.mclean ,(iName+"mclean/F").c_str());
}
vector<PseudoJet> threeHardest(vector<PseudoJet> &iParts, JetDefinition &iJetDef, Selector &iSelector) { 
    // cluster full event (hard + pileup)
    ClusterSequence cs(iParts,iJetDef);   
    vector<PseudoJet> threehardest = iSelector(sorted_by_pt(cs.inclusive_jets()));
    return threehardest;
}
PseudoJet match(PseudoJet &iJet,vector<PseudoJet> &iJets) { 
  for(unsigned int i0 = 0; i0 < iJets.size(); i0++) { 
    double pEta = fabs(iJet.eta()-iJets[i0].eta());
    double pPhi = fabs(iJet.phi() - iJets[i0].phi());
    if(pPhi > 2.*TMath::Pi()-pPhi) pPhi =  2.*TMath::Pi()-pPhi;
    if(sqrt(pEta*pEta+pPhi*pPhi) > 0.3) continue;
    return iJets[i0];
  }
  return PseudoJet();
}
void clear(JetInfo &iJet) { 
  iJet.pt    = -1; 
  iJet.ptraw = -1; 
  iJet.eta   = -1; 
  iJet.phi   = -1; 
  iJet.m     = -1; 
  iJet.mtrim = -1; 
}
int main (int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nev = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  double R = cmdline.value<double>("-R",0.4);
  
  // some definitions
  JetDefinition jet_def(antikt_algorithm,R);         // the jet definition....
  cout << jet_def.description() << endl;             // ....and the output of its description
  Selector selector = SelectorNHardest(3);   // definition of a selector for the three hardest jets
  EventMixer mixer(&cmdline);  
  GridMedianBackgroundEstimator * gmbge = new GridMedianBackgroundEstimator(4.0,0.55);
  FunctionOfPseudoJet<double> * rescaling = new BackgroundRescalingYPolynomial(1.1685397, 0, -0.0246807, 0, 5.94119e-05);
  gmbge -> set_rescaling_class(rescaling); 
  Subtractor sub(gmbge);
 
  // loop over events
  int iev = 0;
  TFile *lFile = new TFile("Output.root","RECREATE");
  TTree *lTree = new TTree("Tree","ATLAS Sucks");
  JetInfo JGen;     setupTree(lTree,JGen    ,"Gen"  );
  JetInfo JPF;      setupTree(lTree,JPF     ,"PF"   );
  JetInfo JPup;     setupTree(lTree,JPup    ,"Puppi");
  JetInfo JCHS;     setupTree(lTree,JCHS    ,"CHS"  );
  JetInfo JCHS2GeV; setupTree(lTree,JCHS2GeV,"CHS2GeV");
  while ( mixer.next_event() && iev < nev ) {
    clear(JGen);
    clear(JPF);
    clear(JPup);
    clear(JCHS);
    clear(JCHS2GeV);
    // increment event number    
    iev++;
    cout << "\nEvent " << iev << endl;
    cout << "nPU = " << mixer.npu() << endl;
    
    // extract particles from event 
    vector<PseudoJet> full_event = mixer.particles() ;
    vector<PseudoJet> hard_event,pf_event, pileup_event,puppi_event,chs_event,chs_event2GeV,neutral_event,chargepu_event,charge_event;     
    SelectorIsHard().sift(full_event, hard_event, pileup_event); 
    //////////////////////////////////////////////////////
    puppiContainer curEvent(hard_event, pileup_event);
    puppi_event     = curEvent.puppiFetch(80);
    pf_event        = curEvent.pfFetch();
    chs_event       = curEvent.pfchsFetch(-1);
    chs_event2GeV   = curEvent.pfchsFetch( 2.);
    charge_event    = curEvent.chargedLVFetch();
    chargepu_event  = curEvent.chargedPUFetch();
    neutral_event   = curEvent.neutralFetch();
    //////////////////////////////////////////////////////
    vector<PseudoJet> genJets     = threeHardest(hard_event   ,jet_def,selector);
    vector<PseudoJet> puppiJets   = threeHardest(puppi_event  ,jet_def,selector);
    vector<PseudoJet> pfJets      = threeHardest(pf_event     ,jet_def,selector);
    vector<PseudoJet> chsJets     = threeHardest(chs_event    ,jet_def,selector);
    vector<PseudoJet> chs2GeVJets = threeHardest(chs_event2GeV,jet_def,selector);
    for(unsigned int i0 = 0; i0 < puppiJets.size(); i0++) {
      PseudoJet puppiJet   = match(genJets[i0],puppiJets);
      PseudoJet pfJet      = match(genJets[i0],pfJets   );
      PseudoJet chsJet     = match(genJets[i0],chsJets  );
      PseudoJet chs2GeVJet = match(genJets[i0],chs2GeVJets);
      setJet(genJets[i0],JGen    ,hard_event   ,gmbge,sub,charge_event,neutral_event,chargepu_event);      
      setJet(pfJet ,     JPF     ,pf_event     ,gmbge,sub,charge_event,neutral_event,chargepu_event);      
      setJet(chsJet,     JCHS    ,chs_event    ,gmbge,sub,charge_event,neutral_event,chargepu_event);
      setJet(chs2GeVJet, JCHS2GeV,chs_event2GeV,gmbge,sub,charge_event,neutral_event,chargepu_event);
      setJet(puppiJet  , JPup    ,puppi_event  ,gmbge,sub,charge_event,neutral_event,chargepu_event);
      lTree->Fill();
    }m
  }
  lFile->cd();
  lTree->Write();
}  
