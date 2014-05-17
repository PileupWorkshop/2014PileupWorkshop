#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/SafeSubtractor.hh"
#include "fastjet/contrib/SoftKiller.hh"
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
void getConstitsForCleansing(vector<PseudoJet> inputs, vector<PseudoJet> &oNeutrals, vector<PseudoJet> &oChargedLV, vector<PseudoJet> &oChargedPU){
    for (unsigned int i = 0; i < inputs.size(); i++){
        if (inputs[i].user_index() <= 1) oNeutrals.push_back(inputs[i]);
        if (inputs[i].user_index() == 2) oChargedLV.push_back(inputs[i]);
        if (inputs[i].user_index() == 3) oChargedPU.push_back(inputs[i]);
    }
}
void setJet(PseudoJet &iJet,JetInfo &iJetI,std::vector<PseudoJet> &iParticles, GridMedianBackgroundEstimator* iGMBE,Subtractor &iSub) { 
  vector<PseudoJet> neutrals,chargedLV,chargedPU;
  getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);
  JetDefinition subjet_def(kt_algorithm,0.2);
  JetCleanser gsn_cleanser(subjet_def,JetCleanser::gaussian_cleansing,JetCleanser::input_nc_separate);
  gsn_cleanser.SetGaussianParameters(0.617,0.62,0.15,0.22);
  PseudoJet     lClean = gsn_cleanser(neutrals,chargedLV,chargedPU);
  
  // define safeAreaSub (PF)
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(4.0)));
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
  Selector rho_range =  SelectorAbsRapMax(4.0);
  ClusterSequenceArea clust_seq_rho(iParticles, jet_def_for_rho, area_def);
  // the two background estimators
  JetMedianBackgroundEstimator bge_rho(rho_range, clust_seq_rho);
  JetMedianBackgroundEstimator bge_rhom(rho_range, clust_seq_rho);
  BackgroundJetPtMDensity m_density;
  bge_rhom.set_jet_density_class(&m_density);
  // declare an area-median subtractor from this
  contrib::SafeAreaSubtractor area_subtractor(&bge_rho, &bge_rhom);

  //iGMBE->set_particles(iParticles);
  PseudoJet lCorr =  area_subtractor(iJet);

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
//vector<PseudoJet> threeHardest(vector<PseudoJet> &iParts, JetDefinition &iJetDef, Selector &iSelector,std::vector<ClusterSequence> &iCSs) { 
  // cluster full event (hard + pileup)
//  vector<PseudoJet> threehardest = iSelector(sorted_by_pt(cs.inclusive_jets()));
//  iCSs.push_back(cs);
//  return threehardest;
//}
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
  iJet.pt      = -1;
  iJet.ptraw   = -1;
  iJet.ptclean = -1;
  iJet.pttrim  = -1;
  iJet.eta     = -1;
  iJet.phi     = -1;
  iJet.m       = -1;
  iJet.mraw    = -1;
  iJet.mclean  = -1;
  iJet.mtrim   = -1;
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
  int lIndex = 0; lTree->Branch("index",&lIndex,"lIndex/F");
  JetInfo JGen;     setupTree(lTree,JGen    ,"Gen"  );
  JetInfo JPF;      setupTree(lTree,JPF     ,"PF"   );
  JetInfo JPup;     setupTree(lTree,JPup    ,"Puppi");
  JetInfo JCHS;     setupTree(lTree,JCHS    ,"CHS"  );
  JetInfo JCHS2GeV; setupTree(lTree,JCHS2GeV,"CHS2GeV");
  JetInfo JSoft;    setupTree(lTree,JSoft   ,"SK"  );
  JetInfo JSoftCHS; setupTree(lTree,JSoftCHS,"SKCHS");
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
    vector<PseudoJet> hard_event;
    vector<PseudoJet> pf_event; 
    vector<PseudoJet> pileup_event;
    vector<PseudoJet> puppi_event;
    vector<PseudoJet> chs_event;
    vector<PseudoJet> chs_event2GeV;
    vector<PseudoJet> neutral_event;
    vector<PseudoJet> chargepu_event;
    vector<PseudoJet> charge_event;     
    SelectorIsHard().sift(full_event, hard_event, pileup_event); 
    //////////////////////////////////////////////////////
    puppiContainer curEvent(hard_event, pileup_event);
    puppi_event     = curEvent.puppiFetch(80);
    pf_event        = curEvent.pfFetch();
    chs_event       = curEvent.pfchsFetch(-1);
    chs_event2GeV   = curEvent.pfchsFetch( 2.);
    //////////////////////////////////////////////////////
    SoftKiller soft_killer   (0.4,0.4);
    SoftKiller soft_killerCHS(4.0,0.5, !SelectorCharged());
    vector<PseudoJet> soft_event    = soft_killer(full_event);
    vector<PseudoJet> softCHS_event = soft_killerCHS(chs_event);

    AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(4.0)));
    ClusterSequenceArea pGen    (hard_event   ,jet_def,area_def);
    ClusterSequenceArea pPup    (puppi_event  ,jet_def,area_def);
    ClusterSequenceArea pPF     (pf_event     ,jet_def,area_def);
    ClusterSequenceArea pCHS    (chs_event    ,jet_def,area_def);
    ClusterSequenceArea pCHS2GeV(chs_event2GeV,jet_def,area_def);
    ClusterSequenceArea pSoft   (soft_event   ,jet_def,area_def);
    ClusterSequenceArea pSoftCHS(softCHS_event,jet_def,area_def);
    vector<PseudoJet> genJets     = selector(sorted_by_pt(pGen    .inclusive_jets())); 
    vector<PseudoJet> puppiJets   = selector(sorted_by_pt(pPup    .inclusive_jets())); 
    vector<PseudoJet> pfJets      = selector(sorted_by_pt(pPF     .inclusive_jets())); 
    vector<PseudoJet> chsJets     = selector(sorted_by_pt(pCHS    .inclusive_jets())); 
    vector<PseudoJet> chs2GeVJets = selector(sorted_by_pt(pCHS2GeV.inclusive_jets())); 
    vector<PseudoJet> softJets    = selector(sorted_by_pt(pSoft   .inclusive_jets())); 
    vector<PseudoJet> softCHSJets = selector(sorted_by_pt(pSoftCHS.inclusive_jets())); 
    for(unsigned int i0 = 0; i0 < genJets.size(); i0++) {
      lIndex = i0;
      PseudoJet puppiJet   = match(genJets[i0],puppiJets);
      PseudoJet pfJet      = match(genJets[i0],pfJets   );
      PseudoJet chsJet     = match(genJets[i0],chsJets  );
      PseudoJet chs2GeVJet = match(genJets[i0],chs2GeVJets);
      PseudoJet softJet    = match(genJets[i0],softJets);
      PseudoJet softCHSJet = match(genJets[i0],softCHSJets);
      setJet(genJets[i0],JGen    ,hard_event   ,gmbge,sub);
      if(pfJet.pt()      != 0) setJet(pfJet ,     JPF     ,pf_event     ,gmbge,sub);
      if(chsJet.pt()     != 0) setJet(chsJet,     JCHS    ,chs_event    ,gmbge,sub);
      if(chs2GeVJet.pt() != 0) setJet(chs2GeVJet, JCHS2GeV,chs_event2GeV,gmbge,sub);
      if(puppiJet.pt()   != 0) setJet(puppiJet  , JPup    ,puppi_event  ,gmbge,sub);
      if(softJet.pt()    != 0) setJet(softJet   , JSoft   ,puppi_event  ,gmbge,sub);
      if(softCHSJet.pt() != 0) setJet(softCHSJet, JSoftCHS,puppi_event  ,gmbge,sub);
      lTree->Fill();
    }
  }
  lFile->cd();
  lTree->Write();
}  
