#include "puppiContainer.hh"
#include "PU14.hh"

#include "Math/ProbFunc.h"
#include "TH2F.h"
#include "TMath.h"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

puppiContainer::puppiContainer(std::vector<PseudoJet> hard_event, std::vector<PseudoJet> pileup_event, bool iDiscretize){
    _allParticles.resize(0);
    _pfParticles.resize(0);
    _genParticles.resize(0);
    _pfchsParticles.resize(0);

    _pfParticles_PU14.resize(0);
    _pfchsParticles_PU14.resize(0);
    
    _chargedLV.resize(0);
    _chargedPU.resize(0);
    _neutrals.resize(0);
    _isPU.resize(0);
    fRandom = new TRandom(0xDEADBEEF);
    // loop over hard_event
    for (unsigned int i = 0; i < hard_event.size(); i++){
        
        if (fabs(hard_event[i].eta()) > 5 || hard_event[i].pt() < 0.1) continue;
        
        PseudoJet curParticle = hard_event[i];
        curParticle.set_user_index(-99);
        int charge_tmp = hard_event[i].user_info<PU14>().charge();


        if (charge_tmp != 0) curParticle.set_user_index(2);
        else curParticle.set_user_index(0);

	if (charge_tmp != 0) hard_event[i].set_user_index(2);
	else hard_event[i].set_user_index(0);

	if (charge_tmp != 0) _chargedLV.push_back(curParticle);
        else _neutrals.push_back(curParticle);


	_genParticles.push_back(curParticle);
	_allParticles       .push_back(curParticle);
	if(charge_tmp != 0) { 
	  _isPU.push_back( 0 );
	  _pfParticles        .push_back(curParticle);
	  _pfchsParticles     .push_back(curParticle);
	  _pfParticles_PU14   .push_back(hard_event[i]);
	  _pfchsParticles_PU14.push_back(hard_event[i]);
	  continue;
	}
	if (charge_tmp == 0 && iDiscretize)  continue;
        _pfParticles.push_back        (curParticle);
	_pfchsParticles.push_back     (curParticle);
	_pfParticles_PU14.push_back   (hard_event[i]);
	_pfchsParticles_PU14.push_back(hard_event[i]);
	_isPU.push_back( 0 );
    }
    // loop over pileup_event
    for (unsigned int i = 0; i < pileup_event.size(); i++){

        if (fabs(pileup_event[i].eta()) > 5 || pileup_event[i].pt() < 0.1) continue;
        
        PseudoJet curParticle = pileup_event[i];
        curParticle.set_user_index(-99);
        int charge_tmp = pileup_event[i].user_info<PU14>().charge();
        
        if (charge_tmp != 0) curParticle.set_user_index(3);
        else curParticle.set_user_index(1);

	if (charge_tmp != 0) pileup_event[i].set_user_index(3);
	else pileup_event[i].set_user_index(1);

        if (charge_tmp != 0) _chargedPU.push_back(curParticle);
        else _neutrals.push_back(curParticle);
	_allParticles       .push_back(curParticle);
	if(charge_tmp  != 0) {
	  _pfParticles        .push_back(curParticle);
	  _pfParticles_PU14   .push_back(pileup_event[i]);
	  continue;
	}
        if (charge_tmp == 0 && iDiscretize)  continue;
	_pfParticles        .push_back(curParticle);
	_pfchsParticles     .push_back(curParticle);
	_pfParticles_PU14   .push_back(pileup_event[i]);
	_pfchsParticles_PU14.push_back(pileup_event[i]);
	_isPU.push_back( 1 );
    }
    if(iDiscretize) discretize(_pfParticles        ,_pfParticles_PU14   ,_allParticles,-1);
    if(iDiscretize) discretize(_pfchsParticles     ,_pfchsParticles_PU14,_allParticles,-1);
}
puppiContainer::~puppiContainer(){}

std::vector<fastjet::PseudoJet> puppiContainer::puppiFetch(int iPU, double iQuant){
    
    std::vector<PseudoJet> particles;

    //Run through all compute mean and RMS
    _vals.resize(0);
    
    double R0 = 0.3;
    //double R1 = 0.3;
    
    // the chi2 2dof version
    getRMSAvg(13,_pfParticles,_chargedLV,_isPU,iQuant,0.5,R0);
    double lMed0=fMed;
    double lRMS0=fRMS;
    int lNEvents  = _vals.size();

    //cout << "pass2" << endl;
    //getRMSAvg(13,_pfParticles,_pfParticles,_isPU,0.5,0.2,R1);
    //double lMed1=fMed;
    //double lRMS1=fRMS;
    
    float wptCutC = 0.5;
    float wptCutF = 1.0;
    
    //a functional form, hard-coded for now
    wptCutC = 0.66667e-2*( (float) iPU ) + 0.6;
    wptCutF = 1.05e-2*( (float) iPU ) + 0.7;

    for(int i0 = 0; i0 < lNEvents; i0++) {
        double pWeight = 1;
        
        pWeight *= compute(0,_vals[i0],lMed0,lRMS0);
        
        if(_pfParticles[i0].user_index() == 2 ) pWeight = 1;
        if(_pfParticles[i0].user_index() == 3 ) pWeight = 0;
        if(_pfParticles[i0].user_index()  < 2 && pWeight*_pfParticles[i0].pt() < wptCutC) continue;
        if(pWeight < 0.1) continue;
        
        PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
        //curjet.set_user_index(_pfParticles[i0].user_index());
	int pdgid   = _pfParticles_PU14[i0].user_info<PU14>().pdg_id();
	int barcode = particles.size();
	int vtx     = _pfParticles_PU14[i0].user_info<PU14>().vertex_number();
	curjet.set_user_info(new PU14(pdgid,barcode,vtx));
        particles.push_back(curjet);
    }
    return particles;
}


///-----------------------------------

void puppiContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,double iQuant,double iPtRMS, double R0) {
    
    ////std::vector<double> lValsPV;
    std::vector<double> lValsPU;
    std::vector<double> lValsPUHEta;
    
    std::vector<double> curVals;
    for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) {
        double pVal = goodVar(iConstits[i0],iParticles,iOpt,R0);
        curVals.push_back(pVal);
        if(iConstits[i0].pt() < iPtRMS) continue;
        if(pVal == 0) continue;
        if(iConstits[i0].user_index()  == 3) lValsPU.push_back(pVal);
        if( fabs(iConstits[i0].eta()) > 2.5   ) {lValsPUHEta.push_back(pVal); }
    }
    
    std::sort (lValsPU.begin(),lValsPU.end());
    std::sort (lValsPUHEta.begin(),lValsPUHEta.end());
    
    double lMedPU = 0; if(lValsPU.size() > 0) lMedPU = lValsPU[int(lValsPU.size()*iQuant+0.5)];
    double lMedPUHEta = 0; if(lValsPUHEta.size() > 0) lMedPUHEta = lValsPUHEta[int(lValsPUHEta.size()*iQuant+0.5)];
    
    // now add back the 0 particles to the vector to compute the RMS
    for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) {
        _vals.push_back(curVals[i0]);
    }
    
    // now compute RMS
    double lRMSPU = 0;
    int lRMSPUctr = 0;
    for(unsigned int i0 = 0 ;i0 < lValsPU.size(); i0++) {
        if ((lValsPU[i0]-lMedPU) > 0) continue;
        lRMSPU += (lValsPU[i0]-lMedPU)*(lValsPU[i0]-lMedPU);
        lRMSPUctr++;
    }
    
    double lRMSPUHEta = 0;
    int lRMSPUHEtactr = 0;
    for(unsigned int i0 = 0 ;i0 < lValsPUHEta.size(); i0++) {
        //if ((lValsPUHEta[i0]-lMedPUHEta) > 0) continue;
        lRMSPUHEta += (lValsPUHEta[i0]-lMedPUHEta)*(lValsPUHEta[i0]-lMedPUHEta);
        lRMSPUHEtactr++;
    }
    
    ////if(lValsPV.size() > 0)  lRMSPV/=lValsPV.size();
    if(lValsPU.size() > 0)  lRMSPU/=lRMSPUctr;
    if(lValsPUHEta.size() > 0)  lRMSPUHEta/=lRMSPUHEtactr;
    
    fMed = lMedPU;
    fRMS = sqrt(lRMSPU);
    fMedHEta = lMedPUHEta;
    fRMSHEta = sqrt(lRMSPUHEta);
    
}

double puppiContainer::compute(int iOpt,double iVal,double iMed,double iRMS) {
    if(iOpt == 1 && iVal < iMed) return 0;
    if(iOpt == 1 && iVal > iMed) return 1;
    double lVal = (iVal-iMed)/iRMS;

    return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal),1.);
}

double puppiContainer::goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt, double R0) {
    double Rsub = R0; //For tracking using 0.2
    //double RLrg = 1.0;
    double lPup = 0;
    lPup = var_within_R(iOpt,iParts,iPart,Rsub);
    if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
    return lPup;
}

double puppiContainer::var_within_R(int iId, const vector<PseudoJet> & particles, const PseudoJet& centre, double R){
    fastjet::Selector sel = fastjet::SelectorCircle(R);
    sel.set_reference(centre);
    vector<PseudoJet> near_particles = sel(particles);
    double var = 0;
    double lSumPt = 0;
    if(iId == 5 || iId == 6) for(unsigned int i=0; i<near_particles.size(); i++) lSumPt += near_particles[i].pt();
    for(unsigned int i=0; i<near_particles.size(); i++){
        double pDEta = near_particles[i].eta()-centre.eta();
        double pDPhi = fabs(near_particles[i].phi()-centre.phi());
        if(pDPhi > 2.*3.14159265-pDPhi) pDPhi =  2.*3.14159265-pDPhi;
        double pDR = sqrt(pDEta*pDEta+pDPhi*pDPhi);
        if(pDR  < 0.001) continue;
        if(pDR  <  0.05) pDR = 0.05;
        if(pDR == 0) continue;
        if(iId == 0) var += 1./pDR/pDR;
        if(iId == 1) var += log(1./pDR)*log(1./pDR);
        if(iId == 2) var += log(1./pDR/pDR)*log(1./pDR/pDR);
        if(iId == 3) var += log(near_particles[i].pt()*near_particles[i].pt()/pDR/pDR)*log(near_particles[i].pt()*near_particles[i].pt()/pDR/pDR);
        if(iId == 4) var += 1./pDR/pDR;;
        if(iId == 5) var += log(near_particles[i].pt()/pDR/lSumPt);
        if(iId == 6) var += near_particles[i].pt()/pDR;
        if(iId == 7) var += log(near_particles[i].pt()/pDR);
        if(iId == 8) var += log(near_particles[i].pt()/pDR)*log(near_particles[i].pt()/pDR);
        if(iId == 9) var += (near_particles[i].pt()/pDR)*(near_particles[i].pt()/pDR);
        if(iId ==10) var += near_particles[i].pt();
        if(iId ==11) var += 1./pDR/pDR;
        if(iId ==12) var += 1./pDR;
        if(iId ==13) var += near_particles[i].pt()/pDR;
        if(iId ==14) var += ((near_particles[i].pt()+centre.pt())*(near_particles[i].pt()+centre.pt())/centre.pt()/near_particles[i].pt()-1)/pDR;
    }
    if((iId == 9 || iId == 4 || iId == 12 || iId == 13 || iId == 14) && var != 0) var = log(var);
    return var;
}

double puppiContainer::pt_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R){
    
    fastjet::Selector sel = fastjet::SelectorCircle(R);
    sel.set_reference(centre);
    vector<PseudoJet> near_particles = sel(particles);
    double answer = 0.0;
    //std::cout << "near particles (pt) = " << near_particles.size() << std::endl;
    
    for(unsigned int i=0; i<near_particles.size(); i++){
        answer += near_particles[i].pt();
    }
    return(answer);
}
std::vector<fastjet::PseudoJet> puppiContainer::pfchsFetch(double iPtCut){
  if(iPtCut < 0) return _pfchsParticles_PU14;
  std::vector<PseudoJet> lParts;
  for(unsigned int i0 = 0; i0 < _pfchsParticles_PU14.size(); i0++) { 
    int charge_tmp = _pfchsParticles_PU14[i0].user_info<PU14>().charge();
    if(charge_tmp && _pfchsParticles_PU14[i0].pt() < iPtCut) continue;
    lParts.push_back(_pfchsParticles_PU14[i0]);
  }
  return lParts;
}

void puppiContainer::discretize(std::vector<fastjet::PseudoJet> &discreteParticles,std::vector<fastjet::PseudoJet> &discreteParticlesPU14,std::vector<fastjet::PseudoJet> &iParticles,double iPtCut) {
  bool lGrid = false;
  //Discretize in 0.1x0.1
  TH2F* h2d      = new TH2F( "h2d"     ,";#eta;#phi;e (GeV)",100, -5,5,63,0,2*TMath::Pi() ); //100 63
  TH2F* h2dPU    = new TH2F( "h2dPU"   ,";#eta;#phi;e (GeV)",100, -5,5,63,0,2*TMath::Pi() ); //100 63
  TH2F* h2dCh    = new TH2F( "h2dCh"   ,";#eta;#phi;e (GeV)",100, -5,5,63,0,2*TMath::Pi() ); //100 63
  std::vector<PseudoJet> lParts; 
  std::vector<PseudoJet> lPartsPU; 
  std::vector<PseudoJet> lPartsCH; 
  std::vector<PseudoJet> lPartsPUCH; 
  for(int i0   = 0; i0 < h2d->GetNbinsX()+1; i0++) { 
    for(int i1 = 0; i1 < h2d->GetNbinsY()+1; i1++) { 
      PseudoJet pPart(0,0,0,0);
      PseudoJet pPartPU(0,0,0,0);
      PseudoJet pPartCH(0,0,0,0);
      PseudoJet pPartPUCH(0,0,0,0);
      lParts.push_back    (pPart);
      lPartsPU.push_back  (pPartPU);
      lPartsCH.push_back  (pPartCH);
      lPartsPUCH.push_back(pPartPUCH);
    }
  }
  for (unsigned int i0 = 0; i0 < iParticles.size(); i0++){
    if (fabs(iParticles[i0].eta()) > 5  ) continue;        
    int curBin = h2d->FindBin(iParticles[i0].eta(),iParticles[i0].phi());
    if (lGrid  && iParticles[i0].user_index() == 0) h2d     ->SetBinContent( curBin, h2d     ->GetBinContent(curBin) + iParticles[i0].e() );        
    if (lGrid  && iParticles[i0].user_index() == 1) h2dPU   ->SetBinContent( curBin, h2dPU   ->GetBinContent(curBin) + iParticles[i0].e() );        
    if (!lGrid && iParticles[i0].user_index() == 0) lParts    [curBin] += iParticles[i0];
    if (!lGrid && iParticles[i0].user_index() == 1) lPartsPU  [curBin] += iParticles[i0];
    if (!lGrid && iParticles[i0].user_index() == 2) lPartsCH  [curBin] += iParticles[i0];
    if (!lGrid && iParticles[i0].user_index() == 3) lPartsPUCH[curBin] += iParticles[i0];
  }
  //For now sumit up
  if(!lGrid) { 
    for(unsigned int i0   = 0; i0 < lParts.size(); i0++) { 
      PseudoJet pPart,pPartPU14;
      if(lParts  [i0].pt() > 0) pPart       += lParts  [i0]; 
      if(lPartsPU[i0].pt() > 0) pPart       += lPartsPU[i0];
      double pEt    = lParts  [i0].pt() + lPartsCH  [i0].pt();
      pEt          += lPartsPU[i0].pt() + lPartsPUCH[i0].pt();
      //pEt           = smear(pEt);
      pEt          -= lPartsCH  [i0].pt();
      pEt          -= lPartsPUCH[i0].pt();
      if(pEt < TMath::Min(2.,0.2*(lPartsCH[i0].pt()))) continue; //PF Like cut
      if(pEt < 0.1) continue;
      pPart    .reset_PtYPhiM(pEt,pPart.eta(),pPart.phi(),pPart.m());
      pPartPU14.reset_PtYPhiM(pEt,pPart.eta(),pPart.phi(),pPart.m());
      int lVtx = (lParts[i0].pt() < lPartsPU[i0].pt() );
      pPart    .set_user_index(lVtx);
      pPartPU14.set_user_info(new PU14(22,i0+5000,lVtx));
      if(pPart.pt() < 0.1   ) continue; 
      if(pPart.pt() < iPtCut) continue; 
      discreteParticles    .push_back(pPart);
      discreteParticlesPU14.push_back(pPartPU14);
    }
  }
  if(lGrid) { 
    for(int i0   = 0; i0 < h2d->GetNbinsX()+1; i0++) { 
      for(int i1 = 0; i1 < h2d->GetNbinsY()+1; i1++) { 
	double pEta   = h2d->GetXaxis()->GetBinCenter(i0);
	double pPhi   = h2d->GetYaxis()->GetBinCenter(i1);
	double pPt    = fabs(h2d  ->GetBinContent(i0,i1)*sin(atan(exp(-pEta))*2.));
	double pPtPU  = fabs(h2dPU->GetBinContent(i0,i1)*sin(atan(exp(-pEta))*2.));
	double pPtCh  = fabs(h2dCh->GetBinContent(i0,i1)*sin(atan(exp(-pEta))*2.));
	double pEt    = pPt + pPtPU + pPtCh;
	//pEt = smear(pEt);
	pEt -= pPtCh;
	PseudoJet pPart,pPartPU14;
	int barcode = discreteParticles.size();
	pPart    .reset_PtYPhiM(pPt+pPtPU,pEta,pPhi);
	pPartPU14.reset_PtYPhiM(pPt+pPtPU,pEta,pPhi);
	int lVtx = (pPt < pPtPU);
	pPart     .set_user_index(lVtx);
	pPartPU14.set_user_info(new PU14(22,barcode,lVtx));
	if(pPart.pt() > 0.2) discreteParticles    .push_back(pPart);
	if(pPart.pt() > 0.2) discreteParticlesPU14.push_back(pPartPU14);
      }
    }
  }
  delete h2d;
  delete h2dPU;
  delete h2dCh;
}
double puppiContainer::smear(double iEt) { 
  double lScale = 1.0;
  if(iEt == 0) return iEt;
  //CMS resolutions because ATLAS sucks
  double lERes = 0.028/sqrt(iEt) + 0.12/iEt + 0.003;
  double lHRes = 1.   /sqrt(iEt) + 0.3/iEt;
  //Take some reasonble E/H ratio
  lERes *= iEt*0.4*lScale;
  lHRes *= iEt*0.6*lScale;
  return (fRandom->Gaus(iEt*0.4,lERes)+  fRandom->Gaus(iEt*0.6,lHRes));
}
