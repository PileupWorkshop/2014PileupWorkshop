#include "fastjet/internal/base.hh"
#include "fastjet/PseudoJet.hh"
#include "TRandom.h"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//......................
class puppiContainer{

public:
    // default ctor
    puppiContainer( std::vector<PseudoJet> hard_event, std::vector<PseudoJet> pileup_event, bool iDiscretize=false);
    ~puppiContainer();
    
    std::vector<fastjet::PseudoJet> puppiFetch(int nPU, double iQuant=0.5);
    std::vector<fastjet::PseudoJet> genFetch(){ return _genParticles; }
    std::vector<fastjet::PseudoJet> pfFetch(){ return _pfParticles_PU14; }
    std::vector<fastjet::PseudoJet> pfchsFetch(double iPtCut);
    std::vector<fastjet::PseudoJet> chargedLVFetch(){ return _chargedLV; }
    std::vector<fastjet::PseudoJet> chargedPUFetch(){ return _chargedPU; }
    std::vector<fastjet::PseudoJet> neutralFetch  (){ return _neutrals; }
protected:

    double var_within_R(int iId, const vector<PseudoJet> & particles, const PseudoJet& centre, double R);
    double pt_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R);
    double goodVar(PseudoJet &iPart,std::vector<PseudoJet> &iParts, int iOpt, double R0);
    double compute(int iOpt,double iVal,double iMed,double iRMS);
    void   getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<int> &iIsPU,double iQuant,double iPtRMS, double R0);
    double smear(double iEt);    
    void discretize(std::vector<fastjet::PseudoJet> &discreteParticles,std::vector<fastjet::PseudoJet> &discreteParticlesPU14,std::vector<fastjet::PseudoJet> &iParticles,double iPtCut);
    
    std::vector<PseudoJet> _allParticles;
    std::vector<PseudoJet> _pfParticles;
    std::vector<PseudoJet> _pfchsParticles;
    std::vector<PseudoJet> _genParticles;

    std::vector<double> _vals;
    double fMed;
    double fRMS;
    double fMedHEta;
    double fRMSHEta;

    std::vector<PseudoJet> _neutrals;
    std::vector<PseudoJet> _chargedLV;
    std::vector<PseudoJet> _chargedPU;
    
    std::vector<int> _isPU;
    std::vector<int> _isCh;
    std::vector<int> _isPFCHS;
  
    std::vector<PseudoJet> _pfParticles_PU14;
    std::vector<PseudoJet> _pfchsParticles_PU14;

    TRandom *fRandom;
};

//FASTJET_END_NAMESPACE
