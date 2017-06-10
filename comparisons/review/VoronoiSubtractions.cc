// a series of simple subtractions basd on the voronoi particle area

#include "VoronoiSubtractions.hh"
#include <algorithm>

using namespace std;

namespace fastjet{
namespace contrib{


//------------------------------------------------------------------------
// class VoronoiNSigma
//------------------------------------------------------------------------
vector<PseudoJet> VoronoiNSigma::result(const vector<PseudoJet> &event) const{
  vector<PseudoJet> subtracted_event;
  
  for (unsigned int i=0;i<event.size();i++){
    double ptsub = _bge->rho(event[i])*_vpa->areas[i]
                 + _N*_bge->sigma(event[i])*sqrt(_vpa->areas[i]);
    if (event[i].pt2() > ptsub*ptsub) {
      PseudoJet j = event[i];
      j.reset_momentum(PtYPhiM(j.pt()-ptsub, j.rap(), j.phi()));
      subtracted_event.push_back(j);
    }
  }

  return subtracted_event;
}

//------------------------------------------------------------------------
// class VoronoiSubAndNSigma
//------------------------------------------------------------------------
vector<PseudoJet> VoronoiSubAndNSigma::result(const vector<PseudoJet> &event) const{
  vector<PseudoJet> subtracted_event;
  
  for (unsigned int i=0;i<event.size();i++){
    PseudoJet j = event[i];
    double ptsub = j.pt() - _bge->rho(j)*_vpa->areas[i];
    if (ptsub > _N*_bge->sigma(j)*sqrt(_vpa->areas[i])) {
      j.reset_momentum(PtYPhiM(ptsub, j.rap(), j.phi()));
      subtracted_event.push_back(j);
    }
  }

  return subtracted_event;
}

//------------------------------------------------------------------------
// class VoronoiSigmaCut
//------------------------------------------------------------------------
vector<PseudoJet> VoronoiSigmaCut::result(const vector<PseudoJet> &event) const{
  vector<PseudoJet> subtracted_event;
  
  for (unsigned int i=0;i<event.size();i++){
    double ptcut = _bge->rho(event[i])*_vpa->areas[i]
                 + _N*_bge->sigma(event[i])*sqrt(_vpa->areas[i]);
    if (event[i].pt2() > ptcut*ptcut) 
      subtracted_event.push_back(event[i]);
  }

  return subtracted_event;
}

//------------------------------------------------------------------------
// class VoronoiSoftKiller2
//
// subtract the softest particles until half of the event is empty
//------------------------------------------------------------------------
vector<PseudoJet> VoronoiSoftKiller2::result(const vector<PseudoJet> &event) const{
  double area_to_subtract = _frac*(twopi*2*_vpa->ymax);

  vector<PseudoJet> subtracted_event=sorted_by_pt(event);
 
  unsigned int i=event.size();
  double area_sub=0.0;
  while (area_sub < area_to_subtract){
    if (i==0){ return vector<PseudoJet>(); }
    i--;
    area_sub += _vpa->areas[i];
  }
  subtracted_event.resize(i);
  
  return subtracted_event;
}

//------------------------------------------------------------------------
// class VoronoiSoftKiller
//
// Sun rho_i A_i for all the particles and subtract it starting from
// the softest particle
//------------------------------------------------------------------------
vector<PseudoJet> VoronoiSoftKiller::result(const vector<PseudoJet> &event) const{
  // compute the default rho A for each cell
  double sub_tot=0.0;
  for (unsigned int i=0; i<event.size(); i++){
    sub_tot += _vpa->areas[i] * _bge->rho(event[i]);
  }

  vector<PseudoJet> subtracted_event=sorted_by_pt(event);
 
  unsigned int i=event.size();
  double pt_sub=0.0;
  while (pt_sub < sub_tot){
    if (i==0){ return vector<PseudoJet>(); }
    i--;
    pt_sub += subtracted_event[i].pt();
  }
  subtracted_event.resize(i);    
  
  return subtracted_event;
}


//------------------------------------------------------------------------

//------------------------------------------------------------------------
vector<PseudoJet> VoronoiProgressiveSubtraction::result(const vector<PseudoJet> &event) const{

  // compute the default rho A for each cell
  vector<pair<double, unsigned int> > subtracted_pt(event.size());
  for (unsigned int i=0; i<event.size(); i++){
    subtracted_pt[i].first  = event[i].pt() - _vpa->areas[i] * _bge->rho(event[i]);
    subtracted_pt[i].second = i;
  }

  // now balance the -ve with the smallest +ve
  //
  //  0. order the subtracted guys
  sort(subtracted_pt.begin(), subtracted_pt.end());

  //  1. sum the -ve part
  unsigned int i=0;
  double left_to_subtract=0.0;
  while (subtracted_pt[i].first<0)
    left_to_subtract -= subtracted_pt[i++].first;

  //  2. subtract in on the +ve
  while (left_to_subtract>0){
    if (left_to_subtract > subtracted_pt[i].first){
      left_to_subtract -= subtracted_pt[i++].first;
    } else {
      subtracted_pt[i].first -= left_to_subtract;
      left_to_subtract = -1.0;
    }
  }

  vector<PseudoJet> subtracted_event;
  while (i< event.size()){
    unsigned int ip = subtracted_pt[i].second;
    PseudoJet s = event[ip];
    s.reset_momentum(PtYPhiM(subtracted_pt[i].first, event[ip].rap(), event[ip].phi()));
    subtracted_event.push_back(s);
    i++;
  }
      
  return subtracted_event;
}

//------------------------------------------------------------------------
// class VoronoiUniformSubtraction
// Subtract rho*voronoi_particle_area from each particle.
//
// The pt is set to 0 if pt<rho*A. And the excess is subtracted
// uniformly on all other cells. More precisely, we do the following:
//
//   0. start with S = rho*Atot (=\sum_i \rho_i A_i). S is the total
//      amount to subtract.
//
//   1. from each positive cell subtract
//         S * (rho_i A_i/[\sum_i \rho_i A_i]).
//      If it is negative set it zero and accumulate the excess in
//      S_new
// 
//   2. iterate 1 with S=S_new
//      [Note: at the 1st iteration, this amounts to subtract rhoi*Ai]
//
//------------------------------------------------------------------------

vector<PseudoJet> VoronoiUniformSubtraction::result(const vector<PseudoJet> &event) const{

  // compute the default rho A for each cell
  vector<double> rhoA(event.size());
  double rhoA_tot=0.0;
  for (unsigned int i=0; i<event.size(); i++){
    rhoA[i] = _vpa->areas[i] * _bge->rho(event[i]);
    rhoA_tot += rhoA[i];
  }

  // start from the full event
  vector<PseudoJet> subtracted = event;

  // now iterate the subtraction
  double tot_to_subtract = rhoA_tot;
  do{
    double still_to_subtract = 0.0;
    double rhoA_tot_next = 0.0;
    bool any_new_negative=false;

    for (unsigned int i=0; i<event.size(); i++){
      if (subtracted[i] == 0) continue;

      double to_subtract = tot_to_subtract * (rhoA[i]/rhoA_tot);

      if (subtracted[i].pt2() > to_subtract*to_subtract){
        subtracted[i].reset_momentum(PtYPhiM(subtracted[i].pt()-to_subtract, subtracted[i].rap(), subtracted[i].phi()));
        rhoA_tot_next += rhoA[i];
      } else {
        still_to_subtract += to_subtract-subtracted[i].pt();
        subtracted[i] = PseudoJet();
        any_new_negative = true;
      }
    }

    // stop if no iteration needed
    if (! any_new_negative) break;

    // prepare the next iteration
    rhoA_tot = rhoA_tot_next;
    tot_to_subtract = still_to_subtract;

  } while (true);

  return (!SelectorIsZero())(subtracted);

}


} // namespace contrib
} // namespace fastjet

