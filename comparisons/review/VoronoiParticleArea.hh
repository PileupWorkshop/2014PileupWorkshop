#ifndef __VORONOI_PARTICLE_AREA_HH__
#define __VORONOI_PARTICLE_AREA_HH__

#include <fastjet/PseudoJet.hh>

/// \class VoronoiParticleAreas
/// associate to each particle the area of its Voronoi cell. 
/// 
/// Periodicity in phi in used and a fixed rapidity cut is used to
/// bound the infinite cells at the edge of the rapidity acceptance
class VoronoiParticleAreas{
public:
  VoronoiParticleAreas(double i_ymax) : ymax(i_ymax){}

  void compute_areas(const std::vector<fastjet::PseudoJet> &event);
  std::vector<double> areas;
  double ymax;
};

#endif
