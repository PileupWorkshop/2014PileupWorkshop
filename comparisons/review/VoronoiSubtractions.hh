// a series of simple subtractions basd on the voronoi particle area

// $Id$
//
// Copyright (c) 2013-, Gregory Soyez and Matteo Cacciari
//
// This code is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FastJet. If not, see <http://www.gnu.org/licenses/>.

#ifndef __VORONOI_SUBTRACTORS_HH__
#define __VORONOI_SUBTRACTORS_HH__

#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include "FunctionOfEvent.hh"
#include "VoronoiParticleArea.hh"

namespace fastjet{
namespace contrib{

//------------------------------------------------------------------------
/// \class VoronoiSubtractionBase
/// a base class to avoid excessive typing (some feature are too gross
/// for the very same purpose)
class VoronoiSubtractionBase : public FunctionOfEvent<std::vector<PseudoJet> >{
public:
  /// ctor with proper initialisation
  ///  \param bge  a pointer to the background estimator used to
  ///              estimate rho
  ///  \param vpa  a pointer to the computed voronoi areas
  VoronoiSubtractionBase(BackgroundEstimatorBase *bge=0, VoronoiParticleAreas *vpa=0) 
    : _bge(bge), _vpa(vpa){}

  /// a description of what this class is doing
  std::string description() const{
    return std::string("One of the many Voronoi subtraction variants");
  }

  // still to be overloaded
  std::vector<PseudoJet> result(const std::vector<PseudoJet> &event) const=0;

protected:
  BackgroundEstimatorBase *_bge; ///< (pointer to) the tool to estimate rho  
  VoronoiParticleAreas *_vpa;     ///< (pointer to) the particle areas
};


//------------------------------------------------------------------------
/// \class VoronoiNSigma
/// subtract rho A + N sigma sqrtA from each particle
class VoronoiNSigma : public VoronoiSubtractionBase{
public:
  VoronoiNSigma(BackgroundEstimatorBase *bge=0, VoronoiParticleAreas *vpa=0, double N=1.0) 
    : VoronoiSubtractionBase(bge, vpa), _N(N){}
  std::vector<PseudoJet> result(const std::vector<PseudoJet> &event) const;
protected:
  double _N;
};

//------------------------------------------------------------------------
/// \class VoronoiSubAndNSigma
/// subtract rho A then remove particles with pt< N sigma sqrtA
class VoronoiSubAndNSigma : public VoronoiSubtractionBase{
public:
  VoronoiSubAndNSigma(BackgroundEstimatorBase *bge=0, VoronoiParticleAreas *vpa=0, double N=1.0) 
    : VoronoiSubtractionBase(bge, vpa), _N(N){}
  std::vector<PseudoJet> result(const std::vector<PseudoJet> &event) const;
protected:
  double _N;
};

//------------------------------------------------------------------------
/// \class VoronoiSigmaCut
/// keep particles with pt > rho A + N sigma sqrtA
class VoronoiSigmaCut : public VoronoiSubtractionBase{
public:
  VoronoiSigmaCut(BackgroundEstimatorBase *bge=0, VoronoiParticleAreas *vpa=0, double N=1.0) 
    : VoronoiSubtractionBase(bge, vpa), _N(N){}
  std::vector<PseudoJet> result(const std::vector<PseudoJet> &event) const;
protected:
  double _N;
};

//------------------------------------------------------------------------
/// \class VoronoiSoftKiller2
/// subtract the softest particles until a fraction 'frac' of the event is empty
class VoronoiSoftKiller2 : public VoronoiSubtractionBase{
public:
  VoronoiSoftKiller2(BackgroundEstimatorBase *bge=0, VoronoiParticleAreas *vpa=0, double frac=0.5) : VoronoiSubtractionBase(bge, vpa), _frac(frac){}
  std::vector<PseudoJet> result(const std::vector<PseudoJet> &event) const;
protected:
  double _frac; ///< the fraction of the enpty event
};


//------------------------------------------------------------------------
/// \class VoronoiSoftKiller
/// Subtract softest particles until the total \sum rhoA_i is reached
///
/// description to be added
class VoronoiSoftKiller : public VoronoiSubtractionBase{
public:
  VoronoiSoftKiller(BackgroundEstimatorBase *bge=0, VoronoiParticleAreas *vpa=0) : VoronoiSubtractionBase(bge, vpa){}
  std::vector<PseudoJet> result(const std::vector<PseudoJet> &event) const;
};


//------------------------------------------------------------------------
/// \class VoronoiProgressiveSubtraction
/// Subtract rho*voronoi_particle_area from each particle.
///
/// description to be added
class VoronoiProgressiveSubtraction : public VoronoiSubtractionBase{
public:
  VoronoiProgressiveSubtraction(BackgroundEstimatorBase *bge=0, VoronoiParticleAreas *vpa=0) : VoronoiSubtractionBase(bge, vpa){}
  std::vector<PseudoJet> result(const std::vector<PseudoJet> &event) const;
};



//------------------------------------------------------------------------
/// \class VoronoiUniformSubtraction
/// Subtract rho*voronoi_particle_area from each particle.
///
/// The pt is set to 0 if pt<rho*A. And the excess is subtracted
/// uniformly on all other cells. More precisely, we do the following:
///
///   0. start with S = rho*Atot (=\sum_i \rho_i A_i). S is the total
///      amount to subtract.
///
///   1. from each positive cell subtract
///         S * (rho_i A_i/[\sum_i \rho_i A_i]).
///      If it is negative set it zero and accumulate the excess in
///      S_new
/// 
///   2. iterate 1 with S=S_new
///      [Note: at the 1st iteration, this amounts to subtract rhoi*Ai]
///
class VoronoiUniformSubtraction : public VoronoiSubtractionBase{
public:
  VoronoiUniformSubtraction(BackgroundEstimatorBase *bge=0, VoronoiParticleAreas *vpa=0) : VoronoiSubtractionBase(bge, vpa){}
  std::vector<PseudoJet> result(const std::vector<PseudoJet> &event) const;
};


} // namespace contrib
} // namespace fastjet

#endif //  __VORONOI_SUBTRACTORS_HH__
