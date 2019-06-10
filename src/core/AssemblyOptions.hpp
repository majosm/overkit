// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ASSEMBLY_OPTIONS_HPP_INCLUDED
#define OVK_CORE_ASSEMBLY_OPTIONS_HPP_INCLUDED

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/AssemblyOptionsBase.h>
#include <ovk/core/Global.hpp>
#include <ovk/core/IDMap.hpp>
#include <ovk/core/IDSet.hpp>

namespace ovk {

enum class occludes {
  NONE = OVK_OCCLUDES_NONE,
  ALL = OVK_OCCLUDES_ALL,
  COARSE = OVK_OCCLUDES_COARSE
};

inline bool ValidOccludes(occludes Occludes) {
  return ovkValidOccludes(ovk_occludes(Occludes));
}

enum class connection_type {
  NONE = OVK_CONNECTION_NONE,
  NEAREST = OVK_CONNECTION_NEAREST,
  LINEAR = OVK_CONNECTION_LINEAR,
  CUBIC = OVK_CONNECTION_CUBIC
};

inline bool ValidConnectionType(connection_type ConnectionType) {
  return ovkValidConnectionType(ovk_connection_type(ConnectionType));
}

class assembly_options {

public:

  assembly_options(int NumDims, const id_set<1> &GridIDs);

  int Dimension() const { return NumDims_; }
  int GridCount() const { return NumGrids_; }

  bool Overlappable(int MGridID, int NGridID) const;
  assembly_options &SetOverlappable(int MGridID, int NGridID, bool Overlappable);

  double OverlapTolerance(int MGridID, int NGridID) const;
  assembly_options &SetOverlapTolerance(int MGridID, int NGridID, double OverlapTolerance);

  double OverlapAccelDepthAdjust(int MGridID) const;
  assembly_options &SetOverlapAccelDepthAdjust(int MGridID, double OverlapAccelDepthAdjust);

  double OverlapAccelResolutionAdjust(int MGridID) const;
  assembly_options &SetOverlapAccelResolutionAdjust(int MGridID, double
    OverlapAccelResolutionAdjust);

  bool InferBoundaries(int GridID) const;
  assembly_options &SetInferBoundaries(int GridID, bool InferBoundaries);

  bool CutBoundaryHoles(int MGridID, int NGridID) const;
  assembly_options &SetCutBoundaryHoles(int MGridID, int NGridID, bool CutBoundaryHoles);

  occludes Occludes(int MGridID, int NGridID) const;
  assembly_options &SetOccludes(int MGridID, int NGridID, occludes Occludes);

  int EdgePadding(int MGridID, int NGridID) const;
  assembly_options &SetEdgePadding(int MGridID, int NGridID, int EdgePadding);

  int EdgeSmoothing(int NGridID) const;
  assembly_options &SetEdgeSmoothing(int NGridID, int EdgeSmoothing);

  connection_type ConnectionType(int MGridID, int NGridID) const;
  assembly_options &SetConnectionType(int MGridID, int NGridID, connection_type ConnectionType);

  int FringeSize(int NGridID) const;
  assembly_options &SetFringeSize(int NGridID, int FringeSize);

  bool MinimizeOverlap(int MGridID, int NGridID) const;
  assembly_options &SetMinimizeOverlap(int MGridID, int NGridID, bool MinimizeOverlap);

private:

  int NumDims_;
  int NumGrids_;
  id_map<2,bool> Overlappable_;
  id_map<2,double> OverlapTolerance_;
  id_map<1,double> OverlapAccelDepthAdjust_;
  id_map<1,double> OverlapAccelResolutionAdjust_;
  id_map<1,bool> InferBoundaries_;
  id_map<2,bool> CutBoundaryHoles_;
  id_map<2,occludes> Occludes_;
  id_map<2,int> EdgePadding_;
  id_map<1,int> EdgeSmoothing_;
  id_map<2,connection_type> ConnectionType_;
  id_map<1,int> FringeSize_;
  id_map<2,bool> MinimizeOverlap_;

  void PrintOptions_();

};

}

#endif
