// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/AssemblyOptions.hpp"

#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/IDMap.hpp"
#include "ovk/core/IDSet.hpp"

namespace ovk {

assembly_options::assembly_options(int NumDims, id_set<1> GridIDs):
  NumDims_(NumDims),
  GridIDs_(std::move(GridIDs))
{
  OVK_DEBUG_ASSERT(NumDims_ == 2 || NumDims_ == 3, "Invalid dimension.");
}

bool assembly_options::Overlappable(int MGridID, int NGridID) const {

  return GetOption_(Overlappable_, MGridID, NGridID, false);

}

assembly_options &assembly_options::SetOverlappable(int MGridID, int NGridID, bool Overlappable) {

  SetOption_(Overlappable_, MGridID, NGridID, Overlappable, false);

  return *this;

}

assembly_options &assembly_options::ResetOverlappable(int MGridID, int NGridID) {

  SetOption_(Overlappable_, MGridID, NGridID, false, false);

  return *this;

}

double assembly_options::OverlapTolerance(int MGridID, int NGridID) const {

  return GetOption_(OverlapTolerance_, MGridID, NGridID, 1.e-12);

}

assembly_options &assembly_options::SetOverlapTolerance(int MGridID, int NGridID, double
  OverlapTolerance) {

  OVK_DEBUG_ASSERT(OverlapTolerance >= 0., "Invalid overlap tolerance value.");

  SetOption_(OverlapTolerance_, MGridID, NGridID, OverlapTolerance, 1.e-12);

  return *this;

}

assembly_options &assembly_options::ResetOverlapTolerance(int MGridID, int NGridID) {

  SetOption_(OverlapTolerance_, MGridID, NGridID, 1.e-12, 1.e-12);

  return *this;

}

double assembly_options::OverlapAccelDepthAdjust(int MGridID) const {

  return GetOption_(OverlapAccelDepthAdjust_, MGridID, 0.);

}

assembly_options &assembly_options::SetOverlapAccelDepthAdjust(int MGridID, double
  OverlapAccelDepthAdjust) {

  SetOption_(OverlapAccelDepthAdjust_, MGridID, OverlapAccelDepthAdjust, 0.);

  return *this;

}

assembly_options &assembly_options::ResetOverlapAccelDepthAdjust(int MGridID) {

  SetOption_(OverlapAccelDepthAdjust_, MGridID, 0., 0.);

  return *this;

}

double assembly_options::OverlapAccelResolutionAdjust(int MGridID) const {

  return GetOption_(OverlapAccelResolutionAdjust_, MGridID, 0.);

}

assembly_options &assembly_options::SetOverlapAccelResolutionAdjust(int MGridID, double
  OverlapAccelResolutionAdjust) {

  SetOption_(OverlapAccelResolutionAdjust_, MGridID, OverlapAccelResolutionAdjust, 0.);

  return *this;

}

assembly_options &assembly_options::ResetOverlapAccelResolutionAdjust(int MGridID) {

  SetOption_(OverlapAccelResolutionAdjust_, MGridID, 0., 0.);

  return *this;

}

bool assembly_options::InferBoundaries(int GridID) const {

  return GetOption_(InferBoundaries_, GridID, false);

}

assembly_options &assembly_options::SetInferBoundaries(int GridID, bool InferBoundaries) {

  SetOption_(InferBoundaries_, GridID, InferBoundaries, false);

  return *this;

}

assembly_options &assembly_options::ResetInferBoundaries(int GridID) {

  SetOption_(InferBoundaries_, GridID, false, false);

  return *this;

}

bool assembly_options::CutBoundaryHoles(int MGridID, int NGridID) const {

  return GetOption_(CutBoundaryHoles_, MGridID, NGridID, false);

}

assembly_options &assembly_options::SetCutBoundaryHoles(int MGridID, int NGridID, bool
  CutBoundaryHoles) {

  SetOption_(CutBoundaryHoles_, MGridID, NGridID, CutBoundaryHoles, false);

  return *this;

}

assembly_options &assembly_options::ResetCutBoundaryHoles(int MGridID, int NGridID) {

  SetOption_(CutBoundaryHoles_, MGridID, NGridID, false, false);

  return *this;

}

occludes assembly_options::Occludes(int MGridID, int NGridID) const {

  return GetOption_(Occludes_, MGridID, NGridID, occludes::NONE);

}

assembly_options &assembly_options::SetOccludes(int MGridID, int NGridID, occludes Occludes) {

  OVK_DEBUG_ASSERT(ValidOccludes(Occludes), "Invalid occludes value.");

  SetOption_(Occludes_, MGridID, NGridID, Occludes, occludes::NONE);

  return *this;

}

assembly_options &assembly_options::ResetOccludes(int MGridID, int NGridID) {

  SetOption_(Occludes_, MGridID, NGridID, occludes::NONE, occludes::NONE);

  return *this;

}

int assembly_options::EdgePadding(int MGridID, int NGridID) const {

  return GetOption_(EdgePadding_, MGridID, NGridID, 0);

}

assembly_options &assembly_options::SetEdgePadding(int MGridID, int NGridID, int EdgePadding) {

  OVK_DEBUG_ASSERT(EdgePadding >= 0, "Invalid edge padding value.");

  SetOption_(EdgePadding_, MGridID, NGridID, EdgePadding, 0);

  return *this;

}

assembly_options &assembly_options::ResetEdgePadding(int MGridID, int NGridID) {

  SetOption_(EdgePadding_, MGridID, NGridID, 0, 0);

  return *this;

}

int assembly_options::EdgeSmoothing(int NGridID) const {

  return GetOption_(EdgeSmoothing_, NGridID, 0);

}

assembly_options &assembly_options::SetEdgeSmoothing(int NGridID, int EdgeSmoothing) {

  OVK_DEBUG_ASSERT(EdgeSmoothing >= 0, "Invalid edge smoothing value.");

  SetOption_(EdgeSmoothing_, NGridID, EdgeSmoothing, 0);

  return *this;

}

assembly_options &assembly_options::ResetEdgeSmoothing(int NGridID) {

  SetOption_(EdgeSmoothing_, NGridID, 0, 0);

  return *this;

}

connection_type assembly_options::ConnectionType(int MGridID, int NGridID) const {

  return GetOption_(ConnectionType_, MGridID, NGridID, connection_type::NONE);

}

assembly_options &assembly_options::SetConnectionType(int MGridID, int NGridID, connection_type
  ConnectionType) {

  OVK_DEBUG_ASSERT(ValidConnectionType(ConnectionType), "Invalid connection type.");

  SetOption_(ConnectionType_, MGridID, NGridID, ConnectionType, connection_type::NONE);

  return *this;

}

assembly_options &assembly_options::ResetConnectionType(int MGridID, int NGridID) {

  SetOption_(ConnectionType_, MGridID, NGridID, connection_type::NONE, connection_type::NONE);

  return *this;

}

int assembly_options::FringeSize(int NGridID) const {

  return GetOption_(FringeSize_, NGridID, 0);

}

assembly_options &assembly_options::SetFringeSize(int NGridID, int FringeSize) {

  OVK_DEBUG_ASSERT(FringeSize >= 0, "Invalid fringe size.");

  SetOption_(FringeSize_, NGridID, FringeSize, 0);

  return *this;

}

assembly_options &assembly_options::ResetFringeSize(int NGridID) {

  SetOption_(FringeSize_, NGridID, 0, 0);

  return *this;

}

bool assembly_options::MinimizeOverlap(int MGridID, int NGridID) const {

  return GetOption_(MinimizeOverlap_, MGridID, NGridID, false);

}

assembly_options &assembly_options::SetMinimizeOverlap(int MGridID, int NGridID, bool
  MinimizeOverlap) {

  SetOption_(MinimizeOverlap_, MGridID, NGridID, MinimizeOverlap, false);

  return *this;

}

assembly_options &assembly_options::ResetMinimizeOverlap(int MGridID, int NGridID) {

  SetOption_(MinimizeOverlap_, MGridID, NGridID, false, false);

  return *this;

}

void assembly_options::PrintOptions_() {

  for (auto &Entry : Overlappable_) {
    std::printf("Overlappable(%i,%i) = %c\n", Entry.Key(0), Entry.Key(1), Entry.Value() ? 'T' :
      'F');
  }

  for (auto &Entry : OverlapTolerance_) {
    std::printf("OverlapTolerance(%i,%i) = %16.8f\n", Entry.Key(0), Entry.Key(1), Entry.Value());
  }

  for (auto &Entry : OverlapAccelDepthAdjust_) {
    std::printf("OverlapAccelDepthAdjust(%i) = %16.8f\n", Entry.Key(0), Entry.Value());
  }

  for (auto &Entry : OverlapAccelResolutionAdjust_) {
    std::printf("OverlapAccelResolutionAdjust(%i) = %16.8f\n", Entry.Key(0), Entry.Value());
  }

  for (auto &Entry : InferBoundaries_) {
    std::printf("InferBoundaries(%i) = %c\n", Entry.Key(0), Entry.Value() ? 'T' : 'F');
  }

  for (auto &Entry : CutBoundaryHoles_) {
    std::printf("CutBoundaryHoles(%i,%i) = %c\n", Entry.Key(0), Entry.Key(1), Entry.Value() ? 'T' :
      'F');
  }

  for (auto &Entry : Occludes_) {
    std::string OccludesString;
    switch (Entry.Value()) {
    case occludes::NONE:
      OccludesString = "NONE";
      break;
    case occludes::ALL:
      OccludesString = "ALL";
      break;
    case occludes::COARSE:
      OccludesString = "COARSE";
      break;
    }
    std::printf("Occludes(%i,%i) = %s\n", Entry.Key(0), Entry.Key(1), OccludesString.c_str());
  }

  for (auto &Entry : EdgePadding_) {
    std::printf("EdgePadding(%i,%i) = %i\n", Entry.Key(0), Entry.Key(1), Entry.Value());
  }

  for (auto &Entry : EdgeSmoothing_) {
    std::printf("EdgeSmoothing(%i) = %i\n", Entry.Key(0), Entry.Value());
  }

  for (auto &Entry : ConnectionType_) {
    std::string ConnectionTypeString;
    switch (Entry.Value()) {
    case connection_type::NONE:
      ConnectionTypeString = "NONE";
      break;
    case connection_type::NEAREST:
      ConnectionTypeString = "NEAREST";
      break;
    case connection_type::LINEAR:
      ConnectionTypeString = "LINEAR";
      break;
    case connection_type::CUBIC:
      ConnectionTypeString = "CUBIC";
      break;
    }
    std::printf("ConnectionType(%i,%i) = %s\n", Entry.Key(0), Entry.Key(1),
      ConnectionTypeString.c_str());
  }

  for (auto &Entry : MinimizeOverlap_) {
    std::printf("MinimizeOverlap(%i,%i) = %c\n", Entry.Key(0), Entry.Key(1), Entry.Value() ? 'T' :
      'F');
  }

  for (auto &Entry : FringeSize_) {
    std::printf("FringeSize(%i) = %i\n", Entry.Key(0), Entry.Value());
  }

}

template <typename T> T assembly_options::GetOption_(const id_map<1,T> &Option, int GridID, T
  DefaultValue) const {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridIDs_.Contains(GridID), "Invalid grid ID.");

  T Value = DefaultValue;

  auto Iter = Option.Find(GridID);
  if (Iter != Option.End()) {
    Value = Iter->Value();
  }

  return Value;

}

template <typename T> T assembly_options::GetOption_(const id_map<2,T> &Option, int MGridID, int
  NGridID, T DefaultValue) const {

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(GridIDs_.Contains(MGridID), "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(GridIDs_.Contains(NGridID), "Invalid N grid ID.");

  T Value = DefaultValue;

  auto Iter = Option.Find(MGridID,NGridID);
  if (Iter != Option.End()) {
    Value = Iter->Value();
  }

  return Value;

}

template <typename T> void assembly_options::SetOption_(id_map<1,T> &Option, int GridID, T Value, T
  DefaultValue) {

  OVK_DEBUG_ASSERT(GridID == ALL_GRIDS || GridIDs_.Contains(GridID), "Invalid grid ID.");

  if (Value != DefaultValue) {
    if (GridID == ALL_GRIDS) {
      for (int ID : GridIDs_) {
        Option.Insert(ID, Value);
      }
    } else {
      Option.Insert(GridID, Value);
    }
  } else {
    if (GridID == ALL_GRIDS) {
      Option.Clear();
    } else {
      Option.Erase(GridID);
    }
  }

}

template <typename T> void assembly_options::SetOption_(id_map<2,T> &Option, int MGridID, int
  NGridID, T Value, T DefaultValue) {

  OVK_DEBUG_ASSERT(MGridID == ALL_GRIDS || GridIDs_.Contains(MGridID), "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID == ALL_GRIDS || GridIDs_.Contains(NGridID), "Invalid N grid ID.");

  if (Value != DefaultValue) {
    if (MGridID == ALL_GRIDS && NGridID == ALL_GRIDS) {
      for (int MID : GridIDs_) {
        for (int NID : GridIDs_) {
          Option.Insert({MID,NID}, Value);
        }
      }
    } else if (MGridID == ALL_GRIDS) {
      for (int MID : GridIDs_) {
        Option.Insert({MID,NGridID}, Value);
      }
    } else if (NGridID == ALL_GRIDS) {
      for (int NID : GridIDs_) {
        Option.Insert({MGridID,NID}, Value);
      }
    } else {
      Option.Insert({MGridID,NGridID}, Value);
    }
  } else {
    if (MGridID == ALL_GRIDS && NGridID == ALL_GRIDS) {
      Option.Clear();
    } else if (MGridID == ALL_GRIDS) {
      for (int MID : GridIDs_) {
        Option.Erase(MID,NGridID);
      }
    } else if (NGridID == ALL_GRIDS) {
      for (int NID : GridIDs_) {
        Option.Erase(MGridID,NID);
      }
    } else {
      Option.Erase(MGridID,NGridID);
    }
  }

}

}
