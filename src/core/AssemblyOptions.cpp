// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/AssemblyOptions.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

#include <map>

namespace ovk {

namespace {

template <typename T> void GetOption(const std::map<int, T> &Option, int GridID, T &Value);
template <typename T> void GetOption(const std::map<int, std::map<int, T>> &Option, int GridID1,
  int GridID2, T &Value);
template <typename T> void SetOption(std::map<int, T> &Option, int GridID, T Value);
template <typename T> void SetOption(std::map<int, std::map<int, T>> &Option, int GridID1,
  int GridID2, T Value);

}

void CreateAssemblyOptions(assembly_options &Options, int NumDims, int NumGrids, int *GridIDs) {

  OVK_DEBUG_ASSERT(NumDims == 2 || NumDims == 3, "Invalid dimension.");
  OVK_DEBUG_ASSERT(NumGrids >= 0, "Invalid grid count.");
  OVK_DEBUG_ASSERT(GridIDs, "Invalid grid IDs pointer.");

  Options.NumDims_ = NumDims;
  Options.NumGrids_ = NumGrids;

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options.Overlappable_[GridIDs[iGrid]][GridIDs[jGrid]] = false;
      Options.OverlapTolerance_[GridIDs[iGrid]][GridIDs[jGrid]] = 1.e-12;
    }
    Options.OverlapAccelDepthAdjust_[GridIDs[iGrid]] = 0.;
    Options.OverlapAccelResolutionAdjust_[GridIDs[iGrid]] = 0.;
  }

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    Options.InferBoundaries_[GridIDs[iGrid]] = false;
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options.CutBoundaryHoles_[GridIDs[iGrid]][GridIDs[jGrid]] = false;
    }
  }

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options.Occludes_[GridIDs[iGrid]][GridIDs[jGrid]] = occludes::NONE;
      Options.EdgePadding_[GridIDs[iGrid]][GridIDs[jGrid]] = 0;
    }
    Options.EdgeSmoothing_[GridIDs[iGrid]] = 0;
  }

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options.ConnectionType_[GridIDs[iGrid]][GridIDs[jGrid]] = connection_type::NONE;
    }
    Options.FringeSize_[GridIDs[iGrid]] = 0;
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options.MinimizeOverlap_[GridIDs[iGrid]][GridIDs[jGrid]] = false;
    }
  }

}

void DestroyAssemblyOptions(assembly_options &Options) {

  Options.Overlappable_.clear();
  Options.OverlapTolerance_.clear();
  Options.OverlapAccelDepthAdjust_.clear();
  Options.OverlapAccelResolutionAdjust_.clear();

  Options.InferBoundaries_.clear();
  Options.CutBoundaryHoles_.clear();

  Options.Occludes_.clear();
  Options.EdgePadding_.clear();
  Options.EdgeSmoothing_.clear();

  Options.ConnectionType_.clear();
  Options.FringeSize_.clear();
  Options.MinimizeOverlap_.clear();

}

void GetAssemblyOptionsDimension(const assembly_options &Options, int &NumDims) {

  NumDims = Options.NumDims_;

}

void GetAssemblyOptionsGridCount(const assembly_options &Options, int &NumGrids) {

  NumGrids = Options.NumGrids_;

}

void GetAssemblyOptionOverlappable(const assembly_options &Options, int OverlappingGridID,
  int OverlappedGridID, bool &Overlappable) {

  OVK_DEBUG_ASSERT(OverlappingGridID >= 0, "Invalid overlapping grid ID.");
  OVK_DEBUG_ASSERT(OverlappedGridID >= 0, "Invalid overlapped grid ID.");

  GetOption(Options.Overlappable_, OverlappingGridID, OverlappedGridID, Overlappable);

}

void SetAssemblyOptionOverlappable(assembly_options &Options, int OverlappingGridID,
  int OverlappedGridID, bool Overlappable) {

  OVK_DEBUG_ASSERT(OverlappingGridID >= 0 || OverlappingGridID == ALL_GRIDS, "Invalid overlapping "
    "grid ID.");
  OVK_DEBUG_ASSERT(OverlappedGridID >= 0 || OverlappedGridID == ALL_GRIDS, "Invalid overlapped "
    "grid ID.");

  SetOption(Options.Overlappable_, OverlappingGridID, OverlappedGridID, Overlappable);

}

void GetAssemblyOptionOverlapTolerance(const assembly_options &Options, int OverlappingGridID,
  int OverlappedGridID, double &OverlapTolerance) {

  OVK_DEBUG_ASSERT(OverlappingGridID >= 0, "Invalid overlapping grid ID.");
  OVK_DEBUG_ASSERT(OverlappedGridID >= 0, "Invalid overlapped grid ID.");

  GetOption(Options.OverlapTolerance_, OverlappingGridID, OverlappedGridID, OverlapTolerance);

}

void SetAssemblyOptionOverlapTolerance(assembly_options &Options, int OverlappingGridID,
  int OverlappedGridID, double OverlapTolerance) {

  OVK_DEBUG_ASSERT(OverlappingGridID >= 0 || OverlappingGridID == ALL_GRIDS, "Invalid overlapping "
    "grid ID.");
  OVK_DEBUG_ASSERT(OverlappedGridID >= 0 || OverlappedGridID == ALL_GRIDS, "Invalid overlapped "
    "grid ID.");
  OVK_DEBUG_ASSERT(OverlapTolerance >= 0., "Invalid overlap tolerance value.");

  GetOption(Options.OverlapTolerance_, OverlappingGridID, OverlappedGridID, OverlapTolerance);

}

void GetAssemblyOptionOverlapAccelDepthAdjust(const assembly_options &Options,
  int OverlappingGridID, double &OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(OverlappingGridID >= 0, "Invalid overlapping grid ID.");

  GetOption(Options.OverlapAccelDepthAdjust_, OverlappingGridID, OverlapAccelDepthAdjust);

}

void SetAssemblyOptionOverlapAccelDepthAdjust(assembly_options &Options,
  int OverlappingGridID, double OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(OverlappingGridID >= 0 || OverlappingGridID == ALL_GRIDS, "Invalid overlapping "
    "grid ID.");

  SetOption(Options.OverlapAccelDepthAdjust_, OverlappingGridID, OverlapAccelDepthAdjust);

}

void GetAssemblyOptionOverlapAccelResolutionAdjust(const assembly_options &Options,
  int OverlappingGridID, double &OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(OverlappingGridID >= 0, "Invalid overlapping grid ID.");

  GetOption(Options.OverlapAccelResolutionAdjust_, OverlappingGridID, OverlapAccelResolutionAdjust);

}

void SetAssemblyOptionOverlapAccelResolutionAdjust(assembly_options &Options,
  int OverlappingGridID, double OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(OverlappingGridID >= 0 || OverlappingGridID == ALL_GRIDS, "Invalid overlapping "
    "grid ID.");

  auto Iter = Options.OverlapAccelResolutionAdjust_.find(OverlappingGridID);
  OVK_DEBUG_ASSERT(Iter != Options.OverlapAccelResolutionAdjust_.end(), "Invalid overlapping "
    "grid ID.");

  SetOption(Options.OverlapAccelResolutionAdjust_, OverlappingGridID, OverlapAccelResolutionAdjust);

}

void GetAssemblyOptionInferBoundaries(const assembly_options &Options, int GridID,
  bool &InferBoundaries) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");

  GetOption(Options.InferBoundaries_, GridID, InferBoundaries);

}

void SetAssemblyOptionInferBoundaries(assembly_options &Options, int GridID,
  bool InferBoundaries) {

  OVK_DEBUG_ASSERT(GridID >= 0 || GridID == ALL_GRIDS, "Invalid grid ID.");

  SetOption(Options.InferBoundaries_, GridID, InferBoundaries);

}

void GetAssemblyOptionCutBoundaryHoles(const assembly_options &Options, int CuttingGridID,
  int CutGridID, bool &CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(CuttingGridID >= 0, "Invalid cutting grid ID.");
  OVK_DEBUG_ASSERT(CutGridID >= 0, "Invalid cut grid ID.");

  GetOption(Options.CutBoundaryHoles_, CuttingGridID, CutGridID, CutBoundaryHoles);

}

void SetAssemblyOptionCutBoundaryHoles(assembly_options &Options, int CuttingGridID,
  int CutGridID, bool CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(CuttingGridID >= 0 || CuttingGridID == ALL_GRIDS, "Invalid cutting grid ID.");
  OVK_DEBUG_ASSERT(CutGridID >= 0 || CutGridID == ALL_GRIDS, "Invalid cut grid ID.");

  SetOption(Options.CutBoundaryHoles_, CuttingGridID, CutGridID, CutBoundaryHoles);

}

void GetAssemblyOptionOccludes(const assembly_options &Options, int OccludingGridID,
  int OccludedGridID, occludes &Occludes) {

  OVK_DEBUG_ASSERT(OccludingGridID >= 0, "Invalid occluding grid ID.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0, "Invalid occluded grid ID.");

  GetOption(Options.Occludes_, OccludingGridID, OccludedGridID, Occludes);

}

void SetAssemblyOptionOccludes(assembly_options &Options, int OccludingGridID,
  int OccludedGridID, occludes Occludes) {

  OVK_DEBUG_ASSERT(OccludingGridID >= 0 || OccludingGridID == ALL_GRIDS, "Invalid occluding "
    "grid ID.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0 || OccludedGridID == ALL_GRIDS, "Invalid occluded "
    "grid ID.");
  OVK_DEBUG_ASSERT(ValidOccludes(Occludes), "Invalid occludes value.");

  SetOption(Options.Occludes_, OccludingGridID, OccludedGridID, Occludes);

}

void GetAssemblyOptionEdgePadding(const assembly_options &Options, int OccludingGridID,
  int OccludedGridID, int &EdgePadding) {

  OVK_DEBUG_ASSERT(OccludingGridID >= 0, "Invalid occluding grid ID.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0, "Invalid occluded grid ID.");

  GetOption(Options.EdgePadding_, OccludingGridID, OccludedGridID, EdgePadding);

}

void SetAssemblyOptionEdgePadding(assembly_options &Options, int OccludingGridID,
  int OccludedGridID, int EdgePadding) {

  OVK_DEBUG_ASSERT(OccludingGridID >= 0 || OccludingGridID == ALL_GRIDS, "Invalid occluding "
    "grid ID.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0 || OccludedGridID == ALL_GRIDS, "Invalid occluded "
    "grid ID.");
  OVK_DEBUG_ASSERT(EdgePadding >= 0, "Invalid edge padding value.");

  SetOption(Options.EdgePadding_, OccludingGridID, OccludedGridID, EdgePadding);

}

void GetAssemblyOptionEdgeSmoothing(const assembly_options &Options, int OccludedGridID,
  int &EdgeSmoothing) {

  OVK_DEBUG_ASSERT(OccludedGridID >= 0, "Invalid occluded grid ID.");

  GetOption(Options.EdgeSmoothing_, OccludedGridID, EdgeSmoothing);

}

void SetAssemblyOptionEdgeSmoothing(assembly_options &Options, int OccludedGridID,
  int EdgeSmoothing) {

  OVK_DEBUG_ASSERT(OccludedGridID >= 0 || OccludedGridID == ALL_GRIDS, "Invalid occluded grid ID.");
  OVK_DEBUG_ASSERT(EdgeSmoothing >= 0, "Invalid edge smoothing value.");

  SetOption(Options.EdgeSmoothing_, OccludedGridID, EdgeSmoothing);

}

void GetAssemblyOptionConnectionType(const assembly_options &Options, int DonorGridID,
  int ReceiverGridID, connection_type &ConnectionType) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");

  GetOption(Options.ConnectionType_, DonorGridID, ReceiverGridID, ConnectionType);

}

void SetAssemblyOptionConnectionType(assembly_options &Options, int DonorGridID,
  int ReceiverGridID, connection_type ConnectionType) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0 || DonorGridID == ALL_GRIDS, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0 || ReceiverGridID == ALL_GRIDS, "Invalid receiver "
    "grid ID.");
  OVK_DEBUG_ASSERT(ValidConnectionType(ConnectionType), "Invalid connection type.");

  SetOption(Options.ConnectionType_, DonorGridID, ReceiverGridID, ConnectionType);

}

void GetAssemblyOptionFringeSize(const assembly_options &Options, int GridID, int &FringeSize) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");

  GetOption(Options.FringeSize_, GridID, FringeSize);

}

void SetAssemblyOptionFringeSize(assembly_options &Options, int GridID, int FringeSize) {

  OVK_DEBUG_ASSERT(GridID >= 0 || GridID == ALL_GRIDS, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(FringeSize >= 0, "Invalid fringe size.");

  SetOption(Options.FringeSize_, GridID, FringeSize);

}

void GetAssemblyOptionMinimizeOverlap(const assembly_options &Options, int DonorGridID,
  int ReceiverGridID, bool &MinimizeOverlap) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");

  GetOption(Options.MinimizeOverlap_, DonorGridID, ReceiverGridID, MinimizeOverlap);

}

void SetAssemblyOptionMinimizeOverlap(assembly_options &Options, int DonorGridID,
  int ReceiverGridID, bool MinimizeOverlap) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0 || DonorGridID == ALL_GRIDS, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0 || ReceiverGridID == ALL_GRIDS, "Invalid receiver "
    "grid ID.");

  SetOption(Options.MinimizeOverlap_, DonorGridID, ReceiverGridID, MinimizeOverlap);

}

namespace {

template <typename T> void GetOption(const std::map<int, T> &Option, int GridID, T &Value) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");

  auto Iter = Option.find(GridID);
  OVK_DEBUG_ASSERT(Iter != Option.end(), "Invalid grid ID.");
  Value = Iter->second;

}

template <typename T> void GetOption(const std::map<int, std::map<int, T>> &Option, int GridID1,
  int GridID2, T &Value) {

  OVK_DEBUG_ASSERT(GridID1 >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridID2 >= 0, "Invalid grid ID.");

  auto RowIter = Option.find(GridID1);
  OVK_DEBUG_ASSERT(RowIter != Option.end(), "Invalid grid ID.");
  auto Iter = RowIter->second.find(GridID2);
  OVK_DEBUG_ASSERT(Iter != RowIter->second.end(), "Invalid grid ID.");
  Value = Iter->second;

}

template <typename T> void SetOption(std::map<int, T> &Option, int GridID, T Value) {

  OVK_DEBUG_ASSERT(GridID >= 0 || GridID == ALL_GRIDS, "Invalid grid ID.");

  if (GridID == ALL_GRIDS) {
    for (auto &Pair : Option) {
      Pair.second = Value;
    }
  } else {
    auto Iter = Option.find(GridID);
    OVK_DEBUG_ASSERT(Iter != Option.end(), "Invalid grid ID.");
    Iter->second = Value;
  }

}

template <typename T> void SetOption(std::map<int, std::map<int, T>> &Option, int GridID1,
  int GridID2, T Value) {

  OVK_DEBUG_ASSERT(GridID1 >= 0 || GridID1 == ALL_GRIDS, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridID2 >= 0 || GridID2 == ALL_GRIDS, "Invalid grid ID.");

  if (GridID1 == ALL_GRIDS && GridID2 == ALL_GRIDS) {
    for (auto &RowPair : Option) {
      for (auto &Pair : RowPair.second) {
        Pair.second = Value;
      }
    }
  } else if (GridID1 == ALL_GRIDS) {
    for (auto &RowPair : Option) {
      auto Iter = RowPair.second.find(GridID2);
      OVK_DEBUG_ASSERT(Iter != RowPair.second.end(), "Invalid grid ID.");
      Iter->second = Value;
    }
  } else if (GridID2 == ALL_GRIDS) {
    auto RowIter = Option.find(GridID1);
    OVK_DEBUG_ASSERT(RowIter != Option.end(), "Invalid grid ID.");
    for (auto &Pair : RowIter->second) {
      Pair.second = Value;
    }
  } else {
    auto RowIter = Option.find(GridID1);
    OVK_DEBUG_ASSERT(RowIter != Option.end(), "Invalid grid ID.");
    auto Iter = RowIter->second.find(GridID2);
    OVK_DEBUG_ASSERT(Iter != RowIter->second.end(), "Invalid grid ID.");
    Iter->second = Value;
  }

}

} 

}
