// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/AssemblyOptions.h"

#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Exchange.h"
#include "ovk/core/Global.h"
#include "ovk/core/Logger.h"
#include "ovk/core/OrderedMap.h"

#include <vector>
#include <map>

struct ovk_assembly_options {
  int num_dims;
  int num_grids;
  std::map<int, std::map<int, bool>> overlappable;
  std::map<int, std::map<int, double>> overlap_tolerance;
  std::map<int, double> overlap_accel_depth_adjust;
  std::map<int, double> overlap_accel_resolution_adjust;
  std::map<int, bool> infer_boundaries;
  std::map<int, std::map<int, bool>> cut_boundary_holes;
  std::map<int, std::map<int, ovk_occludes>> occludes;
  std::map<int, std::map<int, int>> edge_padding;
  std::map<int, int> edge_smoothing;
  std::map<int, std::map<int, ovk_connection_type>> connection_type;
  std::map<int, int> fringe_size;
  std::map<int, std::map<int, bool>> minimize_overlap;
};

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

  OVK_DEBUG_ASSERT(GridID >= 0 || GridID == OVK_ALL_GRIDS, "Invalid grid ID.");

  if (GridID == OVK_ALL_GRIDS) {
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

  OVK_DEBUG_ASSERT(GridID1 >= 0 || GridID1 == OVK_ALL_GRIDS, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridID2 >= 0 || GridID2 == OVK_ALL_GRIDS, "Invalid grid ID.");

  if (GridID1 == OVK_ALL_GRIDS && GridID2 == OVK_ALL_GRIDS) {
    for (auto &RowPair : Option) {
      for (auto &Pair : RowPair.second) {
        Pair.second = Value;
      }
    }
  } else if (GridID1 == OVK_ALL_GRIDS) {
    for (auto &RowPair : Option) {
      auto Iter = RowPair.second.find(GridID2);
      OVK_DEBUG_ASSERT(Iter != RowPair.second.end(), "Invalid grid ID.");
      Iter->second = Value;
    }
  } else if (GridID2 == OVK_ALL_GRIDS) {
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

void ovkCreateAssemblyOptions(ovk_assembly_options **Options_, int NumDims, int NumGrids,
  int *GridIDs) {

  OVK_DEBUG_ASSERT(Options_, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(NumDims == 2 || NumDims == 3, "Invalid dimension.");
  OVK_DEBUG_ASSERT(NumGrids >= 0, "Invalid grid count.");
  OVK_DEBUG_ASSERT(GridIDs, "Invalid grid IDs pointer.");

  *Options_ = new ovk_assembly_options();
  ovk_assembly_options *Options = *Options_;

  Options->num_dims = NumDims;
  Options->num_grids = NumGrids;

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options->overlappable[GridIDs[iGrid]][GridIDs[jGrid]] = false;
      Options->overlap_tolerance[GridIDs[iGrid]][GridIDs[jGrid]] = 1.e-12;
    }
    Options->overlap_accel_depth_adjust[GridIDs[iGrid]] = 0.;
    Options->overlap_accel_resolution_adjust[GridIDs[iGrid]] = 0.;
  }

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    Options->infer_boundaries[GridIDs[iGrid]] = false;
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options->cut_boundary_holes[GridIDs[iGrid]][GridIDs[jGrid]] = false;
    }
  }

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options->occludes[GridIDs[iGrid]][GridIDs[jGrid]] = OVK_OCCLUDES_NONE;
      Options->edge_padding[GridIDs[iGrid]][GridIDs[jGrid]] = 0;
    }
    Options->edge_smoothing[GridIDs[iGrid]] = 0;
  }

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options->connection_type[GridIDs[iGrid]][GridIDs[jGrid]] = OVK_CONNECTION_NONE;
    }
    Options->fringe_size[GridIDs[iGrid]] = 0;
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      Options->minimize_overlap[GridIDs[iGrid]][GridIDs[jGrid]] = false;
    }
  }

}

void ovkDestroyAssemblyOptions(ovk_assembly_options **Options) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(*Options, "Invalid options pointer.");

  delete *Options;
  *Options = NULL;

}

void ovkGetAssemblyOptionsDimension(const ovk_assembly_options *Options, int *NumDims) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Options->num_dims;

}

void ovkGetAssemblyOptionsGridCount(const ovk_assembly_options *Options, int *NumGrids) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(NumGrids, "Invalid num grids pointer.");

  *NumGrids = Options->num_grids;

}

void ovkGetAssemblyOptionOverlappable(const ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, bool *Overlappable) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlappingGridID >= 0, "Invalid overlapping grid ID.");
  OVK_DEBUG_ASSERT(OverlappedGridID >= 0, "Invalid overlapped grid ID.");
  OVK_DEBUG_ASSERT(Overlappable, "Invalid overlappable pointer.");

  GetOption(Options->overlappable, OverlappingGridID, OverlappedGridID, *Overlappable);

}

void ovkSetAssemblyOptionOverlappable(ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, bool Overlappable) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlappingGridID >= 0 || OverlappingGridID == OVK_ALL_GRIDS, "Invalid "
    "overlapping grid ID.");
  OVK_DEBUG_ASSERT(OverlappedGridID >= 0 || OverlappedGridID == OVK_ALL_GRIDS, "Invalid "
    "overlapped grid ID.");

  SetOption(Options->overlappable, OverlappingGridID, OverlappedGridID, Overlappable);

}

void ovkGetAssemblyOptionOverlapTolerance(const ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, double *OverlapTolerance) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlappingGridID >= 0, "Invalid overlapping grid ID.");
  OVK_DEBUG_ASSERT(OverlappedGridID >= 0, "Invalid overlapped grid ID.");
  OVK_DEBUG_ASSERT(OverlapTolerance, "Invalid overlap tolerance pointer.");

  GetOption(Options->overlap_tolerance, OverlappingGridID, OverlappedGridID, *OverlapTolerance);

}

void ovkSetAssemblyOptionOverlapTolerance(ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, double OverlapTolerance) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlappingGridID >= 0 || OverlappingGridID == OVK_ALL_GRIDS, "Invalid "
    "overlapping grid ID.");
  OVK_DEBUG_ASSERT(OverlappedGridID >= 0 || OverlappedGridID == OVK_ALL_GRIDS, "Invalid "
    "overlapped grid ID.");
  OVK_DEBUG_ASSERT(OverlapTolerance >= 0., "Invalid overlap tolerance value.");

  GetOption(Options->overlap_tolerance, OverlappingGridID, OverlappedGridID, OverlapTolerance);

}

void ovkGetAssemblyOptionOverlapAccelDepthAdjust(const ovk_assembly_options *Options,
  int OverlappingGridID, double *OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlappingGridID >= 0, "Invalid overlapping grid ID.");
  OVK_DEBUG_ASSERT(OverlapAccelDepthAdjust, "Invalid overlap accel depth adjust pointer.");

  GetOption(Options->overlap_accel_depth_adjust, OverlappingGridID, *OverlapAccelDepthAdjust);

}

void ovkSetAssemblyOptionOverlapAccelDepthAdjust(ovk_assembly_options *Options,
  int OverlappingGridID, int OverlappedGridID, double OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlappingGridID >= 0 || OverlappingGridID == OVK_ALL_GRIDS, "Invalid "
    "overlapping grid ID.");

  SetOption(Options->overlap_accel_depth_adjust, OverlappingGridID, OverlapAccelDepthAdjust);

}

void ovkGetAssemblyOptionOverlapAccelResolutionAdjust(const ovk_assembly_options *Options,
  int OverlappingGridID, double *OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlappingGridID >= 0, "Invalid overlapping grid ID.");
  OVK_DEBUG_ASSERT(OverlapAccelResolutionAdjust, "Invalid overlap accel resolution adjust pointer.");

  GetOption(Options->overlap_accel_resolution_adjust, OverlappingGridID, *OverlapAccelResolutionAdjust);

}

void ovkSetAssemblyOptionOverlapAccelResolutionAdjust(ovk_assembly_options *Options,
  int OverlappingGridID, double OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlappingGridID >= 0 || OverlappingGridID == OVK_ALL_GRIDS, "Invalid "
    "overlapping grid ID.");

  auto Iter = Options->overlap_accel_resolution_adjust.find(OverlappingGridID);
  OVK_DEBUG_ASSERT(Iter != Options->overlap_accel_resolution_adjust.end(), "Invalid overlapping "
    "grid ID.");

  SetOption(Options->overlap_accel_resolution_adjust, OverlappingGridID, OverlapAccelResolutionAdjust);

}

void ovkGetAssemblyOptionInferBoundaries(const ovk_assembly_options *Options, int GridID,
  bool *InferBoundaries) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(InferBoundaries, "Invalid infer boundaries pointer.");

  GetOption(Options->infer_boundaries, GridID, *InferBoundaries);

}

void ovkSetAssemblyOptionInferBoundaries(ovk_assembly_options *Options, int GridID,
  bool InferBoundaries) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0 || GridID == OVK_ALL_GRIDS, "Invalid grid ID.");

  SetOption(Options->infer_boundaries, GridID, InferBoundaries);

}

void ovkGetAssemblyOptionCutBoundaryHoles(const ovk_assembly_options *Options, int CuttingGridID,
  int CutGridID, bool *CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(CuttingGridID >= 0, "Invalid cutting grid ID.");
  OVK_DEBUG_ASSERT(CutGridID >= 0, "Invalid cut grid ID.");
  OVK_DEBUG_ASSERT(CutBoundaryHoles, "Invalid cut boundary holes pointer.");

  GetOption(Options->cut_boundary_holes, CuttingGridID, CutGridID, *CutBoundaryHoles);

}

void ovkSetAssemblyOptionCutBoundaryHoles(ovk_assembly_options *Options, int CuttingGridID,
  int CutGridID, bool CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(CuttingGridID >= 0 || CuttingGridID == OVK_ALL_GRIDS, "Invalid cutting grid ID.");
  OVK_DEBUG_ASSERT(CutGridID >= 0 || CutGridID == OVK_ALL_GRIDS, "Invalid cut grid ID.");

  SetOption(Options->cut_boundary_holes, CuttingGridID, CutGridID, CutBoundaryHoles);

}

void ovkGetAssemblyOptionOccludes(const ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, ovk_occludes *Occludes) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OccludingGridID >= 0, "Invalid occluding grid ID.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0, "Invalid occluded grid ID.");
  OVK_DEBUG_ASSERT(Occludes, "Invalid occludes pointer.");

  GetOption(Options->occludes, OccludingGridID, OccludedGridID, *Occludes);

}

void ovkSetAssemblyOptionOccludes(ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, ovk_occludes Occludes) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OccludingGridID >= 0 || OccludingGridID == OVK_ALL_GRIDS, "Invalid occluding "
    "grid ID.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0 || OccludedGridID == OVK_ALL_GRIDS, "Invalid occluded "
    "grid ID.");
  OVK_DEBUG_ASSERT(ValidOccludes(Occludes), "Invalid occludes value.");

  SetOption(Options->occludes, OccludingGridID, OccludedGridID, Occludes);

}

void ovkGetAssemblyOptionEdgePadding(const ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, int *EdgePadding) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OccludingGridID >= 0, "Invalid occluding grid ID.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0, "Invalid occluded grid ID.");
  OVK_DEBUG_ASSERT(EdgePadding, "Invalid edge padding pointer.");

  GetOption(Options->edge_padding, OccludingGridID, OccludedGridID, *EdgePadding);

}

void ovkSetAssemblyOptionEdgePadding(ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, int EdgePadding) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OccludingGridID >= 0 || OccludingGridID == OVK_ALL_GRIDS, "Invalid occluding "
    "grid ID.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0 || OccludedGridID == OVK_ALL_GRIDS, "Invalid occluded "
    "grid ID.");
  OVK_DEBUG_ASSERT(EdgePadding >= 0, "Invalid edge padding value.");

  SetOption(Options->edge_padding, OccludingGridID, OccludedGridID, EdgePadding);

}

void ovkGetAssemblyOptionEdgeSmoothing(const ovk_assembly_options *Options, int OccludedGridID,
  int *EdgeSmoothing) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0, "Invalid occluded grid ID.");
  OVK_DEBUG_ASSERT(EdgeSmoothing, "Invalid edge smoothing pointer.");

  GetOption(Options->edge_smoothing, OccludedGridID, *EdgeSmoothing);

}

void ovkSetAssemblyOptionEdgeSmoothing(ovk_assembly_options *Options, int OccludedGridID,
  int EdgeSmoothing) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OccludedGridID >= 0 || OccludedGridID == OVK_ALL_GRIDS, "Invalid occluded grid ID.");
  OVK_DEBUG_ASSERT(EdgeSmoothing >= 0, "Invalid edge smoothing value.");

  SetOption(Options->edge_smoothing, OccludedGridID, EdgeSmoothing);

}

void ovkGetAssemblyOptionConnectionType(const ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, ovk_connection_type *ConnectionType) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ConnectionType, "Invalid connection type pointer.");

  GetOption(Options->connection_type, DonorGridID, ReceiverGridID, *ConnectionType);

}

void ovkSetAssemblyOptionConnectionType(ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, ovk_connection_type ConnectionType) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0 || DonorGridID == OVK_ALL_GRIDS, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0 || ReceiverGridID == OVK_ALL_GRIDS, "Invalid receiver "
    "grid ID.");
  OVK_DEBUG_ASSERT(ValidConnectionType(ConnectionType), "Invalid connection type.");

  SetOption(Options->connection_type, DonorGridID, ReceiverGridID, ConnectionType);

}

void ovkGetAssemblyOptionFringeSize(const ovk_assembly_options *Options, int GridID, int *FringeSize) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(FringeSize, "Invalid fringe size pointer.");

  GetOption(Options->fringe_size, GridID, *FringeSize);

}

void ovkSetAssemblyOptionFringeSize(ovk_assembly_options *Options, int GridID, int FringeSize) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0 || GridID == OVK_ALL_GRIDS, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(FringeSize >= 0, "Invalid fringe size.");

  SetOption(Options->fringe_size, GridID, FringeSize);

}

void ovkGetAssemblyOptionMinimizeOverlap(const ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, bool *MinimizeOverlap) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(MinimizeOverlap, "Invalid minimize overlap pointer.");

  GetOption(Options->minimize_overlap, DonorGridID, ReceiverGridID, *MinimizeOverlap);

}

void ovkSetAssemblyOptionMinimizeOverlap(ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, bool MinimizeOverlap) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0 || DonorGridID == OVK_ALL_GRIDS, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0 || ReceiverGridID == OVK_ALL_GRIDS, "Invalid receiver "
    "grid ID.");

  SetOption(Options->minimize_overlap, DonorGridID, ReceiverGridID, MinimizeOverlap);

}
