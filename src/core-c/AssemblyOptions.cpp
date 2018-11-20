// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/AssemblyOptions.h"

#include "ovk/core-c/Constants.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/AssemblyOptions.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

extern "C" {

void ovkCreateAssemblyOptions(ovk_assembly_options **Options, int NumDims, int NumGrids,
  int *GridIDs) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto OptionsCPPPtr = new ovk::assembly_options();

  ovk::CreateAssemblyOptions(*OptionsCPPPtr, NumDims, NumGrids, GridIDs);

  *Options = reinterpret_cast<ovk_assembly_options *>(OptionsCPPPtr);

}

void ovkDestroyAssemblyOptions(ovk_assembly_options **Options) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(*Options, "Invalid options pointer.");

  auto OptionsCPPPtr = reinterpret_cast<ovk::assembly_options *>(*Options);

  ovk::DestroyAssemblyOptions(*OptionsCPPPtr);

  delete OptionsCPPPtr;

  *Options = NULL;

}

void ovkGetAssemblyOptionsDimension(const ovk_assembly_options *Options, int *NumDims) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionsDimension(OptionsCPP, *NumDims);

}

void ovkGetAssemblyOptionsGridCount(const ovk_assembly_options *Options, int *NumGrids) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(NumGrids, "Invalid num grids pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionsGridCount(OptionsCPP, *NumGrids);

}

void ovkGetAssemblyOptionOverlappable(const ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, bool *Overlappable) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(Overlappable, "Invalid overlappable pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionOverlappable(OptionsCPP, OverlappingGridID, OverlappedGridID,
    *Overlappable);

}

void ovkSetAssemblyOptionOverlappable(ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, bool Overlappable) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionOverlappable(OptionsCPP, OverlappingGridID, OverlappedGridID,
    Overlappable);

}

void ovkGetAssemblyOptionOverlapTolerance(const ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, double *OverlapTolerance) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlapTolerance, "Invalid overlap tolerance pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionOverlapTolerance(OptionsCPP, OverlappingGridID, OverlappedGridID,
    *OverlapTolerance);

}

void ovkSetAssemblyOptionOverlapTolerance(ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, double OverlapTolerance) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionOverlapTolerance(OptionsCPP, OverlappingGridID, OverlappedGridID,
    OverlapTolerance);

}

void ovkGetAssemblyOptionOverlapAccelDepthAdjust(const ovk_assembly_options *Options,
  int OverlappingGridID, double *OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlapAccelDepthAdjust, "Invalid overlap accel depth adjust pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionOverlapAccelDepthAdjust(OptionsCPP, OverlappingGridID,
    *OverlapAccelDepthAdjust);

}

void ovkSetAssemblyOptionOverlapAccelDepthAdjust(ovk_assembly_options *Options,
  int OverlappingGridID, double OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionOverlapAccelDepthAdjust(OptionsCPP, OverlappingGridID,
    OverlapAccelDepthAdjust);

}

void ovkGetAssemblyOptionOverlapAccelResolutionAdjust(const ovk_assembly_options *Options,
  int OverlappingGridID, double *OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlapAccelResolutionAdjust, "Invalid overlap accel resolution adjust pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionOverlapAccelResolutionAdjust(OptionsCPP, OverlappingGridID,
    *OverlapAccelResolutionAdjust);

}

void ovkSetAssemblyOptionOverlapAccelResolutionAdjust(ovk_assembly_options *Options,
  int OverlappingGridID, double OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionOverlapAccelResolutionAdjust(OptionsCPP, OverlappingGridID,
    OverlapAccelResolutionAdjust);

}

void ovkGetAssemblyOptionInferBoundaries(const ovk_assembly_options *Options, int GridID,
  bool *InferBoundaries) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(InferBoundaries, "Invalid infer boundaries pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionInferBoundaries(OptionsCPP, GridID, *InferBoundaries);

}

void ovkSetAssemblyOptionInferBoundaries(ovk_assembly_options *Options, int GridID,
  bool InferBoundaries) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionInferBoundaries(OptionsCPP, GridID, InferBoundaries);

}

void ovkGetAssemblyOptionCutBoundaryHoles(const ovk_assembly_options *Options, int CuttingGridID,
  int CutGridID, bool *CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(CutBoundaryHoles, "Invalid cut boundary holes pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionCutBoundaryHoles(OptionsCPP, CuttingGridID, CutGridID, *CutBoundaryHoles);

}

void ovkSetAssemblyOptionCutBoundaryHoles(ovk_assembly_options *Options, int CuttingGridID,
  int CutGridID, bool CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionCutBoundaryHoles(OptionsCPP, CuttingGridID, CutGridID, CutBoundaryHoles);

}

void ovkGetAssemblyOptionOccludes(const ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, ovk_occludes *Occludes) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(Occludes, "Invalid occludes pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);

  ovk::occludes OccludesCPP;
  ovk::GetAssemblyOptionOccludes(OptionsCPP, OccludingGridID, OccludedGridID, OccludesCPP);

  *Occludes = ovk_occludes(OccludesCPP);

}

void ovkSetAssemblyOptionOccludes(ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, ovk_occludes Occludes) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionOccludes(OptionsCPP, OccludingGridID, OccludedGridID,
    ovk::occludes(Occludes));

}

void ovkGetAssemblyOptionEdgePadding(const ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, int *EdgePadding) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(EdgePadding, "Invalid edge padding pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionEdgePadding(OptionsCPP, OccludingGridID, OccludedGridID, *EdgePadding);

}

void ovkSetAssemblyOptionEdgePadding(ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, int EdgePadding) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionEdgePadding(OptionsCPP, OccludingGridID, OccludedGridID, EdgePadding);

}

void ovkGetAssemblyOptionEdgeSmoothing(const ovk_assembly_options *Options, int OccludedGridID,
  int *EdgeSmoothing) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(EdgeSmoothing, "Invalid edge smoothing pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionEdgeSmoothing(OptionsCPP, OccludedGridID, *EdgeSmoothing);

}

void ovkSetAssemblyOptionEdgeSmoothing(ovk_assembly_options *Options, int OccludedGridID,
  int EdgeSmoothing) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionEdgeSmoothing(OptionsCPP, OccludedGridID, EdgeSmoothing);

}

void ovkGetAssemblyOptionConnectionType(const ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, ovk_connection_type *ConnectionType) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(ConnectionType, "Invalid connection type pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);

  ovk::connection_type ConnectionTypeCPP;
  ovk::GetAssemblyOptionConnectionType(OptionsCPP, DonorGridID, ReceiverGridID, ConnectionTypeCPP);

  *ConnectionType = ovk_connection_type(ConnectionTypeCPP);

}

void ovkSetAssemblyOptionConnectionType(ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, ovk_connection_type ConnectionType) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionConnectionType(OptionsCPP, DonorGridID, ReceiverGridID,
    ovk::connection_type(ConnectionType));

}

void ovkGetAssemblyOptionFringeSize(const ovk_assembly_options *Options, int GridID,
  int *FringeSize) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(FringeSize, "Invalid fringe size pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionFringeSize(OptionsCPP, GridID, *FringeSize);

}

void ovkSetAssemblyOptionFringeSize(ovk_assembly_options *Options, int GridID, int FringeSize) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionFringeSize(OptionsCPP, GridID, FringeSize);

}

void ovkGetAssemblyOptionMinimizeOverlap(const ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, bool *MinimizeOverlap) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(MinimizeOverlap, "Invalid minimize overlap pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::GetAssemblyOptionMinimizeOverlap(OptionsCPP, DonorGridID, ReceiverGridID,
    *MinimizeOverlap);

}

void ovkSetAssemblyOptionMinimizeOverlap(ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, bool MinimizeOverlap) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  ovk::SetAssemblyOptionMinimizeOverlap(OptionsCPP, DonorGridID, ReceiverGridID,
    MinimizeOverlap);

}

}
