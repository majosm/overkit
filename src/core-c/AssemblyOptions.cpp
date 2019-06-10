// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/AssemblyOptions.h"

#include "ovk/core-c/Global.h"
#include "ovk/core/AssemblyOptions.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/IDSet.hpp"

extern "C" {

void ovkCreateAssemblyOptions(ovk_assembly_options **Options, int NumDims, int NumGrids, int
  *GridIDs) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  ovk::id_set<1> GridIDsCPP(GridIDs, GridIDs+NumGrids);

  auto OptionsCPPPtr = new ovk::assembly_options(NumDims, GridIDsCPP);

  *Options = reinterpret_cast<ovk_assembly_options *>(OptionsCPPPtr);

}

void ovkDestroyAssemblyOptions(ovk_assembly_options **Options) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(*Options, "Invalid options pointer.");

  auto OptionsCPPPtr = reinterpret_cast<ovk::assembly_options *>(*Options);

  delete OptionsCPPPtr;

  *Options = NULL;

}

void ovkGetAssemblyOptionsDimension(const ovk_assembly_options *Options, int *NumDims) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *NumDims = OptionsCPP.Dimension();

}

void ovkGetAssemblyOptionsGridCount(const ovk_assembly_options *Options, int *NumGrids) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(NumGrids, "Invalid num grids pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *NumGrids = OptionsCPP.GridCount();

}

void ovkGetAssemblyOptionOverlappable(const ovk_assembly_options *Options, int MGridID, int NGridID,
  bool *Overlappable) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(Overlappable, "Invalid overlappable pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *Overlappable = OptionsCPP.Overlappable(MGridID, NGridID);

}

void ovkSetAssemblyOptionOverlappable(ovk_assembly_options *Options, int MGridID, int NGridID, bool
  Overlappable) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetOverlappable(MGridID, NGridID, Overlappable);

}

void ovkGetAssemblyOptionOverlapTolerance(const ovk_assembly_options *Options, int MGridID, int
  NGridID, double *OverlapTolerance) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlapTolerance, "Invalid overlap tolerance pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *OverlapTolerance = OptionsCPP.OverlapTolerance(MGridID, NGridID);

}

void ovkSetAssemblyOptionOverlapTolerance(ovk_assembly_options *Options, int MGridID, int NGridID,
  double OverlapTolerance) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetOverlapTolerance(MGridID, NGridID, OverlapTolerance);

}

void ovkGetAssemblyOptionOverlapAccelDepthAdjust(const ovk_assembly_options *Options, int MGridID,
  double *OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlapAccelDepthAdjust, "Invalid overlap accel depth adjust pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *OverlapAccelDepthAdjust = OptionsCPP.OverlapAccelDepthAdjust(MGridID);

}

void ovkSetAssemblyOptionOverlapAccelDepthAdjust(ovk_assembly_options *Options, int MGridID, double
  OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetOverlapAccelDepthAdjust(MGridID, OverlapAccelDepthAdjust);

}

void ovkGetAssemblyOptionOverlapAccelResolutionAdjust(const ovk_assembly_options *Options, int
  MGridID, double *OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlapAccelResolutionAdjust, "Invalid overlap accel resolution adjust pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *OverlapAccelResolutionAdjust = OptionsCPP.OverlapAccelResolutionAdjust(MGridID);

}

void ovkSetAssemblyOptionOverlapAccelResolutionAdjust(ovk_assembly_options *Options, int MGridID,
  double OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetOverlapAccelResolutionAdjust(MGridID, OverlapAccelResolutionAdjust);

}

void ovkGetAssemblyOptionInferBoundaries(const ovk_assembly_options *Options, int GridID, bool
  *InferBoundaries) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(InferBoundaries, "Invalid infer boundaries pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *InferBoundaries = OptionsCPP.InferBoundaries(GridID);

}

void ovkSetAssemblyOptionInferBoundaries(ovk_assembly_options *Options, int GridID, bool
  InferBoundaries) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetInferBoundaries(GridID, InferBoundaries);

}

void ovkGetAssemblyOptionCutBoundaryHoles(const ovk_assembly_options *Options, int MGridID, int
  NGridID, bool *CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(CutBoundaryHoles, "Invalid cut boundary holes pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *CutBoundaryHoles = OptionsCPP.CutBoundaryHoles(MGridID, NGridID);

}

void ovkSetAssemblyOptionCutBoundaryHoles(ovk_assembly_options *Options, int MGridID, int NGridID,
  bool CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetCutBoundaryHoles(MGridID, NGridID, CutBoundaryHoles);

}

void ovkGetAssemblyOptionOccludes(const ovk_assembly_options *Options, int MGridID, int NGridID,
  ovk_occludes *Occludes) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(Occludes, "Invalid occludes pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *Occludes = ovk_occludes(OptionsCPP.Occludes(MGridID, NGridID));

}

void ovkSetAssemblyOptionOccludes(ovk_assembly_options *Options, int MGridID, int NGridID,
  ovk_occludes Occludes) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetOccludes(MGridID, NGridID, ovk::occludes(Occludes));

}

void ovkGetAssemblyOptionEdgePadding(const ovk_assembly_options *Options, int MGridID, int NGridID,
  int *EdgePadding) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(EdgePadding, "Invalid edge padding pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *EdgePadding = OptionsCPP.EdgePadding(MGridID, NGridID);

}

void ovkSetAssemblyOptionEdgePadding(ovk_assembly_options *Options, int MGridID, int NGridID, int
  EdgePadding) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetEdgePadding(MGridID, NGridID, EdgePadding);

}

void ovkGetAssemblyOptionEdgeSmoothing(const ovk_assembly_options *Options, int NGridID, int
  *EdgeSmoothing) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(EdgeSmoothing, "Invalid edge smoothing pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *EdgeSmoothing = OptionsCPP.EdgeSmoothing(NGridID);

}

void ovkSetAssemblyOptionEdgeSmoothing(ovk_assembly_options *Options, int NGridID, int
  EdgeSmoothing) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetEdgeSmoothing(NGridID, EdgeSmoothing);

}

void ovkGetAssemblyOptionConnectionType(const ovk_assembly_options *Options, int MGridID, int
  NGridID, ovk_connection_type *ConnectionType) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(ConnectionType, "Invalid connection type pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *ConnectionType = ovk_connection_type(OptionsCPP.ConnectionType(MGridID, NGridID));

}

void ovkSetAssemblyOptionConnectionType(ovk_assembly_options *Options, int MGridID, int NGridID,
  ovk_connection_type ConnectionType) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetConnectionType(MGridID, NGridID, ovk::connection_type(ConnectionType));

}

void ovkGetAssemblyOptionFringeSize(const ovk_assembly_options *Options, int NGridID, int
  *FringeSize) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(FringeSize, "Invalid fringe size pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *FringeSize = OptionsCPP.FringeSize(NGridID);

}

void ovkSetAssemblyOptionFringeSize(ovk_assembly_options *Options, int NGridID, int FringeSize) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetFringeSize(NGridID, FringeSize);

}

void ovkGetAssemblyOptionMinimizeOverlap(const ovk_assembly_options *Options, int MGridID, int
  NGridID, bool *MinimizeOverlap) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(MinimizeOverlap, "Invalid minimize overlap pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  *MinimizeOverlap = OptionsCPP.MinimizeOverlap(MGridID, NGridID);

}

void ovkSetAssemblyOptionMinimizeOverlap(ovk_assembly_options *Options, int MGridID, int NGridID,
  bool MinimizeOverlap) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembly_options *>(Options);
  OptionsCPP.SetMinimizeOverlap(MGridID, NGridID, MinimizeOverlap);

}

}
