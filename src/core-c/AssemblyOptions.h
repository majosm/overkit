// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_ASSEMBLY_OPTIONS_H_INCLUDED
#define OVK_CORE_C_ASSEMBLY_OPTIONS_H_INCLUDED

#include <ovk/core-c/Constants.h>
#include <ovk/core-c/Global.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_assembly_options;
typedef struct ovk_assembly_options ovk_assembly_options;

void ovkCreateAssemblyOptions(ovk_assembly_options **Options, int NumDims, int NumGrids,
  int *GridIDs);
void ovkDestroyAssemblyOptions(ovk_assembly_options **Options);

void ovkGetAssemblyOptionsDimension(const ovk_assembly_options *Options, int *NumDims);
void ovkGetAssemblyOptionsGridCount(const ovk_assembly_options *Options, int *NumGrids);

void ovkGetAssemblyOptionOverlappable(const ovk_assembly_options *Options, int MGridID, int NGridID,
  bool *Overlappable);
void ovkSetAssemblyOptionOverlappable(ovk_assembly_options *Options, int MGridID, int NGridID, bool
  Overlappable);

void ovkGetAssemblyOptionOverlapTolerance(const ovk_assembly_options *Options, int MGridID, int
  NGridID, double *OverlapTolerance);
void ovkSetAssemblyOptionOverlapTolerance(ovk_assembly_options *Options, int MGridID, int NGridID,
  double OverlapTolerance);

void ovkGetAssemblyOptionOverlapAccelDepthAdjust(const ovk_assembly_options *Options, int MGridID,
  double *OverlapAccelDepthAdjust);
void ovkSetAssemblyOptionOverlapAccelDepthAdjust(ovk_assembly_options *Options, int MGridID, double
  OverlapAccelDepthAdjust);

void ovkGetAssemblyOptionOverlapAccelResolutionAdjust(const ovk_assembly_options *Options, int
  MGridID, double *OverlapAccelResolutionAdjust);
void ovkSetAssemblyOptionOverlapAccelResolutionAdjust(ovk_assembly_options *Options, int MGridID,
  double OverlapAccelResolutionAdjust);

void ovkGetAssemblyOptionInferBoundaries(const ovk_assembly_options *Options, int GridID, bool
  *InferBoundaries);
void ovkSetAssemblyOptionInferBoundaries(ovk_assembly_options *Options, int GridID, bool
  InferBoundaries);

void ovkGetAssemblyOptionCutBoundaryHoles(const ovk_assembly_options *Options, int MGridID, int
  NGridID, bool *CutBoundaryHoles);
void ovkSetAssemblyOptionCutBoundaryHoles(ovk_assembly_options *Options, int MGridID, int NGridID,
  bool CutBoundaryHoles);

void ovkGetAssemblyOptionOccludes(const ovk_assembly_options *Options, int MGridID, int NGridID,
  ovk_occludes *Occludes);
void ovkSetAssemblyOptionOccludes(ovk_assembly_options *Options, int MGridID, int NGridID,
  ovk_occludes Occludes);

void ovkGetAssemblyOptionEdgePadding(const ovk_assembly_options *Options, int MGridID, int NGridID,
  int *EdgePadding);
void ovkSetAssemblyOptionEdgePadding(ovk_assembly_options *Options, int MGridID, int NGridID, int
  EdgePadding);

void ovkGetAssemblyOptionEdgeSmoothing(const ovk_assembly_options *Options, int NGridID, int
  *EdgeSmoothing);
void ovkSetAssemblyOptionEdgeSmoothing(ovk_assembly_options *Options, int NGridID, int
  EdgeSmoothing);

void ovkGetAssemblyOptionConnectionType(const ovk_assembly_options *Options, int MGridID, int
  NGridID, ovk_connection_type *ConnectionType);
void ovkSetAssemblyOptionConnectionType(ovk_assembly_options *Options, int MGridID, int NGridID,
  ovk_connection_type ConnectionType);

void ovkGetAssemblyOptionFringeSize(const ovk_assembly_options *Options, int NGridID, int
  *FringeSize);
void ovkSetAssemblyOptionFringeSize(ovk_assembly_options *Options, int NGridID, int FringeSize);

void ovkGetAssemblyOptionMinimizeOverlap(const ovk_assembly_options *Options, int MGridID, int
  NGridID, bool *MinimizeOverlap);
void ovkSetAssemblyOptionMinimizeOverlap(ovk_assembly_options *Options, int MGridID, int NGridID,
  bool MinimizeOverlap);

#ifdef __cplusplus
}
#endif

#endif
