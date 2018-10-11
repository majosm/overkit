// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_ASSEMBLY_OPTIONS_INCLUDED
#define OVK_CORE_PUBLIC_ASSEMBLY_OPTIONS_INCLUDED

#include <ovk/core/ovkGlobal.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_assembly_options;
typedef struct ovk_assembly_options ovk_assembly_options;

void ovkGetAssemblyOptionsDimension(const ovk_assembly_options *Options, int *NumDims);
void ovkGetAssemblyOptionsGridCount(const ovk_assembly_options *Options, int *NumGrids);

void ovkGetAssemblyOptionOverlappable(const ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, bool *Overlappable);
void ovkSetAssemblyOptionOverlappable(ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, bool Overlappable);

void ovkGetAssemblyOptionOverlapTolerance(const ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, double *OverlapTolerance);
void ovkSetAssemblyOptionOverlapTolerance(ovk_assembly_options *Options, int OverlappingGridID,
  int OverlappedGridID, double OverlapTolerance);

void ovkGetAssemblyOptionOverlapAccelDepthAdjust(const ovk_assembly_options *Options,
  int OverlappingGridID, double *OverlapAccelDepthAdjust);
void ovkSetAssemblyOptionOverlapAccelDepthAdjust(ovk_assembly_options *Options,
  int OverlappingGridID, int OverlappedGridID, double OverlapAccelDepthAdjust);

void ovkGetAssemblyOptionOverlapAccelResolutionAdjust(const ovk_assembly_options *Options,
  int OverlappingGridID, double *OverlapAccelResolutionAdjust);
void ovkSetAssemblyOptionOverlapAccelResolutionAdjust(ovk_assembly_options *Options,
  int OverlappingGridID, double OverlapAccelResolutionAdjust);

void ovkGetAssemblyOptionInferBoundaries(const ovk_assembly_options *Options, int GridID,
  bool *InferBoundaries);
void ovkSetAssemblyOptionInferBoundaries(ovk_assembly_options *Options, int GridID,
  bool InferBoundaries);

void ovkGetAssemblyOptionCutBoundaryHoles(const ovk_assembly_options *Options, int CuttingGridID,
  int CutGridID, bool *CutBoundaryHoles);
void ovkSetAssemblyOptionCutBoundaryHoles(ovk_assembly_options *Options, int CuttingGridID,
  int CutGridID, bool CutBoundaryHoles);

void ovkGetAssemblyOptionOccludes(const ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, ovk_occludes *Occludes);
void ovkSetAssemblyOptionOccludes(ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, ovk_occludes Occludes);

void ovkGetAssemblyOptionEdgePadding(const ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, int *EdgePadding);
void ovkSetAssemblyOptionEdgePadding(ovk_assembly_options *Options, int OccludingGridID,
  int OccludedGridID, int EdgePadding);

void ovkGetAssemblyOptionEdgeSmoothing(const ovk_assembly_options *Options, int OccludedGridID,
  int *EdgeSmoothing);
void ovkSetAssemblyOptionEdgeSmoothing(ovk_assembly_options *Options, int OccludedGridID,
  int EdgeSmoothing);

void ovkGetAssemblyOptionConnectionType(const ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, ovk_connection_type *ConnectionType);
void ovkSetAssemblyOptionConnectionType(ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, ovk_connection_type ConnectionType);

void ovkGetAssemblyOptionFringeSize(const ovk_assembly_options *Options, int GridID, int *FringeSize);
void ovkSetAssemblyOptionFringeSize(ovk_assembly_options *Options, int GridID, int FringeSize);

void ovkGetAssemblyOptionMinimizeOverlap(const ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, bool *MinimizeOverlap);
void ovkSetAssemblyOptionMinimizeOverlap(ovk_assembly_options *Options, int DonorGridID,
  int ReceiverGridID, bool MinimizeOverlap);

#ifdef __cplusplus
}
#endif

#endif
