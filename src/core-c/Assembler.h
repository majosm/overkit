// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_ASSEMBLER_H_INCLUDED
#define OVK_CORE_C_ASSEMBLER_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Domain.h>
#include <ovk/core-c/Global.h>
#include <ovk/core/Assembler.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_assembler;
typedef struct ovk_assembler ovk_assembler;

struct ovk_assembler_params;
typedef struct ovk_assembler_params ovk_assembler_params;

struct ovk_assembler_bindings;
typedef struct ovk_assembler_bindings ovk_assembler_bindings;

struct ovk_assembler_options;
typedef struct ovk_assembler_options ovk_assembler_options;

void ovkCreateAssembler(ovk_assembler **Assembler, ovk_shared_context *Context, ovk_assembler_params
  **Params);
void ovkDestroyAssembler(ovk_assembler **Assembler);

void ovkGetAssemblerContextC(const ovk_assembler *Assembler, const ovk_context **Context);
void ovkGetAssemblerContext(ovk_assembler *Assembler, ovk_context **Context);
void ovkGetAssemblerSharedContext(ovk_assembler *Assembler, ovk_shared_context **Context);

bool ovkAssemblerIsBound(const ovk_assembler *Assembler);
void ovkBindAssembler(ovk_assembler *Assembler, ovk_domain *Domain, ovk_assembler_bindings
  **Bindings);
void ovkUnbindAssembler(ovk_assembler *Assembler);

void ovkGetAssemblerDomainC(const ovk_assembler *Assembler, const ovk_domain **Domain);
void ovkGetAssemblerDomain(ovk_assembler *Assembler, ovk_domain **Domain);

void ovkGetAssemblerOptions(const ovk_assembler *Assembler, const ovk_assembler_options **Options);
bool ovkEditingAssemblerOptions(const ovk_assembler *Assembler);
void ovkEditAssemblerOptions(ovk_assembler *Assembler, ovk_assembler_options **Options);
void ovkRestoreAssemblerOptions(ovk_assembler *Assembler, ovk_assembler_options **Options);

void ovkAssemble(ovk_assembler *Assembler);

void ovkCreateAssemblerParams(ovk_assembler_params **Params);
void ovkDestroyAssemblerParams(ovk_assembler_params **Params);
void ovkGetAssemblerParamName(const ovk_assembler_params *Params, char *Name);
void ovkSetAssemblerParamName(ovk_assembler_params *Params, const char *Name);

void ovkCreateAssemblerBindings(ovk_assembler_bindings **Bindings);
void ovkDestroyAssemblerBindings(ovk_assembler_bindings **Bindings);
void ovkGetAssemblerBindingsGeometryComponentID(ovk_assembler_bindings *Bindings, int
  *GeometryComponentID);
void ovkSetAssemblerBindingsGeometryComponentID(ovk_assembler_bindings *Bindings, int
  GeometryComponentID);
void ovkGetAssemblerBindingsStateComponentID(ovk_assembler_bindings *Bindings, int
  *StateComponentID);
void ovkSetAssemblerBindingsStateComponentID(ovk_assembler_bindings *Bindings, int
  StateComponentID);
void ovkGetAssemblerBindingsOverlapComponentID(ovk_assembler_bindings *Bindings, int
  *OverlapComponentID);
void ovkSetAssemblerBindingsOverlapComponentID(ovk_assembler_bindings *Bindings, int
  OverlapComponentID);
void ovkGetAssemblerBindingsConnectivityComponentID(ovk_assembler_bindings *Bindings, int
  *ConnectivityComponentID);
void ovkSetAssemblerBindingsConnectivityComponentID(ovk_assembler_bindings *Bindings, int
  ConnectivityComponentID);

void ovkGetAssemblerOptionOverlappable(const ovk_assembler_options *Options, int MGridID, int
  NGridID, bool *Overlappable);
void ovkSetAssemblerOptionOverlappable(ovk_assembler_options *Options, int MGridID, int NGridID,
  bool Overlappable);
void ovkResetAssemblerOptionOverlappable(ovk_assembler_options *Options, int MGridID, int NGridID);

void ovkGetAssemblerOptionOverlapTolerance(const ovk_assembler_options *Options, int MGridID, int
  NGridID, double *OverlapTolerance);
void ovkSetAssemblerOptionOverlapTolerance(ovk_assembler_options *Options, int MGridID, int NGridID,
  double OverlapTolerance);
void ovkResetAssemblerOptionOverlapTolerance(ovk_assembler_options *Options, int MGridID, int
  NGridID);

void ovkGetAssemblerOptionOverlapAccelDepthAdjust(const ovk_assembler_options *Options, int MGridID,
  double *OverlapAccelDepthAdjust);
void ovkSetAssemblerOptionOverlapAccelDepthAdjust(ovk_assembler_options *Options, int MGridID,
  double OverlapAccelDepthAdjust);
void ovkResetAssemblerOptionOverlapAccelDepthAdjust(ovk_assembler_options *Options, int MGridID);

void ovkGetAssemblerOptionOverlapAccelResolutionAdjust(const ovk_assembler_options *Options, int
  MGridID, double *OverlapAccelResolutionAdjust);
void ovkSetAssemblerOptionOverlapAccelResolutionAdjust(ovk_assembler_options *Options, int MGridID,
  double OverlapAccelResolutionAdjust);
void ovkResetAssemblerOptionOverlapAccelResolutionAdjust(ovk_assembler_options *Options, int
  MGridID);

void ovkGetAssemblerOptionInferBoundaries(const ovk_assembler_options *Options, int GridID, bool
  *InferBoundaries);
void ovkSetAssemblerOptionInferBoundaries(ovk_assembler_options *Options, int GridID, bool
  InferBoundaries);
void ovkResetAssemblerOptionInferBoundaries(ovk_assembler_options *Options, int GridID);

void ovkGetAssemblerOptionCutBoundaryHoles(const ovk_assembler_options *Options, int MGridID, int
  NGridID, bool *CutBoundaryHoles);
void ovkSetAssemblerOptionCutBoundaryHoles(ovk_assembler_options *Options, int MGridID, int NGridID,
  bool CutBoundaryHoles);
void ovkResetAssemblerOptionCutBoundaryHoles(ovk_assembler_options *Options, int MGridID, int
  NGridID);

void ovkGetAssemblerOptionOccludes(const ovk_assembler_options *Options, int MGridID, int NGridID,
  ovk_occludes *Occludes);
void ovkSetAssemblerOptionOccludes(ovk_assembler_options *Options, int MGridID, int NGridID,
  ovk_occludes Occludes);
void ovkResetAssemblerOptionOccludes(ovk_assembler_options *Options, int MGridID, int NGridID);

void ovkGetAssemblerOptionEdgePadding(const ovk_assembler_options *Options, int MGridID, int
  NGridID, int *EdgePadding);
void ovkSetAssemblerOptionEdgePadding(ovk_assembler_options *Options, int MGridID, int NGridID, int
  EdgePadding);
void ovkResetAssemblerOptionEdgePadding(ovk_assembler_options *Options, int MGridID, int NGridID);

void ovkGetAssemblerOptionEdgeSmoothing(const ovk_assembler_options *Options, int NGridID, int
  *EdgeSmoothing);
void ovkSetAssemblerOptionEdgeSmoothing(ovk_assembler_options *Options, int NGridID, int
  EdgeSmoothing);
void ovkResetAssemblerOptionEdgeSmoothing(ovk_assembler_options *Options, int NGridID);

void ovkGetAssemblerOptionConnectionType(const ovk_assembler_options *Options, int MGridID, int
  NGridID, ovk_connection_type *ConnectionType);
void ovkSetAssemblerOptionConnectionType(ovk_assembler_options *Options, int MGridID, int NGridID,
  ovk_connection_type ConnectionType);
void ovkResetAssemblerOptionConnectionType(ovk_assembler_options *Options, int MGridID, int
  NGridID);

void ovkGetAssemblerOptionFringeSize(const ovk_assembler_options *Options, int NGridID, int
  *FringeSize);
void ovkSetAssemblerOptionFringeSize(ovk_assembler_options *Options, int NGridID, int FringeSize);
void ovkResetAssemblerOptionFringeSize(ovk_assembler_options *Options, int NGridID);

void ovkGetAssemblerOptionMinimizeOverlap(const ovk_assembler_options *Options, int MGridID, int
  NGridID, bool *MinimizeOverlap);
void ovkSetAssemblerOptionMinimizeOverlap(ovk_assembler_options *Options, int MGridID, int NGridID,
  bool MinimizeOverlap);
void ovkResetAssemblerOptionMinimizeOverlap(ovk_assembler_options *Options, int MGridID, int
  NGridID);

void ovkGetAssemblerOptionDisjointConnections(const ovk_assembler_options *Options, int MGridID, int
  NGridID, bool *DisjointConnections);
void ovkSetAssemblerOptionDisjointConnections(ovk_assembler_options *Options, int MGridID, int
  NGridID, bool DisjointConnections);
void ovkResetAssemblerOptionDisjointConnections(ovk_assembler_options *Options, int MGridID, int
  NGridID);

#ifdef __cplusplus
}
#endif

#endif
