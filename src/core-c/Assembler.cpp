// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Assembler.h"

#include "ovk/core-c/Domain.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Assembler.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <cstring>
#include <string>
#include <utility>

extern "C" {

void ovkCreateAssembler(ovk_assembler **Assembler, ovk_shared_context *Context, ovk_assembler_params
  **Params) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ContextCPP = *reinterpret_cast<std::shared_ptr<ovk::context> *>(Context);

  ovk::assembler::params *ParamsCPPPtr = nullptr;
  if (Params && *Params) {
    ParamsCPPPtr = reinterpret_cast<ovk::assembler::params *>(*Params);
  }

  ovk::assembler *AssemblerCPPPtr;
  if (ParamsCPPPtr) {
    AssemblerCPPPtr = new ovk::assembler(ovk::CreateAssembler(ContextCPP, std::move(
      *ParamsCPPPtr)));
  } else {
    AssemblerCPPPtr = new ovk::assembler(ovk::CreateAssembler(ContextCPP));
  }

  *Assembler = reinterpret_cast<ovk_assembler *>(AssemblerCPPPtr);

  if (ParamsCPPPtr) {
    delete ParamsCPPPtr;
    *Params = nullptr;
  }

}

void ovkDestroyAssembler(ovk_assembler **Assembler) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(*Assembler, "Invalid assembler pointer.");

  auto AssemblerCPPPtr = reinterpret_cast<ovk::assembler *>(*Assembler);

  delete AssemblerCPPPtr;

  *Assembler = nullptr;

}

void ovkGetAssemblerContextC(const ovk_assembler *Assembler, const ovk_context **Context) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &AssemblerCPP = *reinterpret_cast<const ovk::assembler *>(Assembler);
  *Context = reinterpret_cast<const ovk_context *>(&AssemblerCPP.Context());

}

void ovkGetAssemblerContext(ovk_assembler *Assembler, ovk_context **Context) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &AssemblerCPP = *reinterpret_cast<ovk::assembler *>(Assembler);
  *Context = reinterpret_cast<ovk_context *>(&AssemblerCPP.Context());

}

void ovkGetAssemblerSharedContext(ovk_assembler *Assembler, ovk_shared_context **Context) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &AssemblerCPP = *reinterpret_cast<ovk::assembler *>(Assembler);
  auto &ContextCPP = AssemblerCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

bool ovkAssemblerIsBound(const ovk_assembler *Assembler) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");

  auto &AssemblerCPP = *reinterpret_cast<const ovk::assembler *>(Assembler);
  return AssemblerCPP.Bound();

}

void ovkBindAssembler(ovk_assembler *Assembler, ovk_domain *Domain, ovk_assembler_bindings
  **Bindings) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");
  OVK_DEBUG_ASSERT(*Bindings, "Invalid bindings pointer.");

  auto &AssemblerCPP = *reinterpret_cast<ovk::assembler *>(Assembler);
  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  auto BindingsCPPPtr = reinterpret_cast<ovk::assembler::bindings *>(*Bindings);

  AssemblerCPP.Bind(DomainCPP, std::move(*BindingsCPPPtr));

  delete BindingsCPPPtr;

  *Bindings = nullptr;

}

void ovkUnbindAssembler(ovk_assembler *Assembler) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");

  auto &AssemblerCPP = *reinterpret_cast<ovk::assembler *>(Assembler);
  AssemblerCPP.Unbind();

}

void ovkGetAssemblerDomainC(const ovk_assembler *Assembler, const ovk_domain **Domain) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &AssemblerCPP = *reinterpret_cast<const ovk::assembler *>(Assembler);
  *Domain = reinterpret_cast<const ovk_domain *>(&AssemblerCPP.Domain());

}

void ovkGetAssemblerDomain(ovk_assembler *Assembler, ovk_domain **Domain) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &AssemblerCPP = *reinterpret_cast<ovk::assembler *>(Assembler);
  *Domain = reinterpret_cast<ovk_domain *>(&AssemblerCPP.Domain());

}

void ovkGetAssemblerOptions(const ovk_assembler *Assembler, const ovk_assembler_options **Options) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &AssemblerCPP = *reinterpret_cast<const ovk::assembler *>(Assembler);
  *Options = reinterpret_cast<const ovk_assembler_options *>(&AssemblerCPP.Options());

}

bool ovkEditingAssemblerOptions(const ovk_assembler *Assembler) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");

  auto &AssemblerCPP = *reinterpret_cast<const ovk::assembler *>(Assembler);
  return AssemblerCPP.EditingOptions();

}

void ovkEditAssemblerOptions(ovk_assembler *Assembler, ovk_assembler_options **Options) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &AssemblerCPP = *reinterpret_cast<ovk::assembler *>(Assembler);

  ovk::edit_handle<ovk::assembler::options> EditHandle = AssemblerCPP.EditOptions();
  auto &OptionsCPP = *EditHandle.Release();

  *Options = reinterpret_cast<ovk_assembler_options *>(&OptionsCPP);

}

void ovkRestoreAssemblerOptions(ovk_assembler *Assembler, ovk_assembler_options **Options) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");
  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &AssemblerCPP = *reinterpret_cast<ovk::assembler *>(Assembler);
  AssemblerCPP.RestoreOptions();

  *Options = nullptr;

}

void ovkAssemble(ovk_assembler *Assembler) {

  OVK_DEBUG_ASSERT(Assembler, "Invalid assembler pointer.");

  auto &AssemblerCPP = *reinterpret_cast<ovk::assembler *>(Assembler);
  AssemblerCPP.Assemble();

}

void ovkCreateAssemblerParams(ovk_assembler_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::assembler::params();

  *Params = reinterpret_cast<ovk_assembler_params *>(ParamsCPPPtr);

}

void ovkDestroyAssemblerParams(ovk_assembler_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::assembler::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetAssemblerParamName(const ovk_assembler_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::assembler::params *>(Params);
  std::strcpy(Name, ParamsCPP.Name().c_str());

}

void ovkSetAssemblerParamName(ovk_assembler_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::assembler::params *>(Params);
  ParamsCPP.SetName(Name);

}

void ovkCreateAssemblerBindings(ovk_assembler_bindings **Bindings) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");

  auto BindingsCPPPtr = new ovk::assembler::bindings();

  *Bindings = reinterpret_cast<ovk_assembler_bindings *>(BindingsCPPPtr);

}

void ovkDestroyAssemblerBindings(ovk_assembler_bindings **Bindings) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");
  OVK_DEBUG_ASSERT(*Bindings, "Invalid bindings pointer.");

  auto BindingsCPPPtr = reinterpret_cast<ovk::assembler::bindings *>(*Bindings);

  delete BindingsCPPPtr;

  *Bindings = nullptr;

}

void ovkGetAssemblerBindingsGeometryComponentID(ovk_assembler_bindings *Bindings, int
  *GeometryComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");
  OVK_DEBUG_ASSERT(GeometryComponentID, "Invalid geometry component ID pointer.");

  auto &BindingsCPP = *reinterpret_cast<const ovk::assembler::bindings *>(Bindings);
  *GeometryComponentID = BindingsCPP.GeometryComponentID();

}

void ovkSetAssemblerBindingsGeometryComponentID(ovk_assembler_bindings *Bindings, int
  GeometryComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");

  auto &BindingsCPP = *reinterpret_cast<ovk::assembler::bindings *>(Bindings);
  BindingsCPP.SetGeometryComponentID(GeometryComponentID);

}

void ovkGetAssemblerBindingsStateComponentID(ovk_assembler_bindings *Bindings, int
  *StateComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");
  OVK_DEBUG_ASSERT(StateComponentID, "Invalid state component ID pointer.");

  auto &BindingsCPP = *reinterpret_cast<const ovk::assembler::bindings *>(Bindings);
  *StateComponentID = BindingsCPP.StateComponentID();

}

void ovkSetAssemblerBindingsStateComponentID(ovk_assembler_bindings *Bindings, int
  StateComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");

  auto &BindingsCPP = *reinterpret_cast<ovk::assembler::bindings *>(Bindings);
  BindingsCPP.SetStateComponentID(StateComponentID);

}

void ovkGetAssemblerBindingsOverlapComponentID(ovk_assembler_bindings *Bindings, int
  *OverlapComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");
  OVK_DEBUG_ASSERT(OverlapComponentID, "Invalid overlap component ID pointer.");

  auto &BindingsCPP = *reinterpret_cast<const ovk::assembler::bindings *>(Bindings);
  *OverlapComponentID = BindingsCPP.OverlapComponentID();

}

void ovkSetAssemblerBindingsOverlapComponentID(ovk_assembler_bindings *Bindings, int
  OverlapComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");

  auto &BindingsCPP = *reinterpret_cast<ovk::assembler::bindings *>(Bindings);
  BindingsCPP.SetOverlapComponentID(OverlapComponentID);

}

void ovkGetAssemblerBindingsConnectivityComponentID(ovk_assembler_bindings *Bindings, int
  *ConnectivityComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");
  OVK_DEBUG_ASSERT(ConnectivityComponentID, "Invalid connectivity component ID pointer.");

  auto &BindingsCPP = *reinterpret_cast<const ovk::assembler::bindings *>(Bindings);
  *ConnectivityComponentID = BindingsCPP.ConnectivityComponentID();

}

void ovkSetAssemblerBindingsConnectivityComponentID(ovk_assembler_bindings *Bindings, int
  ConnectivityComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");

  auto &BindingsCPP = *reinterpret_cast<ovk::assembler::bindings *>(Bindings);
  BindingsCPP.SetConnectivityComponentID(ConnectivityComponentID);

}

void ovkGetAssemblerOptionOverlappable(const ovk_assembler_options *Options, int MGridID, int
  NGridID, bool *Overlappable) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(Overlappable, "Invalid overlappable pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *Overlappable = OptionsCPP.Overlappable(MGridID, NGridID);

}

void ovkSetAssemblerOptionOverlappable(ovk_assembler_options *Options, int MGridID, int NGridID,
  bool Overlappable) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetOverlappable(MGridID, NGridID, Overlappable);

}

void ovkResetAssemblerOptionOverlappable(ovk_assembler_options *Options, int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetOverlappable(MGridID, NGridID);

}

void ovkGetAssemblerOptionOverlapTolerance(const ovk_assembler_options *Options, int MGridID, int
  NGridID, double *OverlapTolerance) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlapTolerance, "Invalid overlap tolerance pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *OverlapTolerance = OptionsCPP.OverlapTolerance(MGridID, NGridID);

}

void ovkSetAssemblerOptionOverlapTolerance(ovk_assembler_options *Options, int MGridID, int NGridID,
  double OverlapTolerance) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetOverlapTolerance(MGridID, NGridID, OverlapTolerance);

}

void ovkResetAssemblerOptionOverlapTolerance(ovk_assembler_options *Options, int MGridID, int
  NGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetOverlapTolerance(MGridID, NGridID);

}

void ovkGetAssemblerOptionOverlapAccelDepthAdjust(const ovk_assembler_options *Options, int MGridID,
  double *OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlapAccelDepthAdjust, "Invalid overlap accel depth adjust pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *OverlapAccelDepthAdjust = OptionsCPP.OverlapAccelDepthAdjust(MGridID);

}

void ovkSetAssemblerOptionOverlapAccelDepthAdjust(ovk_assembler_options *Options, int MGridID,
  double OverlapAccelDepthAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetOverlapAccelDepthAdjust(MGridID, OverlapAccelDepthAdjust);

}

void ovkResetAssemblerOptionOverlapAccelDepthAdjust(ovk_assembler_options *Options, int MGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetOverlapAccelDepthAdjust(MGridID);

}

void ovkGetAssemblerOptionOverlapAccelResolutionAdjust(const ovk_assembler_options *Options, int
  MGridID, double *OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(OverlapAccelResolutionAdjust, "Invalid overlap accel resolution adjust pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *OverlapAccelResolutionAdjust = OptionsCPP.OverlapAccelResolutionAdjust(MGridID);

}

void ovkSetAssemblerOptionOverlapAccelResolutionAdjust(ovk_assembler_options *Options, int MGridID,
  double OverlapAccelResolutionAdjust) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetOverlapAccelResolutionAdjust(MGridID, OverlapAccelResolutionAdjust);

}

void ovkResetAssemblerOptionOverlapAccelResolutionAdjust(ovk_assembler_options *Options, int
  MGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetOverlapAccelResolutionAdjust(MGridID);

}

void ovkGetAssemblerOptionInferBoundaries(const ovk_assembler_options *Options, int GridID, bool
  *InferBoundaries) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(InferBoundaries, "Invalid infer boundaries pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *InferBoundaries = OptionsCPP.InferBoundaries(GridID);

}

void ovkSetAssemblerOptionInferBoundaries(ovk_assembler_options *Options, int GridID, bool
  InferBoundaries) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetInferBoundaries(GridID, InferBoundaries);

}

void ovkResetAssemblerOptionInferBoundaries(ovk_assembler_options *Options, int GridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetInferBoundaries(GridID);

}

void ovkGetAssemblerOptionCutBoundaryHoles(const ovk_assembler_options *Options, int MGridID, int
  NGridID, bool *CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(CutBoundaryHoles, "Invalid cut boundary holes pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *CutBoundaryHoles = OptionsCPP.CutBoundaryHoles(MGridID, NGridID);

}

void ovkSetAssemblerOptionCutBoundaryHoles(ovk_assembler_options *Options, int MGridID, int NGridID,
  bool CutBoundaryHoles) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetCutBoundaryHoles(MGridID, NGridID, CutBoundaryHoles);

}

void ovkResetAssemblerOptionCutBoundaryHoles(ovk_assembler_options *Options, int MGridID, int
  NGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetCutBoundaryHoles(MGridID, NGridID);

}

void ovkGetAssemblerOptionOccludes(const ovk_assembler_options *Options, int MGridID, int NGridID,
  ovk_occludes *Occludes) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(Occludes, "Invalid occludes pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *Occludes = ovk_occludes(OptionsCPP.Occludes(MGridID, NGridID));

}

void ovkSetAssemblerOptionOccludes(ovk_assembler_options *Options, int MGridID, int NGridID,
  ovk_occludes Occludes) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetOccludes(MGridID, NGridID, ovk::occludes(Occludes));

}

void ovkResetAssemblerOptionOccludes(ovk_assembler_options *Options, int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetOccludes(MGridID, NGridID);

}

void ovkGetAssemblerOptionEdgePadding(const ovk_assembler_options *Options, int MGridID, int NGridID,
  int *EdgePadding) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(EdgePadding, "Invalid edge padding pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *EdgePadding = OptionsCPP.EdgePadding(MGridID, NGridID);

}

void ovkSetAssemblerOptionEdgePadding(ovk_assembler_options *Options, int MGridID, int NGridID, int
  EdgePadding) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetEdgePadding(MGridID, NGridID, EdgePadding);

}

void ovkResetAssemblerOptionEdgePadding(ovk_assembler_options *Options, int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetEdgePadding(MGridID, NGridID);

}

void ovkGetAssemblerOptionEdgeSmoothing(const ovk_assembler_options *Options, int NGridID, int
  *EdgeSmoothing) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(EdgeSmoothing, "Invalid edge smoothing pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *EdgeSmoothing = OptionsCPP.EdgeSmoothing(NGridID);

}

void ovkSetAssemblerOptionEdgeSmoothing(ovk_assembler_options *Options, int NGridID, int
  EdgeSmoothing) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetEdgeSmoothing(NGridID, EdgeSmoothing);

}

void ovkResetAssemblerOptionEdgeSmoothing(ovk_assembler_options *Options, int NGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetEdgeSmoothing(NGridID);

}

void ovkGetAssemblerOptionConnectionType(const ovk_assembler_options *Options, int MGridID, int
  NGridID, ovk_connection_type *ConnectionType) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(ConnectionType, "Invalid connection type pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *ConnectionType = ovk_connection_type(OptionsCPP.ConnectionType(MGridID, NGridID));

}

void ovkSetAssemblerOptionConnectionType(ovk_assembler_options *Options, int MGridID, int NGridID,
  ovk_connection_type ConnectionType) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetConnectionType(MGridID, NGridID, ovk::connection_type(ConnectionType));

}

void ovkResetAssemblerOptionConnectionType(ovk_assembler_options *Options, int MGridID, int NGridID)
  {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetConnectionType(MGridID, NGridID);

}

void ovkGetAssemblerOptionFringeSize(const ovk_assembler_options *Options, int NGridID, int
  *FringeSize) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(FringeSize, "Invalid fringe size pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *FringeSize = OptionsCPP.FringeSize(NGridID);

}

void ovkSetAssemblerOptionFringeSize(ovk_assembler_options *Options, int NGridID, int FringeSize) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetFringeSize(NGridID, FringeSize);

}

void ovkResetAssemblerOptionFringeSize(ovk_assembler_options *Options, int NGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetFringeSize(NGridID);

}

void ovkGetAssemblerOptionMinimizeOverlap(const ovk_assembler_options *Options, int MGridID, int
  NGridID, bool *MinimizeOverlap) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");
  OVK_DEBUG_ASSERT(MinimizeOverlap, "Invalid minimize overlap pointer.");

  auto &OptionsCPP = *reinterpret_cast<const ovk::assembler::options *>(Options);
  *MinimizeOverlap = OptionsCPP.MinimizeOverlap(MGridID, NGridID);

}

void ovkSetAssemblerOptionMinimizeOverlap(ovk_assembler_options *Options, int MGridID, int NGridID,
  bool MinimizeOverlap) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.SetMinimizeOverlap(MGridID, NGridID, MinimizeOverlap);

}

void ovkResetAssemblerOptionMinimizeOverlap(ovk_assembler_options *Options, int MGridID, int
  NGridID) {

  OVK_DEBUG_ASSERT(Options, "Invalid options pointer.");

  auto &OptionsCPP = *reinterpret_cast<ovk::assembler::options *>(Options);
  OptionsCPP.ResetMinimizeOverlap(MGridID, NGridID);

}

}
