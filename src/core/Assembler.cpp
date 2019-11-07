// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Assembler.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/ConnectivityComponent.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/ElemSet.hpp"
#include "ovk/core/Event.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/GeometryComponent.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/OverlapComponent.hpp"
#include "ovk/core/Set.hpp"
#include "ovk/core/StateComponent.hpp"

#include <mpi.h>

#include <memory>
#include <string>
#include <utility>

namespace ovk {

assembler::assembler(std::shared_ptr<context> &&Context, params &&Params):
  Context_(std::move(Context)),
  Name_(std::move(*Params.Name_))
{}

assembler::~assembler() noexcept {

  if (Domain_) {
    Unbind();
  }

}

assembler assembler::internal_Create(std::shared_ptr<context> &&Context, params &&Params) {

  return {std::move(Context), std::move(Params)};

}

assembler CreateAssembler(std::shared_ptr<context> Context, assembler::params Params) {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return assembler::internal_Create(std::move(Context), std::move(Params));

}

const context &assembler::Context() const {

  return *Context_;

}

context &assembler::Context() {

  return *Context_;

}

const std::shared_ptr<context> &assembler::SharedContext() const {

  return Context_;

}

bool assembler::Bound() const {

  return static_cast<bool>(Domain_);

}

void assembler::Bind(domain &Domain, bindings Bindings) {

  OVK_DEBUG_ASSERT(!Domain_, "Assembler is already bound to a domain.");

  Domain_ = Domain.GetFloatingRef();

  floating_ref<assembler> FloatingRef = FloatingRefGenerator_.Generate(*this);

  GridEventListener_ = Domain.AddGridEventListener([FloatingRef](int GridID, grid_event_flags Flags,
    bool LastInSequence) {
    assembler &Assembler = *FloatingRef;
    Assembler.OnGridEvent_(GridID, Flags, LastInSequence);
  });

  ComponentEventListener_ = Domain.AddComponentEventListener([FloatingRef](int ComponentID,
    component_event_flags Flags) {
    assembler &Assembler = *FloatingRef;
    Assembler.OnComponentEvent_(ComponentID, Flags);
  });

  GeometryComponentID_ = Bindings.GeometryComponentID_;
  OVK_DEBUG_ASSERT(GeometryComponentID_ >= 0, "Invalid geometry component ID.");
  OVK_DEBUG_ASSERT(Domain.ComponentExists(GeometryComponentID_), "Component %i does not exist.",
    GeometryComponentID_);

  auto &GeometryComponent = Domain.Component<geometry_component>(GeometryComponentID_);
  GeometryEventListener_ = GeometryComponent.AddGeometryEventListener([FloatingRef](int GridID,
    geometry_event_flags Flags, bool LastInSequence) {
    assembler &Assembler = *FloatingRef;
    Assembler.OnGeometryEvent_(GridID, Flags, LastInSequence);
  });

  StateComponentID_ = Bindings.StateComponentID_;
  OVK_DEBUG_ASSERT(StateComponentID_ >= 0, "Invalid state component ID.");
  OVK_DEBUG_ASSERT(Domain.ComponentExists(StateComponentID_), "Component %i does not exist.",
    StateComponentID_);

  auto &StateComponent = Domain.Component<state_component>(StateComponentID_);
  StateEventListener_ = StateComponent.AddStateEventListener([FloatingRef](int GridID,
    state_event_flags Flags, bool LastInSequence) {
    assembler &Assembler = *FloatingRef;
    Assembler.OnStateEvent_(GridID, Flags, LastInSequence);
  });

  OverlapComponentID_ = Bindings.OverlapComponentID_;
  OVK_DEBUG_ASSERT(OverlapComponentID_ >= 0, "Invalid overlap component ID.");
  OVK_DEBUG_ASSERT(Domain.ComponentExists(OverlapComponentID_), "Component %i does not exist.",
    OverlapComponentID_);

  auto &OverlapComponent = Domain.Component<overlap_component>(OverlapComponentID_);
  OverlapEventListener_ = OverlapComponent.AddOverlapEventListener([FloatingRef](const elem<int,2>
    &OverlapID, overlap_event_flags Flags, bool LastInSequence) {
    assembler &Assembler = *FloatingRef;
    Assembler.OnOverlapEvent_(OverlapID, Flags, LastInSequence);
  });

  ConnectivityComponentID_ = Bindings.ConnectivityComponentID_;
  OVK_DEBUG_ASSERT(ConnectivityComponentID_ >= 0, "Invalid connectivity component ID.");
  OVK_DEBUG_ASSERT(Domain.ComponentExists(ConnectivityComponentID_), "Component %i does not exist.",
    ConnectivityComponentID_);

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);
  ConnectivityEventListener_ = ConnectivityComponent.AddConnectivityEventListener([FloatingRef](
    const elem<int,2> &ConnectivityID, connectivity_event_flags Flags, bool LastInSequence) {
    assembler &Assembler = *FloatingRef;
    Assembler.OnConnectivityEvent_(ConnectivityID, Flags, LastInSequence);
  });

  MPI_Barrier(Domain.Comm());

  AssemblyData_ = assembly_data(Domain.Dimension(), Domain.Comm());

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.Comm().Rank() == 0, "Bound assembler %s to domain %s.", *Name_,
    Domain.Name());

  for (int GridID : Domain.GridIDs()) {
    UpdateManifest_.AddGridsToOptions.Insert(GridID);
  }

  Update_();

}

void assembler::Unbind() {

  OVK_DEBUG_ASSERT(Domain_, "Assembler is not bound to a domain.");

  const domain &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  Options_ = options();

  ConnectivityEventListener_.Reset();
  ConnectivityComponentID_ = -1;

  OverlapEventListener_.Reset();
  OverlapComponentID_ = -1;

  StateEventListener_.Reset();
  StateComponentID_ = -1;

  GeometryEventListener_.Reset();
  GeometryComponentID_ = -1;

  ComponentEventListener_.Reset();
  GridEventListener_.Reset();
  Domain_.Reset();

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.Comm().Rank() == 0, "Unbound assembler %s.", *Name_);

}

const domain &assembler::Domain() const {

  OVK_DEBUG_ASSERT(Domain_, "Assembler is not bound to a domain.");

  return *Domain_;

}

domain &assembler::Domain() {

  OVK_DEBUG_ASSERT(Domain_, "Assembler is not bound to a domain.");

  return *Domain_;

}

bool assembler::EditingOptions() const {

  return OptionsEditor_.Active();

}

edit_handle<assembler::options> assembler::EditOptions() {

  const domain &Domain = *Domain_;

  if (!OptionsEditor_.Active()) {
    MPI_Barrier(Domain.Comm());
    OnOptionsStartEdit_();
    floating_ref<assembler> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      assembler &Assembler = *FloatingRef;
      Assembler.OnOptionsEndEdit_();
      MPI_Barrier(Assembler.Domain().Comm());
    };
    OptionsEditor_.Activate(std::move(DeactivateFunc));
  }

  return OptionsEditor_.Edit(Options_);

}

void assembler::RestoreOptions() {

  if (OVK_DEBUG) {
    OVK_DEBUG_ASSERT(OptionsEditor_.Active(), "Unable to restore options for assembler %s; not "
      "currently being edited.", *Name_);
  }

  OptionsEditor_.Restore();

}

void assembler::OnOptionsStartEdit_() {

  CachedOptions_ = Options_;

}

void assembler::OnOptionsEndEdit_() {

  const domain &Domain = *Domain_;

  // No self-intersections, etc.
  for (int GridID : Domain.GridIDs()) {
    elem<int,2> SelfPair = {GridID,GridID};
    Options_.ResetOverlappable(SelfPair);
    Options_.ResetOverlapTolerance(SelfPair);
    Options_.ResetCutBoundaryHoles(SelfPair);
    Options_.ResetOccludes(SelfPair);
    Options_.ResetEdgePadding(SelfPair);
    Options_.ResetConnectionType(SelfPair);
    Options_.ResetMinimizeOverlap(SelfPair);
  }

  // TODO: Make this more fine-grained

  const set<int> &GridIDs = Domain.GridIDs();
  elem_set<int,2> GridIDPairs;

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      if (MGridID != NGridID) {
        GridIDPairs.Insert({MGridID,NGridID});
      }
    }
  }

  AssemblyManifest_.DetectOverlap = GridIDPairs;
  AssemblyManifest_.InferBoundaries = GridIDs;
  AssemblyManifest_.CutBoundaryHoles = GridIDPairs;
  AssemblyManifest_.ComputeOcclusion = GridIDPairs;
  AssemblyManifest_.ApplyPadding = GridIDPairs;
  AssemblyManifest_.ApplySmoothing = GridIDs;
  AssemblyManifest_.MinimizeOverlap = GridIDPairs;
  AssemblyManifest_.GenerateConnectivity = GridIDPairs;

  CachedOptions_ = options();

}

void assembler::OnGridEvent_(int GridID, grid_event_flags Flags, bool LastInSequence) {

  auto FlagsMatchAny = [&Flags](grid_event_flags Mask) -> bool {
    return (Flags & Mask) != grid_event_flags::NONE;
  };

  bool Create = FlagsMatchAny(grid_event_flags::CREATE);
  bool Destroy = FlagsMatchAny(grid_event_flags::DESTROY);

  if (Create) {
    UpdateManifest_.AddGridsToOptions.Insert(GridID);
  }

  if (Destroy) {
    UpdateManifest_.RemoveGridsFromOptions.Insert(GridID);
    UpdateManifest_.RemoveAssemblyManifestEntries.Insert(GridID);
  }

  if (LastInSequence) {
    Update_();
  }

}

void assembler::OnComponentEvent_(int ComponentID, component_event_flags Flags) {

  auto FlagsMatchAny = [&Flags](component_event_flags Mask) -> bool {
    return (Flags & Mask) != component_event_flags::NONE;
  };

  bool DestroyGeometryComponent = ComponentID == GeometryComponentID_ &&
    FlagsMatchAny(component_event_flags::DESTROY);
  bool DestroyStateComponent = ComponentID == StateComponentID_ &&
    FlagsMatchAny(component_event_flags::DESTROY);
  bool DestroyOverlapComponent = ComponentID == OverlapComponentID_ &&
    FlagsMatchAny(component_event_flags::DESTROY);
  bool DestroyConnectivityComponent = ComponentID == ConnectivityComponentID_ &&
    FlagsMatchAny(component_event_flags::DESTROY);

  if (DestroyGeometryComponent || DestroyStateComponent || DestroyOverlapComponent ||
    DestroyConnectivityComponent) {
    Unbind();
  }

}

void assembler::OnGeometryEvent_(int GridID, geometry_event_flags Flags, bool LastInSequence) {

  // TODO: Make this more fine-grained

  const domain &Domain = *Domain_;

  const set<int> &GridIDs = Domain.GridIDs();
  elem_set<int,2> GridIDPairs;

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      if (MGridID != NGridID) {
        GridIDPairs.Insert({MGridID,NGridID});
      }
    }
  }

  AssemblyManifest_.DetectOverlap = GridIDPairs;
  AssemblyManifest_.InferBoundaries = GridIDs;
  AssemblyManifest_.CutBoundaryHoles = GridIDPairs;
  AssemblyManifest_.ComputeOcclusion = GridIDPairs;
  AssemblyManifest_.ApplyPadding = GridIDPairs;
  AssemblyManifest_.ApplySmoothing = GridIDs;
  AssemblyManifest_.MinimizeOverlap = GridIDPairs;
  AssemblyManifest_.GenerateConnectivity = GridIDPairs;

}

void assembler::OnStateEvent_(int GridID, state_event_flags Flags, bool LastInSequence) {

  // TODO: Make this more fine-grained

  const domain &Domain = *Domain_;

  const set<int> &GridIDs = Domain.GridIDs();
  elem_set<int,2> GridIDPairs;

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      if (MGridID != NGridID) {
        GridIDPairs.Insert({MGridID,NGridID});
      }
    }
  }

  AssemblyManifest_.DetectOverlap = GridIDPairs;
  AssemblyManifest_.InferBoundaries = GridIDs;
  AssemblyManifest_.CutBoundaryHoles = GridIDPairs;
  AssemblyManifest_.ComputeOcclusion = GridIDPairs;
  AssemblyManifest_.ApplyPadding = GridIDPairs;
  AssemblyManifest_.ApplySmoothing = GridIDs;
  AssemblyManifest_.MinimizeOverlap = GridIDPairs;
  AssemblyManifest_.GenerateConnectivity = GridIDPairs;

}

void assembler::OnOverlapEvent_(const elem<int,2> &OverlapID, overlap_event_flags Flags, bool
  LastInSequence) {

  // TODO: Make this more fine-grained

  const domain &Domain = *Domain_;

  const set<int> &GridIDs = Domain.GridIDs();
  elem_set<int,2> GridIDPairs;

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      if (MGridID != NGridID) {
        GridIDPairs.Insert({MGridID,NGridID});
      }
    }
  }

  AssemblyManifest_.InferBoundaries = GridIDs;
  AssemblyManifest_.CutBoundaryHoles = GridIDPairs;
  AssemblyManifest_.ComputeOcclusion = GridIDPairs;
  AssemblyManifest_.ApplyPadding = GridIDPairs;
  AssemblyManifest_.ApplySmoothing = GridIDs;
  AssemblyManifest_.MinimizeOverlap = GridIDPairs;
  AssemblyManifest_.GenerateConnectivity = GridIDPairs;

}

void assembler::OnConnectivityEvent_(const elem<int,2> &ConnectivityID, connectivity_event_flags
  Flags, bool LastInSequence) {

  // Nothing to do

}

void assembler::Update_() {

  const domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  Logger.LogStatus(Domain.Comm().Rank() == 0, "Updating assembler %s...", *Name_);
  auto Level1 = Logger.IncreaseStatusLevelAndIndent();

  AddGridsToOptions_();
  RemoveGridsFromOptions_();
  RemoveAssemblyManifestEntries_();

  UpdateManifest_.AddGridsToOptions.Clear();
  UpdateManifest_.RemoveGridsFromOptions.Clear();
  UpdateManifest_.RemoveAssemblyManifestEntries.Clear();

  Level1.Reset();
  Logger.LogStatus(Domain.Comm().Rank() == 0, "Done updating assembler %s.", *Name_);

}

void assembler::AddGridsToOptions_() {

  if (UpdateManifest_.AddGridsToOptions.Empty()) return;

  Options_.AddGrids(UpdateManifest_.AddGridsToOptions);

}

void assembler::RemoveGridsFromOptions_() {

  if (UpdateManifest_.RemoveGridsFromOptions.Empty()) return;

  Options_.RemoveGrids(UpdateManifest_.RemoveGridsFromOptions);

}


void assembler::RemoveAssemblyManifestEntries_() {

  if (UpdateManifest_.RemoveAssemblyManifestEntries.Empty()) return;

  auto MatchesGridToRemove = [&](int GridID) -> bool {
    return UpdateManifest_.RemoveAssemblyManifestEntries.Contains(GridID);
  };

  auto MatchesGridToRemovePair = [&](const elem<int,2> &IDPair) -> bool {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    return UpdateManifest_.RemoveAssemblyManifestEntries.Contains(MGridID) ||
      UpdateManifest_.RemoveAssemblyManifestEntries.Contains(NGridID);
  };

  AssemblyManifest_.DetectOverlap.EraseIf(MatchesGridToRemovePair);
  AssemblyManifest_.InferBoundaries.EraseIf(MatchesGridToRemove);
  AssemblyManifest_.CutBoundaryHoles.EraseIf(MatchesGridToRemovePair);
  AssemblyManifest_.ComputeOcclusion.EraseIf(MatchesGridToRemovePair);
  AssemblyManifest_.ApplyPadding.EraseIf(MatchesGridToRemovePair);
  AssemblyManifest_.ApplySmoothing.EraseIf(MatchesGridToRemove);
  AssemblyManifest_.MinimizeOverlap.EraseIf(MatchesGridToRemovePair);
  AssemblyManifest_.GenerateConnectivity.EraseIf(MatchesGridToRemovePair);

}

assembler::params &assembler::params::SetName(std::string Name) {

  Name_ = std::move(Name);

  return *this;

}

assembler::bindings &assembler::bindings::SetGeometryComponentID(int GeometryComponentID) {

  OVK_DEBUG_ASSERT(GeometryComponentID >= 0, "Invalid geometry component ID.");

  GeometryComponentID_ = GeometryComponentID;

  return *this;

}

assembler::bindings &assembler::bindings::SetStateComponentID(int StateComponentID) {

  OVK_DEBUG_ASSERT(StateComponentID >= 0, "Invalid state component ID.");

  StateComponentID_ = StateComponentID;

  return *this;

}

assembler::bindings &assembler::bindings::SetOverlapComponentID(int OverlapComponentID) {

  OVK_DEBUG_ASSERT(OverlapComponentID >= 0, "Invalid overlap component ID.");

  OverlapComponentID_ = OverlapComponentID;

  return *this;

}

assembler::bindings &assembler::bindings::SetConnectivityComponentID(int ConnectivityComponentID) {

  OVK_DEBUG_ASSERT(ConnectivityComponentID >= 0, "Invalid connectivity component ID.");

  ConnectivityComponentID_ = ConnectivityComponentID;

  return *this;

}

}
