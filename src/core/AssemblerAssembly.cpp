// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Assembler.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayOps.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/ConnectivityComponent.hpp"
#include "ovk/core/ConnectivityM.hpp"
#include "ovk/core/ConnectivityN.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/DistributedField.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/Event.hpp"
#include "ovk/core/Field.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Geometry.hpp"
#include "ovk/core/GeometryComponent.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/OverlapComponent.hpp"
#include "ovk/core/OverlapM.hpp"
#include "ovk/core/OverlapN.hpp"
#include "ovk/core/Partition.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Set.hpp"
#include "ovk/core/State.hpp"
#include "ovk/core/StateComponent.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <memory>
#include <string>
#include <utility>

namespace ovk {

namespace {

void GenerateActiveMask(const grid &Grid, const distributed_field<state_flags> &Flags,
  distributed_field<bool> &ActiveMask);
void GenerateCellActiveMask(const grid &Grid, const distributed_field<state_flags> &Flags,
  distributed_field<bool> &CellActiveMask);
void GenerateDomainBoundaryMask(const grid &Grid, const distributed_field<state_flags> &Flags,
  distributed_field<bool> &DomainBoundaryMask);
void GenerateInternalBoundaryMask(const grid &Grid, const distributed_field<state_flags> &Flags,
  distributed_field<bool> &InternalBoundaryMask);

}

void assembler::Assemble() {

  const domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Assembling domain %s using assembler %s...",
    Domain.Name(), *Name_);

  InitializeAssembly_();

  AssemblyManifest_.DetectOverlap.Clear();
  AssemblyManifest_.InferBoundaries.Clear();
  AssemblyManifest_.CutBoundaryHoles.Clear();
  AssemblyManifest_.ComputeOcclusion.Clear();
  AssemblyManifest_.ApplyPadding.Clear();
  AssemblyManifest_.ApplySmoothing.Clear();
  AssemblyManifest_.MinimizeOverlap.Clear();
  AssemblyManifest_.GenerateConnectivity.Clear();

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done assembling domain %s using assembler %s.",
    Domain.Name(), *Name_);

}

void assembler::InitializeAssembly_() {

  const domain &Domain = *Domain_;
  int NumDims = Domain.Dimension();
  auto &GeometryComponent = Domain.Component<geometry_component>(GeometryComponentID_);
  auto &StateComponent = Domain.Component<state_component>(StateComponentID_);
  assembly_data &AssemblyData = *AssemblyData_;

  for (int GridID : Domain.GridIDs()) {
    OVK_DEBUG_ASSERT(GeometryComponent.GeometryExists(GridID), "No geometry data for grid %s.",
      Domain.GridInfo(GridID).Name());
    OVK_DEBUG_ASSERT(StateComponent.StateExists(GridID), "No state data for grid %s.",
      Domain.GridInfo(GridID).Name());
  }

  if (OVK_DEBUG) {
    ValidateOptions_();
  }

  range VertexOffsetRange = MakeEmptyRange(NumDims);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    VertexOffsetRange.Begin(iDim) = 0;
    VertexOffsetRange.End(iDim) = 2;
  }

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    auto &Flags = StateComponent.State(GridID).Flags();
    core::partition_pool PartitionPool(Context_, Grid.Comm(), Grid.Partition().NeighborRanks());
    PartitionPool.Insert(Grid.SharedPartition());
    PartitionPool.Insert(Grid.SharedCellPartition());
    local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData.Insert(GridID,
      std::move(PartitionPool));
    GenerateActiveMask(Grid, Flags, GridAuxData.ActiveMask);
    GenerateCellActiveMask(Grid, Flags, GridAuxData.CellActiveMask);
    GenerateDomainBoundaryMask(Grid, Flags, GridAuxData.DomainBoundaryMask);
    GenerateInternalBoundaryMask(Grid, Flags, GridAuxData.InternalBoundaryMask);
  }

}

void assembler::ValidateOptions_() {

  domain &Domain = *Domain_;

  for (int MGridID : Domain.GridIDs()) {
    for (int NGridID : Domain.GridIDs()) {
      const std::string &MGridName = Domain.GridInfo(MGridID).Name();
      const std::string &NGridName = Domain.GridInfo(NGridID).Name();
      if (Options_.CutBoundaryHoles({MGridID,NGridID})) {
        OVK_DEBUG_ASSERT(Options_.Overlappable({MGridID,NGridID}), "Grid %s being boundary-hole-"
          "cut by grid %s requires %s to be overlappable by %s.", NGridName, MGridName, NGridName,
          MGridName);
        OVK_DEBUG_ASSERT(Options_.Overlappable({NGridID,MGridID}), "Grid %s being boundary-hole-"
          "cut by grid %s requires %s to be overlappable by %s.", NGridName, MGridName, MGridName,
          NGridName);
      }
    }
  }

  // This is incomplete; add rest

}

namespace {

void GenerateActiveMask(const grid &Grid, const distributed_field<state_flags> &Flags,
  distributed_field<bool> &ActiveMask) {

  long long NumExtended = Grid.ExtendedRange().Count();

  ActiveMask.Assign(Grid.SharedPartition());

  for (long long l = 0; l < NumExtended; ++l) {
    ActiveMask[l] = (Flags[l] & state_flags::ACTIVE) != state_flags::NONE;
  }

}

void GenerateCellActiveMask(const grid &Grid, const distributed_field<state_flags> &Flags,
  distributed_field<bool> &CellActiveMask) {

  int NumDims = Grid.Dimension();
  const range &CellLocalRange = Grid.CellLocalRange();

  CellActiveMask.Assign(Grid.SharedCellPartition());

  for (int k = CellLocalRange.Begin(2); k < CellLocalRange.End(2); ++k) {
    for (int j = CellLocalRange.Begin(1); j < CellLocalRange.End(1); ++j) {
      for (int i = CellLocalRange.Begin(0); i < CellLocalRange.End(0); ++i) {
        tuple<int> Cell = {i,j,k};
        CellActiveMask(Cell) = true;
        range NeighborRange;
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          NeighborRange.Begin(iDim) = Cell(iDim);
          NeighborRange.End(iDim) = Cell(iDim)+2;
        }
        for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
          NeighborRange.Begin(iDim) = 0;
          NeighborRange.End(iDim) = 1;
        }
        for (int o = NeighborRange.Begin(2); o < NeighborRange.End(2); ++o) {
          for (int n = NeighborRange.Begin(1); n < NeighborRange.End(1); ++n) {
            for (int m = NeighborRange.Begin(0); m < NeighborRange.End(0); ++m) {
              tuple<int> Point = {m,n,o};
              CellActiveMask(Cell) = CellActiveMask(Cell) && (Flags(Point) & state_flags::ACTIVE)
                != state_flags::NONE;
            }
          }
        }
      }
    }
  }

  CellActiveMask.Exchange();

}

void GenerateDomainBoundaryMask(const grid &Grid, const distributed_field<state_flags> &Flags,
  distributed_field<bool> &DomainBoundaryMask) {

  long long NumExtended = Grid.ExtendedRange().Count();

  DomainBoundaryMask.Assign(Grid.SharedPartition());

  auto MatchesAll = [](state_flags Flags, state_flags Mask) -> bool {
    return (Flags & Mask) == Mask;
  };

  for (long long l = 0; l < NumExtended; ++l) {
    DomainBoundaryMask[l] = MatchesAll(Flags[l], state_flags::ACTIVE |
      state_flags::DOMAIN_BOUNDARY);
  }

}

void GenerateInternalBoundaryMask(const grid &Grid, const distributed_field<state_flags> &Flags,
  distributed_field<bool> &InternalBoundaryMask) {

  long long NumExtended = Grid.ExtendedRange().Count();

  InternalBoundaryMask.Assign(Grid.SharedPartition());

  auto MatchesAll = [](state_flags Flags, state_flags Mask) -> bool {
    return (Flags & Mask) == Mask;
  };

  for (long long l = 0; l < NumExtended; ++l) {
    InternalBoundaryMask[l] = MatchesAll(Flags[l], state_flags::ACTIVE |
      state_flags::INTERNAL_BOUNDARY);
  }

}

}

}
