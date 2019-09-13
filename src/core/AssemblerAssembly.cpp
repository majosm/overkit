// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Assembler.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayOps.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Collect.hpp"
#include "ovk/core/CollectMap.hpp"
#include "ovk/core/ConnectivityComponent.hpp"
#include "ovk/core/ConnectivityM.hpp"
#include "ovk/core/ConnectivityN.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Disperse.hpp"
#include "ovk/core/DisperseMap.hpp"
#include "ovk/core/DistributedField.hpp"
#include "ovk/core/DistributedFieldOps.hpp"
#include "ovk/core/DistributedRegionHash.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/ElemMap.hpp"
#include "ovk/core/ElemSet.hpp"
#include "ovk/core/Event.hpp"
#include "ovk/core/Field.hpp"
#include "ovk/core/FieldOps.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Geometry.hpp"
#include "ovk/core/GeometryComponent.hpp"
#include "ovk/core/GeometryOps.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Math.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/OverlapComponent.hpp"
#include "ovk/core/OverlapM.hpp"
#include "ovk/core/OverlapN.hpp"
#include "ovk/core/Partition.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Recv.hpp"
#include "ovk/core/RecvMap.hpp"
#include "ovk/core/ScalarOps.hpp"
#include "ovk/core/Send.hpp"
#include "ovk/core/SendMap.hpp"
#include "ovk/core/Set.hpp"
#include "ovk/core/State.hpp"
#include "ovk/core/StateComponent.hpp"
#include "ovk/core/TextProcessing.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>

namespace ovk {

namespace {

constexpr long long NO_CELL = std::numeric_limits<long long>::min();

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
  DetectOverlap_();
  InferBoundaries_();
  CutBoundaryHoles_();
  LocateOuterFringe_();
  DetectOccluded_();
  MinimizeOverlap_();
  GenerateConnectivityData_();

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

assembler::assembly_data::assembly_data(int NumDims, comm_view Comm):
  BoundingBoxHash(NumDims, Comm)
{}

void assembler::InitializeAssembly_() {

  domain &Domain = *Domain_;
  auto &GeometryComponent = Domain.Component<geometry_component>(GeometryComponentID_);
  assembly_data &AssemblyData = *AssemblyData_;

  auto StateComponentEditHandle = Domain.EditComponent<state_component>(StateComponentID_);
  state_component &StateComponent = *StateComponentEditHandle;

  for (int GridID : Domain.GridIDs()) {
    OVK_DEBUG_ASSERT(GeometryComponent.GeometryExists(GridID), "No geometry data for grid %s.",
      Domain.GridInfo(GridID).Name());
    OVK_DEBUG_ASSERT(StateComponent.StateExists(GridID), "No state data for grid %s.",
      Domain.GridInfo(GridID).Name());
  }

  if (OVK_DEBUG) {
    ValidateOptions_();
  }

  constexpr state_flags AllAssemblyFlags =
    state_flags::OVERLAPPED |
    state_flags::INFERRED_DOMAIN_BOUNDARY |
    state_flags::BOUNDARY_HOLE |
    state_flags::OCCLUDED |
    state_flags::FRINGE |
    state_flags::OUTER_FRINGE |
    state_flags::INNER_FRINGE |
    state_flags::OVERLAP_MINIMIZED |
    state_flags::RECEIVER |
    state_flags::ORPHAN;

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    {
      auto StateEditHandle = StateComponent.EditState(GridID);
      auto FlagsEditHandle = StateEditHandle->EditFlags();
      distributed_field<state_flags> &Flags = *FlagsEditHandle;
      for (long long l = 0; l < NumExtended; ++l) {
        if ((Flags[l] & (state_flags::BOUNDARY_HOLE | state_flags::OVERLAP_MINIMIZED)) !=
          state_flags::NONE) {
          Flags[l] = Flags[l] | state_flags::ACTIVE;
        }
        if ((Flags[l] & state_flags::INFERRED_DOMAIN_BOUNDARY) != state_flags::NONE) {
          Flags[l] = Flags[l] & ~state_flags::DOMAIN_BOUNDARY;
        }
        Flags[l] = Flags[l] & ~AllAssemblyFlags;
      }
    }
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

void assembler::DetectOverlap_() {

  domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Detecting overlap between grids...");

  int NumDims = Domain.Dimension();
  auto &GeometryComponent = Domain.Component<geometry_component>(GeometryComponentID_);
  assembly_data &AssemblyData = *AssemblyData_;

  Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Generating distributed bounding box hash...");

  range VertexOffsetRange = MakeEmptyRange(NumDims);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    VertexOffsetRange.Begin(iDim) = 0;
    VertexOffsetRange.End(iDim) = 2;
  }

  // Range consisting of all cells having vertices in local range
  auto MakeCellCoverRange = [](const cart &CellCart, const range &CellLocalRange) -> range {
    range CellCoverRange = MakeEmptyRange(CellCart.Dimension());
    for (int iDim = 0; iDim < CellCart.Dimension(); ++iDim) {
      if (CellLocalRange.Begin(iDim) > CellCart.Range().Begin(iDim) || (CellCart.Periodic(iDim) &&
        CellLocalRange.End(iDim) != CellCart.Range().End(iDim))) {
        CellCoverRange.Begin(iDim) = CellLocalRange.Begin(iDim)-1;
      } else {
        CellCoverRange.Begin(iDim) = CellLocalRange.Begin(iDim);
      }
      CellCoverRange.End(iDim) = CellLocalRange.End(iDim);
    }
    return CellCoverRange;
  };

  array<box> LocalGridBounds;
  LocalGridBounds.Reserve(Domain.LocalGridCount());

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    range CellCoverRange = MakeCellCoverRange(Grid.CellCart(), Grid.CellLocalRange());
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(GridID);
    const distributed_field<bool> &CellActiveMask = GridAuxData.CellActiveMask;
    const geometry &Geometry = GeometryComponent.Geometry(GridID);
    auto &Coords = Geometry.Coords();
    box &Bounds = LocalGridBounds.Append();
    Bounds = MakeEmptyBox(NumDims);
    for (int k = CellCoverRange.Begin(2); k < CellCoverRange.End(2); ++k) {
      for (int j = CellCoverRange.Begin(1); j < CellCoverRange.End(1); ++j) {
        for (int i = CellCoverRange.Begin(0); i < CellCoverRange.End(0); ++i) {
          tuple<int> Cell = {i,j,k};
          if (!CellActiveMask(Cell)) continue;
          for (int o = VertexOffsetRange.Begin(2); o < VertexOffsetRange.End(2); ++o) {
            for (int n = VertexOffsetRange.Begin(1); n < VertexOffsetRange.End(1); ++n) {
              for (int m = VertexOffsetRange.Begin(0); m < VertexOffsetRange.End(0); ++m) {
                tuple<int> Vertex = {Cell(0)+m,Cell(1)+n,Cell(2)+o};
                tuple<double> VertexCoords = {
                  Coords(0)(Vertex),
                  Coords(1)(Vertex),
                  Coords(2)(Vertex)
                };
                Bounds = ExtendBox(Bounds, VertexCoords);
              }
            }
          }
        }
      }
    }
  }

  array<int> LocalGridIDs(Domain.LocalGridIDs());

  bounding_box_hash &BoundingBoxHash = AssemblyData.BoundingBoxHash;
  BoundingBoxHash = bounding_box_hash(NumDims, Domain.Comm(), Domain.LocalGridCount(),
    LocalGridBounds, LocalGridIDs);

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done generating distributed bounding box hash.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Mapping local coordinates into hash bins...");
  }

  map<int,field<int>> LocalPointOverlappingBinIndices;
  set<int> UniqueOverlappingBinIndices;

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    const range &LocalRange = Grid.LocalRange();
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(GridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    const geometry &Geometry = GeometryComponent.Geometry(GridID);
    auto &Coords = Geometry.Coords();
    field<int> &BinIndices = LocalPointOverlappingBinIndices.Insert(GridID);
    BinIndices.Resize(LocalRange, -1);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          if (!ActiveMask(Point)) continue;
          tuple<double> PointCoords = {
            Coords(0)(Point),
            Coords(1)(Point),
            Coords(2)(Point)
          };
          tuple<int> BinLoc = BoundingBoxHash.MapPointToBin(PointCoords);
          int BinIndex = BoundingBoxHash.BinIndexer().ToIndex(BinLoc);
          BinIndices(Point) = BinIndex;
          UniqueOverlappingBinIndices.Insert(BinIndex);
        }
      }
    }
  }

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done mapping local coordinates into hash bins.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Retrieving remote hash bins...");
  }

  map<int,bounding_box_hash_bin> Bins;

  for (int BinIndex : UniqueOverlappingBinIndices) {
    Bins.Insert(BinIndex);
  }

  BoundingBoxHash.RetrieveBins(Bins);

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done retrieving remote hash bins.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Establishing communication between potentially-"
      "overlapping ranks...");
  }

  map<int,map<int,set<int>>> OverlappingMGridIDsAndRanksForLocalNGrid;

  for (int NGridID : Domain.LocalGridIDs()) {
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    const geometry &Geometry = GeometryComponent.Geometry(NGridID);
    auto &Coords = Geometry.Coords();
    field<int> &BinIndices = LocalPointOverlappingBinIndices(NGridID);
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid.Insert(NGridID);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          int BinIndex = BinIndices(Point);
          if (BinIndex < 0) continue;
          bounding_box_hash_bin &Bin = Bins(BinIndex);
          tuple<double> PointCoords = {
            Coords(0)(Point),
            Coords(1)(Point),
            Coords(2)(Point)
          };
          for (auto &Region : Bin.Regions()) {
            int MGridID = Region.Tag;
            if (Options_.Overlappable({MGridID,NGridID}) && Region.Extents.Contains(PointCoords)) {
              MGridIDsAndRanks.Fetch(MGridID).Insert(Region.Rank);
            }
          }
        }
      }
    }
  }

  set<int> RemoteMRanksSet;

  for (auto &NEntry : OverlappingMGridIDsAndRanksForLocalNGrid) {
    auto &MGridIDsAndRanks = NEntry.Value();
    for (auto &MEntry : MGridIDsAndRanks) {
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        if (Rank != Domain.Comm().Rank()) {
          RemoteMRanksSet.Insert(Rank);
        }
      }
    }
  }

  array<int> RemoteMRanks(RemoteMRanksSet);

  array<int> RemoteNRanks = core::DynamicHandshake(Domain.Comm(), RemoteMRanks);

  array<MPI_Request> MPIRequests;
  MPIRequests.Reserve(RemoteMRanks.Count() + RemoteNRanks.Count());

  map<int,int> NumGridIDPairsFromRank;
  NumGridIDPairsFromRank.Reserve(RemoteNRanks.Count());

  for (int Rank : RemoteNRanks) {
    int &NumPairs = NumGridIDPairsFromRank.Insert(Rank, 0);
    MPI_Irecv(&NumPairs, 1, MPI_INT, Rank, 0, Domain.Comm(), &MPIRequests.Append());
  }

  map<int,int> NumGridIDPairsToRank;
  NumGridIDPairsToRank.Reserve(RemoteMRanks.Count());

  for (auto &NEntry : OverlappingMGridIDsAndRanksForLocalNGrid) {
    auto &MGridIDsAndRanks = NEntry.Value();
    for (auto &MEntry : MGridIDsAndRanks) {
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        ++NumGridIDPairsToRank.Fetch(Rank, 0);
      }
    }
  }

  for (int Rank : RemoteMRanks) {
    int &NumPairs = NumGridIDPairsToRank(Rank);
    MPI_Isend(&NumPairs, 1, MPI_INT, Rank, 0, Domain.Comm(), &MPIRequests.Append());
  }

  MPI_Waitall(MPIRequests.Count(), MPIRequests.Data(), MPI_STATUSES_IGNORE);

  MPIRequests.Clear();

  map<int,array<int>> GridIDPairRecvData;
  GridIDPairRecvData.Reserve(RemoteNRanks.Count());

  for (int Rank : RemoteNRanks) {
    int NumPairs = NumGridIDPairsFromRank(Rank);
    array<int> &Pairs = GridIDPairRecvData.Insert(Rank);
    Pairs.Resize({2*NumPairs});
    MPI_Irecv(Pairs.Data(), 2*NumPairs, MPI_INT, Rank, 0, Domain.Comm(), &MPIRequests.Append());
  }

  map<int,array<int>> GridIDPairSendData;
  GridIDPairSendData.Reserve(RemoteMRanks.Count());

  for (int Rank : RemoteMRanks) {
    int NumPairs = NumGridIDPairsToRank(Rank);
    array<int> &Pairs = GridIDPairSendData.Insert(Rank);
    Pairs.Reserve(2*NumPairs);
  }

  for (auto &NEntry : OverlappingMGridIDsAndRanksForLocalNGrid) {
    int NGridID = NEntry.Key();
    auto &MGridIDsAndRanks = NEntry.Value();
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        if (Rank != Domain.Comm().Rank()) {
          array<int> &Pairs = GridIDPairSendData(Rank);
          Pairs.Append(MGridID);
          Pairs.Append(NGridID);
        }
      }
    }
  }

  for (int Rank : RemoteMRanks) {
    int NumPairs = NumGridIDPairsToRank(Rank);
    array<int> &Pairs = GridIDPairSendData(Rank);
    MPI_Isend(Pairs.Data(), 2*NumPairs, MPI_INT, Rank, 0, Domain.Comm(), &MPIRequests.Append());
  }

  MPI_Waitall(MPIRequests.Count(), MPIRequests.Data(), MPI_STATUSES_IGNORE);

  MPIRequests.Clear();

  map<int,map<int,set<int>>> OverlappingNGridIDsAndRanksForLocalMGrid;

  for (int MGridID : Domain.LocalGridIDs()) {
    OverlappingNGridIDsAndRanksForLocalMGrid.Insert(MGridID);
  }

  for (auto &NEntry : OverlappingMGridIDsAndRanksForLocalNGrid) {
    int NGridID = NEntry.Key();
    auto &MGridIDsAndRanks = NEntry.Value();
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        if (Rank == Domain.Comm().Rank()) {
          OverlappingNGridIDsAndRanksForLocalMGrid(MGridID).Fetch(NGridID).Insert(Rank);
        }
      }
    }
  }

  for (int Rank : RemoteNRanks) {
    int NumPairs = NumGridIDPairsFromRank(Rank);
    array<int> &Pairs = GridIDPairRecvData(Rank);
    for (int iPair = 0; iPair < NumPairs; ++iPair) {
      int MGridID = Pairs(2*iPair);
      int NGridID = Pairs(2*iPair+1);
      OverlappingNGridIDsAndRanksForLocalMGrid(MGridID).Fetch(NGridID).Insert(Rank);
    }
  }

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done establishing communication between "
      "potentially-overlapping ranks.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Transferring coordinate data...");
  }

  elem_set<int,2> MGridDataSends;

  for (int MGridID : Domain.LocalGridIDs()) {
    auto &NGridIDsAndRanks = OverlappingNGridIDsAndRanksForLocalMGrid(MGridID);
    for (auto &NEntry : NGridIDsAndRanks) {
      const set<int> &NGridRanks = NEntry.Value();
      for (int Rank : NGridRanks) {
        if (Rank != Domain.Comm().Rank()) {
          MGridDataSends.Insert({MGridID,Rank});
        }
      }
    }
  }

  elem_set<int,2> MGridDataRecvs;

  for (int NGridID : Domain.LocalGridIDs()) {
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        if (Rank != Domain.Comm().Rank()) {
          MGridDataRecvs.Insert({MGridID,Rank});
        }
      }
    }
  }

  struct partition_data {
    range ExtendedRange;
    range CellLocalRange;
    range CellExtendedRange;
    range CellCoverRange;
  };

  elem_map<int,2,partition_data> MGridPartitionData;

  for (int NGridID : Domain.LocalGridIDs()) {
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        MGridPartitionData.Insert({MGridID,Rank});
      }
    }
  }

  MPIRequests.Reserve(6*(MGridDataSends.Count() + MGridDataRecvs.Count()));

  for (auto &MGridIDAndRankPair : MGridDataRecvs) {
    int Rank = MGridIDAndRankPair(1);
    partition_data &Data = MGridPartitionData(MGridIDAndRankPair);
    MPI_Irecv(Data.ExtendedRange.Begin().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Irecv(Data.ExtendedRange.End().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Irecv(Data.CellLocalRange.Begin().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Irecv(Data.CellLocalRange.End().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Irecv(Data.CellExtendedRange.Begin().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Irecv(Data.CellExtendedRange.End().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
  }

  for (auto &MGridIDAndSendToRankPair : MGridDataSends) {
    int MGridID = MGridIDAndSendToRankPair(0);
    int Rank = MGridIDAndSendToRankPair(1);
    const grid &MGrid = Domain.Grid(MGridID);
    MPI_Isend(MGrid.ExtendedRange().Begin().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Isend(MGrid.ExtendedRange().End().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Isend(MGrid.CellLocalRange().Begin().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Isend(MGrid.CellLocalRange().End().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Isend(MGrid.CellExtendedRange().Begin().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Isend(MGrid.CellExtendedRange().End().Data(), MAX_DIMS, MPI_INT, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
  }

  for (int NGridID : Domain.LocalGridIDs()) {
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        if (Rank == Domain.Comm().Rank()) {
          const grid &MGrid = Domain.Grid(MGridID);
          partition_data &Data = MGridPartitionData({MGridID,Rank});
          Data.ExtendedRange = MGrid.ExtendedRange();
          Data.CellLocalRange = MGrid.CellLocalRange();
          Data.CellExtendedRange = MGrid.CellExtendedRange();
        }
      }
    }
  }

  MPI_Waitall(MPIRequests.Count(), MPIRequests.Data(), MPI_STATUSES_IGNORE);

  MPIRequests.Clear();

  for (int NGridID : Domain.LocalGridIDs()) {
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const cart &CellCart = MGridInfo.CellCart();
      for (int Rank : MGridRanks) {
        partition_data &Data = MGridPartitionData({MGridID,Rank});
        Data.CellCoverRange = MakeCellCoverRange(CellCart, Data.CellLocalRange);
      }
    }
  }

  struct cell_coord_data {
    array<field<double>> Coords;
    geometry_type GeometryType;
    field<bool> CellActiveMask;
  };

  elem_map<int,2,cell_coord_data> CellCoordData;

  for (int NGridID : Domain.LocalGridIDs()) {
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        const partition_data &PartitionData = MGridPartitionData({MGridID,Rank});
        cell_coord_data &CoordData = CellCoordData.Insert({MGridID,Rank});
        CoordData.Coords.Resize({MAX_DIMS});
        CoordData.Coords(0).Resize(PartitionData.ExtendedRange);
        CoordData.Coords(1).Resize(PartitionData.ExtendedRange);
        CoordData.Coords(2).Resize(PartitionData.ExtendedRange);
        CoordData.CellActiveMask.Resize(PartitionData.CellExtendedRange);
      }
    }
  }

  MPIRequests.Reserve(5*(MGridDataSends.Count() + MGridDataRecvs.Count()));

  for (auto &MGridIDAndRankPair : MGridDataRecvs) {
    int Rank = MGridIDAndRankPair(1);
    const partition_data &PartitionData = MGridPartitionData(MGridIDAndRankPair);
    cell_coord_data &CoordData = CellCoordData(MGridIDAndRankPair);
    long long NumExtended = PartitionData.ExtendedRange.Count();
    long long NumCellExtended = PartitionData.CellExtendedRange.Count();
    MPI_Irecv(CoordData.Coords(0).Data(), NumExtended, MPI_DOUBLE, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Irecv(CoordData.Coords(1).Data(), NumExtended, MPI_DOUBLE, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Irecv(CoordData.Coords(2).Data(), NumExtended, MPI_DOUBLE, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Irecv(&CoordData.GeometryType, 1, core::GetMPIDataType<geometry_type>(), Rank, 0,
      Domain.Comm(), &MPIRequests.Append());
    MPI_Irecv(CoordData.CellActiveMask.Data(), NumCellExtended, MPI_C_BOOL, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
  }

  // Need to store geometry type values somewhere for Isend buffer
  map<int,geometry_type> GeometryTypeStorage;

  for (int MGridID : Domain.LocalGridIDs()) {
    const geometry &Geometry = GeometryComponent.Geometry(MGridID);
    GeometryTypeStorage.Insert(MGridID, Geometry.Type());
  }

  for (auto &MGridIDAndSendToRankPair : MGridDataSends) {
    int MGridID = MGridIDAndSendToRankPair(0);
    int Rank = MGridIDAndSendToRankPair(1);
    const grid &MGrid = Domain.Grid(MGridID);
    const geometry &Geometry = GeometryComponent.Geometry(MGridID);
    const array<distributed_field<double>> &Coords = Geometry.Coords();
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(MGridID);
    const distributed_field<bool> &CellActiveMask = GridAuxData.CellActiveMask;
    long long NumExtended = MGrid.ExtendedRange().Count();
    long long NumCellExtended = MGrid.CellExtendedRange().Count();
    MPI_Isend(Coords(0).Data(), NumExtended, MPI_DOUBLE, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Isend(Coords(1).Data(), NumExtended, MPI_DOUBLE, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Isend(Coords(2).Data(), NumExtended, MPI_DOUBLE, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
    MPI_Isend(&GeometryTypeStorage(MGridID), 1, core::GetMPIDataType<geometry_type>(), Rank, 0,
      Domain.Comm(), &MPIRequests.Append());
    MPI_Isend(CellActiveMask.Data(), NumCellExtended, MPI_C_BOOL, Rank, 0, Domain.Comm(),
      &MPIRequests.Append());
  }

  for (int NGridID : Domain.LocalGridIDs()) {
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        if (Rank == Domain.Comm().Rank()) {
          const geometry &Geometry = GeometryComponent.Geometry(MGridID);
          const array<distributed_field<double>> &Coords = Geometry.Coords();
          const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(MGridID);
          const distributed_field<bool> &CellActiveMask = GridAuxData.CellActiveMask;
          cell_coord_data &Data = CellCoordData({MGridID,Rank});
          Data.GeometryType = Geometry.Type();
          Data.Coords(0).Fill(Coords(0));
          Data.Coords(1).Fill(Coords(1));
          Data.Coords(2).Fill(Coords(2));
          Data.CellActiveMask.Fill(CellActiveMask);
        }
      }
    }
  }

  MPI_Waitall(MPIRequests.Count(), MPIRequests.Data(), MPI_STATUSES_IGNORE);

  MPIRequests.Clear();

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done transferring coordinate data.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Searching for overlapping cells...");
  }

  // Brute force for now
  auto FindOverlappingCell = [NumDims](const range &CellRange, const array<field<double>> &Coords,
    geometry_type GeometryType, const field<bool> &CellActiveMask, double Tolerance, const
    tuple<double> &PointCoords) -> optional<tuple<int>> {
    for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
      for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
        for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
          tuple<int> Cell = {i,j,k};
          if (!CellActiveMask(Cell)) continue;
          if (core::OverlapsCell(NumDims, Coords, GeometryType, Tolerance, Cell, PointCoords)) {
            return Cell;
          }
        }
      }
    }
    return {};
  };

  struct overlapping_cell_data {
    bool Allocated = false;
    field_indexer Indexer;
    field<long long> Cells;
    overlapping_cell_data() = default;
    overlapping_cell_data(const range &MGridCellGlobalRange, const range &LocalRange):
      Allocated(true),
      Indexer(MGridCellGlobalRange),
      Cells(LocalRange, NO_CELL)
    {}
  };

  elem_map<int,2,overlapping_cell_data> OverlappingCellData;

  for (int NGridID : Domain.LocalGridIDs()) {
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(NGridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    const geometry &Geometry = GeometryComponent.Geometry(NGridID);
    auto &Coords = Geometry.Coords();
    field<int> &BinIndices = LocalPointOverlappingBinIndices(NGridID);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          if (!ActiveMask(Point)) continue;
          int BinIndex = BinIndices(Point);
          if (BinIndex < 0) continue;
          bounding_box_hash_bin &Bin = Bins(BinIndex);
          tuple<double> PointCoords = {
            Coords(0)(Point),
            Coords(1)(Point),
            Coords(2)(Point)
          };
          for (auto &Region : Bin.Regions()) {
            int MGridID = Region.Tag;
            if (!Region.Extents.Contains(PointCoords) || !Options_.Overlappable({MGridID,NGridID}))
              continue;
            elem<int,2> IDPair = {MGridID,NGridID};
            const partition_data &PartitionData = MGridPartitionData({MGridID,Region.Rank});
            const cell_coord_data &CoordData = CellCoordData({MGridID,Region.Rank});
            auto MaybeCell = FindOverlappingCell(PartitionData.CellLocalRange, CoordData.Coords,
              CoordData.GeometryType, CoordData.CellActiveMask, Options_.OverlapTolerance(IDPair),
              PointCoords);
            if (MaybeCell) {
              const tuple<int> &Cell = *MaybeCell;
              overlapping_cell_data &CellData = OverlappingCellData.Fetch(IDPair);
              if (!CellData.Allocated) {
                const range &MGridCellGlobalRange = Domain.GridInfo(MGridID).CellGlobalRange();
                CellData = overlapping_cell_data(MGridCellGlobalRange, LocalRange);
              }
              CellData.Cells(Point) = CellData.Indexer.ToIndex(Cell);
            }
          }
        }
      }
    }
  }

  map<int,elem_map<int,2,long long>> NumOverlappingFromMGridAndRankForLocalNGrid;
  map<int,elem_map<int,2,long long>> NumOverlappingFromNGridAndRankForLocalMGrid;

  int NumSends = 0;
  for (int NGridID : Domain.LocalGridIDs()) {
    auto &NumFromMGridAndRank = NumOverlappingFromMGridAndRankForLocalNGrid.Insert(NGridID);
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        NumFromMGridAndRank.Insert({MGridID,Rank}, 0);
        ++NumSends;
      }
    }
  }

  int NumRecvs = 0;
  for (int MGridID : Domain.LocalGridIDs()) {
    auto &NumFromNGridAndRank = NumOverlappingFromNGridAndRankForLocalMGrid.Insert(MGridID);
    auto &NGridIDsAndRanks = OverlappingNGridIDsAndRanksForLocalMGrid(MGridID);
    for (auto &NEntry : NGridIDsAndRanks) {
      int NGridID = NEntry.Key();
      const set<int> &NGridRanks = NEntry.Value();
      for (int Rank : NGridRanks) {
        NumFromNGridAndRank.Insert({NGridID,Rank}, 0);
        ++NumRecvs;
      }
    }
  }

  MPIRequests.Reserve(NumSends + NumRecvs);

  for (int MGridID : Domain.LocalGridIDs()) {
    auto &NGridIDsAndRanks = OverlappingNGridIDsAndRanksForLocalMGrid(MGridID);
    auto &NumFromNGridAndRank = NumOverlappingFromNGridAndRankForLocalMGrid(MGridID);
    for (auto &NEntry : NGridIDsAndRanks) {
      int NGridID = NEntry.Key();
      set<int> &NGridRanks = NEntry.Value();
      for (int Rank : NGridRanks) {
        long long &NumOverlapping = NumFromNGridAndRank({NGridID,Rank});
        MPI_Irecv(&NumOverlapping, 1, MPI_LONG_LONG, Rank, MGridID, Domain.Comm(),
          &MPIRequests.Append());
      }
    }
  }

  for (int NGridID : Domain.LocalGridIDs()) {
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    auto &NumFromMGridAndRank = NumOverlappingFromMGridAndRankForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const cart &MGridCart = MGridInfo.Cart();
      set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        long long &NumOverlapping = NumFromMGridAndRank({MGridID,Rank});
        const partition_data &PartitionData = MGridPartitionData({MGridID,Rank});
        auto Iter = OverlappingCellData.Find({MGridID,NGridID});
        if (Iter != OverlappingCellData.End()) {
          const overlapping_cell_data &CellData = Iter->Value();
          for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
            for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
              for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
                tuple<int> Point = {i,j,k};
                if (CellData.Cells(Point) != NO_CELL) {
                  tuple<int> Cell = CellData.Indexer.ToTuple(CellData.Cells(Point));
                  if (MGridCart.MapToRange(PartitionData.CellCoverRange, Cell).Present()) {
                    ++NumOverlapping;
                  }
                }
              }
            }
          }
        }
        MPI_Isend(&NumOverlapping, 1, MPI_LONG_LONG, Rank, MGridID, Domain.Comm(),
          &MPIRequests.Append());
      }
    }
  }

  MPI_Waitall(MPIRequests.Count(), MPIRequests.Data(), MPI_STATUSES_IGNORE);

  MPIRequests.Clear();

  for (int NGridID : Domain.LocalGridIDs()) {
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    auto &NumFromMGridAndRank = NumOverlappingFromMGridAndRankForLocalNGrid(NGridID);
    NumFromMGridAndRank.EraseIf([](const typename elem_map<int,2,long long>::entry &Entry) -> bool {
      return Entry.Value() == 0;
    });
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      set<int> &MGridRanks = MEntry.Value();
      MGridRanks.EraseIf([&](int Rank) -> bool {
        return !NumFromMGridAndRank.Contains({MGridID,Rank});
      });
    }
    MGridIDsAndRanks.EraseIf([](const typename map<int,set<int>>::entry &MEntry) -> bool {
      return MEntry.Value().Empty();
    });
  }

  for (int MGridID : Domain.LocalGridIDs()) {
    auto &NGridIDsAndRanks = OverlappingNGridIDsAndRanksForLocalMGrid(MGridID);
    auto &NumFromNGridAndRank = NumOverlappingFromNGridAndRankForLocalMGrid(MGridID);
    NumFromNGridAndRank.EraseIf([](const typename elem_map<int,2,long long>::entry &Entry) -> bool {
      return Entry.Value() == 0;
    });
    for (auto &NEntry : NGridIDsAndRanks) {
      int NGridID = NEntry.Key();
      set<int> &NGridRanks = NEntry.Value();
      NGridRanks.EraseIf([&](int Rank) -> bool {
        return !NumFromNGridAndRank.Contains({NGridID,Rank});
      });
    }
    NGridIDsAndRanks.EraseIf([](const typename map<int,set<int>>::entry &NEntry) -> bool {
      return NEntry.Value().Empty();
    });
  }

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    elem_map<int,2,long long> NumOverlappedByMGridForLocalNGrid;
    for (int NGridID : Domain.LocalGridIDs()) {
      const grid &NGrid = Domain.Grid(NGridID);
      auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
      for (int MGridID : Domain.GridIDs()) {
        elem<int,2> OverlapID = {MGridID,NGridID};
        if (!Options_.Overlappable(OverlapID)) continue;
        long long &NumOverlapped = NumOverlappedByMGridForLocalNGrid.Insert(OverlapID, 0);
        if (MGridIDsAndRanks.Contains(MGridID)) {
          const overlapping_cell_data &CellData = OverlappingCellData(OverlapID);
          for (long long l = 0; l < CellData.Cells.Count(); ++l) {
            if (CellData.Cells[l] != NO_CELL) {
              ++NumOverlapped;
            }
          }
        }
        MPI_Allreduce(MPI_IN_PLACE, &NumOverlapped, 1, MPI_LONG_LONG, MPI_SUM, NGrid.Comm());
      }
    }
    for (int MGridID : Domain.GridIDs()) {
      for (int NGridID : Domain.GridIDs()) {
        elem<int,2> OverlapID = {MGridID,NGridID};
        if (Options_.Overlappable(OverlapID) && Domain.GridIsLocal(NGridID)) {
          const grid &NGrid = Domain.Grid(NGridID);
          long long NumOverlapped = NumOverlappedByMGridForLocalNGrid(OverlapID);
          if (NumOverlapped > 0) {
            const grid_info &MGridInfo = Domain.GridInfo(MGridID);
            std::string NumOverlappedString = core::FormatNumber(NumOverlapped, "points", "point");
            Logger.LogDebug(NGrid.Comm().Rank() == 0, 3, "Detected %s overlapped by grid %s on grid "
              "%s.", NumOverlappedString, MGridInfo.Name(), NGrid.Name());
          }
        }
        MPI_Barrier(Domain.Comm());
      }
    }
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done searching for overlapping cells.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Creating and filling overlap data "
      "structures...");
  }

  elem_set<int,2> OverlappingGridIDs;

  for (int NGridID : Domain.LocalGridIDs()) {
    const grid &NGrid = Domain.Grid(NGridID);
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      OverlappingGridIDs.Insert({MGridID,NGridID});
    }
    for (int MGridID : Domain.GridIDs()) {
      elem<int,2> OverlapID = {MGridID,NGridID};
      if (Options_.Overlappable(OverlapID)) {
        int Overlaps = OverlappingGridIDs.Contains(OverlapID);
        if (NGrid.Comm().Rank() > 0) {
          MPI_Reduce(&Overlaps, nullptr, 1, MPI_INT, MPI_MAX, 0, NGrid.Comm());
        } else {
          MPI_Reduce(MPI_IN_PLACE, &Overlaps, 1, MPI_INT, MPI_MAX, 0, NGrid.Comm());
          if (Overlaps) {
            OverlappingGridIDs.Insert(OverlapID);
          }
        }
      }
    }
  }

  for (int NGridID : Domain.GridIDs()) {
    bool IsNGridRoot = false;
    int NGridRootRank;
    if (Domain.GridIsLocal(NGridID)) {
      const grid &NGrid = Domain.Grid(NGridID);
      IsNGridRoot = NGrid.Comm().Rank() == 0;
      if (IsNGridRoot) {
        NGridRootRank = Domain.Comm().Rank();
      }
    }
    core::BroadcastAnySource(&NGridRootRank, 1, MPI_INT, IsNGridRoot, Domain.Comm());
    for (int MGridID : Domain.GridIDs()) {
      elem<int,2> OverlapID = {MGridID,NGridID};
      if (Options_.Overlappable(OverlapID)) {
        int Overlaps;
        if (IsNGridRoot) Overlaps = OverlappingGridIDs.Contains(OverlapID);
        MPI_Bcast(&Overlaps, 1, MPI_INT, NGridRootRank, Domain.Comm());
        if (Overlaps) {
          OverlappingGridIDs.Insert(OverlapID);
        }
      }
    }
  }

  auto OverlapComponentEditHandle = Domain.EditComponent<overlap_component>(OverlapComponentID_);
  overlap_component &OverlapComponent = *OverlapComponentEditHandle;

  OverlapComponent.ClearOverlaps();
  OverlapComponent.CreateOverlaps(OverlappingGridIDs);

  struct overlap_m_data {
    long long NumOverlapping;
    array<int,2> Cells;
    array<double,2> Coords;
    array<int,2> Destinations;
    explicit overlap_m_data(long long NumOverlapping_):
      NumOverlapping(NumOverlapping_),
      Cells({{MAX_DIMS,NumOverlapping_}}),
      Coords({{MAX_DIMS,NumOverlapping_}}),
      Destinations({{MAX_DIMS,NumOverlapping_}})
    {}
  };

  map<int,elem_map<int,2,overlap_m_data>> OverlapMSendDataForLocalNGrid;
  map<int,elem_map<int,2,overlap_m_data>> OverlapMRecvDataForLocalMGrid;
  elem_map<int,2,overlap_m_data> OverlapMLocalToLocalData;

  for (int NGridID : Domain.LocalGridIDs()) {
    auto &OverlapMSendData = OverlapMSendDataForLocalNGrid.Insert(NGridID);
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    auto &NumFromMGridAndRank = NumOverlappingFromMGridAndRankForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        long long NumOverlapping = NumFromMGridAndRank({MGridID,Rank});
        if (Rank != Domain.Comm().Rank()) {
          OverlapMSendData.Insert({MGridID,Rank}, NumOverlapping);
        } else {
          OverlapMLocalToLocalData.Insert({MGridID,NGridID}, NumOverlapping);
        }
      }
    }
  }

  for (int MGridID : Domain.LocalGridIDs()) {
    auto &OverlapMRecvData = OverlapMRecvDataForLocalMGrid.Insert(MGridID);
    auto &NGridIDsAndRanks = OverlappingNGridIDsAndRanksForLocalMGrid(MGridID);
    auto &NumFromNGridAndRank = NumOverlappingFromNGridAndRankForLocalMGrid(MGridID);
    for (auto &NEntry : NGridIDsAndRanks) {
      int NGridID = NEntry.Key();
      const set<int> &NGridRanks = NEntry.Value();
      for (int Rank : NGridRanks) {
        if (Rank != Domain.Comm().Rank()) {
          long long NumOverlapping = NumFromNGridAndRank({NGridID,Rank});
          OverlapMRecvData.Insert({NGridID,Rank}, NumOverlapping);
        }
      }
    }
  }

  NumSends = 0;
  for (int NGridID : Domain.LocalGridIDs()) {
    NumSends += OverlapMSendDataForLocalNGrid(NGridID).Count();
  }

  NumRecvs = 0;
  for (int MGridID : Domain.LocalGridIDs()) {
    NumRecvs += OverlapMRecvDataForLocalMGrid(MGridID).Count();
  }

  MPIRequests.Reserve(3*(NumSends + NumRecvs));

  for (int MGridID : Domain.LocalGridIDs()) {
    auto &OverlapMRecvData = OverlapMRecvDataForLocalMGrid(MGridID);
    for (auto &Entry : OverlapMRecvData) {
      int Rank = Entry.Key()(1);
      overlap_m_data &OverlapMData = Entry.Value();
      long long NumOverlapping = OverlapMData.NumOverlapping;
      MPI_Irecv(OverlapMData.Cells.Data(), MAX_DIMS*NumOverlapping, MPI_INT, Rank, MGridID,
        Domain.Comm(), &MPIRequests.Append());
      MPI_Irecv(OverlapMData.Coords.Data(), MAX_DIMS*NumOverlapping, MPI_DOUBLE, Rank, MGridID,
        Domain.Comm(), &MPIRequests.Append());
      MPI_Irecv(OverlapMData.Destinations.Data(), MAX_DIMS*NumOverlapping, MPI_INT, Rank, MGridID,
        Domain.Comm(), &MPIRequests.Append());
    }
  }

  for (int NGridID : Domain.LocalGridIDs()) {
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    const geometry &Geometry = GeometryComponent.Geometry(NGridID);
    auto &Coords = Geometry.Coords();
    auto &OverlapMSendData = OverlapMSendDataForLocalNGrid(NGridID);
    for (auto &Entry : OverlapMSendData) {
      int MGridID = Entry.Key()(0);
      int Rank = Entry.Key()(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const cart &MGridCart = MGridInfo.Cart();
      overlap_m_data &OverlapMData = Entry.Value();
      const partition_data &PartitionData = MGridPartitionData({MGridID,Rank});
      const cell_coord_data &CoordData = CellCoordData({MGridID,Rank});
      const overlapping_cell_data &CellData = OverlappingCellData({MGridID,NGridID});
      long long iOverlapping = 0;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (CellData.Cells(Point) != NO_CELL) {
              tuple<int> Cell = CellData.Indexer.ToTuple(CellData.Cells(Point));
              auto MaybeCoverCell = MGridCart.MapToRange(PartitionData.CellCoverRange, Cell);
              if (MaybeCoverCell) {
                tuple<int> MappedCell;
                auto MaybeOwnCell = MGridCart.MapToRange(PartitionData.CellLocalRange, Cell);
                if (MaybeOwnCell) {
                  MappedCell = *MaybeOwnCell;
                } else {
                  MappedCell = *MaybeCoverCell;
                }
                tuple<double> PointCoords = {
                  Coords(0)(Point),
                  Coords(1)(Point),
                  Coords(2)(Point)
                };
                OverlapMData.Cells(0,iOverlapping) = MappedCell(0);
                OverlapMData.Cells(1,iOverlapping) = MappedCell(1);
                OverlapMData.Cells(2,iOverlapping) = MappedCell(2);
                auto MaybeLocalCoords = core::CoordsInCell(NumDims, CoordData.Coords,
                  CoordData.GeometryType, MappedCell, PointCoords);
                if (MaybeLocalCoords) {
                  const tuple<double> &LocalCoords = *MaybeLocalCoords;
                  OverlapMData.Coords(0,iOverlapping) = LocalCoords(0);
                  OverlapMData.Coords(1,iOverlapping) = LocalCoords(1);
                  OverlapMData.Coords(2,iOverlapping) = LocalCoords(2);
                } else {
                  Logger.LogWarning(true, "Failed to compute local coordinates of point "
                    "(%i,%i,%i) of grid %s inside cell (%i,%i,%i) of grid %s.", Point(0),
                    Point(1), Point(2), NGrid.Name(), MappedCell(0), MappedCell(1),
                    MappedCell(2), Domain.GridInfo(MGridID).Name());
                }
                OverlapMData.Destinations(0,iOverlapping) = Point(0);
                OverlapMData.Destinations(1,iOverlapping) = Point(1);
                OverlapMData.Destinations(2,iOverlapping) = Point(2);
                ++iOverlapping;
              }
            }
          }
        }
      }
      long long NumOverlapping = OverlapMData.NumOverlapping;
      MPI_Isend(OverlapMData.Cells.Data(), MAX_DIMS*NumOverlapping, MPI_INT, Rank, MGridID,
        Domain.Comm(), &MPIRequests.Append());
      MPI_Isend(OverlapMData.Coords.Data(), MAX_DIMS*NumOverlapping, MPI_DOUBLE, Rank, MGridID,
        Domain.Comm(), &MPIRequests.Append());
      MPI_Isend(OverlapMData.Destinations.Data(), MAX_DIMS*NumOverlapping, MPI_INT, Rank, MGridID,
        Domain.Comm(), &MPIRequests.Append());
    }
  }

  for (int NGridID : Domain.LocalGridIDs()) {
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    const geometry &Geometry = GeometryComponent.Geometry(NGridID);
    auto &Coords = Geometry.Coords();
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const cart &MGridCart = MGridInfo.Cart();
      const set<int> &MGridRanks = MEntry.Value();
      for (int Rank : MGridRanks) {
        if (Rank != Domain.Comm().Rank()) continue;
        const partition_data &PartitionData = MGridPartitionData({MGridID,Rank});
        const cell_coord_data &CoordData = CellCoordData({MGridID,Rank});
        const overlapping_cell_data &CellData = OverlappingCellData({MGridID,NGridID});
        overlap_m_data &OverlapMData = OverlapMLocalToLocalData({MGridID,NGridID});
        long long iOverlapping = 0;
        for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
          for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
            for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
              tuple<int> Point = {i,j,k};
              if (CellData.Cells(Point) != NO_CELL) {
                tuple<int> Cell = CellData.Indexer.ToTuple(CellData.Cells(Point));
                auto MaybeCoverCell = MGridCart.MapToRange(PartitionData.CellCoverRange, Cell);
                if (MaybeCoverCell) {
                  tuple<int> MappedCell;
                  auto MaybeOwnCell = MGridCart.MapToRange(PartitionData.CellLocalRange, Cell);
                  if (MaybeOwnCell) {
                    MappedCell = *MaybeOwnCell;
                  } else {
                    MappedCell = *MaybeCoverCell;
                  }
                  tuple<double> PointCoords = {
                    Coords(0)(Point),
                    Coords(1)(Point),
                    Coords(2)(Point)
                  };
                  OverlapMData.Cells(0,iOverlapping) = MappedCell(0);
                  OverlapMData.Cells(1,iOverlapping) = MappedCell(1);
                  OverlapMData.Cells(2,iOverlapping) = MappedCell(2);
                  auto MaybeLocalCoords = core::CoordsInCell(NumDims, CoordData.Coords,
                    CoordData.GeometryType, MappedCell, PointCoords);
                  if (MaybeLocalCoords) {
                    const tuple<double> &LocalCoords = *MaybeLocalCoords;
                    OverlapMData.Coords(0,iOverlapping) = LocalCoords(0);
                    OverlapMData.Coords(1,iOverlapping) = LocalCoords(1);
                    OverlapMData.Coords(2,iOverlapping) = LocalCoords(2);
                  } else {
                    Logger.LogWarning(true, "Failed to compute local coordinates of point "
                      "(%i,%i,%i) of grid %s inside cell (%i,%i,%i) of grid %s.", Point(0),
                      Point(1), Point(2), NGrid.Name(), MappedCell(0), MappedCell(1),
                      MappedCell(2), Domain.GridInfo(MGridID).Name());
                  }
                  OverlapMData.Destinations(0,iOverlapping) = Point(0);
                  OverlapMData.Destinations(1,iOverlapping) = Point(1);
                  OverlapMData.Destinations(2,iOverlapping) = Point(2);
                  ++iOverlapping;
                }
              }
            }
          }
        }
      }
    }
  }

  MPI_Waitall(MPIRequests.Count(), MPIRequests.Data(), MPI_STATUSES_IGNORE);

  MPIRequests.Clear();

  struct overlap_m_edit {
    edit_handle<overlap_m> Overlap;
    long long NumOverlapping = 0;
    edit_handle<array<int,2>> Cells;
    edit_handle<array<double,2>> Coords;
    edit_handle<array<int,2>> Destinations;
    edit_handle<array<int>> DestinationRanks;
  };

  struct overlap_n_edit {
    edit_handle<overlap_n> Overlap;
    long long NumOverlapping = 0;
    edit_handle<array<int,2>> Points;
    edit_handle<array<int,2>> Sources;
    edit_handle<array<int>> SourceRanks;
  };

  elem_map<int,2,overlap_m_edit> OverlapMEdits;
  elem_map<int,2,overlap_n_edit> OverlapNEdits;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    overlap_m_edit &Edit = OverlapMEdits.Insert(OverlapID);
    Edit.Overlap = OverlapComponent.EditOverlapM(OverlapID);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    overlap_n_edit &Edit = OverlapNEdits.Insert(OverlapID);
    Edit.Overlap = OverlapComponent.EditOverlapN(OverlapID);
  }

  for (int MGridID : Domain.LocalGridIDs()) {
    auto &NGridIDsAndRanks = OverlappingNGridIDsAndRanksForLocalMGrid(MGridID);
    auto &NumFromNGridAndRank = NumOverlappingFromNGridAndRankForLocalMGrid(MGridID);
    for (auto &NEntry : NGridIDsAndRanks) {
      int NGridID = NEntry.Key();
      const set<int> &NGridRanks = NEntry.Value();
      for (int Rank : NGridRanks) {
        long long NumOverlapping = NumFromNGridAndRank({NGridID,Rank});
        OverlapMEdits({MGridID,NGridID}).NumOverlapping += NumOverlapping;
      }
    }
  }

  for (int NGridID : Domain.LocalGridIDs()) {
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const overlapping_cell_data &CellData = OverlappingCellData({MGridID,NGridID});
      long long NumOverlapping = 0;
      for (long long l = 0; l < CellData.Cells.Count(); ++l) {
        if (CellData.Cells[l] != NO_CELL) {
          ++NumOverlapping;
        }
      }
      OverlapNEdits({MGridID,NGridID}).NumOverlapping = NumOverlapping;
    }
  }

  for (auto &Entry : OverlapMEdits) {
    overlap_m_edit &Edit = Entry.Value();
    Edit.Overlap->Resize(Edit.NumOverlapping);
    Edit.Cells = Edit.Overlap->EditCells();
    Edit.Coords = Edit.Overlap->EditCoords();
    Edit.Destinations = Edit.Overlap->EditDestinations();
    Edit.DestinationRanks = Edit.Overlap->EditDestinationRanks();
  }

  for (auto &Entry : OverlapNEdits) {
    overlap_n_edit &Edit = Entry.Value();
    Edit.Overlap->Resize(Edit.NumOverlapping);
    Edit.Points = Edit.Overlap->EditPoints();
    Edit.Sources = Edit.Overlap->EditSources();
    Edit.SourceRanks = Edit.Overlap->EditSourceRanks();
  }

  for (int MGridID : Domain.LocalGridIDs()) {
    const grid &MGrid = Domain.Grid(MGridID);
    const range &CellLocalRange = MGrid.CellLocalRange();
    auto &NGridIDsAndRanks = OverlappingNGridIDsAndRanksForLocalMGrid(MGridID);
    auto &OverlapMRecvData = OverlapMRecvDataForLocalMGrid(MGridID);
    for (auto &NEntry : NGridIDsAndRanks) {
      int NGridID = NEntry.Key();
      const set<int> &NGridRanks = NEntry.Value();
      field_indexer NGridGlobalIndexer(Domain.GridInfo(NGridID).GlobalRange());
      overlap_m_edit &Edit = OverlapMEdits({MGridID,NGridID});
      long long NumOverlapping = Edit.Overlap->Count();
      // Want to have the same order as overlap N data
      array<long long> DestinationPointIndices;
      DestinationPointIndices.Reserve(NumOverlapping);
      for (int Rank : NGridRanks) {
        overlap_m_data *OverlapMData;
        if (Rank != Domain.Comm().Rank()) {
          OverlapMData = &OverlapMRecvData({NGridID,Rank});
        } else {
          OverlapMData = &OverlapMLocalToLocalData({MGridID,NGridID});
        }
        for (long long iFromRank = 0; iFromRank < OverlapMData->NumOverlapping; ++iFromRank) {
          tuple<int> DestinationPoint = {
            OverlapMData->Destinations(0,iFromRank),
            OverlapMData->Destinations(1,iFromRank),
            OverlapMData->Destinations(2,iFromRank)
          };
          DestinationPointIndices.Append(NGridGlobalIndexer.ToIndex(DestinationPoint));
        }
      }
      array<long long> ROrder = ArrayOrder(DestinationPointIndices);
      // Need to use order on LHS since RHS is not one contiguous array
      array<long long> LOrder({NumOverlapping});
      for (long long iOverlapping = 0; iOverlapping < NumOverlapping; ++iOverlapping) {
        LOrder(ROrder(iOverlapping)) = iOverlapping;
      }
      long long iOverlapping = 0;
      for (int Rank : NGridRanks) {
        overlap_m_data *OverlapMData;
        if (Rank != Domain.Comm().Rank()) {
          OverlapMData = &OverlapMRecvData({NGridID,Rank});
        } else {
          OverlapMData = &OverlapMLocalToLocalData({MGridID,NGridID});
        }
        for (long long iFromRank = 0; iFromRank < OverlapMData->NumOverlapping; ++iFromRank) {
          long long iOrder = LOrder(iOverlapping);
          tuple<int> Cell = {
            OverlapMData->Cells(0,iFromRank),
            OverlapMData->Cells(1,iFromRank),
            OverlapMData->Cells(2,iFromRank)
          };
          (*Edit.Cells)(0,iOrder) = Cell(0);
          (*Edit.Cells)(1,iOrder) = Cell(1);
          (*Edit.Cells)(2,iOrder) = Cell(2);
          (*Edit.Coords)(0,iOrder) = OverlapMData->Coords(0,iFromRank);
          (*Edit.Coords)(1,iOrder) = OverlapMData->Coords(1,iFromRank);
          (*Edit.Coords)(2,iOrder) = OverlapMData->Coords(2,iFromRank);
          (*Edit.Destinations)(0,iOrder) = OverlapMData->Destinations(0,iFromRank);
          (*Edit.Destinations)(1,iOrder) = OverlapMData->Destinations(1,iFromRank);
          (*Edit.Destinations)(2,iOrder) = OverlapMData->Destinations(2,iFromRank);
          (*Edit.DestinationRanks)(iOrder) = CellLocalRange.Contains(Cell) ? Rank : -1;
          ++iOverlapping;
        }
      }
    }
  }

  for (int NGridID : Domain.LocalGridIDs()) {
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    auto &MGridIDsAndRanks = OverlappingMGridIDsAndRanksForLocalNGrid(NGridID);
    for (auto &MEntry : MGridIDsAndRanks) {
      int MGridID = MEntry.Key();
      const set<int> &MGridRanks = MEntry.Value();
      overlap_n_edit &Edit = OverlapNEdits({MGridID,NGridID});
      const overlapping_cell_data &CellData = OverlappingCellData({MGridID,NGridID});
      long long iOverlapping = 0;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (CellData.Cells(Point) != NO_CELL) {
              tuple<int> Cell = CellData.Indexer.ToTuple(CellData.Cells(Point));
              int SourceRank = -1;
              for (int Rank : MGridRanks) {
                const partition_data &PartitionData = MGridPartitionData({MGridID,Rank});
                if (PartitionData.CellCoverRange.Contains(Cell)) {
                  SourceRank = Rank;
                  break;
                }
              }
              (*Edit.Points)(0,iOverlapping) = Point(0);
              (*Edit.Points)(1,iOverlapping) = Point(1);
              (*Edit.Points)(2,iOverlapping) = Point(2);
              (*Edit.Sources)(0,iOverlapping) = Cell(0);
              (*Edit.Sources)(1,iOverlapping) = Cell(1);
              (*Edit.Sources)(2,iOverlapping) = Cell(2);
              (*Edit.SourceRanks)(iOverlapping) = SourceRank;
              ++iOverlapping;
            }
          }
        }
      }
    }
  }

  OverlapMEdits.Clear();
  OverlapNEdits.Clear();

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done creating and filling overlap data "
      "structures.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Setting up overlap exchanges...");
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    AssemblyData.LocalOverlapMAuxData.Insert(OverlapID);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    AssemblyData.LocalOverlapNAuxData.Insert(OverlapID);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    local_overlap_m_aux_data &OverlapMAuxData = AssemblyData.LocalOverlapMAuxData(OverlapID);
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    const array<int,2> &Cells = OverlapM.Cells();
    array<int,3> CellExtents({{2,MAX_DIMS,OverlapM.Count()}});
    for (long long iOverlapping = 0; iOverlapping < OverlapM.Count(); ++iOverlapping) {
      for (int iDim = 0; iDim < NumDims; ++iDim) {
        CellExtents(0,iDim,iOverlapping) = Cells(iDim,iOverlapping);
        CellExtents(1,iDim,iOverlapping) = Cells(iDim,iOverlapping)+2;
      }
      for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
        CellExtents(0,iDim,iOverlapping) = 0;
        CellExtents(1,iDim,iOverlapping) = 1;
      }
    }
    OverlapMAuxData.CollectMap = core::collect_map(MGrid.Partition(), std::move(CellExtents));
    OverlapMAuxData.SendMap = core::send_map(OverlapM.DestinationRanks());
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    OverlapNAuxData.RecvMap = core::recv_map(OverlapN.SourceRanks());
    OverlapNAuxData.DisperseMap = core::disperse_map(OverlapN.Points());
  }

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done setting up overlap exchanges.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Creating auxiliary overlap data...");
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const array<int,2> &Points = OverlapN.Points();
    OverlapMask.Assign(NGrid.SharedPartition(), false);
    for (long long iOverlapping = 0; iOverlapping < OverlapN.Count(); ++iOverlapping) {
      tuple<int> Point = {
        Points(0,iOverlapping),
        Points(1,iOverlapping),
        Points(2,iOverlapping)
      };
      OverlapMask(Point) = true;
    }
    OverlapMask.Exchange();
  }

  struct exchange_m {
    floating_ref_generator FloatingRefGenerator;
    array<double,3> InterpCoefs;
    core::collect Collect;
    core::send Send;
    array<double> SendBuffer;
  };

  struct exchange_n {
    core::recv Recv;
  };

  elem_map<int,2,exchange_m> ExchangeMs;
  elem_map<int,2,exchange_n> ExchangeNs;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    const local_overlap_m_aux_data &OverlapMAuxData = AssemblyData.LocalOverlapMAuxData(OverlapID);
    exchange_m &ExchangeM = ExchangeMs.Insert(OverlapID);
    array<double,3> &InterpCoefs = ExchangeM.InterpCoefs;
    InterpCoefs.Resize({{MAX_DIMS,2,OverlapM.Count()}});
    for (long long iOverlapping = 0; iOverlapping < OverlapM.Count(); ++iOverlapping) {
      for (int iDim = 0; iDim < NumDims; ++iDim) {
        elem<double,2> Coefs = core::LagrangeInterpLinear(OverlapM.Coords()(iDim,iOverlapping));
        InterpCoefs(iDim,0,iOverlapping) = Coefs(0);
        InterpCoefs(iDim,1,iOverlapping) = Coefs(1);
      }
      for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
        InterpCoefs(iDim,0,iOverlapping) = 1.;
        InterpCoefs(iDim,1,iOverlapping) = 0.;
      }
    }
    auto InterpCoefsRef = ExchangeM.FloatingRefGenerator.Generate(InterpCoefs);
    ExchangeM.Collect = core::CreateCollectInterp(Context_, MGrid.Comm(), MGrid.Cart(),
      MGrid.LocalRange(), OverlapMAuxData.CollectMap, data_type::DOUBLE, 1, MGrid.ExtendedRange(),
      array_layout::COLUMN_MAJOR, InterpCoefsRef);
    ExchangeM.Send = core::CreateSend(Context_, Domain.Comm(), OverlapMAuxData.SendMap,
      data_type::DOUBLE, 1, 0);
    ExchangeM.SendBuffer.Resize({OverlapM.Count()});
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    exchange_n &ExchangeN = ExchangeNs.Insert(OverlapID);
    ExchangeN.Recv = core::CreateRecv(Context_, Domain.Comm(), OverlapNAuxData.RecvMap,
      data_type::DOUBLE, 1, 0);
    OverlapNAuxData.Volumes.Resize({OverlapN.Count()});
  }

  array<request> Requests;
  Requests.Reserve(OverlapComponent.LocalOverlapMCount() + OverlapComponent.LocalOverlapNCount());

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::recv &Recv = ExchangeN.Recv;
    request &Request = Requests.Append();
    double *VolumesData = OverlapNAuxData.Volumes.Data();
    Request = Recv.Recv(&VolumesData);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const geometry &Geometry = GeometryComponent.Geometry(MGridID);
    const distributed_field<double> &Volumes = Geometry.Volumes();
    exchange_m &ExchangeM = ExchangeMs(OverlapID);
    core::collect &Collect = ExchangeM.Collect;
    core::send &Send = ExchangeM.Send;
    const double *VolumesData = Volumes.Data();
    double *SendBufferData = ExchangeM.SendBuffer.Data();
    Collect.Collect(&VolumesData, &SendBufferData);
    request &Request = Requests.Append();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  ExchangeMs.Clear();
  ExchangeNs.Clear();

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done creating auxiliary overlap data.");
  }

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Done detecting overlap between grids.");

}

void assembler::InferBoundaries_() {

  domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Inferring non-overlapping boundaries...");

  auto &OverlapComponent = Domain.Component<overlap_component>(OverlapComponentID_);
  assembly_data &AssemblyData = *AssemblyData_;

  auto StateComponentEditHandle = Domain.EditComponent<state_component>(StateComponentID_);
  state_component &StateComponent = *StateComponentEditHandle;

  map<int,long long> NumInferredForGrid;
  if (Logger.LoggingDebug()) {
    for (int GridID : Domain.LocalGridIDs()) {
      NumInferredForGrid.Insert(GridID, 0);
    }
  }

  for (int GridID : Domain.LocalGridIDs()) {
    if (!Options_.InferBoundaries(GridID)) continue;
    const grid &Grid = Domain.Grid(GridID);
    const range &LocalRange = Grid.LocalRange();
    long long NumExtended = Grid.ExtendedRange().Count();
    local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(GridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    distributed_field<bool> &DomainBoundaryMask = GridAuxData.DomainBoundaryMask;
    auto StateEditHandle = StateComponent.EditState(GridID);
    auto FlagsEditHandle = StateEditHandle->EditFlags();
    distributed_field<state_flags> &Flags = *FlagsEditHandle;
    distributed_field<bool> InferredBoundaryMask(Grid.SharedPartition());
    core::DetectEdge(ActiveMask, core::edge_type::INNER, core::mask_bc::FALSE, false,
      InferredBoundaryMask);
    for (long long l = 0; l < NumExtended; ++l) {
      InferredBoundaryMask[l] = InferredBoundaryMask[l] && !DomainBoundaryMask[l];
    }
    for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
      if (OverlapID(1) == GridID) {
        const field<bool> &OverlapMask = OverlapComponent.OverlapN(OverlapID).Mask();
        for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
          for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
            for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
              InferredBoundaryMask(i,j,k) = InferredBoundaryMask(i,j,k) && !OverlapMask(i,j,k);
            }
          }
        }
      }
    }
    InferredBoundaryMask.Exchange();
    for (long long l = 0; l < NumExtended; ++l) {
      if (InferredBoundaryMask[l]) {
        Flags[l] |= state_flags::DOMAIN_BOUNDARY | state_flags::INFERRED_DOMAIN_BOUNDARY;
      }
    }
    for (long long l = 0; l < NumExtended; ++l) {
      DomainBoundaryMask[l] = DomainBoundaryMask[l] || InferredBoundaryMask[l];
    }
    if (Logger.LoggingDebug()) {
      NumInferredForGrid.Insert(GridID, core::CountDistributedMask(InferredBoundaryMask));
    }
  }

  if (Logger.LoggingDebug()) {
    for (int GridID : Domain.GridIDs()) {
      if (Domain.GridIsLocal(GridID)) {
        const grid &Grid = Domain.Grid(GridID);
        long long NumInferred = NumInferredForGrid(GridID);
        if (NumInferred > 0) {
          std::string NumInferredString = core::FormatNumber(NumInferred, "points", "point");
          Logger.LogDebug(Grid.Comm().Rank() == 0, 2, "%s marked as boundaries on grid %s.",
            NumInferredString, Grid.Name());
        }
      }
      MPI_Barrier(Domain.Comm());
    }
  }

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Done inferring non-overlapping boundaries.");

}

void assembler::CutBoundaryHoles_() {

  domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Cutting boundary holes...");

  int NumDims = Domain.Dimension();
  assembly_data &AssemblyData = *AssemblyData_;

  auto &OverlapComponent = Domain.Component<overlap_component>(OverlapComponentID_);

  auto StateComponentEditHandle = Domain.EditComponent<state_component>(StateComponentID_);
  state_component &StateComponent = *StateComponentEditHandle;

  elem_set<int,2> LocalCutMPairIDs;
  elem_set<int,2> LocalCutNPairIDs;
  set<int> LocalCutMGridIDs;
  set<int> LocalCutNGridIDs;

  for (int MGridID : Domain.GridIDs()) {
    for (int NGridID : Domain.GridIDs()) {
      elem<int,2> IDPair = {MGridID,NGridID};
      if (Options_.CutBoundaryHoles(IDPair)) {
        if (Domain.GridIsLocal(MGridID)) {
          LocalCutMPairIDs.Insert(IDPair);
          LocalCutMGridIDs.Insert(MGridID);
        }
        if (Domain.GridIsLocal(NGridID)) {
          LocalCutNPairIDs.Insert(IDPair);
          LocalCutNGridIDs.Insert(NGridID);
        }
      }
    }
  }

  // Exchanging in reverse, from points to cells

  struct reverse_exchange_m {
    core::recv_map RecvMap;
    core::recv Recv;
    array<bool> RecvBuffer;
  };

  struct reverse_exchange_n {
    core::send_map SendMap;
    core::send Send;
    array<bool> SendBuffer;
  };

  elem_map<int,2,reverse_exchange_m> ReverseExchangeMs;
  elem_map<int,2,reverse_exchange_n> ReverseExchangeNs;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (!LocalCutNPairIDs.Contains({OverlapID(1),OverlapID(0)})) continue;
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    reverse_exchange_m &ExchangeM = ReverseExchangeMs.Insert(OverlapID);
    ExchangeM.RecvMap = core::recv_map(OverlapM.DestinationRanks());
    ExchangeM.Recv = core::CreateRecv(Context_, Domain.Comm(), ExchangeM.RecvMap, data_type::BOOL,
      1, 0);
    ExchangeM.RecvBuffer.Resize({OverlapM.Count()});
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (!LocalCutMPairIDs.Contains({OverlapID(1),OverlapID(0)})) continue;
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    reverse_exchange_n &ExchangeN = ReverseExchangeNs.Insert(OverlapID);
    ExchangeN.SendMap = core::send_map(OverlapN.SourceRanks());
    ExchangeN.Send = core::CreateSend(Context_, Domain.Comm(), ExchangeN.SendMap, data_type::BOOL,
      1, 0);
    ExchangeN.SendBuffer.Resize({OverlapN.Count()});
  }

  array<request> Requests;

  Requests.Reserve(OverlapComponent.LocalOverlapMCount() + OverlapComponent.LocalOverlapNCount());

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (!LocalCutNPairIDs.Contains({OverlapID(1),OverlapID(0)})) continue;
    reverse_exchange_m &ExchangeM = ReverseExchangeMs(OverlapID);
    core::recv &Recv = ExchangeM.Recv;
    request &Request = Requests.Append();
    bool *RecvBufferData = ExchangeM.RecvBuffer.Data();
    Request = Recv.Recv(&RecvBufferData);
  }

  elem_map<int,2,distributed_field<bool>> OverlapEdgeMasks;

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (!LocalCutMPairIDs.Contains({OverlapID(1),OverlapID(0)})) continue;
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const array<int,2> &Points = OverlapN.Points();
    const local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    const distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    distributed_field<bool> &OverlapEdgeMask = OverlapEdgeMasks.Insert(OverlapID);
    core::DetectEdge(OverlapMask, core::edge_type::INNER, core::mask_bc::FALSE, false,
      OverlapEdgeMask);
    reverse_exchange_n &ExchangeN = ReverseExchangeNs(OverlapID);
    core::send &Send = ExchangeN.Send;
    array<bool> &SendBuffer = ExchangeN.SendBuffer;
    for (long long iOverlapped = 0; iOverlapped < OverlapN.Count(); ++iOverlapped) {
      tuple<int> Point = {
        Points(0,iOverlapped),
        Points(1,iOverlapped),
        Points(2,iOverlapped)
      };
      SendBuffer(iOverlapped) = OverlapEdgeMask(Point);
    }
    request &Request = Requests.Append();
    bool *SendBufferData = SendBuffer.Data();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  OverlapEdgeMasks.Clear();

  elem_map<int,2,distributed_field<bool>> CoverMasks;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (!LocalCutNPairIDs.Contains({OverlapID(1),OverlapID(0)})) continue;
    int MGridID = OverlapID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const range &LocalRange = MGrid.LocalRange();
    const range &CellExtendedRange = MGrid.CellExtendedRange();
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    const array<int,2> &Cells = OverlapM.Cells();
    reverse_exchange_m &ExchangeM = ReverseExchangeMs(OverlapID);
    const array<bool> &RecvBuffer = ExchangeM.RecvBuffer;
    distributed_field<bool> CellCoverMask(MGrid.SharedPartition(), false);
    for (long long iOverlapping = 0; iOverlapping < OverlapM.Count(); ++iOverlapping) {
      tuple<int> Cell = {
        Cells(0,iOverlapping),
        Cells(1,iOverlapping),
        Cells(2,iOverlapping)
      };
      CellCoverMask(Cell) = CellCoverMask(Cell) || RecvBuffer(iOverlapping);
    }
    CellCoverMask.Exchange();
    distributed_field<bool> &CoverMask = CoverMasks.Insert(OverlapID, MGrid.SharedPartition(),
      false);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          range NeighborCellRange;
          for (int iDim = 0; iDim < NumDims; ++iDim) {
            NeighborCellRange.Begin(iDim) = Point(iDim)-1;
            NeighborCellRange.End(iDim) = Point(iDim)+1;
          }
          for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
            NeighborCellRange.Begin(iDim) = 0;
            NeighborCellRange.End(iDim) = 1;
          }
          NeighborCellRange = IntersectRanges(CellExtendedRange, NeighborCellRange);
          for (int o = NeighborCellRange.Begin(2); o < NeighborCellRange.End(2); ++o) {
            for (int n = NeighborCellRange.Begin(1); n < NeighborCellRange.End(1); ++n) {
              for (int m = NeighborCellRange.Begin(0); m < NeighborCellRange.End(0); ++m) {
                tuple<int> NeighborCell = {m,n,o};
                CoverMask(Point) = CoverMask(Point) || CellCoverMask(NeighborCell);
              }
            }
          }
        }
      }
    }
    CoverMask.Exchange();
  }

  elem_map<int,2,int> NumDilates;

  auto GlobalAny = [](const distributed_field<bool> &Mask) -> bool {
    int Any = ArrayAny(Mask);
    MPI_Allreduce(MPI_IN_PLACE, &Any, 1, MPI_INT, MPI_LOR, Mask.Comm());
    return Any;
  };

  for (auto &OverlapID : LocalCutNPairIDs) {
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    long long NumExtended = NGrid.ExtendedRange().Count();
    const local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    const distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    distributed_field<bool> UncoveredEdgeMask;
    core::DetectEdge(OverlapMask, core::edge_type::OUTER, core::mask_bc::MIRROR, false,
      UncoveredEdgeMask);
    distributed_field<bool> &CoverMask = CoverMasks({NGridID,MGridID});
    for (long long l = 0; l < NumExtended; ++l) {
      UncoveredEdgeMask[l] = UncoveredEdgeMask[l] && !CoverMask[l];
    }
    int &d = NumDilates.Insert(OverlapID, 0);
    while (GlobalAny(UncoveredEdgeMask)) {
      core::DilateMask(CoverMask, 1, core::mask_bc::FALSE);
      for (long long l = 0; l < NumExtended; ++l) {
        UncoveredEdgeMask[l] = UncoveredEdgeMask[l] && !CoverMask[l];
      }
      ++d;
    }
  }

  Requests.Reserve(OverlapComponent.LocalOverlapMCount() + OverlapComponent.LocalOverlapNCount());

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (!LocalCutNPairIDs.Contains({OverlapID(1),OverlapID(0)})) continue;
    reverse_exchange_m &ExchangeM = ReverseExchangeMs(OverlapID);
    core::recv &Recv = ExchangeM.Recv;
    request &Request = Requests.Append();
    bool *RecvBufferData = ExchangeM.RecvBuffer.Data();
    Request = Recv.Recv(&RecvBufferData);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (!LocalCutMPairIDs.Contains({OverlapID(1),OverlapID(0)})) continue;
    int NGridID = OverlapID(1);
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(NGridID);
    const distributed_field<bool> &BoundaryMask = GridAuxData.DomainBoundaryMask;
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const array<int,2> &Points = OverlapN.Points();
    reverse_exchange_n &ExchangeN = ReverseExchangeNs(OverlapID);
    core::send &Send = ExchangeN.Send;
    array<bool> &SendBuffer = ExchangeN.SendBuffer;
    for (long long iOverlapped = 0; iOverlapped < OverlapN.Count(); ++iOverlapped) {
      tuple<int> Point = {
        Points(0,iOverlapped),
        Points(1,iOverlapped),
        Points(2,iOverlapped)
      };
      SendBuffer(iOverlapped) = BoundaryMask(Point);
    }
    request &Request = Requests.Append();
    bool *SendBufferData = SendBuffer.Data();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (!LocalCutNPairIDs.Contains({OverlapID(1),OverlapID(0)})) continue;
    int MGridID = OverlapID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const range &LocalRange = MGrid.LocalRange();
    const range &CellExtendedRange = MGrid.CellExtendedRange();
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    const array<int,2> &Cells = OverlapM.Cells();
    reverse_exchange_m &ExchangeM = ReverseExchangeMs(OverlapID);
    const array<bool> &RecvBuffer = ExchangeM.RecvBuffer;
    distributed_field<bool> CellCoverMask(MGrid.SharedPartition(), false);
    for (long long iOverlapping = 0; iOverlapping < OverlapM.Count(); ++iOverlapping) {
      tuple<int> Cell = {
        Cells(0,iOverlapping),
        Cells(1,iOverlapping),
        Cells(2,iOverlapping)
      };
      CellCoverMask(Cell) = CellCoverMask(Cell) || RecvBuffer(iOverlapping);
    }
    CellCoverMask.Exchange();
    distributed_field<bool> &CoverMask = CoverMasks(OverlapID);
    CoverMask.Fill(false);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          range NeighborCellRange;
          for (int iDim = 0; iDim < NumDims; ++iDim) {
            NeighborCellRange.Begin(iDim) = Point(iDim)-1;
            NeighborCellRange.End(iDim) = Point(iDim)+1;
          }
          for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
            NeighborCellRange.Begin(iDim) = 0;
            NeighborCellRange.End(iDim) = 1;
          }
          NeighborCellRange = IntersectRanges(CellExtendedRange, NeighborCellRange);
          for (int o = NeighborCellRange.Begin(2); o < NeighborCellRange.End(2); ++o) {
            for (int n = NeighborCellRange.Begin(1); n < NeighborCellRange.End(1); ++n) {
              for (int m = NeighborCellRange.Begin(0); m < NeighborCellRange.End(0); ++m) {
                tuple<int> NeighborCell = {m,n,o};
                CoverMask(Point) = CoverMask(Point) || CellCoverMask(NeighborCell);
              }
            }
          }
        }
      }
    }
    CoverMask.Exchange();
  }

  ReverseExchangeMs.Clear();
  ReverseExchangeNs.Clear();

  elem_map<int,2,distributed_field<bool>> &ProjectedBoundaryMasks = AssemblyData
    .ProjectedBoundaryMasks;

  for (auto &OverlapID : LocalCutNPairIDs) {
    if (!Options_.CutBoundaryHoles(OverlapID)) continue;
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    long long NumExtended = NGrid.ExtendedRange().Count();
    const local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    const distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    distributed_field<bool> &CoverMask = CoverMasks({NGridID,MGridID});
    core::DilateMask(CoverMask, NumDilates(OverlapID), core::mask_bc::FALSE);
    for (long long l = 0; l < NumExtended; ++l) {
      CoverMask[l] = CoverMask[l] && !OverlapMask[l];
    }
    distributed_field<bool> &ProjectedBoundaryMask = ProjectedBoundaryMasks.Insert(OverlapID);
    core::DetectEdge(CoverMask, core::edge_type::OUTER, core::mask_bc::FALSE, false,
      ProjectedBoundaryMask);
    for (long long l = 0; l < NumExtended; ++l) {
      ProjectedBoundaryMask[l] = ProjectedBoundaryMask[l] && OverlapMask[l];
    }
  }

  NumDilates.Clear();

  map<int,distributed_field<bool>> BoundaryMasks;
  map<int,distributed_field<bool>> InteriorMasks;

  for (int GridID : LocalCutNGridIDs) {
    const grid &Grid = Domain.Grid(GridID);
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(GridID);
    const distributed_field<bool> &OwnBoundaryMask = GridAuxData.DomainBoundaryMask;
    distributed_field<bool> &BoundaryMask = BoundaryMasks.Insert(GridID);
    BoundaryMask = OwnBoundaryMask;
    distributed_field<bool> &InteriorMask = InteriorMasks.Insert(GridID);
    InteriorMask.Assign(Grid.SharedPartition(), false);
  }

  for (auto &OverlapID : LocalCutNPairIDs) {
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    long long NumExtended = NGrid.ExtendedRange().Count();
    const local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    const distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    distributed_field<bool> &BoundaryMask = BoundaryMasks(NGridID);
    const distributed_field<bool> &ProjectedBoundaryMask = ProjectedBoundaryMasks(OverlapID);
    for (long long l = 0; l < NumExtended; ++l) {
      BoundaryMask[l] = BoundaryMask[l] || ProjectedBoundaryMask[l];
    }
    distributed_field<bool> ProjectedInteriorEdgeMask;
    core::DetectEdge(ProjectedBoundaryMask, core::edge_type::OUTER, core::mask_bc::FALSE, false,
      ProjectedInteriorEdgeMask);
    for (long long l = 0; l < NumExtended; ++l) {
      ProjectedInteriorEdgeMask[l] = ProjectedInteriorEdgeMask[l] && OverlapMask[l];
    }
    distributed_field<bool> &InteriorMask = InteriorMasks(NGridID);
    for (long long l = 0; l < NumExtended; ++l) {
      InteriorMask[l] = InteriorMask[l] || ProjectedInteriorEdgeMask[l];
    }
  }

  map<int,distributed_field<bool>> BoundaryHoleMasks;

  for (int GridID : LocalCutNGridIDs) {
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(GridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    const distributed_field<bool> &OwnBoundaryMask = GridAuxData.DomainBoundaryMask;
    distributed_field<bool> &BoundaryMask = BoundaryMasks(GridID);
    distributed_field<bool> &InteriorMask = InteriorMasks(GridID);
    distributed_field<bool> &BoundaryHoleMask = BoundaryHoleMasks.Insert(GridID,
      Grid.SharedPartition());
    long long NumInterior = core::CountDistributedMask(InteriorMask);
    if (NumInterior > 0) {
      core::FloodMask(InteriorMask, BoundaryMask);
      distributed_field<bool> ActiveEdgeMask;
      core::DetectEdge(ActiveMask, core::edge_type::INNER, core::mask_bc::FALSE, false,
        ActiveEdgeMask);
      distributed_field<bool> InteriorEdgeMask;
      core::DetectEdge(InteriorMask, core::edge_type::OUTER, core::mask_bc::FALSE, false,
        InteriorEdgeMask);
      for (long long l = 0; l < NumExtended; ++l) {
        if (OwnBoundaryMask[l] && ActiveEdgeMask[l] && !InteriorEdgeMask[l]) {
          BoundaryMask[l] = false;
        }
      }
      for (long long l = 0; l < NumExtended; ++l) {
        BoundaryHoleMask[l] = ActiveMask[l] && !InteriorMask[l] && !BoundaryMask[l];
      }
    } else {
      BoundaryHoleMask.Fill(false);
    }
  }

  map<int,long long> NumRemovedForGrid;

  for (int GridID : LocalCutNGridIDs) {
    long long &NumRemoved = NumRemovedForGrid.Insert(GridID);
    NumRemoved = core::CountDistributedMask(BoundaryHoleMasks(GridID));
  }

  for (int GridID : LocalCutNGridIDs) {
    if (NumRemovedForGrid(GridID) == 0) continue;
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(GridID);
    {
      auto StateEditHandle = StateComponent.EditState(GridID);
      auto FlagsEditHandle = StateEditHandle->EditFlags();
      distributed_field<state_flags> &Flags = *FlagsEditHandle;
      const distributed_field<bool> &BoundaryHoleMask = BoundaryHoleMasks(GridID);
      for (long long l = 0; l < NumExtended; ++l) {
        if (BoundaryHoleMask[l]) {
          Flags[l] = (Flags[l] & ~state_flags::ACTIVE) | state_flags::BOUNDARY_HOLE;
        }
      }
    }
    auto &Flags = StateComponent.State(GridID).Flags();
    GenerateActiveMask(Grid, Flags, GridAuxData.ActiveMask);
    GenerateCellActiveMask(Grid, Flags, GridAuxData.CellActiveMask);
    GenerateDomainBoundaryMask(Grid, Flags, GridAuxData.DomainBoundaryMask);
    GenerateInternalBoundaryMask(Grid, Flags, GridAuxData.InternalBoundaryMask);
  }

  for (auto &OverlapID : LocalCutNPairIDs) {
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    long long NumExtended = NGrid.ExtendedRange().Count();
    local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(NGridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    for (long long l = 0; l < NumExtended; ++l) {
      OverlapMask[l] = OverlapMask[l] && ActiveMask[l];
    }
  }

  struct exchange_m {
    core::collect Collect;
    core::send Send;
    array<bool> SendBuffer;
  };

  struct exchange_n {
    core::recv Recv;
    array<bool> RecvBuffer;
    core::disperse Disperse;
  };

  elem_map<int,2,exchange_m> ExchangeMs;
  elem_map<int,2,exchange_n> ExchangeNs;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    const local_overlap_m_aux_data &OverlapMAuxData = AssemblyData.LocalOverlapMAuxData(OverlapID);
    exchange_m &ExchangeM = ExchangeMs.Insert(OverlapID);
    ExchangeM.Collect = core::CreateCollectAll(Context_, MGrid.Comm(), MGrid.Cart(),
      MGrid.LocalRange(), OverlapMAuxData.CollectMap, data_type::BOOL, 1, MGrid.ExtendedRange(),
      array_layout::COLUMN_MAJOR);
    ExchangeM.Send = core::CreateSend(Context_, Domain.Comm(), OverlapMAuxData.SendMap,
      data_type::BOOL, 1, 0);
    ExchangeM.SendBuffer.Resize({OverlapM.Count()});
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    exchange_n &ExchangeN = ExchangeNs.Insert(OverlapID);
    ExchangeN.Recv = core::CreateRecv(Context_, Domain.Comm(), OverlapNAuxData.RecvMap,
      data_type::BOOL, 1, 0);
    ExchangeN.RecvBuffer.Resize({OverlapN.Count()});
    ExchangeN.Disperse = core::CreateDisperseOverwrite(Context_, OverlapNAuxData.DisperseMap,
      data_type::BOOL, 1, NGrid.ExtendedRange(), array_layout::COLUMN_MAJOR);
  }

  Requests.Reserve(OverlapComponent.LocalOverlapMCount() + OverlapComponent.LocalOverlapNCount());

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::recv &Recv = ExchangeN.Recv;
    request &Request = Requests.Append();
    bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    Request = Recv.Recv(&RecvBufferData);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(MGridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    exchange_m &ExchangeM = ExchangeMs(OverlapID);
    core::collect &Collect = ExchangeM.Collect;
    core::send &Send = ExchangeM.Send;
    const bool *ActiveMaskData = ActiveMask.Data();
    bool *SendBufferData = ExchangeM.SendBuffer.Data();
    Collect.Collect(&ActiveMaskData, &SendBufferData);
    request &Request = Requests.Append();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::disperse &Disperse = ExchangeN.Disperse;
    const bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    bool *OverlapMaskData = OverlapMask.Data();
    Disperse.Disperse(&RecvBufferData, &OverlapMaskData);
    OverlapMask.Exchange();
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingDebug()) {
    for (int GridID : Domain.GridIDs()) {
      if (LocalCutNGridIDs.Contains(GridID)) {
        const grid &Grid = Domain.Grid(GridID);
        long long NumRemoved = NumRemovedForGrid(GridID);
        if (NumRemoved > 0) {
          std::string NumRemovedString = core::FormatNumber(NumRemoved, "points", "point");
          Logger.LogDebug(Grid.Comm().Rank() == 0, 2, "%s removed from grid %s.",
            NumRemovedString, Grid.Name());
        }
      }
      MPI_Barrier(Domain.Comm());
    }
    Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Done cutting boundary holes.");
  }

}

void assembler::LocateOuterFringe_() {

  domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Locating outer fringe points...");

  assembly_data &AssemblyData = *AssemblyData_;

  auto StateComponentEditHandle = Domain.EditComponent<state_component>(StateComponentID_);
  state_component &StateComponent = *StateComponentEditHandle;

  map<int,distributed_field<bool>> &OuterFringeMasks = AssemblyData.OuterFringeMasks;

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    OuterFringeMasks.Insert(GridID, Grid.SharedPartition(), false);
  }

  for (int GridID : Domain.LocalGridIDs()) {
    if (Options_.FringeSize(GridID) == 0) continue;
    const grid &Grid = Domain.Grid(GridID);
    const range &ExtendedRange = Grid.ExtendedRange();
    long long NumExtended = ExtendedRange.Count();
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(GridID);
    const core::partition_pool &PartitionPool = GridAuxData.PartitionPool;
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    const distributed_field<bool> &DomainBoundaryMask = GridAuxData.DomainBoundaryMask;
    const distributed_field<bool> &InternalBoundaryMask = GridAuxData.InternalBoundaryMask;
    distributed_field<bool> BoundaryMask(Grid.SharedPartition());
    for (long long l = 0; l < NumExtended; ++l) {
      BoundaryMask[l] = DomainBoundaryMask[l] || InternalBoundaryMask[l];
    }
    distributed_field<bool> BoundaryEdgeMask;
    core::DetectEdge(BoundaryMask, core::edge_type::OUTER, core::mask_bc::FALSE, true,
      BoundaryEdgeMask, &PartitionPool);
    distributed_field<bool> NonBoundaryMask(Grid.SharedPartition());
    for (long long l = 0; l < NumExtended; ++l) {
      NonBoundaryMask[l] = ActiveMask[l] && !BoundaryMask[l];
    }
    distributed_field<bool> NonBoundaryEdgeMask;
    core::DetectEdge(NonBoundaryMask, core::edge_type::OUTER, core::mask_bc::FALSE, true,
      NonBoundaryEdgeMask, &PartitionPool);
    distributed_field<bool> CoverMask;
    core::DetectEdge(ActiveMask, core::edge_type::OUTER, core::mask_bc::FALSE, true,
      CoverMask, &PartitionPool);
    for (long long l = 0; l < CoverMask.Extents().Count(); ++l) {
      CoverMask[l] = CoverMask[l] && (NonBoundaryEdgeMask[l] || !BoundaryEdgeMask[l]);
    }
    core::DilateMask(CoverMask, Options_.FringeSize(GridID), core::mask_bc::FALSE);
    distributed_field<bool> &OuterFringeMask = OuterFringeMasks(GridID);
    for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
      for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
        for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          OuterFringeMask(Point) = ActiveMask(Point) && CoverMask(Point);
        }
      }
    }
  }

  map<int,long long> NumOuterFringeForGrid;

  for (int GridID : Domain.LocalGridIDs()) {
    long long &NumOuterFringe = NumOuterFringeForGrid.Insert(GridID);
    NumOuterFringe = core::CountDistributedMask(OuterFringeMasks(GridID));
  }

  for (int GridID : Domain.LocalGridIDs()) {
    if (Options_.FringeSize(GridID) == 0) continue;
    if (NumOuterFringeForGrid(GridID) == 0) continue;
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    const distributed_field<bool> &OuterFringeMask = OuterFringeMasks(GridID);
    auto StateEditHandle = StateComponent.EditState(GridID);
    auto FlagsEditHandle = StateEditHandle->EditFlags();
    distributed_field<state_flags> &Flags = *FlagsEditHandle;
    for (long long l = 0; l < NumExtended; ++l) {
      if (OuterFringeMask[l]) {
        Flags[l] = Flags[l] | (state_flags::FRINGE | state_flags::OUTER_FRINGE);
      }
    }
  }

  if (Logger.LoggingDebug()) {
    for (int GridID : Domain.GridIDs()) {
      if (Domain.GridIsLocal(GridID)) {
        const grid &Grid = Domain.Grid(GridID);
        long long NumOuterFringe = NumOuterFringeForGrid(GridID);
        if (NumOuterFringe > 0) {
          std::string NumOuterFringeString = core::FormatNumber(NumOuterFringe,
            "outer fringe points", "outer fringe point");
          Logger.LogDebug(Grid.Comm().Rank() == 0, 2, "%s on grid %s.", NumOuterFringeString,
            Grid.Name());
        }
      }
      MPI_Barrier(Domain.Comm());
    }
  }

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Done locating outer fringe points.");

}

void assembler::DetectOccluded_() {

  domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Detecting occluded points...");

  assembly_data &AssemblyData = *AssemblyData_;

  auto &GeometryComponent = Domain.Component<geometry_component>(GeometryComponentID_);

  auto StateComponentEditHandle = Domain.EditComponent<state_component>(StateComponentID_);
  state_component &StateComponent = *StateComponentEditHandle;

  auto &OverlapComponent = Domain.Component<overlap_component>(OverlapComponentID_);

  Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Computing pairwise occlusion...");

  elem_map<int,2,distributed_field<bool>> &PairwiseOcclusionMasks = AssemblyData
    .PairwiseOcclusionMasks;

  constexpr double TOLERANCE = 1.e-10;

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    const geometry &Geometry = GeometryComponent.Geometry(NGridID);
    const distributed_field<double> &Volumes = Geometry.Volumes();
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const field<bool> &BaseOverlapMask = OverlapN.Mask();
    const local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    const distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    const array<double> &OverlapVolumes = OverlapNAuxData.Volumes;
    distributed_field<bool> &PairwiseOcclusionMask = PairwiseOcclusionMasks.Insert(OverlapID);
    switch (Options_.Occludes(OverlapID)) {
    case occludes::COARSE: {
      PairwiseOcclusionMask.Assign(NGrid.SharedPartition(), false);
      long long iOverlapping = 0;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (BaseOverlapMask(Point)) {
              PairwiseOcclusionMask(Point) = OverlapMask(Point) && Volumes(Point) > (1.+TOLERANCE) *
                OverlapVolumes(iOverlapping);
              ++iOverlapping;
            }
          }
        }
      }
      PairwiseOcclusionMask.Exchange();
      break;
    }
    case occludes::ALL:
      PairwiseOcclusionMask.Assign(OverlapMask);
      break;
    case occludes::NONE:
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      break;
    }
  }

  // Exclude points that are overlapped by occluded points

  struct exchange_m {
    core::collect Collect;
    core::send Send;
    array<bool> SendBuffer;
  };

  struct exchange_n {
    core::recv Recv;
    array<bool> RecvBuffer;
    core::disperse Disperse;
  };

  elem_map<int,2,exchange_m> ExchangeMs;
  elem_map<int,2,exchange_n> ExchangeNs;

  auto MutuallyOccludes = [&](int MGridID, int NGridID) -> bool {
    return Options_.Occludes({MGridID,NGridID}) == occludes::COARSE &&
      Options_.Occludes({NGridID,MGridID}) == occludes::COARSE;
  };

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    const local_overlap_m_aux_data &OverlapMAuxData = AssemblyData.LocalOverlapMAuxData(OverlapID);
    exchange_m &ExchangeM = ExchangeMs.Insert(OverlapID);
    ExchangeM.Collect = core::CreateCollectNone(Context_, MGrid.Comm(), MGrid.Cart(),
      MGrid.LocalRange(), OverlapMAuxData.CollectMap, data_type::BOOL, 1, MGrid.ExtendedRange(),
      array_layout::COLUMN_MAJOR);
    ExchangeM.Send = core::CreateSend(Context_, Domain.Comm(), OverlapMAuxData.SendMap,
      data_type::BOOL, 1, 0);
    ExchangeM.SendBuffer.Resize({OverlapM.Count()});
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    exchange_n &ExchangeN = ExchangeNs.Insert(OverlapID);
    ExchangeN.Recv = core::CreateRecv(Context_, Domain.Comm(), OverlapNAuxData.RecvMap,
      data_type::BOOL, 1, 0);
    ExchangeN.RecvBuffer.Resize({OverlapN.Count()});
    ExchangeN.Disperse = core::CreateDisperseOverwrite(Context_, OverlapNAuxData.DisperseMap,
      data_type::BOOL, 1, NGrid.ExtendedRange(), array_layout::COLUMN_MAJOR);
  }

  array<request> Requests;
  Requests.Reserve(OverlapComponent.LocalOverlapMCount() + OverlapComponent.LocalOverlapNCount());

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    if (!MutuallyOccludes(MGridID, NGridID)) continue;
    if (NGridID < MGridID) continue;
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::recv &Recv = ExchangeN.Recv;
    request &Request = Requests.Append();
    bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    Request = Recv.Recv(&RecvBufferData);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    if (!MutuallyOccludes(MGridID, NGridID)) continue;
    if (NGridID < MGridID) continue;
    exchange_m &ExchangeM = ExchangeMs(OverlapID);
    core::collect &Collect = ExchangeM.Collect;
    core::send &Send = ExchangeM.Send;
    const bool *PairwiseOcclusionMaskData = PairwiseOcclusionMasks({NGridID,MGridID}).Data();
    bool *SendBufferData = ExchangeM.SendBuffer.Data();
    Collect.Collect(&PairwiseOcclusionMaskData, &SendBufferData);
    request &Request = Requests.Append();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    if (!MutuallyOccludes(MGridID, NGridID)) continue;
    if (NGridID < MGridID) continue;
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::disperse &Disperse = ExchangeN.Disperse;
    distributed_field<bool> &PairwiseOcclusionMask = PairwiseOcclusionMasks(OverlapID);
    const bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    bool *PairwiseOcclusionMaskData = PairwiseOcclusionMask.Data();
    Disperse.Disperse(&RecvBufferData, &PairwiseOcclusionMaskData);
    PairwiseOcclusionMask.Exchange();
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    if (!MutuallyOccludes(MGridID, NGridID)) continue;
    if (NGridID > MGridID) continue;
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::recv &Recv = ExchangeN.Recv;
    request &Request = Requests.Append();
    bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    Request = Recv.Recv(&RecvBufferData);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    if (!MutuallyOccludes(MGridID, NGridID)) continue;
    if (NGridID > MGridID) continue;
    exchange_m &ExchangeM = ExchangeMs(OverlapID);
    core::collect &Collect = ExchangeM.Collect;
    core::send &Send = ExchangeM.Send;
    const bool *PairwiseOcclusionMaskData = PairwiseOcclusionMasks({NGridID,MGridID}).Data();
    bool *SendBufferData = ExchangeM.SendBuffer.Data();
    Collect.Collect(&PairwiseOcclusionMaskData, &SendBufferData);
    request &Request = Requests.Append();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    if (!MutuallyOccludes(MGridID, NGridID)) continue;
    if (NGridID > MGridID) continue;
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::disperse &Disperse = ExchangeN.Disperse;
    distributed_field<bool> &PairwiseOcclusionMask = PairwiseOcclusionMasks(OverlapID);
    const bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    bool *PairwiseOcclusionMaskData = PairwiseOcclusionMask.Data();
    Disperse.Disperse(&RecvBufferData, &PairwiseOcclusionMaskData);
    PairwiseOcclusionMask.Exchange();
  }

  if (Logger.LoggingDebug()) {
    elem_map<int,2,long long> NumOccludedForGridPair;
    for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
      if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
      long long &NumOccluded = NumOccludedForGridPair.Insert(OverlapID);
      NumOccluded = core::CountDistributedMask(PairwiseOcclusionMasks(OverlapID));
    }
    MPI_Barrier(Domain.Comm());
    for (auto &OverlapID : OverlapComponent.OverlapIDs()) {
      int MGridID = OverlapID(0);
      int NGridID = OverlapID(1);
      if (Options_.Occludes(OverlapID) != occludes::NONE && Domain.GridIsLocal(NGridID)) {
        const grid &NGrid = Domain.Grid(NGridID);
        long long NumOccluded = NumOccludedForGridPair(OverlapID);
        if (NumOccluded > 0) {
          const grid_info &MGridInfo = Domain.GridInfo(MGridID);
          std::string NumOccludedString = core::FormatNumber(NumOccluded, "points", "point");
          Logger.LogDebug(NGrid.Comm().Rank() == 0, 3, "%s occluded by grid %s on grid %s.",
            NumOccludedString, MGridInfo.Name(), NGrid.Name());
        }
      }
      MPI_Barrier(Domain.Comm());
    }
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done computing pairwise occlusion.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Applying edge padding and smoothing...");
  }

  elem_map<int,2,distributed_field<bool>> DisallowMasks;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    const grid &MGrid = Domain.Grid(MGridID);
    long long NumExtended = MGrid.ExtendedRange().Count();
    const distributed_field<state_flags> &Flags = StateComponent.State(MGridID).Flags();
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(MGridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    distributed_field<bool> &DisallowMask = DisallowMasks.Insert(OverlapID,
      MGrid.SharedPartition());
    for (long long l = 0; l < NumExtended; ++l) {
      DisallowMask[l] = (Flags[l] & state_flags::OUTER_FRINGE) != state_flags::NONE;
    }
    if (Options_.Occludes({NGridID,MGridID}) != occludes::NONE) {
      const distributed_field<bool> &PairwiseOcclusionMask = PairwiseOcclusionMasks({NGridID,
        MGridID});
      // Not sure if the "|| !ActiveMask[l]" should be applied unconditionally? Below is how it is
      // in serial Overkit
      for (long long l = 0; l < NumExtended; ++l) {
        DisallowMask[l] = DisallowMask[l] || PairwiseOcclusionMask[l] || !ActiveMask[l];
      }
    }
    core::DilateMask(DisallowMask, Options_.EdgePadding(OverlapID), core::mask_bc::MIRROR);
  }

  Requests.Reserve(OverlapComponent.LocalOverlapMCount() + OverlapComponent.LocalOverlapNCount());

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::recv &Recv = ExchangeN.Recv;
    request &Request = Requests.Append();
    bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    Request = Recv.Recv(&RecvBufferData);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
    exchange_m &ExchangeM = ExchangeMs(OverlapID);
    core::collect &Collect = ExchangeM.Collect;
    core::send &Send = ExchangeM.Send;
    const bool *DisallowMaskData = DisallowMasks(OverlapID).Data();
    bool *SendBufferData = ExchangeM.SendBuffer.Data();
    Collect.Collect(&DisallowMaskData, &SendBufferData);
    request &Request = Requests.Append();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  DisallowMasks.Clear();

  elem_map<int,2,distributed_field<bool>> PaddingMasks;
  map<int,distributed_field<bool>> BaseOcclusionMasks;
  map<int,distributed_field<bool>> &OcclusionMasks = AssemblyData.OcclusionMasks;

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    BaseOcclusionMasks.Insert(GridID, Grid.SharedPartition(), false);
    OcclusionMasks.Insert(GridID, Grid.SharedPartition(), false);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    long long NumExtended = NGrid.ExtendedRange().Count();
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::disperse &Disperse = ExchangeN.Disperse;
    distributed_field<bool> &PairwiseOcclusionMask = PairwiseOcclusionMasks(OverlapID);
    const bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    distributed_field<bool> AllowMask(NGrid.SharedPartition());
    bool *AllowMaskData = AllowMask.Data();
    Disperse.Disperse(&RecvBufferData, &AllowMaskData);
    AllowMask.Exchange();
    distributed_field<bool> &PaddingMask = PaddingMasks.Insert(OverlapID, NGrid.SharedPartition());
    for (long long l = 0; l < NumExtended; ++l) {
      PaddingMask[l] = !AllowMask[l] && PairwiseOcclusionMask[l];
    }
    distributed_field<bool> &BaseOcclusionMask = BaseOcclusionMasks(NGridID);
    distributed_field<bool> &OcclusionMask = OcclusionMasks(NGridID);
    for (long long l = 0; l < NumExtended; ++l) {
      BaseOcclusionMask[l] = BaseOcclusionMask[l] || PairwiseOcclusionMask[l];
    }
    for (long long l = 0; l < NumExtended; ++l) {
      OcclusionMask[l] = OcclusionMask[l] || (PairwiseOcclusionMask[l] && !PaddingMask[l]);
    }
  }

  for (int GridID : Domain.LocalGridIDs()) {
    if (Options_.EdgeSmoothing(GridID) == 0) continue;
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    const distributed_field<bool> &BaseOcclusionMask = BaseOcclusionMasks(GridID);
    distributed_field<bool> &OcclusionMask = OcclusionMasks(GridID);
    core::DilateMask(OcclusionMask, Options_.EdgeSmoothing(GridID), core::mask_bc::MIRROR);
    core::ErodeMask(OcclusionMask, Options_.EdgeSmoothing(GridID), core::mask_bc::MIRROR);
    for (long long l = 0; l < NumExtended; ++l) {
      OcclusionMask[l] = OcclusionMask[l] && BaseOcclusionMask[l];
    }
  }

  BaseOcclusionMasks.Clear();

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    long long NumExtended = NGrid.ExtendedRange().Count();
    distributed_field<bool> &PaddingMask = PaddingMasks(OverlapID);
    distributed_field<bool> &PairwiseOcclusionMask = PairwiseOcclusionMasks(OverlapID);
    const distributed_field<bool> &OcclusionMask = OcclusionMasks(NGridID);
    for (long long l = 0; l < NumExtended; ++l) {
      PaddingMask[l] = PaddingMask[l] && !OcclusionMask[l];
    }
    for (long long l = 0; l < NumExtended; ++l) {
      PairwiseOcclusionMask[l] = PairwiseOcclusionMask[l] && !PaddingMask[l];
    }
  }

  if (Logger.LoggingDebug()) {
    elem_map<int,2,long long> NumPaddedForGridPair;
    for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
      if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
      long long &NumPadded = NumPaddedForGridPair.Insert(OverlapID);
      NumPadded = core::CountDistributedMask(PaddingMasks(OverlapID));
    }
    MPI_Barrier(Domain.Comm());
    for (auto &OverlapID : OverlapComponent.OverlapIDs()) {
      if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
      int MGridID = OverlapID(0);
      int NGridID = OverlapID(1);
      if (Domain.GridIsLocal(NGridID)) {
        const grid &NGrid = Domain.Grid(NGridID);
        long long NumPadded = NumPaddedForGridPair(OverlapID);
        if (NumPadded > 0) {
          const grid_info &MGridInfo = Domain.GridInfo(MGridID);
          std::string NumPaddedString = core::FormatNumber(NumPadded, "points", "point");
          Logger.LogDebug(NGrid.Comm().Rank() == 0, 3, "%s marked as not occluded by grid "
            "%s on grid %s.", NumPaddedString, MGridInfo.Name(), NGrid.Name());
        }
      }
      MPI_Barrier(Domain.Comm());
    }
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done applying edge padding and smoothing.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Accumulating occlusion...");
  }

  PaddingMasks.Clear();

  for (int GridID : Domain.LocalGridIDs()) {
    distributed_field<bool> &OcclusionMask = OcclusionMasks(GridID);
    OcclusionMask.Fill(false);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    long long NumExtended = NGrid.ExtendedRange().Count();
    distributed_field<bool> &PairwiseOcclusionMask = PairwiseOcclusionMasks(OverlapID);
    distributed_field<bool> &OcclusionMask = OcclusionMasks(NGridID);
    for (long long l = 0; l < NumExtended; ++l) {
      OcclusionMask[l] = OcclusionMask[l] || PairwiseOcclusionMask[l];
    }
  }

  map<int,long long> NumOccludedForGrid;

  for (int GridID : Domain.LocalGridIDs()) {
    long long &NumOccluded = NumOccludedForGrid.Insert(GridID);
    NumOccluded = core::CountDistributedMask(OcclusionMasks(GridID));
  }

  for (int GridID : Domain.LocalGridIDs()) {
    if (NumOccludedForGrid(GridID) == 0) continue;
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    auto StateEditHandle = StateComponent.EditState(GridID);
    auto FlagsEditHandle = StateEditHandle->EditFlags();
    distributed_field<state_flags> &Flags = *FlagsEditHandle;
    const distributed_field<bool> &OcclusionMask = OcclusionMasks(GridID);
    for (long long l = 0; l < NumExtended; ++l) {
      if (OcclusionMask[l]) {
        Flags[l] |= state_flags::OCCLUDED;
      }
    }
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingDebug()) {
    for (int GridID : Domain.GridIDs()) {
      if (Domain.GridIsLocal(GridID)) {
        const grid &Grid = Domain.Grid(GridID);
        long long NumOccluded = NumOccludedForGrid(GridID);
        if (NumOccluded > 0) {
          std::string NumOccludedString = core::FormatNumber(NumOccluded, "occluded points",
            "occluded point");
          Logger.LogDebug(Grid.Comm().Rank() == 0, 3, "%s on grid %s.", NumOccludedString,
            Grid.Name());
        }
      }
      MPI_Barrier(Domain.Comm());
    }
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done accumulating occlusion.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Done detecting occluded points.");
  }

}

void assembler::MinimizeOverlap_() {

  domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Minimizing overlap...");

  assembly_data &AssemblyData = *AssemblyData_;

  auto StateComponentEditHandle = Domain.EditComponent<state_component>(StateComponentID_);
  state_component &StateComponent = *StateComponentEditHandle;

  auto &OverlapComponent = Domain.Component<overlap_component>(OverlapComponentID_);

  const elem_map<int,2,distributed_field<bool>> &PairwiseOcclusionMasks = AssemblyData
    .PairwiseOcclusionMasks;
  const map<int,distributed_field<bool>> &OcclusionMasks = AssemblyData.OcclusionMasks;

  map<int,distributed_field<bool>> &OverlapMinimizationMasks = AssemblyData
    .OverlapMinimizationMasks;
  map<int,distributed_field<bool>> &InnerFringeMasks = AssemblyData.InnerFringeMasks;

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    OverlapMinimizationMasks.Insert(GridID, Grid.SharedPartition(), false);
    InnerFringeMasks.Insert(GridID, Grid.SharedPartition(), false);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (Options_.Occludes(OverlapID) == occludes::NONE) continue;
    if (!Options_.MinimizeOverlap(OverlapID)) continue;
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    long long NumExtended = NGrid.ExtendedRange().Count();
    const distributed_field<bool> &PairwiseOcclusionMask = PairwiseOcclusionMasks(OverlapID);
    distributed_field<bool> &OverlapMinimizationMask = OverlapMinimizationMasks(NGridID);
    for (long long l = 0; l < NumExtended; ++l) {
      OverlapMinimizationMask[l] = OverlapMinimizationMask[l] || PairwiseOcclusionMask[l];
    }
  }

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(GridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    distributed_field<bool> &OverlapMinimizationMask = OverlapMinimizationMasks(GridID);
    distributed_field<bool> &InnerFringeMask = InnerFringeMasks(GridID);
    if (Options_.FringeSize(GridID) > 0) {
      const distributed_field<bool> &OcclusionMask = OcclusionMasks(GridID);
      distributed_field<bool> RemovableMask(Grid.SharedPartition());
      for (long long l = 0; l < NumExtended; ++l) {
        RemovableMask[l] = OcclusionMask[l] || !ActiveMask[l];
      }
      core::ErodeMask(RemovableMask, Options_.FringeSize(GridID), core::mask_bc::TRUE);
      for (long long l = 0; l < NumExtended; ++l) {
        OverlapMinimizationMask[l] = OverlapMinimizationMask[l] && (RemovableMask[l] &&
          ActiveMask[l]);
      }
      InnerFringeMask.Fill(OverlapMinimizationMask);
      core::DilateMask(InnerFringeMask, Options_.FringeSize(GridID), core::mask_bc::FALSE);
      for (long long l = 0; l < NumExtended; ++l) {
        InnerFringeMask[l] = InnerFringeMask[l] && (ActiveMask[l] && !OverlapMinimizationMask[l]);
      }
    } else {
      // Serial Overkit has this, but I'm not sure why...
//       for (long long l = 0; l < NumExtended; ++l) {
//         OverlapMinimizationMask[l] = OverlapMinimizationMask[l] && ActiveMask[l];
//       }
    }
  }

  map<int,long long> NumRemovedForGrid;
  map<int,long long> NumInnerFringeForGrid;

  for (int GridID : Domain.LocalGridIDs()) {
    long long &NumRemoved = NumRemovedForGrid.Insert(GridID);
    NumRemoved = core::CountDistributedMask(OverlapMinimizationMasks(GridID));
    long long &NumInnerFringe = NumInnerFringeForGrid.Insert(GridID);
    NumInnerFringe = core::CountDistributedMask(InnerFringeMasks(GridID));
  }

  for (int GridID : Domain.LocalGridIDs()) {
    if (NumRemovedForGrid(GridID) == 0 && NumInnerFringeForGrid(GridID) == 0) continue;
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(GridID);
    {
      auto StateEditHandle = StateComponent.EditState(GridID);
      auto FlagsEditHandle = StateEditHandle->EditFlags();
      distributed_field<state_flags> &Flags = *FlagsEditHandle;
      const distributed_field<bool> &OverlapMinimizationMask = OverlapMinimizationMasks(GridID);
      const distributed_field<bool> &InnerFringeMask = InnerFringeMasks(GridID);
      for (long long l = 0; l < NumExtended; ++l) {
        if (OverlapMinimizationMask[l]) {
          Flags[l] = (Flags[l] & ~(state_flags::ACTIVE | state_flags::OUTER_FRINGE)) |
            state_flags::OVERLAP_MINIMIZED;
        } else if (InnerFringeMask[l]) {
          Flags[l] = Flags[l] | (state_flags::FRINGE | state_flags::INNER_FRINGE);
        }
      }
    }
    if (NumRemovedForGrid(GridID) > 0) {
      auto &Flags = StateComponent.State(GridID).Flags();
      GenerateActiveMask(Grid, Flags, GridAuxData.ActiveMask);
      GenerateCellActiveMask(Grid, Flags, GridAuxData.CellActiveMask);
      GenerateDomainBoundaryMask(Grid, Flags, GridAuxData.DomainBoundaryMask);
      GenerateInternalBoundaryMask(Grid, Flags, GridAuxData.InternalBoundaryMask);
      distributed_field<bool> &OuterFringeMask = AssemblyData.OuterFringeMasks(GridID);
      for (long long l = 0; l < NumExtended; ++l) {
        OuterFringeMask[l] = OuterFringeMask[l] && (Flags[l] & state_flags::ACTIVE) !=
          state_flags::NONE;
      }
    }
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    long long NumExtended = NGrid.ExtendedRange().Count();
    local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(NGridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    for (long long l = 0; l < NumExtended; ++l) {
      OverlapMask[l] = OverlapMask[l] && ActiveMask[l];
    }
  }

  struct exchange_m {
    core::collect Collect;
    core::send Send;
    array<bool> SendBuffer;
  };

  struct exchange_n {
    core::recv Recv;
    array<bool> RecvBuffer;
    core::disperse Disperse;
  };

  elem_map<int,2,exchange_m> ExchangeMs;
  elem_map<int,2,exchange_n> ExchangeNs;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    const local_overlap_m_aux_data &OverlapMAuxData = AssemblyData.LocalOverlapMAuxData(OverlapID);
    exchange_m &ExchangeM = ExchangeMs.Insert(OverlapID);
    ExchangeM.Collect = core::CreateCollectAll(Context_, MGrid.Comm(), MGrid.Cart(),
      MGrid.LocalRange(), OverlapMAuxData.CollectMap, data_type::BOOL, 1, MGrid.ExtendedRange(),
      array_layout::COLUMN_MAJOR);
    ExchangeM.Send = core::CreateSend(Context_, Domain.Comm(), OverlapMAuxData.SendMap,
      data_type::BOOL, 1, 0);
    ExchangeM.SendBuffer.Resize({OverlapM.Count()});
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    exchange_n &ExchangeN = ExchangeNs.Insert(OverlapID);
    ExchangeN.Recv = core::CreateRecv(Context_, Domain.Comm(), OverlapNAuxData.RecvMap,
      data_type::BOOL, 1, 0);
    ExchangeN.RecvBuffer.Resize({OverlapN.Count()});
    ExchangeN.Disperse = core::CreateDisperseOverwrite(Context_, OverlapNAuxData.DisperseMap,
      data_type::BOOL, 1, NGrid.ExtendedRange(), array_layout::COLUMN_MAJOR);
  }

  array<request> Requests;
  Requests.Reserve(OverlapComponent.LocalOverlapMCount() + OverlapComponent.LocalOverlapNCount());

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::recv &Recv = ExchangeN.Recv;
    request &Request = Requests.Append();
    bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    Request = Recv.Recv(&RecvBufferData);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const local_grid_aux_data &GridAuxData = AssemblyData.LocalGridAuxData(MGridID);
    const distributed_field<bool> &ActiveMask = GridAuxData.ActiveMask;
    exchange_m &ExchangeM = ExchangeMs(OverlapID);
    core::collect &Collect = ExchangeM.Collect;
    core::send &Send = ExchangeM.Send;
    const bool *ActiveMaskData = ActiveMask.Data();
    bool *SendBufferData = ExchangeM.SendBuffer.Data();
    Collect.Collect(&ActiveMaskData, &SendBufferData);
    request &Request = Requests.Append();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    local_overlap_n_aux_data &OverlapNAuxData = AssemblyData.LocalOverlapNAuxData(OverlapID);
    distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::disperse &Disperse = ExchangeN.Disperse;
    const bool *RecvBufferData = ExchangeN.RecvBuffer.Data();
    bool *OverlapMaskData = OverlapMask.Data();
    Disperse.Disperse(&RecvBufferData, &OverlapMaskData);
    OverlapMask.Exchange();
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingDebug()) {
    for (int GridID : Domain.GridIDs()) {
      if (Domain.GridIsLocal(GridID)) {
        const grid &Grid = Domain.Grid(GridID);
        long long NumRemoved = NumRemovedForGrid(GridID);
        if (NumRemoved > 0) {
          std::string NumRemovedString = core::FormatNumber(NumRemoved, "points", "point");
          Logger.LogDebug(Grid.Comm().Rank() == 0, 3, "%s removed from grid %s.", NumRemovedString,
            Grid.Name());
        }
      }
      MPI_Barrier(Domain.Comm());
    }
    Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Done minimizing overlap.");
  }

}

void assembler::GenerateConnectivityData_() {

  domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Generating connectivity data...");

  int NumDims = Domain.Dimension();
  const assembly_data &AssemblyData = *AssemblyData_;

  const map<int,local_grid_aux_data> &LocalGridAuxData = AssemblyData.LocalGridAuxData;
  const elem_map<int,2,local_overlap_m_aux_data> &LocalOverlapMAuxData = AssemblyData
    .LocalOverlapMAuxData;
  const elem_map<int,2,local_overlap_n_aux_data> &LocalOverlapNAuxData = AssemblyData
    .LocalOverlapNAuxData;

  auto StateComponentEditHandle = Domain.EditComponent<state_component>(StateComponentID_);
  state_component &StateComponent = *StateComponentEditHandle;

  auto &OverlapComponent = Domain.Component<overlap_component>(OverlapComponentID_);

  Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Locating receiver points...");

  const map<int,distributed_field<bool>> &OcclusionMasks = AssemblyData.OcclusionMasks;
  const map<int,distributed_field<bool>> &OuterFringeMasks = AssemblyData.OuterFringeMasks;
  const map<int,distributed_field<bool>> &InnerFringeMasks = AssemblyData.InnerFringeMasks;

  map<int,distributed_field<bool>> ReceiverMasks;

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    const distributed_field<bool> &ActiveMask = LocalGridAuxData(GridID).ActiveMask;
    const distributed_field<bool> &OcclusionMask = OcclusionMasks(GridID);
    const distributed_field<bool> &OuterFringeMask = OuterFringeMasks(GridID);
    const distributed_field<bool> &InnerFringeMask = InnerFringeMasks(GridID);
    distributed_field<bool> &ReceiverMask = ReceiverMasks.Insert(GridID, Grid.SharedPartition());
    for (long long l = 0; l < NumExtended; ++l) {
      ReceiverMask[l] = OuterFringeMask[l] || InnerFringeMask[l] || (OcclusionMask[l] &&
        ActiveMask[l]);
    }
  }

  map<int,long long> NumReceiversForGrid;

  for (int GridID : Domain.LocalGridIDs()) {
    long long &NumReceivers = NumReceiversForGrid.Insert(GridID);
    NumReceivers = core::CountDistributedMask(ReceiverMasks(GridID));
  }

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    for (int GridID : Domain.GridIDs()) {
      if (Domain.GridIsLocal(GridID)) {
        const grid &Grid = Domain.Grid(GridID);
        long long NumReceivers = NumReceiversForGrid(GridID);
        if (NumReceivers > 0) {
          std::string NumReceiversString = core::FormatNumber(NumReceivers, "receiver points",
            "receiver point");
          Logger.LogDebug(Grid.Comm().Rank() == 0, 3, "%s on grid %s.", NumReceiversString,
            Grid.Name());
        }
      }
      MPI_Barrier(Domain.Comm());
    }
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done locating receiver points.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Choosing donors...");
  }

  map<int,int> MaxReceiverDistances;

  for (int GridID : Domain.LocalGridIDs()) {
    MaxReceiverDistances.Insert(GridID, 0);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    int &MaxDistance = MaxReceiverDistances(MGridID);
    MaxDistance = Max(MaxDistance, Options_.EdgePadding(OverlapID));
  }

  map<int,distributed_field<int>> ReceiverDistancesForGrid;

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    const distributed_field<bool> &ReceiverMask = ReceiverMasks(GridID);
    int MaxDistance = MaxReceiverDistances(GridID);
    distributed_field<int> &ReceiverDistances = ReceiverDistancesForGrid.Insert(GridID,
      Grid.SharedPartition(), MaxDistance);
    distributed_field<bool> CoverMask;
    core::DetectEdge(ReceiverMask, core::edge_type::INNER, core::mask_bc::MIRROR, false,
      CoverMask);
    for (int Distance = 0; Distance < MaxDistance; ++Distance) {
      for (long long l = 0; l < NumExtended; ++l) {
        if (ReceiverDistances[l] == MaxDistance && CoverMask[l]) {
          ReceiverDistances[l] = ReceiverMask[l] ? -Distance : Distance;
        }
      }
      core::DilateMask(CoverMask, 1, core::mask_bc::MIRROR);
    }
  }

  struct exchange_m {
    core::collect Collect;
    core::send Send;
    array<int> SendBuffer;
  };

  struct exchange_n {
    core::recv Recv;
  };

  elem_map<int,2,exchange_m> ExchangeMs;
  elem_map<int,2,exchange_n> ExchangeNs;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    const local_overlap_m_aux_data &OverlapMAuxData = LocalOverlapMAuxData(OverlapID);
    exchange_m &ExchangeM = ExchangeMs.Insert(OverlapID);
    ExchangeM.Collect = core::CreateCollectMin(Context_, MGrid.Comm(), MGrid.Cart(),
      MGrid.LocalRange(), OverlapMAuxData.CollectMap, data_type::INT, 1, MGrid.ExtendedRange(),
      array_layout::COLUMN_MAJOR);
    ExchangeM.Send = core::CreateSend(Context_, Domain.Comm(), OverlapMAuxData.SendMap,
      data_type::INT, 1, 0);
    ExchangeM.SendBuffer.Resize({OverlapM.Count()});
  }

  elem_map<int,2,array<int>> OverlapReceiverDistances;

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const local_overlap_n_aux_data &OverlapNAuxData = LocalOverlapNAuxData(OverlapID);
    exchange_n &ExchangeN = ExchangeNs.Insert(OverlapID);
    ExchangeN.Recv = core::CreateRecv(Context_, Domain.Comm(), OverlapNAuxData.RecvMap,
      data_type::INT, 1, 0);
    array<int> &ReceiverDistances = OverlapReceiverDistances.Insert(OverlapID);
    ReceiverDistances.Resize({OverlapN.Count()});
  }

  array<request> Requests;
  Requests.Reserve(OverlapComponent.LocalOverlapMCount() + OverlapComponent.LocalOverlapNCount());

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    exchange_n &ExchangeN = ExchangeNs(OverlapID);
    core::recv &Recv = ExchangeN.Recv;
    request &Request = Requests.Append();
    array<int> &ReceiverDistances = OverlapReceiverDistances(OverlapID);
    int *ReceiverDistancesData = ReceiverDistances.Data();
    Request = Recv.Recv(&ReceiverDistancesData);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    int MGridID = OverlapID(0);
    const distributed_field<int> &ReceiverDistances = ReceiverDistancesForGrid(MGridID);
    exchange_m &ExchangeM = ExchangeMs(OverlapID);
    core::collect &Collect = ExchangeM.Collect;
    core::send &Send = ExchangeM.Send;
    const int *ReceiverDistancesData = ReceiverDistances.Data();
    int *SendBufferData = ExchangeM.SendBuffer.Data();
    Collect.Collect(&ReceiverDistancesData, &SendBufferData);
    request &Request = Requests.Append();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  ExchangeMs.Clear();
  ExchangeNs.Clear();

  map<int,field<int>> DonorGridIDsForLocalGrid;
  map<int,field<double>> DonorNormalizedDistancesForLocalGrid;
  map<int,field<double>> DonorVolumesForLocalGrid;

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    DonorGridIDsForLocalGrid.Insert(GridID, Grid.LocalRange(), -1);
    DonorNormalizedDistancesForLocalGrid.Insert(GridID, Grid.LocalRange());
    DonorVolumesForLocalGrid.Insert(GridID, Grid.LocalRange());
  }

  constexpr double TOLERANCE = 1.e-12;

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (Options_.ConnectionType(OverlapID) == connection_type::NONE) continue;
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const field<bool> &BaseOverlapMask = OverlapN.Mask();
    const local_overlap_n_aux_data &OverlapNAuxData = LocalOverlapNAuxData(OverlapID);
    const distributed_field<bool> &OverlapMask = OverlapNAuxData.OverlapMask;
    const array<double> &OverlapVolumes = OverlapNAuxData.Volumes;
    const distributed_field<bool> &ReceiverMask = ReceiverMasks(NGridID);
    const array<int> &ReceiverDistances = OverlapReceiverDistances(OverlapID);
    field<int> &DonorGridIDs = DonorGridIDsForLocalGrid(NGridID);
    field<double> &DonorNormalizedDistances = DonorNormalizedDistancesForLocalGrid(NGridID);
    field<double> &DonorVolumes = DonorVolumesForLocalGrid(NGridID);
    long long iOverlapping = 0;
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          if (BaseOverlapMask(Point)) {
            if (ReceiverMask(Point) && OverlapMask(Point)) {
              double ReceiverDistance = double(ReceiverDistances(iOverlapping));
              double MaxDistance = double(Max(Options_.EdgePadding(OverlapID),1));
              double NormalizedDistance = Min(ReceiverDistance/MaxDistance, 1.);
              double Volume = OverlapVolumes(iOverlapping);
              if (DonorGridIDs(Point) < 0) {
                DonorGridIDs(Point) = MGridID;
                DonorNormalizedDistances(Point) = NormalizedDistance;
                DonorVolumes(Point) = Volume;
              } else {
                bool BetterDonor;
                if (std::abs(NormalizedDistance-DonorNormalizedDistances(Point)) > TOLERANCE) {
                  BetterDonor = NormalizedDistance > DonorNormalizedDistances(Point);
                } else {
                  BetterDonor = Volume < DonorVolumes(Point);
                }
                if (BetterDonor) {
                  DonorGridIDs(Point) = MGridID;
                  DonorNormalizedDistances(Point) = NormalizedDistance;
                  DonorVolumes(Point) = Volume;
                }
              }
            }
            ++iOverlapping;
          }
        }
      }
    }
  }

  map<int,distributed_field<bool>> OrphanMasks;
  map<int,long long> NumOrphansForGrid;

  for (int GridID : Domain.LocalGridIDs()) {
    const grid &Grid = Domain.Grid(GridID);
    const range &LocalRange = Grid.LocalRange();
    const distributed_field<bool> &ReceiverMask = ReceiverMasks(GridID);
    const field<int> &DonorGridIDs = DonorGridIDsForLocalGrid(GridID);
    distributed_field<bool> &OrphanMask = OrphanMasks.Insert(GridID, Grid.SharedPartition());
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          OrphanMask(Point) = ReceiverMask(Point) && DonorGridIDs(Point) < 0;
          if (OrphanMask(Point)) {
            Logger.LogWarning(true, "Could not find suitable donor for point (%i,%i,%i) of grid "
              "%s.", Point(0), Point(1), Point(2), Grid.Name());
          }
        }
      }
    }
    OrphanMask.Exchange();
    long long &NumOrphans = NumOrphansForGrid.Insert(GridID);
    NumOrphans = core::CountDistributedMask(OrphanMask);
  }

  for (int GridID : Domain.LocalGridIDs()) {
    if (NumReceiversForGrid(GridID) == 0 && NumOrphansForGrid(GridID) == 0) continue;
    const grid &Grid = Domain.Grid(GridID);
    long long NumExtended = Grid.ExtendedRange().Count();
    const distributed_field<bool> &ReceiverMask = ReceiverMasks(GridID);
    const distributed_field<bool> &OrphanMask = OrphanMasks(GridID);
    auto StateEditHandle = StateComponent.EditState(GridID);
    auto FlagsEditHandle = StateEditHandle->EditFlags();
    distributed_field<state_flags> &Flags = *FlagsEditHandle;
    for (long long l = 0; l < NumExtended; ++l) {
      if (ReceiverMask[l]) {
        Flags[l] = Flags[l] | state_flags::RECEIVER;
      }
      if (OrphanMask[l]) {
        Flags[l] = Flags[l] | state_flags::ORPHAN;
      }
    }
  }

  // Exchanging in reverse, from points to cells

  struct reverse_exchange_m {
    core::recv_map RecvMap;
    core::recv Recv;
  };

  struct reverse_exchange_n {
    core::send_map SendMap;
    core::send Send;
    array<bool> SendBuffer;
  };

  elem_map<int,2,reverse_exchange_m> ReverseExchangeMs;
  elem_map<int,2,reverse_exchange_n> ReverseExchangeNs;

  elem_map<int,2,array<bool>> OverlappingCellDonates;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (Options_.ConnectionType(OverlapID) == connection_type::NONE) continue;
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    reverse_exchange_m &ExchangeM = ReverseExchangeMs.Insert(OverlapID);
    ExchangeM.RecvMap = core::recv_map(OverlapM.DestinationRanks());
    ExchangeM.Recv = core::CreateRecv(Context_, Domain.Comm(), ExchangeM.RecvMap, data_type::BOOL,
      1, 0);
    array<bool> &Donates = OverlappingCellDonates.Insert(OverlapID);
    Donates.Resize({OverlapM.Count()});
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (Options_.ConnectionType(OverlapID) == connection_type::NONE) continue;
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    reverse_exchange_n &ExchangeN = ReverseExchangeNs.Insert(OverlapID);
    ExchangeN.SendMap = core::send_map(OverlapN.SourceRanks());
    ExchangeN.Send = core::CreateSend(Context_, Domain.Comm(), ExchangeN.SendMap, data_type::BOOL,
      1, 0);
    ExchangeN.SendBuffer.Resize({OverlapN.Count()});
  }

  Requests.Reserve(OverlapComponent.LocalOverlapMCount() + OverlapComponent.LocalOverlapNCount());

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (Options_.ConnectionType(OverlapID) == connection_type::NONE) continue;
    reverse_exchange_m &ExchangeM = ReverseExchangeMs(OverlapID);
    core::recv &Recv = ExchangeM.Recv;
    array<bool> &Donates = OverlappingCellDonates(OverlapID);
    request &Request = Requests.Append();
    bool *DonatesData = Donates.Data();
    Request = Recv.Recv(&DonatesData);
  }

  for (auto &OverlapID : OverlapComponent.LocalOverlapNIDs()) {
    if (Options_.ConnectionType(OverlapID) == connection_type::NONE) continue;
    int MGridID = OverlapID(0);
    int NGridID = OverlapID(1);
    const field<int> &DonorGridIDs = DonorGridIDsForLocalGrid(NGridID);
    const overlap_n &OverlapN = OverlapComponent.OverlapN(OverlapID);
    const array<int,2> &Points = OverlapN.Points();
    reverse_exchange_n &ExchangeN = ReverseExchangeNs(OverlapID);
    core::send &Send = ExchangeN.Send;
    array<bool> &SendBuffer = ExchangeN.SendBuffer;
    for (long long iOverlapped = 0; iOverlapped < OverlapN.Count(); ++iOverlapped) {
      tuple<int> Point = {
        Points(0,iOverlapped),
        Points(1,iOverlapped),
        Points(2,iOverlapped)
      };
      SendBuffer(iOverlapped) = DonorGridIDs(Point) == MGridID;
    }
    request &Request = Requests.Append();
    bool *SendBufferData = SendBuffer.Data();
    Request = Send.Send(&SendBufferData);
  }

  WaitAll(Requests);
  Requests.Clear();

  ReverseExchangeMs.Clear();
  ReverseExchangeNs.Clear();

  if (Logger.LoggingDebug()) {
    MPI_Barrier(Domain.Comm());
    for (int GridID : Domain.GridIDs()) {
      if (Domain.GridIsLocal(GridID)) {
        const grid &Grid = Domain.Grid(GridID);
        long long NumReceivers = NumReceiversForGrid(GridID);
        long long NumOrphans = NumOrphansForGrid(GridID);
        if (NumReceivers > 0) {
          std::string NumOrphansString = core::FormatNumber(NumOrphans, "orphans", "orphan");
          Logger.LogDebug(Grid.Comm().Rank() == 0, 3, "%s on grid %s.", NumOrphansString,
            Grid.Name());
        }
      }
      MPI_Barrier(Domain.Comm());
    }
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done choosing donors.");
    Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Creating and filling connectivity data "
      "structures...");
  }

  elem_map<int,2,long long> NumLocalDonorsForGridPair;

  elem_set<int,2> ConnectedGridIDs;

  for (auto &OverlapID : OverlapComponent.LocalOverlapMIDs()) {
    if (Options_.ConnectionType(OverlapID) == connection_type::NONE) continue;
    int MGridID = OverlapID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const overlap_m &OverlapM = OverlapComponent.OverlapM(OverlapID);
    const array<bool> &Donates = OverlappingCellDonates(OverlapID);
    long long &NumLocalDonors = NumLocalDonorsForGridPair.Insert(OverlapID, 0);
    for (long long iOverlapping = 0; iOverlapping < OverlapM.Count(); ++iOverlapping) {
      if (Donates(iOverlapping)) {
        ++NumLocalDonors;
      }
    }
    if (MGrid.Comm().Rank() > 0) {
      MPI_Reduce(&NumLocalDonors, nullptr, 1, MPI_LONG_LONG, MPI_MAX, 0, MGrid.Comm());
    } else {
      long long MaxLocalDonors = NumLocalDonors;
      MPI_Reduce(MPI_IN_PLACE, &MaxLocalDonors, 1, MPI_LONG_LONG, MPI_MAX, 0, MGrid.Comm());
      if (MaxLocalDonors > 0) {
        ConnectedGridIDs.Insert(OverlapID);
      }
    }
  }

  for (int MGridID : Domain.GridIDs()) {
    bool IsMGridRoot = false;
    int MGridRootRank;
    if (Domain.GridIsLocal(MGridID)) {
      const grid &MGrid = Domain.Grid(MGridID);
      IsMGridRoot = MGrid.Comm().Rank() == 0;
      if (IsMGridRoot) {
        MGridRootRank = Domain.Comm().Rank();
      }
    }
    core::BroadcastAnySource(&MGridRootRank, 1, MPI_INT, IsMGridRoot, Domain.Comm());
    for (int NGridID : Domain.GridIDs()) {
      elem<int,2> OverlapID = {MGridID,NGridID};
      if (Options_.ConnectionType(OverlapID) != connection_type::NONE) {
        int Connected;
        if (IsMGridRoot) Connected = ConnectedGridIDs.Contains(OverlapID);
        MPI_Bcast(&Connected, 1, MPI_INT, MGridRootRank, Domain.Comm());
        if (Connected) {
          ConnectedGridIDs.Insert(OverlapID);
        }
      }
    }
  }

  auto ConnectivityComponentEditHandle = Domain.EditComponent<connectivity_component>(
    ConnectivityComponentID_);
  connectivity_component &ConnectivityComponent = *ConnectivityComponentEditHandle;

  ConnectivityComponent.ClearConnectivities();
  ConnectivityComponent.CreateConnectivities(ConnectedGridIDs);

  struct connectivity_m_data {
    long long NumDonors = 0;
    array<int,3> Extents;
    array<double,2> Coords;
    array<int,2> Destinations;
    array<int> DestinationRanks;
    array<long long> Order;
  };

  elem_map<int,2,connectivity_m_data> ConnectivityMDataForLocalDonors;

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityMIDs()) {
    const overlap_m &OverlapM = OverlapComponent.OverlapM(ConnectivityID);
    const array<int,2> &OverlapCells = OverlapM.Cells();
    const array<double,2> &OverlapCoords = OverlapM.Coords();
    const array<int,2> &OverlapDestinations = OverlapM.Destinations();
    const array<int> &OverlapDestinationRanks = OverlapM.DestinationRanks();
    const array<bool> &Donates = OverlappingCellDonates(ConnectivityID);
    long long NumLocalDonors = NumLocalDonorsForGridPair(ConnectivityID);
    connectivity_m_data &ConnectivityMData = ConnectivityMDataForLocalDonors.Insert(ConnectivityID);
    ConnectivityMData.NumDonors = NumLocalDonors;
    ConnectivityMData.Extents.Resize({{2,MAX_DIMS,NumLocalDonors}});
    ConnectivityMData.Coords.Resize({{MAX_DIMS,NumLocalDonors}});
    ConnectivityMData.Destinations.Resize({{MAX_DIMS,NumLocalDonors}});
    ConnectivityMData.DestinationRanks.Resize({NumLocalDonors});
    switch (Options_.ConnectionType(ConnectivityID)) {
    case connection_type::NEAREST:
    case connection_type::LINEAR: {
      long long iDonor = 0;
      for (long long iOverlapping = 0; iOverlapping < OverlapM.Count(); ++iOverlapping) {
        if (!Donates(iOverlapping)) continue;
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          ConnectivityMData.Extents(0,iDim,iDonor) = OverlapCells(iDim,iOverlapping);
          ConnectivityMData.Extents(1,iDim,iDonor) = OverlapCells(iDim,iOverlapping)+2;
        }
        for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
          ConnectivityMData.Extents(0,iDim,iDonor) = 0;
          ConnectivityMData.Extents(1,iDim,iDonor) = 1;
        }
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          ConnectivityMData.Coords(iDim,iDonor) = OverlapCoords(iDim,iOverlapping);
          ConnectivityMData.Destinations(iDim,iDonor) = OverlapDestinations(iDim,iOverlapping);
        }
        ConnectivityMData.DestinationRanks(iDonor) = OverlapDestinationRanks(iOverlapping);
        ++iDonor;
      }
      break;
    }
    case connection_type::CUBIC:
      OVK_DEBUG_ASSERT(false, "Cubic interpolation not yet implemented.");
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      break;
    }
  }

  elem_map<int,2,array<connectivity_m_data>> ConnectivityMSends;
  elem_map<int,2,array<connectivity_m_data>> ConnectivityMRecvs;

  int NumSends = 0;
  int NumRecvs = 0;

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityMIDs()) {
    int MGridID = ConnectivityID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const map<int,partition::neighbor_info> &Neighbors = MGrid.Partition().Neighbors();
    array<connectivity_m_data> &Sends = ConnectivityMSends.Insert(ConnectivityID);
    array<connectivity_m_data> &Recvs = ConnectivityMRecvs.Insert(ConnectivityID);
    Sends.Resize({Neighbors.Count()});
    Recvs.Resize({Neighbors.Count()});
    NumSends += Sends.Count();
    NumRecvs += Recvs.Count();
  }

  array<MPI_Request> MPIRequests;
  MPIRequests.Reserve(NumSends + NumRecvs);

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityMIDs()) {
    int MGridID = ConnectivityID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const map<int,partition::neighbor_info> &Neighbors = MGrid.Partition().Neighbors();
    array<connectivity_m_data> &Recvs = ConnectivityMRecvs(ConnectivityID);
    for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
      connectivity_m_data &Recv = Recvs(iNeighbor);
      MPI_Irecv(&Recv.NumDonors, 1, MPI_LONG_LONG, Neighbors[iNeighbor].Key(), 0, MGrid.Comm(),
        &MPIRequests.Append());
    }
  }

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityMIDs()) {
    int MGridID = ConnectivityID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const cart &Cart = MGrid.Cart();
    const map<int,partition::neighbor_info> &Neighbors = MGrid.Partition().Neighbors();
    const connectivity_m_data &ConnectivityMData = ConnectivityMDataForLocalDonors(ConnectivityID);
    array<connectivity_m_data> &Sends = ConnectivityMSends(ConnectivityID);
    for (long long iDonor = 0; iDonor < ConnectivityMData.NumDonors; ++iDonor) {
      range DonorRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorRange.Begin(iDim) = ConnectivityMData.Extents(0,iDim,iDonor);
        DonorRange.End(iDim) = ConnectivityMData.Extents(1,iDim,iDonor);
      }
      for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
        const partition::neighbor_info &Neighbor = Neighbors[iNeighbor].Value();
        auto MaybeMappedDonorRange = Cart.MapToRange(Neighbor.LocalRange, DonorRange);
        if (MaybeMappedDonorRange) {
          ++Sends(iNeighbor).NumDonors;
        }
      }
    }
    for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
      connectivity_m_data &Send = Sends(iNeighbor);
      MPI_Isend(&Send.NumDonors, 1, MPI_LONG_LONG, Neighbors[iNeighbor].Key(), 0, MGrid.Comm(),
        &MPIRequests.Append());
    }
  }

  MPI_Waitall(MPIRequests.Count(), MPIRequests.Data(), MPI_STATUSES_IGNORE);
  MPIRequests.Clear();

  MPIRequests.Reserve(4*(NumSends + NumRecvs));

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityMIDs()) {
    int MGridID = ConnectivityID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const map<int,partition::neighbor_info> &Neighbors = MGrid.Partition().Neighbors();
    array<connectivity_m_data> &Recvs = ConnectivityMRecvs(ConnectivityID);
    for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
      connectivity_m_data &Recv = Recvs(iNeighbor);
      long long NumDonors = Recv.NumDonors;
      Recv.Extents.Resize({{2,MAX_DIMS,NumDonors}});
      Recv.Coords.Resize({{MAX_DIMS,NumDonors}});
      Recv.Destinations.Resize({{MAX_DIMS,NumDonors}});
      Recv.DestinationRanks.Resize({NumDonors});
      MPI_Irecv(Recv.Extents.Data(), 2*MAX_DIMS*NumDonors, MPI_INT, Neighbors[iNeighbor].Key(), 0,
        MGrid.Comm(), &MPIRequests.Append());
      MPI_Irecv(Recv.Coords.Data(), MAX_DIMS*NumDonors, MPI_DOUBLE, Neighbors[iNeighbor].Key(), 0,
        MGrid.Comm(), &MPIRequests.Append());
      MPI_Irecv(Recv.Destinations.Data(), MAX_DIMS*NumDonors, MPI_INT, Neighbors[iNeighbor].Key(),
        0, MGrid.Comm(), &MPIRequests.Append());
      MPI_Irecv(Recv.DestinationRanks.Data(), NumDonors, MPI_INT, Neighbors[iNeighbor].Key(), 0,
        MGrid.Comm(), &MPIRequests.Append());
    }
  }

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityMIDs()) {
    int MGridID = ConnectivityID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const cart &Cart = MGrid.Cart();
    const map<int,partition::neighbor_info> &Neighbors = MGrid.Partition().Neighbors();
    const connectivity_m_data &ConnectivityMData = ConnectivityMDataForLocalDonors(
      ConnectivityID);
    array<connectivity_m_data> &Sends = ConnectivityMSends(ConnectivityID);
    for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
      connectivity_m_data &Send = Sends(iNeighbor);
      long long NumDonors = Send.NumDonors;
      Send.Extents.Resize({{2,MAX_DIMS,NumDonors}});
      Send.Coords.Resize({{MAX_DIMS,NumDonors}});
      Send.Destinations.Resize({{MAX_DIMS,NumDonors}});
      Send.DestinationRanks.Resize({NumDonors});
      Send.NumDonors = 0;
    }
    for (long long iDonor = 0; iDonor < ConnectivityMData.NumDonors; ++iDonor) {
      range DonorRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorRange.Begin(iDim) = ConnectivityMData.Extents(0,iDim,iDonor);
        DonorRange.End(iDim) = ConnectivityMData.Extents(1,iDim,iDonor);
      }
      for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
        const partition::neighbor_info &Neighbor = Neighbors[iNeighbor].Value();
        auto MaybeMappedDonorRange = Cart.MapToRange(Neighbor.LocalRange, DonorRange);
        if (MaybeMappedDonorRange) {
          const range &MappedDonorRange = *MaybeMappedDonorRange;
          connectivity_m_data &Send = Sends(iNeighbor);
          long long &iSendEntry = Send.NumDonors;
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            Send.Extents(0,iDim,iSendEntry) = MappedDonorRange.Begin(iDim);
            Send.Extents(1,iDim,iSendEntry) = MappedDonorRange.End(iDim);
            Send.Coords(iDim,iSendEntry) = ConnectivityMData.Coords(iDim,iDonor);
            Send.Destinations(iDim,iSendEntry) = ConnectivityMData.Destinations(iDim,iDonor);
          }
          Send.DestinationRanks(iSendEntry) = ConnectivityMData.DestinationRanks(iDonor);
          ++iSendEntry;
        }
      }
    }
    for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
      connectivity_m_data &Send = Sends(iNeighbor);
      long long NumDonors = Send.NumDonors;
      MPI_Isend(Send.Extents.Data(), 2*MAX_DIMS*NumDonors, MPI_INT, Neighbors[iNeighbor].Key(), 0,
        MGrid.Comm(), &MPIRequests.Append());
      MPI_Isend(Send.Coords.Data(), MAX_DIMS*NumDonors, MPI_DOUBLE, Neighbors[iNeighbor].Key(), 0,
        MGrid.Comm(), &MPIRequests.Append());
      MPI_Isend(Send.Destinations.Data(), MAX_DIMS*NumDonors, MPI_INT, Neighbors[iNeighbor].Key(),
        0, MGrid.Comm(), &MPIRequests.Append());
      MPI_Isend(Send.DestinationRanks.Data(), NumDonors, MPI_INT, Neighbors[iNeighbor].Key(), 0,
        MGrid.Comm(), &MPIRequests.Append());
    }
  }

  MPI_Waitall(MPIRequests.Count(), MPIRequests.Data(), MPI_STATUSES_IGNORE);
  MPIRequests.Clear();

  // Not sure if initiating all of the edits up front is necessary anymore...

  struct connectivity_m_edit {
    edit_handle<connectivity_m> Connectivity;
    long long NumConnections = 0;
    edit_handle<array<int,3>> Extents;
    edit_handle<array<double,2>> Coords;
    edit_handle<array<double,3>> InterpCoefs;
    edit_handle<array<int,2>> Destinations;
    edit_handle<array<int>> DestinationRanks;
  };

  struct connectivity_n_edit {
    edit_handle<connectivity_n> Connectivity;
    long long NumConnections = 0;
    edit_handle<array<int,2>> Points;
    edit_handle<array<int,2>> Sources;
    edit_handle<array<int>> SourceRanks;
  };

  elem_map<int,2,connectivity_m_edit> ConnectivityMEdits;
  elem_map<int,2,connectivity_n_edit> ConnectivityNEdits;

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityMIDs()) {
    int MGridID = ConnectivityID(0);
    const grid &MGrid = Domain.Grid(MGridID);
    const map<int,partition::neighbor_info> &Neighbors = MGrid.Partition().Neighbors();
    connectivity_m_data &LocalDonors = ConnectivityMDataForLocalDonors(ConnectivityID);
    array<connectivity_m_data> &Recvs = ConnectivityMRecvs(ConnectivityID);
    long long NumDonors = LocalDonors.NumDonors;
    for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
      connectivity_m_data &Recv = Recvs(iNeighbor);
      NumDonors += Recv.NumDonors;
    }
    int MaxSize;
    switch (Options_.ConnectionType(ConnectivityID)) {
    case connection_type::NEAREST:
    case connection_type::LINEAR:
      MaxSize = 2;
      break;
    case connection_type::CUBIC:
      MaxSize = 4;
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      MaxSize = 2;
      break;
    }
    connectivity_m_edit &Edit = ConnectivityMEdits.Insert(ConnectivityID);
    Edit.Connectivity = ConnectivityComponent.EditConnectivityM(ConnectivityID);
    Edit.NumConnections = NumDonors;
    Edit.Connectivity->Resize(Edit.NumConnections, MaxSize);
    Edit.Extents = Edit.Connectivity->EditExtents();
    Edit.Coords = Edit.Connectivity->EditCoords();
    Edit.InterpCoefs = Edit.Connectivity->EditInterpCoefs();
    Edit.Destinations = Edit.Connectivity->EditDestinations();
    Edit.DestinationRanks = Edit.Connectivity->EditDestinationRanks();
  }

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityNIDs()) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    const distributed_field<bool> &ReceiverMask = ReceiverMasks(NGridID);
    const field<int> &DonorGridIDs = DonorGridIDsForLocalGrid(NGridID);
    long long NumReceivers = 0;
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          if (ReceiverMask(Point) && DonorGridIDs(Point) == MGridID) {
            ++NumReceivers;
          }
        }
      }
    }
    connectivity_n_edit &Edit = ConnectivityNEdits.Insert(ConnectivityID);
    Edit.Connectivity = ConnectivityComponent.EditConnectivityN(ConnectivityID);
    Edit.NumConnections = NumReceivers;
    Edit.Connectivity->Resize(Edit.NumConnections);
    Edit.Points = Edit.Connectivity->EditPoints();
    Edit.Sources = Edit.Connectivity->EditSources();
    Edit.SourceRanks = Edit.Connectivity->EditSourceRanks();
  }

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityMIDs()) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    const grid &MGrid = Domain.Grid(MGridID);
    const map<int,partition::neighbor_info> &Neighbors = MGrid.Partition().Neighbors();
    field_indexer NGridGlobalIndexer(Domain.GridInfo(NGridID).GlobalRange());
    connectivity_m_data &LocalDonors = ConnectivityMDataForLocalDonors(ConnectivityID);
    array<connectivity_m_data> &Recvs = ConnectivityMRecvs(ConnectivityID);
    connectivity_m_edit &Edit = ConnectivityMEdits(ConnectivityID);
    array<long long> DestinationPointIndices({Edit.NumConnections});
    LocalDonors.Order.Resize({LocalDonors.NumDonors});
    long long iDonor = 0;
    for (long long iLocalDonor = 0; iLocalDonor < LocalDonors.NumDonors; ++iLocalDonor) {
      tuple<int> DestinationPoint = {
        LocalDonors.Destinations(0,iLocalDonor),
        LocalDonors.Destinations(1,iLocalDonor),
        LocalDonors.Destinations(2,iLocalDonor)
      };
      DestinationPointIndices(iDonor) = NGridGlobalIndexer.ToIndex(DestinationPoint);
      LocalDonors.Order(iLocalDonor) = iDonor;
      ++iDonor;
    }
    for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
      connectivity_m_data &Recv = Recvs(iNeighbor);
      Recv.Order.Resize({Recv.NumDonors});
      for (long long iRecvDonor = 0; iRecvDonor < Recv.NumDonors; ++iRecvDonor) {
        tuple<int> DestinationPoint = {
          Recv.Destinations(0,iRecvDonor),
          Recv.Destinations(1,iRecvDonor),
          Recv.Destinations(2,iRecvDonor)
        };
        DestinationPointIndices(iDonor) = NGridGlobalIndexer.ToIndex(DestinationPoint);
        Recv.Order(iRecvDonor) = iDonor;
        ++iDonor;
      }
    }
    array<long long> ROrder = ArrayOrder(DestinationPointIndices);
    // Need to use order on LHS since RHS is not one contiguous array
    array<long long> LOrder({Edit.NumConnections});
    for (long long iDonor = 0; iDonor < Edit.NumConnections; ++iDonor) {
      LOrder(ROrder(iDonor)) = iDonor;
    }
    iDonor = 0;
    for (long long iLocalDonor = 0; iLocalDonor < LocalDonors.NumDonors; ++iLocalDonor) {
      long long iOrder = LOrder(iDonor);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        (*Edit.Extents)(0,iDim,iOrder) = LocalDonors.Extents(0,iDim,iLocalDonor);
        (*Edit.Extents)(1,iDim,iOrder) = LocalDonors.Extents(1,iDim,iLocalDonor);
        (*Edit.Coords)(iDim,iOrder) = LocalDonors.Coords(iDim,iLocalDonor);
        (*Edit.Destinations)(iDim,iOrder) = LocalDonors.Destinations(iDim,iLocalDonor);
      }
      // Let exchanger detect source/destination ranks for now -- revisit later (rank containing
      // lower corner may not be the same as overlap source rank; receiver grid needs to be made
      // aware of this)
//       (*Edit.DestinationRanks)(iOrder) = LocalDonors.DestinationRanks(iLocalDonor);
      ++iDonor;
    }
    for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
      connectivity_m_data &Recv = Recvs(iNeighbor);
      for (long long iRecvDonor = 0; iRecvDonor < Recv.NumDonors; ++iRecvDonor) {
        long long iOrder = LOrder(iDonor);
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          (*Edit.Extents)(0,iDim,iOrder) = Recv.Extents(0,iDim,iRecvDonor);
          (*Edit.Extents)(1,iDim,iOrder) = Recv.Extents(1,iDim,iRecvDonor);
          (*Edit.Coords)(iDim,iOrder) = Recv.Coords(iDim,iRecvDonor);
          (*Edit.Destinations)(iDim,iOrder) = Recv.Destinations(iDim,iRecvDonor);
        }
//         (*Edit.DestinationRanks)(iOrder) = Recv.DestinationRanks(iRecvDonor);
        ++iDonor;
      }
    }
    switch (Options_.ConnectionType(ConnectivityID)) {
    case connection_type::NEAREST:
      for (iDonor = 0; iDonor < Edit.NumConnections; ++iDonor) {
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          (*Edit.InterpCoefs)(iDim,0,iDonor) = (*Edit.Coords)(iDim,iDonor) <= 0.5 ? 1. : 0.;
          (*Edit.InterpCoefs)(iDim,1,iDonor) = (*Edit.Coords)(iDim,iDonor) <= 0.5 ? 0. : 1.;
        }
      }
      break;
    case connection_type::LINEAR:
      for (iDonor = 0; iDonor < Edit.NumConnections; ++iDonor) {
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          elem<double,2> Coefs = core::LagrangeInterpLinear((*Edit.Coords)(iDim,iDonor));
          (*Edit.InterpCoefs)(iDim,0,iDonor) = Coefs(0);
          (*Edit.InterpCoefs)(iDim,1,iDonor) = Coefs(1);
        }
      }
      break;
    case connection_type::CUBIC:
      for (iDonor = 0; iDonor < Edit.NumConnections; ++iDonor) {
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          elem<double,4> Coefs = core::LagrangeInterpCubic((*Edit.Coords)(iDim,iDonor));
          (*Edit.InterpCoefs)(iDim,0,iDonor) = Coefs(0);
          (*Edit.InterpCoefs)(iDim,1,iDonor) = Coefs(1);
          (*Edit.InterpCoefs)(iDim,2,iDonor) = Coefs(2);
          (*Edit.InterpCoefs)(iDim,3,iDonor) = Coefs(3);
        }
      }
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      break;
    }
  }

  for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityNIDs()) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    const grid &NGrid = Domain.Grid(NGridID);
    const range &LocalRange = NGrid.LocalRange();
    const overlap_n &OverlapN = OverlapComponent.OverlapN(ConnectivityID);
    const field<bool> &BaseOverlapMask = OverlapN.Mask();
    const array<int,2> &OverlapSources = OverlapN.Sources();
//     const array<int> &OverlapSourceRanks = OverlapN.SourceRanks();
    const distributed_field<bool> &ReceiverMask = ReceiverMasks(NGridID);
    const field<int> &DonorGridIDs = DonorGridIDsForLocalGrid(NGridID);
    connectivity_n_edit &Edit = ConnectivityNEdits(ConnectivityID);
    long long iOverlapping = 0;
    long long iReceiver = 0;
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          if (BaseOverlapMask(Point)) {
            if (ReceiverMask(Point) && DonorGridIDs(Point) == MGridID) {
              for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
                (*Edit.Points)(iDim,iReceiver) = Point(iDim);
                (*Edit.Sources)(iDim,iReceiver) = OverlapSources(iDim,iOverlapping);
              }
//               (*Edit.SourceRanks)(iReceiver) = OverlapSourceRanks(iOverlapping);
              ++iReceiver;
            }
            ++iOverlapping;
          }
        }
      }
    }
  }

  ConnectivityMEdits.Clear();
  ConnectivityNEdits.Clear();

  MPI_Barrier(Domain.Comm());

  Logger.LogDebug(Domain.Comm().Rank() == 0, 2, "Done creating and filling connectivity data "
    "structures.");
  Logger.LogDebug(Domain.Comm().Rank() == 0, 1, "Done generating connectivity data.");

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
