// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/OverlapAccel.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Box.hpp"
#include "ovk/core/Field.hpp"
#include "ovk/core/FieldOps.hpp"
#include "ovk/core/GeometryManipulator.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Optional.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Tuple.hpp"

#include <limits>
#include <memory>
#include <utility>

namespace ovk {
namespace core {

namespace {
struct compute_cell_bounds {
  template <typename T> void operator()(const T &Manipulator, int NumDims, const range &CellRange,
    const array<field_view<const double>> &Coords, field_view<const bool> CellMask, double
    MaxTolerance, field_view<box> CellBounds) const {
    tuple<double> ScaleFactor = MakeUniformTuple<double>(NumDims, 1.+2.*MaxTolerance, 1.);
    for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
      for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
        for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
          tuple<int> Cell = {i,j,k};
          if (CellMask(Cell)) {
            CellBounds(Cell) = ScaleBox(Manipulator.CellBounds(Coords, Cell), ScaleFactor);
          } else {
            CellBounds(Cell) = MakeEmptyBox(NumDims);
          }
        }
      }
    }
  }
};

struct compute_cell_volumes {
  template <typename T> void operator()(const T &Manipulator, int NumDims, const range &CellRange,
    const array<field_view<const double>> &Coords, field_view<const bool> CellMask,
    field_view<double> CellVolumes) const {
    for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
      for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
        for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
          tuple<int> Cell = {i,j,k};
          if (CellMask(Cell)) {
            CellVolumes(Cell) = Manipulator.CellVolume(Coords, Cell);
          } else {
            CellVolumes(Cell) = 0.;
          }
        }
      }
    }
  }
};
}

overlap_accel::overlap_accel(geometry_type GeometryType, int NumDims):
  GeometryType_(GeometryType),
  NumDims_(NumDims),
  GeometryManipulator_(GeometryType, NumDims),
  CellRange_(MakeEmptyRange(NumDims)),
  CellIndexer_(CellRange_),
  Coords_({MAX_DIMS}),
  Bounds_(MakeEmptyBox(NumDims))
{}

overlap_accel::overlap_accel(geometry_type GeometryType, int NumDims, const range &CellRange,
  array_view<const field_view<const double>> Coords, field_view<const bool> CellMask, double
  MaxTolerance, long long NumCellsLeaf, double MaxNodeUnoccupiedVolume, double
  MaxNodeCellVolumeVariation, double BinScale):
  GeometryType_(GeometryType),
  NumDims_(NumDims),
  GeometryManipulator_(GeometryType, NumDims),
  CellRange_(CellRange),
  CellIndexer_(CellRange),
  Coords_({MAX_DIMS})
{

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Coords_(iDim) = Coords(iDim);
  }

  field<box> CellBounds(CellRange);

  GeometryManipulator_.Apply(compute_cell_bounds(), NumDims, CellRange, Coords_, CellMask,
    MaxTolerance, CellBounds);

  long long NumCells = CellRange.Count();

  Bounds_ = MakeEmptyBox(NumDims);
  for (long long iCell = 0; iCell < NumCells; ++iCell) {
    Bounds_ = UnionBoxes(Bounds_, CellBounds[iCell]);
  }

  if (!Bounds_.Empty()) {

    field<double> CellVolumes(CellRange);

    GeometryManipulator_.Apply(compute_cell_volumes(), NumDims, CellRange, Coords_, CellMask,
      CellVolumes);

    long long NumContainedCells = 0;
    for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
      for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
        for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
          if (CellMask(i,j,k)) ++NumContainedCells;
        }
      }
    }

    array<long long> CellIndices({NumContainedCells});

    long long iContainedCell = 0;
    for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
      for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
        for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
          tuple<int> Cell = {i,j,k};
          if (CellMask(Cell)) {
            CellIndices(iContainedCell) = CellIndexer_.ToIndex(Cell);
            ++iContainedCell;
          }
        }
      }
    }

    int MaxDepth = 0;
    while ((NumContainedCells/NumCellsLeaf >> MaxDepth) != 0) {
      MaxDepth += 1;
    }

    Root_ = CreateNode_(Bounds_, CellBounds, CellVolumes, CellIndices, 0, MaxDepth, NumCellsLeaf,
      MaxNodeUnoccupiedVolume, MaxNodeCellVolumeVariation, BinScale);

  }

}

overlap_accel::node overlap_accel::CreateNode_(const box &AccelBounds, const field<box> &CellBounds,
  const field<double> &CellVolumes, const array<long long> &ContainedCellIndices, int Depth, int
  MaxDepth, long long NumCellsLeaf, double MaxUnoccupiedVolume, double MaxCellVolumeVariation,
  double BinScale) const {

  constexpr long long MAX_HASH_SIZE = 1L << 26;

  long long NumContainedCells = ContainedCellIndices.Count();

  box NodeBounds = MakeEmptyBox(NumDims_);
  for (long long iCell : ContainedCellIndices) {
    NodeBounds = UnionBoxes(NodeBounds, CellBounds[iCell]);
  }

  bool LeafNode;

  if (Depth == MaxDepth || NumContainedCells <= NumCellsLeaf) {

    LeafNode = true;

  } else {

    double TotalCellVolume = 0.;
    for (long long iCell : ContainedCellIndices) {
      TotalCellVolume += CellVolumes[iCell];
    }
    double MeanCellVolume = TotalCellVolume/double(NumContainedCells);

    double CellVolumeVariation = 0.;
    for (long long iCell : ContainedCellIndices) {
      double Deviation = CellVolumes[iCell] - MeanCellVolume;
      CellVolumeVariation += Deviation*Deviation;
    }
    CellVolumeVariation = std::sqrt(CellVolumeVariation/double(NumContainedCells-1))/MeanCellVolume;

    auto BoxVolume = [&](const box &Box) -> double {
      double Volume = 1.;
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        Volume *= Box.Size(iDim);
      }
      return Volume;
    };

    double UnoccupiedVolume = Max(BoxVolume(NodeBounds)-TotalCellVolume, 0.);

    bool OccupiedEnough = UnoccupiedVolume <= MaxUnoccupiedVolume*BoxVolume(AccelBounds);
    bool UniformEnough = CellVolumeVariation <= MaxCellVolumeVariation;
    bool NotTooBig = NumContainedCells <= MAX_HASH_SIZE;

    LeafNode = OccupiedEnough and UniformEnough and NotTooBig;

  }

  node Node;

  if (!LeafNode) {

    Node.SplitDim = 1;
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      if (NodeBounds.Size(iDim) > NodeBounds.Size(Node.SplitDim)) {
        Node.SplitDim = iDim;
      }
    }

    Node.Split = 0.;
    for (long long iCell : ContainedCellIndices) {
      const box &Bounds = CellBounds[iCell];
      Node.Split += 0.5*(Bounds.Begin(Node.SplitDim) + Bounds.End(Node.SplitDim));
    }
    Node.Split /= Max(double(NumContainedCells), 1.);

    long long NumLeftChildCells = 0;
    long long NumRightChildCells = 0;

    for (long long iCell : ContainedCellIndices) {
      const box &Bounds = CellBounds[iCell];
      if (Bounds.Begin(Node.SplitDim) <= Node.Split) ++NumLeftChildCells;
      if (Bounds.End(Node.SplitDim) >= Node.Split) ++NumRightChildCells;
    }

    array<long long> LeftChildIndices, RightChildIndices;
    LeftChildIndices.Reserve(NumLeftChildCells);
    RightChildIndices.Reserve(NumRightChildCells);

    for (long long iCell : ContainedCellIndices) {
      const box &Bounds = CellBounds[iCell];
      if (Bounds.Begin(Node.SplitDim) <= Node.Split) LeftChildIndices.Append(iCell);
      if (Bounds.End(Node.SplitDim) >= Node.Split) RightChildIndices.Append(iCell);
    }

    Node.LeftChild.reset(new node());
    Node.RightChild.reset(new node());

    *Node.LeftChild = CreateNode_(AccelBounds, CellBounds, CellVolumes, LeftChildIndices,
      Depth+1, MaxDepth, NumCellsLeaf, MaxUnoccupiedVolume, MaxCellVolumeVariation, BinScale);
    *Node.RightChild = CreateNode_(AccelBounds, CellBounds, CellVolumes, RightChildIndices,
      Depth+1, MaxDepth, NumCellsLeaf, MaxUnoccupiedVolume, MaxCellVolumeVariation, BinScale);

  } else {

    tuple<double> MeanCellBoundsSize = {0.,0.,0.};
    for (long long iCell : ContainedCellIndices) {
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        MeanCellBoundsSize(iDim) += CellBounds[iCell].Size(iDim);
      }
    }
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      MeanCellBoundsSize(iDim) /= Max(double(NumContainedCells), 1.);
    }

    tuple<double> DesiredBinSize;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      DesiredBinSize(iDim) = BinScale * MeanCellBoundsSize(iDim);
    }

    tuple<int> NumBins = {1,1,1};
    if (NumContainedCells > 1) {
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        NumBins(iDim) = Max(int(std::ceil(NodeBounds.Size(iDim)/DesiredBinSize(iDim))), 1);
      }
    }

    array<box> ContainedCellBounds({NumContainedCells});

    for (long long iContainedCell = 0; iContainedCell < NumContainedCells; ++iContainedCell) {
      long long iCell = ContainedCellIndices(iContainedCell);
      ContainedCellBounds(iContainedCell) = CellBounds[iCell];
    }

    Node.CellIndices = ContainedCellIndices;
    Node.Hash.reset(new bounding_box_hash(NumDims_, NumBins, ContainedCellBounds));

  }

  return Node;

}

optional<tuple<int>> overlap_accel::FindCell(const tuple<double> &PointCoords, double Tolerance)
  const {

  optional<tuple<int>> MaybeCell;

  if (Bounds_.Contains(PointCoords)) {
    MaybeCell = FindCellInNode_(*Root_, PointCoords, Tolerance);
  }

  return MaybeCell;

}

namespace {
struct find_cell_in_bin {
  template <typename T> void operator()(const T &Manipulator, int NumDims, const array<field_view<
    const double>> &Coords, const array_view<const long long> &BinCells, const array<long long>
    &NodeContainedCells, const field_indexer &CellIndexer, const tuple<double> &PointCoords, double
    Tolerance, optional<tuple<int>> &MaybeCell) const {

    // Try to find overlap without tolerance first
    for (long long iContainedCell : BinCells) {
      long long iCell = NodeContainedCells(iContainedCell);
      tuple<int> Cell = CellIndexer.ToTuple(iCell);
      if (Manipulator.OverlapsCell(Coords, 0., Cell, PointCoords)) {
        MaybeCell = Cell;
        break;
      }
    }

    // If that fails, try to project onto edge of nearest cell
    if (!MaybeCell) {
      double BestDisplacementSq = std::numeric_limits<double>::max();
      for (long long iContainedCell : BinCells) {
        long long iCell = NodeContainedCells(iContainedCell);
        tuple<int> Cell = CellIndexer.ToTuple(iCell);
        if (!Manipulator.OverlapsCell(Coords, Tolerance, Cell, PointCoords)) continue;
        auto MaybeCellCoords = Manipulator.CoordsInCell(Coords, Cell, PointCoords);
        if (!MaybeCellCoords) continue;
        const tuple<double> &CellCoords = *MaybeCellCoords;
        double DisplacementSq = 0.;
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          double ClampOffset = 0.;
          if (CellCoords(iDim) > 1.) {
            ClampOffset = 1.-CellCoords(iDim);
          } else if (CellCoords(iDim) < 0.) {
            ClampOffset = -CellCoords(iDim);
          }
          DisplacementSq += ClampOffset*ClampOffset;
        }
        if (DisplacementSq < BestDisplacementSq) {
          MaybeCell = Cell;
          BestDisplacementSq = DisplacementSq;
        }
      }
    }

  }
};
}

optional<tuple<int>> overlap_accel::FindCellInNode_(const node &Node, const tuple<double>
  &PointCoords, double Tolerance) const {

  optional<tuple<int>> MaybeCell;

  bool LeafNode = Node.Hash != nullptr;

  if (!LeafNode) {

    if (PointCoords(Node.SplitDim) <= Node.Split) {
      MaybeCell = FindCellInNode_(*Node.LeftChild, PointCoords, Tolerance);
    } else {
      MaybeCell = FindCellInNode_(*Node.RightChild, PointCoords, Tolerance);
    }

  } else {

    long long iBin = Node.Hash->MapToBin(PointCoords);

    if (iBin >= 0) {
      array_view<const long long> BinCells = Node.Hash->RetrieveBin(iBin);
      GeometryManipulator_.Apply(find_cell_in_bin(), NumDims_, Coords_, BinCells, Node.CellIndices,
        CellIndexer_, PointCoords, Tolerance, MaybeCell);
    }

  }

  return MaybeCell;

}

}}
