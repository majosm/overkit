// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_OVERLAP_ACCEL_HPP_INCLUDED
#define OVK_CORE_OVERLAP_ACCEL_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Box.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/GeometryManipulator.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/RegionHash.hpp>
#include <ovk/core/Tuple.hpp>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

class overlap_accel {

public:

  overlap_accel(geometry_type GeometryType, int NumDims, const range &CellRange, array_view<const
    field_view<const double>> Coords, field_view<const bool> CellMask, double MaxTolerance,
    long long NumCellsLeaf, double MaxNodeUnoccupiedVolume, double MaxNodeCellVolumeVariation,
    double BinScale);

  optional<tuple<int>> FindCell(const tuple<double> &PointCoords, double Tolerance) const;

private:

  using bounding_box_hash = region_hash<double>;

  struct node {
    int SplitDim;
    double Split;
    std::unique_ptr<node> LeftChild;
    std::unique_ptr<node> RightChild;
    array<long long> CellIndices;
    std::unique_ptr<bounding_box_hash> Hash;
  };

  geometry_type GeometryType_;
  int NumDims_;

  geometry_manipulator GeometryManipulator_;

  range CellRange_;
  field_indexer CellIndexer_;

  array<field_view<const double>> Coords_;

  box Bounds_;

  optional<node> Root_;

  node CreateNode_(const box &AccelBounds, const field<box> &CellBounds, const field<double>
    &CellVolumes, const array<long long> &ContainedCellIndices, int Depth, int MaxDepth, long long
    NumCellsLeaf, double MaxUnoccupiedVolume, double MaxCellVolumeVariation, double BinScale) const;

  optional<tuple<int>> FindCellInNode_(const node &Node, const tuple<double> &PointCoords, double
    Tolerance) const;

};

}}

#endif
