// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISTRIBUTED_REGION_HASH_HPP_INCLUDED
#define OVK_CORE_DISTRIBUTED_REGION_HASH_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Box.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Misc.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <cmath>
#include <memory>
#include <numeric>
#include <type_traits>

namespace ovk {
namespace core {

namespace distributed_region_internal {
template <typename CoordType> struct region_traits;
template <> struct region_traits<int> {
  using extents_type = range;
  static range MakeEmptyExtents(int NumDims) {
    return MakeEmptyRange(NumDims);
  }
  static range UnionExtents(const range &Left, const range &Right) {
    return UnionRanges(Left, Right);
  }
  static range IntersectExtents(const range &Left, const range &Right) {
    return IntersectRanges(Left, Right);
  }
  static tuple<int> GetExtentsLowerCorner(const range &Extents) {
    return Extents.Begin();
  }
  static tuple<int> GetExtentsUpperCorner(const range &Extents) {
    tuple<int> UpperCorner;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      UpperCorner(iDim) = Extents.End(iDim)-1;
    }
    return UpperCorner;
  }
};
template <> struct region_traits<double> {
  using extents_type = box;
  static box MakeEmptyExtents(int NumDims) {
    return MakeEmptyBox(NumDims);
  }
  static box UnionExtents(const box &Left, const box &Right) {
    return UnionBoxes(Left, Right);
  }
  static box IntersectExtents(const box &Left, const box &Right) {
    return IntersectBoxes(Left, Right);
  }
  static tuple<double> GetExtentsLowerCorner(const box &Extents) {
    return Extents.Begin();
  }
  static tuple<double> GetExtentsUpperCorner(const box &Extents) {
    return Extents.End();
  }
};
}

template <typename CoordType> struct distributed_region_data {
  using traits = distributed_region_internal::region_traits<CoordType>;
  using extents_type = typename traits::extents_type;
  int Rank;
  extents_type Extents;
  int Tag;
};

template <typename CoordType> class distributed_region_hash;

template <typename CoordType> class distributed_region_hash_bin {

public:

  static_assert(std::is_same<CoordType, int>::value || std::is_same<CoordType, double>::value,
    "Coord type must be int or double.");

  using coord_type = CoordType;
  using traits = distributed_region_internal::region_traits<CoordType>;
  using extents_type = typename traits::extents_type;

  using region = distributed_region_data<CoordType>;

  distributed_region_hash_bin():
    Index_(-1),
    Extents_(traits::MakeEmptyExtents(2))
  {}

  explicit operator bool() const { return Index_ >= 0; }

  const extents_type &Extents() const { return Extents_; }

  const array<region> &Regions() const { return Regions_; }

  const region &Region(int iRegion) const { return Regions_(iRegion); }
  region &Region(int iRegion) { return Regions_(iRegion); }

private:

  int Index_;

  extents_type Extents_;

  array<region> Regions_;

  distributed_region_hash_bin(int Index, const extents_type &Extents, int NumRegions):
    Index_(Index),
    Extents_(Extents),
    Regions_({NumRegions})
  {}

  friend class distributed_region_hash<coord_type>;

};

template <typename CoordType> class distributed_region_hash {

public:

  static_assert(std::is_same<CoordType, int>::value || std::is_same<CoordType, double>::value,
    "Coord type must be int or double.");

  using coord_type = CoordType;
  using traits = distributed_region_internal::region_traits<CoordType>;
  using extents_type = typename traits::extents_type;

  using bin = distributed_region_hash_bin<CoordType>;

  using region = distributed_region_data<CoordType>;

  distributed_region_hash(int NumDims, comm_view Comm);
  distributed_region_hash(int NumDims, comm_view Comm, int NumLocalRegions, array_view<const
    extents_type> LocalRegionExtents, array_view<const int> LocalRegionTags);

  distributed_region_hash(const distributed_region_hash &Other) = delete;
  distributed_region_hash(distributed_region_hash &&Other) noexcept = default;

  distributed_region_hash &operator=(const distributed_region_hash &Other) = delete;
  distributed_region_hash &operator=(distributed_region_hash &&Other) noexcept = default;

  tuple<int> MapPointToBin(const tuple<coord_type> &Point) const;
  range MapExtentsToBinRange(const extents_type &Extents) const;

  const range_indexer<int> &BinIndexer() const { return BinIndexer_; }

  void RetrieveBins(map<int,bin> &Bins) const;

  bool HasBin() const { return static_cast<bool>(Bin_); }

  const bin &Bin() const { return Bin_; }
  bin &Bin() { return Bin_; }

private:

  int NumDims_;

  comm_view Comm_;

  extents_type GlobalExtents_;

  range BinRange_;
  range_indexer<int> BinIndexer_;
  tuple<coord_type> BinSize_;

  bin Bin_;

  static tuple<int> BinDecomp_(int NumDims, const extents_type &GlobalExtents, int MaxBins);

  static tuple<int> GetBinSize_(const range &GlobalExtents, const tuple<int> &NumBins);
  static tuple<double> GetBinSize_(const box &GlobalExtents, const tuple<int> &NumBins);

  static tuple<int> MapToUniformCell_(int NumDims, const tuple<int> &Origin, const tuple<int>
    &CellSize, const tuple<int> &Point);
  static tuple<int> MapToUniformCell_(int NumDims, const tuple<double> &Origin, const tuple<double>
    &CellSize, const tuple<double> &Point);

};

}}

#include <ovk/core/DistributedRegionHash.inl>

#endif
