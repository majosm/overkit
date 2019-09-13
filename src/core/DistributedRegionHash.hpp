// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISTRIBUTED_REGION_HASH_HPP_INCLUDED
#define OVK_CORE_DISTRIBUTED_REGION_HASH_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Box.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

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

template <typename CoordType> struct distributed_region_hash_retrieved_bins {
  using region = distributed_region_data<CoordType>;
  array<region> Regions;
  map<int,interval<long long>> BinRegionIndicesIntervals;
  array<int> BinRegionIndices;
};

template <typename CoordType> class distributed_region_hash {

public:

  static_assert(std::is_same<CoordType, int>::value || std::is_same<CoordType, double>::value,
    "Coord type must be int or double.");

  using coord_type = CoordType;
  using traits = distributed_region_internal::region_traits<coord_type>;
  using extents_type = typename traits::extents_type;

  using region = distributed_region_data<coord_type>;
  using retrieved_bins = distributed_region_hash_retrieved_bins<coord_type>;

  distributed_region_hash(int NumDims, comm_view Comm);
  distributed_region_hash(int NumDims, comm_view Comm, int NumLocalRegions, array_view<const
    extents_type> LocalRegionExtents, array_view<const int> LocalRegionTags);

  // Can't define these here due to issues with GCC < 6.3 and Intel < 17
//   distributed_region_hash(const distributed_region_hash &Other) = delete;
//   distributed_region_hash(distributed_region_hash &&Other) noexcept = default;

//   distributed_region_hash &operator=(const distributed_region_hash &Other) = delete;
//   distributed_region_hash &operator=(distributed_region_hash &&Other) noexcept = default;

  elem<int,2> MapToBin(const tuple<coord_type> &Point) const;

  map<int,retrieved_bins> RetrieveBins(array_view<const elem<int,2>> BinIDs) const;

private:

  int NumDims_;

  comm_view Comm_;

  extents_type GlobalExtents_;

  range ProcRange_;
  range_indexer<int> ProcIndexer_;
  tuple<coord_type> ProcSize_;

  array<int> ProcToBinMultipliers_;

  array<region> Regions_;
  range BinRange_;
  field<int> NumRegionsPerBin_;
  field<long long> BinRegionIndicesStarts_;
  array<int> BinRegionIndices_;

  static tuple<int> BinDecomp_(int NumDims, const extents_type &GlobalExtents, int MaxBins);

  static tuple<int> GetBinSize_(const range &GlobalExtents, const tuple<int> &NumBins);
  static tuple<double> GetBinSize_(const box &GlobalExtents, const tuple<int> &NumBins);

  static tuple<int> MapToUniformCell_(int NumDims, const tuple<int> &Origin, const tuple<int>
    &CellSize, const tuple<int> &Point);
  static tuple<int> MapToUniformCell_(int NumDims, const tuple<double> &Origin, const tuple<double>
    &CellSize, const tuple<double> &Point);

};

extern template class distributed_region_hash<int>;
extern template class distributed_region_hash<double>;

}}

#endif
