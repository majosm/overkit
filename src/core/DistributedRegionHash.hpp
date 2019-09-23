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
#include <ovk/core/RegionTraits.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <type_traits>

namespace ovk {
namespace core {

template <typename CoordType> struct distributed_region_data {
  using region_type = typename std::conditional<std::is_same<CoordType, double>::value, box,
    range>::type;
  region_type Region;
  int Rank;
  int Tag;
};

template <typename RegionType> struct distributed_region_hash_retrieved_bins {
  using region_data = distributed_region_data<RegionType>;
  array<region_data> RegionData;
  map<int,interval<long long>> BinRegionIndicesIntervals;
  array<int> BinRegionIndices;
};

template <typename CoordType> class distributed_region_hash {

public:

  static_assert(std::is_same<CoordType, int>::value || std::is_same<CoordType, double>::value,
    "Coord type must be int or double.");

  using coord_type = CoordType;
  using region_type = typename std::conditional<std::is_same<coord_type, double>::value, box,
    range>::type;
  using traits = region_traits<region_type>;

  using region_data = distributed_region_data<coord_type>;
  using retrieved_bins = distributed_region_hash_retrieved_bins<coord_type>;

  distributed_region_hash(int NumDims, comm_view Comm);
  distributed_region_hash(int NumDims, comm_view Comm, int NumLocalRegions, array_view<const
    region_type> LocalRegionExtents, array_view<const int> LocalRegionTags);

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

  region_type GlobalExtents_;

  range ProcRange_;
  range_indexer<int> ProcIndexer_;
  tuple<coord_type> ProcSize_;

  array<int> ProcToBinMultipliers_;

  array<region_data> RegionData_;
  range BinRange_;
  field<int> NumRegionsPerBin_;
  field<long long> BinRegionIndicesStarts_;
  array<int> BinRegionIndices_;

  static tuple<int> BinDecomp_(int NumDims, const region_type &GlobalExtents, int MaxBins);

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
