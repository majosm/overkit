// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISTRIBUTED_REGION_HASH_HPP_INCLUDED
#define OVK_CORE_DISTRIBUTED_REGION_HASH_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Box.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/CommunicationOps.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/HashableRegionTraits.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Math.hpp>
#include <ovk/core/MPISerializableTraits.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <cmath>
#include <cstring>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

template <typename RegionType> class distributed_region_hash;

template <typename RegionType> class distributed_region_data {

public:

  using region_type = RegionType;

  const region_type &Region() const { return Region_; }

  int Rank() const { return Rank_; }

private:

  region_type Region_;
  int Rank_;

  friend class distributed_region_hash<RegionType>;

};

template <typename RegionType> class distributed_region_hash_retrieved_bins {

public:

  using region_data = distributed_region_data<RegionType>;

  const region_data &RegionData(int iRegion) const { return RegionData_(iRegion); }

  array_view<const int> BinRegionIndices(int iBin) const {
    const interval<long long> &Interval = BinRegionIndicesIntervals_(iBin);
    return {BinRegionIndices_.Data(Interval.Begin()), Interval};
  }

private:

  array<region_data> RegionData_;
  map<int,interval<long long>> BinRegionIndicesIntervals_;
  array<int> BinRegionIndices_;

  friend class distributed_region_hash<RegionType>;

};

template <typename RegionType> class distributed_region_hash {

public:

  static_assert(IsHashableRegion<RegionType>(), "Invalid region type (not hashable).");
  static_assert(IsMPISerializable<RegionType>(), "Invalid region type (not MPI-serializable).");

  using region_type = RegionType;
  using region_traits = hashable_region_traits<region_type>;
  using coord_type = typename region_traits::coord_type;
  using mpi_traits = mpi_serializable_traits<region_type>;

  static_assert(std::is_same<coord_type, int>::value || std::is_same<coord_type, double>::value,
    "Coord type must be int or double.");

  using extents_type = interval<coord_type,MAX_DIMS>;

  using region_data = distributed_region_data<region_type>;
  using retrieved_bins = distributed_region_hash_retrieved_bins<region_type>;

  distributed_region_hash(int NumDims, comm_view Comm);
  distributed_region_hash(int NumDims, comm_view Comm, array_view<const region_type> LocalRegions);

  distributed_region_hash(const distributed_region_hash &Other) = delete;
  distributed_region_hash(distributed_region_hash &&Other) noexcept = default;

  distributed_region_hash &operator=(const distributed_region_hash &Other) = delete;
  distributed_region_hash &operator=(distributed_region_hash &&Other) noexcept = default;

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

  array<region_data> RegionData_;
  range BinRange_;
  field<int> NumRegionsPerBin_;
  field<long long> BinRegionIndicesStarts_;
  array<int> BinRegionIndices_;

  template <hashable_region_maps_to MapsTo> struct maps_to_tag {};

  template <typename IndexerType> set<typename IndexerType::index_type> MapToBins_(const range
    &BinRange, const IndexerType &BinIndexer, const tuple<coord_type> &LowerCorner, const
    tuple<coord_type> &BinSize, const region_type &Region, maps_to_tag<
    hashable_region_maps_to::SET>) const;

  template <typename IndexerType> set<typename IndexerType::index_type> MapToBins_(const range
    &BinRange, const IndexerType &BinIndexer, const tuple<coord_type> &LowerCorner, const
    tuple<coord_type> &BinSize, const region_type &Region, maps_to_tag<
    hashable_region_maps_to::RANGE>) const;

  template <typename T> struct coord_type_tag {};

  static interval<int,MAX_DIMS> MakeEmptyExtents_(int NumDims, coord_type_tag<int>);
  static interval<double,MAX_DIMS> MakeEmptyExtents_(int NumDims, coord_type_tag<double>);

  static interval<int,MAX_DIMS> UnionExtents_(const interval<int,MAX_DIMS> &Left, const
    interval<int,MAX_DIMS> &Right);
  static interval<double,MAX_DIMS> UnionExtents_(const interval<double,MAX_DIMS> &Left, const
    interval<double,MAX_DIMS> &Right);

  static tuple<int> BinDecomp_(int NumDims, const extents_type &Extents, int MaxBins);

  static tuple<int> GetBinSize_(const interval<int,MAX_DIMS> &Extents, const tuple<int> &NumBins);
  static tuple<double> GetBinSize_(const interval<double,MAX_DIMS> &Extents, const tuple<int>
    &NumBins);

};

}}

#include <ovk/core/DistributedRegionHash.inl>

#endif
