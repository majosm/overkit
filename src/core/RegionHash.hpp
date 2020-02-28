// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_REGION_HASH_HPP_INCLUDED
#define OVK_CORE_REGION_HASH_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Box.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/HashableRegionTraits.hpp>
#include <ovk/core/Math.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/Tuple.hpp>

#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

template <typename RegionType> class region_hash {

public:

  static_assert(IsHashableRegion<RegionType>(), "Invalid region type (not hashable).");

  using region_type = RegionType;
  using region_traits = hashable_region_traits<region_type>;
  using coord_type = typename region_traits::coord_type;

  static_assert(std::is_same<coord_type, int>::value || std::is_same<coord_type, double>::value,
    "Coord type must be int or double.");

  using extents_type = interval<coord_type,MAX_DIMS>;

  explicit region_hash(int NumDims);
  region_hash(int NumDims, const tuple<int> &NumBins, array_view<const region_type> Regions);

  long long MapToBin(const tuple<coord_type> &Point) const;

  array_view<const long long> RetrieveBin(long long iBin) const;

  int Dimension() const { return NumDims_; }

  const range &BinRange() const { return BinRange_; }

  const extents_type &Extents() const { return Extents_; }

private:

  int NumDims_;
  range BinRange_;
  field_indexer BinIndexer_;
  extents_type Extents_;
  tuple<double> BinSize_;
  array<long long> BinRegionIndicesStarts_;
  array<long long> BinRegionIndices_;

  void MapToBins_(const region_type &Region, range &Bins) const;
  void MapToBins_(const region_type &Region, set<long long> &Bins) const;

  void AccumulateBinRegionCounts_(const range &Bins, field<long long> &NumRegionsInBin) const;
  void AccumulateBinRegionCounts_(const set<long long> &Bins, field<long long> &NumRegionsInBin)
    const;

  void AddToBins_(int iRegion, const range &Bins, const array<long long> &BinRegionIndicesStarts,
    array<long long> &BinRegionIndices, field<long long> &NumRegionsAddedToBin) const;
  void AddToBins_(int iRegion, const set<long long> &Bins, const array<long long>
    &BinRegionIndicesStarts, array<long long> &BinRegionIndices, field<long long>
    &NumRegionsAddedToBin) const ;

  template <typename T> struct coord_type_tag {};

  interval<int,MAX_DIMS> MakeEmptyExtents_(int NumDims, coord_type_tag<int>);
  interval<double,MAX_DIMS> MakeEmptyExtents_(int NumDims, coord_type_tag<double>);

  interval<int,MAX_DIMS> UnionExtents_(const interval<int,MAX_DIMS> &Left, const
    interval<int,MAX_DIMS> &Right);
  interval<double,MAX_DIMS> UnionExtents_(const interval<double,MAX_DIMS> &Left, const
    interval<double,MAX_DIMS> &Right);

  static tuple<int> GetBinSize_(const interval<int,MAX_DIMS> &Extents, const tuple<int> &NumBins);
  static tuple<double> GetBinSize_(const interval<double,MAX_DIMS> &Extents, const tuple<int>
    &NumBins);

};

}}

#include <ovk/core/RegionHash.inl>

#endif
