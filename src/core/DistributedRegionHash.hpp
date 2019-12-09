// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISTRIBUTED_REGION_HASH_HPP_INCLUDED
#define OVK_CORE_DISTRIBUTED_REGION_HASH_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Box.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/RegionTraits.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <cstring>
#include <type_traits>

namespace ovk {
namespace core {

template <typename CoordType> class distributed_region_hash;

template <typename CoordType> class distributed_region_data {

public:

  using region_type = typename std::conditional<std::is_same<CoordType, double>::value, box,
    range>::type;

  const region_type &Region() const { return Region_; }

  int Rank() const { return Rank_; }

  template <typename T> const T &AuxData() const {
    return *reinterpret_cast<const T *>(AuxData_.Data());
  }

private:

  region_type Region_;
  int Rank_;
  array<byte> AuxData_;

  friend class distributed_region_hash<CoordType>;

};

template <typename CoordType> class distributed_region_hash_retrieved_bins {

public:

  using region_data = distributed_region_data<CoordType>;

  const region_data &RegionData(int iRegion) const { return RegionData_(iRegion); }

  array_view<const int> BinRegionIndices(int iBin) const {
    const interval<long long> &Interval = BinRegionIndicesIntervals_(iBin);
    return {BinRegionIndices_.Data(Interval.Begin()), Interval};
  }

private:

  array<region_data> RegionData_;
  map<int,interval<long long>> BinRegionIndicesIntervals_;
  array<int> BinRegionIndices_;

  friend class distributed_region_hash<CoordType>;

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
    region_type> LocalRegions);
  // is_trivially_copyable not yet supported on Intel 17
//   template <typename ArrayType, OVK_FUNCDECL_REQUIRES(std::is_trivially_copyable<array_value_type<
//     ArrayType>>::value && std::is_trivially_destructible<array_value_type<ArrayType>>::value)>
  template <typename ArrayType, OVK_FUNCDECL_REQUIRES(std::is_trivially_destructible<
    array_value_type<ArrayType>>::value)> distributed_region_hash(int NumDims, comm_view Comm, int
    NumLocalRegions, array_view<const region_type> LocalRegions, const ArrayType
    &LocalRegionAuxData, MPI_Datatype AuxDataMPIType=GetMPIDataType<array_value_type<ArrayType>>());

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

  long long AuxDataNumBytes_;
  MPI_Datatype AuxDataMPIType_;

  array<region_data> RegionData_;
  range BinRange_;
  field<int> NumRegionsPerBin_;
  field<long long> BinRegionIndicesStarts_;
  array<int> BinRegionIndices_;

  distributed_region_hash(int NumDims, comm_view Comm, int NumLocalRegions, array_view<const
    region_type> LocalRegions, const array<const byte *> &LocalRegionAuxData, long long
    AuxDataNumBytes, MPI_Datatype AuxDataMPIType);

  template <typename ArrayType> static array<const byte *> GetBytePointers_(const ArrayType &Array);

  static tuple<int> BinDecomp_(int NumDims, const region_type &GlobalExtents, int MaxBins);

  static tuple<int> GetBinSize_(const range &GlobalExtents, const tuple<int> &NumBins);
  static tuple<double> GetBinSize_(const box &GlobalExtents, const tuple<int> &NumBins);

  static tuple<int> MapToUniformCell_(int NumDims, const tuple<int> &Origin, const tuple<int>
    &CellSize, const tuple<int> &Point);
  static tuple<int> MapToUniformCell_(int NumDims, const tuple<double> &Origin, const tuple<double>
    &CellSize, const tuple<double> &Point);

};

// is_trivially_copyable not yet supported on Intel 17
// template <typename CoordType> template <typename ArrayType, OVK_FUNCDEF_REQUIRES(
//   std::is_trivially_copyable<array_value_type<ArrayType>>::value && std::is_trivially_destructible<
//   array_value_type<ArrayType>>::value)> distributed_region_hash<CoordType>::distributed_region_hash(
template <typename CoordType> template <typename ArrayType, OVK_FUNCDEF_REQUIRES(
  std::is_trivially_destructible<array_value_type<ArrayType>>::value)> distributed_region_hash<
  CoordType>::distributed_region_hash(int NumDims, comm_view Comm, int NumLocalRegions,
  array_view<const region_type> LocalRegions, const ArrayType &LocalRegionAuxData, MPI_Datatype
  AuxDataMPIType):
  distributed_region_hash(NumDims, Comm, NumLocalRegions, LocalRegions, GetBytePointers_(
    LocalRegionAuxData), sizeof(array_value_type<ArrayType>), AuxDataMPIType)
{}

template <typename CoordType> template <typename ArrayType> array<const byte *>
  distributed_region_hash<CoordType>::GetBytePointers_(const ArrayType &Array) {

  array<const byte *> BytePointers({Array.Count()});

  for (int iValue = 0; iValue < Array.Count(); ++iValue) {
    BytePointers(iValue) = reinterpret_cast<const byte *>(ArrayData(Array)+iValue);
  }

  return BytePointers;

}

extern template class distributed_region_hash<int>;
extern template class distributed_region_hash<double>;

}}

#endif
