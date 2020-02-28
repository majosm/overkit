// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename RegionType> region_hash<RegionType>::region_hash(int NumDims):
  NumDims_(NumDims),
  BinRange_(MakeEmptyRange(NumDims)),
  BinIndexer_(BinRange_),
  Extents_(MakeEmptyExtents_(NumDims, coord_type_tag<coord_type>())),
  BinSize_(MakeUniformTuple<coord_type>(NumDims, coord_type(0))),
  BinRegionIndices_({1}, 0)
{}

template <typename RegionType> region_hash<RegionType>::region_hash(int NumDims, const tuple<int>
  &NumBins, array_view<const region_type> Regions):
  NumDims_(NumDims),
  BinRange_({NumBins}),
  BinIndexer_(BinRange_)
{

  Extents_ = MakeEmptyExtents_(NumDims, coord_type_tag<coord_type>());
  for (auto &Region : Regions) {
    Extents_ = UnionExtents_(Extents_, region_traits::ComputeExtents(NumDims, Region));
  }

  BinSize_ = GetBinSize_(Extents_, NumBins);

  using mapped_bins_type = typename std::conditional<region_traits::MapsTo()
    == hashable_region_maps_to::RANGE, range, set<long long>>::type;

  array<mapped_bins_type> RegionOverlappedBins({Regions.Count()});

  for (int iRegion = 0; iRegion < Regions.Count(); ++iRegion) {
    MapToBins_(Regions(iRegion), RegionOverlappedBins(iRegion));
  }

  field<long long> NumRegionsInBin(BinRange_, 0);

  for (int iRegion = 0; iRegion < Regions.Count(); ++iRegion) {
    AccumulateBinRegionCounts_(RegionOverlappedBins(iRegion), NumRegionsInBin);
  }

  BinRegionIndicesStarts_.Resize({BinRange_.Count()+1});

  long long NumBinRegions = 0;
  for (int k = BinRange_.Begin(2); k < BinRange_.End(2); ++k) {
    for (int j = BinRange_.Begin(1); j < BinRange_.End(1); ++j) {
      for (int i = BinRange_.Begin(0); i < BinRange_.End(0); ++i) {
        tuple<int> BinLoc = {i,j,k};
        long long iBin = BinIndexer_.ToIndex(BinLoc);
        BinRegionIndicesStarts_(iBin) = NumBinRegions;
        NumBinRegions += NumRegionsInBin[iBin];
      }
    }
  }
  BinRegionIndicesStarts_((BinRange_.Count()-1)+1) = NumBinRegions;

  BinRegionIndices_.Resize({NumBinRegions});

  // Reset for filling in data
  NumRegionsInBin.Fill(0);

  for (long long iRegion = 0; iRegion < Regions.Count(); ++iRegion) {
    AddToBins_(iRegion, RegionOverlappedBins(iRegion), BinRegionIndicesStarts_, BinRegionIndices_,
      NumRegionsInBin);
  }

}

template <typename RegionType> long long region_hash<RegionType>::MapToBin(const tuple<coord_type>
  &Point) const {

  long long iBin = -1;

  if (Extents_.Contains(Point)) {
    tuple<int> BinLoc = ClampToRange(BinRange_, MapToUniformGridCell(NumDims_, Extents_.Begin(),
      BinSize_, Point));
    iBin = BinIndexer_.ToIndex(BinLoc);
  }

  return iBin;

}

template <typename RegionType> array_view<const long long> region_hash<RegionType>::RetrieveBin(long
  long iBin) const {

  long long BinStart = BinRegionIndicesStarts_(iBin);
  long long BinEnd = BinRegionIndicesStarts_(iBin+1);

  return {BinRegionIndices_.Data(BinStart), {0,BinEnd-BinStart}};

}

template <typename RegionType> void region_hash<RegionType>::MapToBins_(const region_type &Region, 
  range &Bins) const {

  Bins = region_traits::MapToBins(NumDims_, BinRange_, Extents_.Begin(), BinSize_, Region);

}

template <typename RegionType> void region_hash<RegionType>::MapToBins_(const region_type &Region, 
  set<long long> &Bins) const {

  Bins = region_traits::MapToBins(NumDims_, BinRange_, BinIndexer_, Extents_.Begin(), BinSize_,
    Region);

}

template <typename RegionType> void region_hash<RegionType>::AccumulateBinRegionCounts_(const range
  &Bins, field<long long> &NumRegionsInBin) const {

  for (int k = Bins.Begin(2); k < Bins.End(2); ++k) {
    for (int j = Bins.Begin(1); j < Bins.End(1); ++j) {
      for (int i = Bins.Begin(0); i < Bins.End(0); ++i) {
        ++NumRegionsInBin(i,j,k);
      }
    }
  }

}

template <typename RegionType> void region_hash<RegionType>::AccumulateBinRegionCounts_(const
  set<long long> &Bins, field<long long> &NumRegionsInBin) const {

  for (long long iBin : Bins) {
    ++NumRegionsInBin[iBin];
  }

}

template <typename RegionType> void region_hash<RegionType>::AddToBins_(int iRegion, const
  range &Bins, const array<long long> &BinRegionIndicesStarts, array<long long> &BinRegionIndices,
  field<long long> &NumRegionsAddedToBin) const {

  for (int k = Bins.Begin(2); k < Bins.End(2); ++k) {
    for (int j = Bins.Begin(1); j < Bins.End(1); ++j) {
      for (int i = Bins.Begin(0); i < Bins.End(0); ++i) {
        long long iBin = BinIndexer_.ToIndex(i,j,k);
        BinRegionIndices(BinRegionIndicesStarts(iBin)+NumRegionsAddedToBin[iBin]) = iRegion;
        ++NumRegionsAddedToBin[iBin];
      }
    }
  }

}

template <typename RegionType> void region_hash<RegionType>::AddToBins_(int iRegion, const
  set<long long> &Bins, const array<long long> &BinRegionIndicesStarts, array<long long>
  &BinRegionIndices, field<long long> &NumRegionsAddedToBin) const {

  for (long long iBin : Bins) {
    BinRegionIndices(BinRegionIndicesStarts(iBin)+NumRegionsAddedToBin[iBin]) = iRegion;
    ++NumRegionsAddedToBin[iBin];
  }

}

template <typename RegionType> interval<int,MAX_DIMS> region_hash<RegionType>::MakeEmptyExtents_(
  int NumDims, coord_type_tag<int>) {

  return MakeEmptyRange(NumDims);

}

template <typename RegionType> interval<double,MAX_DIMS> region_hash<RegionType>::MakeEmptyExtents_(
  int NumDims, coord_type_tag<double>) {

  return MakeEmptyBox(NumDims);

}

template <typename RegionType> interval<int,MAX_DIMS> region_hash<RegionType>::UnionExtents_(const
  interval<int,MAX_DIMS> &Left, const interval<int,MAX_DIMS> &Right) {

  return UnionRanges(Left, Right);

}

template <typename RegionType> interval<double,MAX_DIMS> region_hash<RegionType>::UnionExtents_(
  const interval<double,MAX_DIMS> &Left, const interval<double,MAX_DIMS> &Right) {

  return UnionBoxes(Left, Right);

}

template <typename RegionType> tuple<int> region_hash<RegionType>::GetBinSize_(const
  interval<int,MAX_DIMS> &Extents, const tuple<int> &NumBins) {

  tuple<int> BinSize;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinSize(iDim) = (Extents.Size(iDim)+NumBins(iDim)-1)/NumBins(iDim);
  }

  return BinSize;

}

template <typename RegionType> tuple<double> region_hash<RegionType>::GetBinSize_(const
  interval<double,MAX_DIMS> &Extents, const tuple<int> &NumBins) {

  tuple<double> BinSize;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinSize(iDim) = Extents.Size(iDim)/double(NumBins(iDim));
  }

  return BinSize;

}

}}
