// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
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

  array<set<int>> RegionOverlappedBins({Regions.Count()});

  for (int iRegion = 0; iRegion < Regions.Count(); ++iRegion) {
    elem_set<int,MAX_DIMS> BinLocs = region_traits::MapToBins(NumDims_, BinRange_, Extents_.Begin(),
      BinSize_, Regions(iRegion));
    set<int> &Bins = RegionOverlappedBins(iRegion);
    Bins.Reserve(BinLocs.Count());
    for (auto &BinLoc : BinLocs) {
      Bins.Insert(BinIndexer_.ToIndex(BinLoc));
    }
  }

  field<long long> NumRegionsInBin(BinRange_, 0);

  for (int iRegion = 0; iRegion < Regions.Count(); ++iRegion) {
    const set<int> &Bins = RegionOverlappedBins(iRegion);
    for (int iBin : Bins) {
      tuple<int> BinLoc = BinIndexer_.ToTuple(iBin);
      ++NumRegionsInBin(BinLoc);
    }
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
    const set<int> &Bins = RegionOverlappedBins(iRegion);
    for (int iBin : Bins) {
      BinRegionIndices_(BinRegionIndicesStarts_(iBin)+NumRegionsInBin[iBin]) = iRegion;
      ++NumRegionsInBin[iBin];
    }
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
