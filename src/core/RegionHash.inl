// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename CoordType> region_hash<CoordType>::region_hash(int NumDims):
  NumDims_(NumDims),
  BinRange_(MakeEmptyRange(NumDims)),
  BinIndexer_(BinRange_),
  Extents_(traits::MakeEmtpyRegion(NumDims)),
  BinSize_(MakeUniformTuple<coord_type>(NumDims, coord_type(0))),
  BinRegionIndices_({1}, 0)
{}

template <typename CoordType> region_hash<CoordType>::region_hash(int NumDims, const tuple<int>
  &NumBins, array_view<const region_type> Regions):
  NumDims_(NumDims),
  BinRange_({NumBins}),
  BinIndexer_(BinRange_)
{

  Extents_ = traits::MakeEmptyRegion(NumDims);
  for (auto &Region : Regions) {
    Extents_ = traits::UnionRegions(Extents_, Region);
  }

  BinSize_ = GetBinSize_(Extents_, NumBins);

  field<long long> NumRegionsInBin(BinRange_, 0);

  for (auto &Region : Regions) {
    tuple<coord_type> LowerCorner = traits::GetRegionLowerCorner(Region);
    tuple<coord_type> UpperCorner = traits::GetRegionUpperCorner(Region);
    tuple<int> BinLocLower = ClampToRange(BinRange_, MapToUniformCell_(NumDims_,
      Extents_.Begin(), BinSize_, LowerCorner));
    tuple<int> BinLocUpper = ClampToRange(BinRange_, MapToUniformCell_(NumDims_,
      Extents_.Begin(), BinSize_, UpperCorner));
    range OverlappedBinRange = MakeEmptyRange(NumDims_);
    OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocLower);
    OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocUpper);
    for (int k = OverlappedBinRange.Begin(2); k < OverlappedBinRange.End(2); ++k) {
      for (int j = OverlappedBinRange.Begin(1); j < OverlappedBinRange.End(1); ++j) {
        for (int i = OverlappedBinRange.Begin(0); i < OverlappedBinRange.End(0); ++i) {
          tuple<int> BinLoc = {i,j,k};
          ++NumRegionsInBin(BinLoc);
        }
      }
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
    const region_type &Region = Regions(iRegion);
    tuple<coord_type> LowerCorner = traits::GetRegionLowerCorner(Region);
    tuple<coord_type> UpperCorner = traits::GetRegionUpperCorner(Region);
    tuple<int> BinLocLower = ClampToRange(BinRange_, MapToUniformCell_(NumDims_,
      Extents_.Begin(), BinSize_, LowerCorner));
    tuple<int> BinLocUpper = ClampToRange(BinRange_, MapToUniformCell_(NumDims_,
      Extents_.Begin(), BinSize_, UpperCorner));
    range OverlappedBinRange = MakeEmptyRange(NumDims_);
    OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocLower);
    OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocUpper);
    for (int k = OverlappedBinRange.Begin(2); k < OverlappedBinRange.End(2); ++k) {
      for (int j = OverlappedBinRange.Begin(1); j < OverlappedBinRange.End(1); ++j) {
        for (int i = OverlappedBinRange.Begin(0); i < OverlappedBinRange.End(0); ++i) {
          tuple<int> BinLoc = {i,j,k};
          long long iBin = BinIndexer_.ToIndex(BinLoc);
          BinRegionIndices_(BinRegionIndicesStarts_(iBin)+NumRegionsInBin[iBin]) = iRegion;
          ++NumRegionsInBin[iBin];
        }
      }
    }
  }

}

template <typename CoordType> long long region_hash<CoordType>::MapToBin(const tuple<coord_type>
  &Point) const {

  long long iBin = -1;

  if (Extents_.Contains(Point)) {
    tuple<int> BinLoc = ClampToRange(BinRange_, MapToUniformCell_(NumDims_, Extents_.Begin(),
      BinSize_, Point));
    iBin = BinIndexer_.ToIndex(BinLoc);
  }

  return iBin;

}

template <typename CoordType> array_view<const long long> region_hash<CoordType>::RetrieveBin(long
  long iBin) const {

  long long BinStart = BinRegionIndicesStarts_(iBin);
  long long BinEnd = BinRegionIndicesStarts_(iBin+1);

  return {BinRegionIndices_.Data(BinStart), {0,BinEnd-BinStart}};

}

template <typename CoordType> tuple<int> region_hash<CoordType>::GetBinSize_(const range &Extents,
  const tuple<int> &NumBins) {

  tuple<int> BinSize;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinSize(iDim) = (Extents.Size(iDim)+NumBins(iDim)-1)/NumBins(iDim);
  }

  return BinSize;

}

template <typename CoordType> tuple<double> region_hash<CoordType>::GetBinSize_(const box &Extents,
  const tuple<int> &NumBins) {

  tuple<double> BinSize;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinSize(iDim) = Extents.Size(iDim)/double(NumBins(iDim));
  }

  return BinSize;

}

template <typename CoordType> tuple<int> region_hash<CoordType>::MapToUniformCell_(int NumDims,
  const tuple<int> &Origin, const tuple<int> &CellSize, const tuple<int> &Point) {

  tuple<int> Cell = MakeUniformTuple<int>(NumDims, 0);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    int Offset = Point(iDim) - Origin(iDim);
    // Division rounding down
    Cell(iDim) = Offset/CellSize(iDim) - (Offset % CellSize(iDim) < 0);
  }

  return Cell;

}

template <typename CoordType> tuple<int> region_hash<CoordType>::MapToUniformCell_(int NumDims,
  const tuple<double> &Origin, const tuple<double> &CellSize, const tuple<double> &Point) {

  tuple<int> Cell = MakeUniformTuple<int>(NumDims, 0);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    double Offset = Point(iDim) - Origin(iDim);
    Cell(iDim) = int(std::floor(Offset/CellSize(iDim)));
  }

  return Cell;

}

}}
