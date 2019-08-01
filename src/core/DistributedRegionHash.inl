// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename CoordType> distributed_region_hash<CoordType>::distributed_region_hash(int
  NumDims, comm_view Comm):
  NumDims_(NumDims),
  Comm_(Comm),
  GlobalExtents_(traits::MakeEmptyExtents(NumDims)),
  BinSize_(MakeUniformTuple<int>(NumDims, 0, 1))
{}

template <typename CoordType> distributed_region_hash<CoordType>::distributed_region_hash(int
  NumDims, comm_view Comm, int NumLocalRegions, array_view<const extents_type> LocalRegionExtents,
  array_view<const int> LocalRegionTags):
  NumDims_(NumDims),
  Comm_(Comm)
{

  MPI_Datatype MPICoordType = GetMPIDataType<coord_type>();

  GlobalExtents_ = traits::MakeEmptyExtents(NumDims);
  for (auto &RegionExtents : LocalRegionExtents) {
    GlobalExtents_ = traits::UnionExtents(GlobalExtents_, RegionExtents);
  }

  MPI_Allreduce(MPI_IN_PLACE, GlobalExtents_.Begin().Data(), NumDims_, MPICoordType, MPI_MIN, Comm_);
  MPI_Allreduce(MPI_IN_PLACE, GlobalExtents_.End().Data(), NumDims_, MPICoordType, MPI_MAX, Comm_);

  tuple<int> NumBins = BinDecomp_(NumDims_, GlobalExtents_, Comm_.Size());

  BinRange_ = range(NumBins);
  BinIndexer_ = range_indexer<int>(BinRange_);

  BinSize_ = GetBinSize_(GlobalExtents_, NumBins);

  array<set<int>> OverlappedBinIndices({NumLocalRegions});
  set<int> AllOverlappedBinIndices;

  for (int iRegion = 0; iRegion < NumLocalRegions; ++iRegion) {
    range OverlappedBinRange = MakeEmptyRange(NumDims_);
    tuple<coord_type> LowerCorner = traits::GetExtentsLowerCorner(LocalRegionExtents(iRegion));
    tuple<coord_type> UpperCorner = traits::GetExtentsUpperCorner(LocalRegionExtents(iRegion));
    tuple<int> BinLocLower = ClampToRange(BinRange_, MapToUniformCell_(NumDims,
      GlobalExtents_.Begin(), BinSize_, LowerCorner));
    tuple<int> BinLocUpper = ClampToRange(BinRange_, MapToUniformCell_(NumDims,
      GlobalExtents_.Begin(), BinSize_, UpperCorner));
    OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocLower);
    OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocUpper);
    for (int k = OverlappedBinRange.Begin(2); k < OverlappedBinRange.End(2); ++k) {
      for (int j = OverlappedBinRange.Begin(1); j < OverlappedBinRange.End(1); ++j) {
        for (int i = OverlappedBinRange.Begin(0); i < OverlappedBinRange.End(0); ++i) {
          int BinIndex = BinIndexer_.ToIndex(i,j,k);
          OverlappedBinIndices(iRegion).Insert(BinIndex);
          AllOverlappedBinIndices.Insert(BinIndex);
        }
      }
    }
  }

  array<int> SendToRanks(AllOverlappedBinIndices);

  AllOverlappedBinIndices.Clear();

  array<int> RecvFromRanks = core::DynamicHandshake(Comm_, SendToRanks);

  map<int,int> NumRegionsToRank;
  map<int,int> NumRegionsFromRank;

  for (auto &IndexSet : OverlappedBinIndices) {
    for (int BinIndex : IndexSet) {
      ++NumRegionsToRank.Fetch(BinIndex, 0);
    }
  }

  OverlappedBinIndices.Clear();

  for (int Rank : RecvFromRanks) {
    NumRegionsFromRank.Insert(Rank, 0);
  }

  array<MPI_Request> Requests;

  Requests.Reserve(SendToRanks.Count() + RecvFromRanks.Count());

  for (int Rank : RecvFromRanks) {
    MPI_Irecv(&NumRegionsFromRank(Rank), 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  for (int Rank : SendToRanks) {
    MPI_Isend(&NumRegionsToRank(Rank), 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  int NumBinRegions = 0;
  for (auto &Entry : NumRegionsFromRank) {
    NumBinRegions += Entry.Value();
  }

  int NumSends = 0;
  for (auto &Entry : NumRegionsToRank) {
    NumSends += 3*Entry.Value();
  }

  int NumRecvs = 3*NumBinRegions;

  Requests.Reserve(NumSends + NumRecvs);

  if (NumBinRegions > 0) {

    int BinIndex = Comm_.Rank();

    tuple<int> BinLoc = BinIndexer_.ToTuple(BinIndex);

    extents_type BinExtents = traits::MakeEmptyExtents(NumDims_);
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      BinExtents.Begin(iDim) = GlobalExtents_.Begin(iDim)+BinSize_(iDim)*coord_type(BinLoc(iDim));
      BinExtents.End(iDim) = GlobalExtents_.Begin(iDim)+BinSize_(iDim)*coord_type(BinLoc(iDim)+1);
    }
    BinExtents = traits::IntersectExtents(BinExtents, GlobalExtents_);

    Bin_ = bin(BinIndex, BinExtents, NumBinRegions);

    int iRegion = 0;
    for (int Rank : RecvFromRanks) {
      for (int iRegionFromRank = 0; iRegionFromRank < NumRegionsFromRank(Rank); ++iRegionFromRank) {
        region &Region = Bin_.Region(iRegion);
        Region.Rank = Rank;
        MPI_Irecv(Region.Extents.Begin().Data(), MAX_DIMS, MPICoordType, Rank, 0, Comm_,
          &Requests.Append());
        MPI_Irecv(Region.Extents.End().Data(), MAX_DIMS, MPICoordType, Rank, 0, Comm_,
          &Requests.Append());
        MPI_Irecv(&Region.Tag, 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
        ++iRegion;
      }
    }

  }

  for (int Rank : SendToRanks) {
    for (int iRegion = 0; iRegion < NumLocalRegions; ++iRegion) {
      MPI_Isend(LocalRegionExtents(iRegion).Begin().Data(), MAX_DIMS, MPICoordType, Rank, 0, Comm_,
        &Requests.Append());
      MPI_Isend(LocalRegionExtents(iRegion).End().Data(), MAX_DIMS, MPICoordType, Rank, 0, Comm_,
        &Requests.Append());
      MPI_Isend(&LocalRegionTags(iRegion), 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
    }
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

}

template <typename CoordType> tuple<int> distributed_region_hash<CoordType>::MapPointToBin(const
  tuple<coord_type> &Point) const {

  tuple<int> BinLoc = MapToUniformCell_(NumDims_, GlobalExtents_.Begin(), BinSize_, Point);
  BinLoc = ClampToRange(BinRange_, BinLoc);

  return BinLoc;

}

template <typename CoordType> range distributed_region_hash<CoordType>::MapExtentsToBinRange(const
  extents_type &Extents) const {

  tuple<int> BinLocLower = MapPointToBin(traits::GetExtentsLowerCorner(Extents));
  tuple<int> BinLocUpper = MapPointToBin(traits::GetExtentsUpperCorner(Extents));

  range BinRange;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinRange.Begin(iDim) = BinLocLower(iDim);
    BinRange.End(iDim) = BinLocUpper(iDim)+1;
  }

  return BinRange;

}

template <typename CoordType> void distributed_region_hash<CoordType>::RetrieveBins(map<int,bin>
  &Bins) const {

  MPI_Datatype MPICoordType = GetMPIDataType<coord_type>();

  array<int> RecvFromRanks;
  for (auto &Entry : Bins) {
    int BinIndex = Entry.Key();
    bin &Bin = Entry.Value();
    if (!Bin) {
      RecvFromRanks.Append(BinIndex);
    }
  }

  array<int> SendToRanks = core::DynamicHandshake(Comm_, RecvFromRanks);

  array<MPI_Request> Requests;

  Requests.Reserve(SendToRanks.Count() + RecvFromRanks.Count());

  map<int,int> NumRegionsInBin;
  NumRegionsInBin.Reserve(RecvFromRanks.Count());

  for (int Rank : RecvFromRanks) {
    NumRegionsInBin.Insert(Rank);
    MPI_Irecv(&NumRegionsInBin(Rank), 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  int NumBinRegions = 0;
  if (Bin_) NumBinRegions = Bin_.Regions().Count();

  for (int Rank : SendToRanks) {
    MPI_Isend(&NumBinRegions, 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  for (int Rank : RecvFromRanks) {
    tuple<int> BinLoc = BinIndexer_.ToTuple(Rank);
    extents_type BinExtents = traits::MakeEmptyExtents(NumDims_);
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      BinExtents.Begin(iDim) = GlobalExtents_.Begin(iDim)+BinSize_(iDim)*coord_type(BinLoc(iDim));
      BinExtents.End(iDim) = GlobalExtents_.Begin(iDim)+BinSize_(iDim)*coord_type(BinLoc(iDim)+1);
    }
    BinExtents = traits::IntersectExtents(BinExtents, GlobalExtents_);
    Bins(Rank) = bin(Rank, BinExtents, NumRegionsInBin(Rank));
  }

  struct bin_region_data {
    int NumRegions = 0;
    array<int> Ranks;
    array<coord_type,3> ExtentsValues;
    array<int> Tags;
    bin_region_data() = default;
    explicit bin_region_data(int NumRegions_):
      NumRegions(NumRegions_),
      Ranks({NumRegions}),
      ExtentsValues({{NumRegions,2,MAX_DIMS}}),
      Tags({NumRegions})
    {}
  };

  bin_region_data BinRegionData;
  if (Bin_) {
    BinRegionData = bin_region_data(Bin_.Regions().Count());
    for (int iRegion = 0; iRegion < Bin_.Regions().Count(); ++iRegion) {
      const region &Region = Bin_.Region(iRegion);
      BinRegionData.Ranks(iRegion) = Region.Rank;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        BinRegionData.ExtentsValues(iRegion,0,iDim) = Region.Extents.Begin(iDim);
        BinRegionData.ExtentsValues(iRegion,1,iDim) = Region.Extents.End(iDim);
      }
      BinRegionData.Tags(iRegion) = Region.Tag;
    }
  }

  map<int,bin_region_data> RetrievedBinRegionData;
  for (int Rank : RecvFromRanks) {
    RetrievedBinRegionData.Insert(Rank, bin_region_data(NumRegionsInBin(Rank)));
  }

  Requests.Reserve(3*SendToRanks.Count() + 3*RecvFromRanks.Count());

  for (int Rank : RecvFromRanks) {
    bin_region_data &RetrievedData = RetrievedBinRegionData(Rank);
    MPI_Irecv(RetrievedData.Ranks.Data(), RetrievedData.NumRegions, MPI_INT, Rank, 0, Comm_,
      &Requests.Append());
    MPI_Irecv(RetrievedData.ExtentsValues.Data(), 2*MAX_DIMS*RetrievedData.NumRegions, MPICoordType,
      Rank, 0, Comm_, &Requests.Append());
    MPI_Irecv(RetrievedData.Tags.Data(), RetrievedData.NumRegions, MPI_INT, Rank, 0, Comm_,
      &Requests.Append());
  }

  for (int Rank : SendToRanks) {
    MPI_Isend(BinRegionData.Ranks.Data(), BinRegionData.NumRegions, MPI_INT, Rank, 0, Comm_,
      &Requests.Append());
    MPI_Isend(BinRegionData.ExtentsValues.Data(), 2*MAX_DIMS*BinRegionData.NumRegions, MPICoordType,
      Rank, 0, Comm_, &Requests.Append());
    MPI_Isend(BinRegionData.Tags.Data(), BinRegionData.NumRegions, MPI_INT, Rank, 0, Comm_,
      &Requests.Append());
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  for (int Rank : RecvFromRanks) {
    bin &Bin = Bins(Rank);
    const bin_region_data &RetrievedData = RetrievedBinRegionData(Rank);
    for (int iRegion = 0; iRegion < Bin.Regions().Count(); ++iRegion) {
      region &Region = Bin.Region(iRegion);
      Region.Rank = RetrievedData.Ranks(iRegion);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Region.Extents.Begin(iDim) = RetrievedData.ExtentsValues(iRegion,0,iDim);
        Region.Extents.End(iDim) = RetrievedData.ExtentsValues(iRegion,1,iDim);
      }
      Region.Tag = RetrievedData.Tags(iRegion);
    }
  }

}

template <typename CoordType> tuple<int> distributed_region_hash<CoordType>::BinDecomp_(int NumDims,
  const extents_type &GlobalExtents, int MaxBins) {

  tuple<int> NumBins = {1,1,1};

  if (NumDims == 1) {

    NumBins(0) = MaxBins;

  } else {

    tuple<double> Length;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      Length(iDim) = double(GlobalExtents.Size(iDim));
    }

    double Volume = 1.;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      Volume *= Length(iDim);
    }

    int iMinLengthDim = 0;
    double MinLength = std::numeric_limits<double>::max();
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (Length(iDim) <= MinLength) {
        iMinLengthDim = iDim;
        MinLength = Length(iDim);
      }
    }

    double Base = std::pow(double(MaxBins)/Volume, 1./double(NumDims));

    NumBins(iMinLengthDim) = Max(int(Length(iMinLengthDim)*Base),1);

    extents_type GlobalExtentsReduced = traits::MakeEmptyExtents(NumDims-1);

    int iReducedDim;

    iReducedDim = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        GlobalExtentsReduced.Begin(iReducedDim) = GlobalExtents.Begin(iDim);
        GlobalExtentsReduced.End(iReducedDim) = GlobalExtents.End(iDim);
        ++iReducedDim;
      }
    }

    int MaxBinsReduced = MaxBins/NumBins(iMinLengthDim);

    tuple<int> NumBinsReduced = BinDecomp_(NumDims-1, GlobalExtentsReduced, MaxBinsReduced);

    iReducedDim = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        NumBins(iDim) = NumBinsReduced(iReducedDim);
        ++iReducedDim;
      }
    }

  }

  return NumBins;

}

template <typename CoordType> tuple<int> distributed_region_hash<CoordType>::GetBinSize_(const range
  &GlobalExtents, const tuple<int> &NumBins) {

  tuple<int> BinSize;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinSize(iDim) = (GlobalExtents.Size(iDim)+NumBins(iDim)-1)/NumBins(iDim);
  }

  return BinSize;

}

template <typename CoordType> tuple<double> distributed_region_hash<CoordType>::GetBinSize_(const
  box &GlobalExtents, const tuple<int> &NumBins) {

  tuple<double> BinSize;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinSize(iDim) = GlobalExtents.Size(iDim)/double(NumBins(iDim));
  }

  return BinSize;

}

template <typename CoordType> tuple<int> distributed_region_hash<CoordType>::MapToUniformCell_(int
  NumDims, const tuple<int> &Origin, const tuple<int> &CellSize, const tuple<int> &Point) {

  tuple<int> Cell = MakeUniformTuple<int>(NumDims, 0);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    int Offset = Point(iDim) - Origin(iDim);
    // Division rounding down
    Cell(iDim) = Offset/CellSize(iDim) - (Offset % CellSize(iDim) < 0);
  }

  return Cell;

}

template <typename CoordType> tuple<int> distributed_region_hash<CoordType>::MapToUniformCell_(int
  NumDims, const tuple<double> &Origin, const tuple<double> &CellSize, const tuple<double> &Point) {

  tuple<int> Cell = MakeUniformTuple<int>(NumDims, 0);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    double Offset = Point(iDim) - Origin(iDim);
    Cell(iDim) = int(std::floor(Offset/CellSize(iDim)));
  }

  return Cell;

}

}}
