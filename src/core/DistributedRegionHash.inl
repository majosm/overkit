// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename RegionType> distributed_region_hash<RegionType>::distributed_region_hash(int
  NumDims, comm_view Comm):
  NumDims_(NumDims),
  Comm_(Comm),
  GlobalExtents_(MakeEmptyExtents_(NumDims, coord_type_tag<coord_type>())),
  ProcRange_(MakeEmptyRange(NumDims)),
  ProcIndexer_(ProcRange_),
  ProcSize_(MakeUniformTuple<coord_type>(NumDims, coord_type(0), coord_type(1))),
  BinRange_(MakeEmptyRange(NumDims))
{}

template <typename RegionType> distributed_region_hash<RegionType>::distributed_region_hash(int
  NumDims, comm_view Comm, array_view<const region_type> LocalRegions):
  NumDims_(NumDims),
  Comm_(Comm)
{

  auto CoordMPIType = mpi_serializable_traits<coord_type>::CreateMPIType();

  GlobalExtents_ = MakeEmptyExtents_(NumDims, coord_type_tag<coord_type>());
  for (auto &Region : LocalRegions) {
    GlobalExtents_ = UnionExtents_(GlobalExtents_, region_traits::ComputeExtents(NumDims, Region));
  }

  MPI_Allreduce(MPI_IN_PLACE, GlobalExtents_.Begin().Data(), NumDims_, CoordMPIType, MPI_MIN, Comm_);
  MPI_Allreduce(MPI_IN_PLACE, GlobalExtents_.End().Data(), NumDims_, CoordMPIType, MPI_MAX, Comm_);

  tuple<int> NumProcs = BinDecomp_(NumDims_, GlobalExtents_, Comm_.Size());

  ProcRange_ = range(NumProcs);
  ProcIndexer_ = range_indexer<int>(ProcRange_);

  ProcSize_ = GetBinSize_(GlobalExtents_, NumProcs);

  array<set<int>> LocalRegionOverlappedProcs({LocalRegions.Count()});

  for (int iRegion = 0; iRegion < LocalRegions.Count(); ++iRegion) {
    LocalRegionOverlappedProcs(iRegion) = MapToBins_(ProcRange_, ProcIndexer_,
      GlobalExtents_.Begin(), ProcSize_, LocalRegions(iRegion), 
      maps_to_tag<region_traits::MapsTo()>());
  }

  set<int> SendToRanks;

  for (auto &Ranks : LocalRegionOverlappedProcs) {
    for (int Rank : Ranks) {
      SendToRanks.Insert(Rank);
    }
  }

  array<int> RecvFromRanks = core::DynamicHandshake(Comm_, SendToRanks);

  map<int,int> NumRegionsToRank;
  map<int,int> NumRegionsFromRank;

  for (auto &Ranks : LocalRegionOverlappedProcs) {
    for (int Rank : Ranks) {
      ++NumRegionsToRank.Fetch(Rank, 0);
    }
  }

  for (int Rank : RecvFromRanks) {
    NumRegionsFromRank.Insert(Rank);
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

  int NumRecvs = 0;
  for (auto &Entry : NumRegionsFromRank) {
    NumRecvs += Entry.Value();
  }

  int NumSends = 0;
  for (auto &Entry : NumRegionsToRank) {
    NumSends += Entry.Value();
  }

  Requests.Reserve(NumSends + NumRecvs);

  using mpi_region_type = typename mpi_traits::packed_type;
  auto RegionMPIType = mpi_traits::CreateMPIType();
  MPI_Type_commit(&RegionMPIType.Get());

  array<mpi_region_type> SendRegions({LocalRegions.Count()});
  array<mpi_region_type> RecvRegions({NumRecvs});

  int iRecv = 0;
  for (int Rank : RecvFromRanks) {
    int NumRegions = NumRegionsFromRank(Rank);
    for (int iRegionFromRank = 0; iRegionFromRank < NumRegions; ++iRegionFromRank) {
      mpi_region_type &RecvRegion = RecvRegions(iRecv);
      MPI_Irecv(&RecvRegion, 1, RegionMPIType, Rank, 0, Comm_, &Requests.Append());
      ++iRecv;
    }
  }

  for (int iRegion = 0; iRegion < LocalRegions.Count(); ++iRegion) {
    auto &Ranks = LocalRegionOverlappedProcs(iRegion);
    mpi_region_type &SendRegion = SendRegions(iRegion);
    SendRegion = mpi_traits::Pack(LocalRegions(iRegion));
    for (int Rank : Ranks) {
      MPI_Isend(&SendRegion, 1, RegionMPIType, Rank, 0, Comm_, &Requests.Append());
    }
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  RegionData_.Resize({NumRecvs});

  iRecv = 0;
  for (int Rank : RecvFromRanks) {
    int NumRegions = NumRegionsFromRank(Rank);
    for (int iRegionFromRank = 0; iRegionFromRank < NumRegions; ++iRegionFromRank) {
      region_data &Data = RegionData_(iRecv);
      const mpi_region_type &RecvRegion = RecvRegions(iRecv);
      Data.Region_ = mpi_traits::Unpack(RecvRegion);
      Data.Rank_ = Rank;
      ++iRecv;
    }
  }

  int ProcToBinMultiplier = 1;

  if (RegionData_.Count() > 0) {

    double AvgBinRegionLength = 0.;
    for (auto &Data : RegionData_) {
      extents_type Extents = region_traits::ComputeExtents(NumDims, Data.Region_);
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        AvgBinRegionLength += double(Extents.Size(iDim));
      }
    }
    AvgBinRegionLength /= double(NumDims_*RegionData_.Count());

    double AvgProcLength = 0.;
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      AvgProcLength += double(ProcSize_(iDim));
    }
    AvgProcLength /= double(NumDims_);

    ProcToBinMultiplier = Max(int(2.*AvgProcLength/AvgBinRegionLength),1);

  }

  ProcToBinMultipliers_.Resize({Comm_.Size()});

  MPI_Allgather(&ProcToBinMultiplier, 1, MPI_INT, ProcToBinMultipliers_.Data(), 1, MPI_INT,
    Comm_);

  array<map<int,set<int>>> LocalRegionOverlappedBins({LocalRegions.Count()});

  Requests.Reserve(NumSends + NumRecvs);

  array<map<int,int>> NumLocalRegionOverlappedBins({LocalRegions.Count()});
  array<int> NumProcRegionOverlappedBins({RegionData_.Count()});

  int iNextRegion = 0;
  for (int Rank : RecvFromRanks) {
    int NumRegions = NumRegionsFromRank(Rank);
    for (int iRegionFromRank = 0; iRegionFromRank < NumRegions; ++iRegionFromRank) {
      int &NumBins = NumProcRegionOverlappedBins(iNextRegion);
      MPI_Irecv(&NumBins, 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
      ++iNextRegion;
    }
  }

  for (int iRegion = 0; iRegion < LocalRegions.Count(); ++iRegion) {
    const set<int> &OverlappedProcs = LocalRegionOverlappedProcs(iRegion);
    auto &BinsForProc = LocalRegionOverlappedBins(iRegion);
    auto &NumBinsForProc = NumLocalRegionOverlappedBins(iRegion);
    for (int iProc : OverlappedProcs) {
      tuple<int> ProcLoc = ProcIndexer_.ToTuple(iProc);
      range ProcBinRange = MakeEmptyRange(NumDims_);
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        ProcBinRange.Begin(iDim) = 0;
        ProcBinRange.End(iDim) = ProcToBinMultipliers_(iProc);
      }
      extents_type ProcExtents = MakeEmptyExtents_(NumDims_, coord_type_tag<coord_type>());
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        ProcExtents.Begin(iDim) = GlobalExtents_.Begin(iDim) + coord_type(ProcLoc(iDim))*
          ProcSize_(iDim);
        ProcExtents.End(iDim) = Min(GlobalExtents_.Begin(iDim) + coord_type(ProcLoc(iDim)+1)*
          ProcSize_(iDim), GlobalExtents_.End(iDim));
      }
      range_indexer_c<int> ProcBinIndexer(ProcBinRange);
      tuple<coord_type> ProcBinSize = GetBinSize_(ProcExtents, ProcBinRange.Size());
      set<int> &Bins = BinsForProc.Insert(iProc);
      Bins = MapToBins_(ProcBinRange, ProcBinIndexer, ProcExtents.Begin(),
        ProcBinSize, LocalRegions(iRegion), maps_to_tag<region_traits::MapsTo()>());
      NumBinsForProc.Insert(iProc, int(Bins.Count()));
    }
    for (auto &Entry : NumBinsForProc) {
      int Rank = Entry.Key();
      const int &NumBins = Entry.Value();
      MPI_Isend(&NumBins, 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
    }
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  array<array<int>> ProcRegionOverlappedBins({RegionData_.Count()});

  iNextRegion = 0;
  for (int Rank : RecvFromRanks) {
    int NumRegions = NumRegionsFromRank(Rank);
    for (int iRegionFromRank = 0; iRegionFromRank < NumRegions; ++iRegionFromRank) {
      int NumBins = NumProcRegionOverlappedBins(iNextRegion);
      array<int> &Bins = ProcRegionOverlappedBins(iNextRegion);
      Bins.Resize({NumBins});
      MPI_Irecv(Bins.Data(), NumBins, MPI_INT, Rank, 0, Comm_, &Requests.Append());
      ++iNextRegion;
    }
  }

  for (int iRegion = 0; iRegion < LocalRegions.Count(); ++iRegion) {
    auto &BinsForProc = LocalRegionOverlappedBins(iRegion);
    for (auto &Entry : BinsForProc) {
      int Rank = Entry.Key();
      const set<int> &Bins = Entry.Value();
      MPI_Isend(Bins.Data(), Bins.Count(), MPI_INT, Rank, 0, Comm_, &Requests.Append());
    }
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  if (RegionData_.Count() > 0) {

    BinRange_ = MakeEmptyRange(NumDims_);
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      BinRange_.Begin(iDim) = 0;
      BinRange_.End(iDim) = ProcToBinMultiplier;
    }

    range_indexer_c<int> BinIndexer(BinRange_);

    NumRegionsPerBin_.Resize(BinRange_, 0);

    for (auto &Bins : ProcRegionOverlappedBins) {
      for (int iBin : Bins) {
        tuple<int> BinLoc = BinIndexer.ToTuple(iBin);
        ++NumRegionsPerBin_(BinLoc);
      }
    }

    BinRegionIndicesStarts_.Resize(BinRange_);

    long long TotalBinRegionIndices = 0;
    for (int k = BinRange_.Begin(2); k < BinRange_.End(2); ++k) {
      for (int j = BinRange_.Begin(1); j < BinRange_.End(1); ++j) {
        for (int i = BinRange_.Begin(0); i < BinRange_.End(0); ++i) {
          tuple<int> BinLoc = {i,j,k};
          BinRegionIndicesStarts_(BinLoc) = TotalBinRegionIndices;
          TotalBinRegionIndices += NumRegionsPerBin_(BinLoc);
        }
      }
    }

    BinRegionIndices_.Resize({TotalBinRegionIndices});

    // Reset for filling in indices
    NumRegionsPerBin_.Fill(0);

    for (int iRegion = 0; iRegion < RegionData_.Count(); ++iRegion) {
      auto &Bins = ProcRegionOverlappedBins(iRegion);
      for (int iBin : Bins) {
        tuple<int> BinLoc = BinIndexer.ToTuple(iBin);
        long long iBinRegionIndex = BinRegionIndicesStarts_(BinLoc) + NumRegionsPerBin_(BinLoc);
        BinRegionIndices_(iBinRegionIndex) = iRegion;
        ++NumRegionsPerBin_(BinLoc);
      }
    }

  }

}

template <typename RegionType> elem<int,2> distributed_region_hash<RegionType>::MapToBin(const
  tuple<coord_type> &Point) const {

  tuple<int> ProcLoc = ClampToRange(ProcRange_, MapToUniformGridCell(NumDims_,
    GlobalExtents_.Begin(), ProcSize_, Point));

  int Rank = ProcIndexer_.ToIndex(ProcLoc);

  range ProcBinRange = MakeEmptyRange(NumDims_);
  for (int iDim = 0; iDim < NumDims_; ++iDim) {
    ProcBinRange.Begin(iDim) = 0;
    ProcBinRange.End(iDim) = ProcToBinMultipliers_(Rank);
  }

  range_indexer_c<int> ProcBinIndexer(ProcBinRange);

  extents_type ProcExtents = MakeEmptyExtents_(NumDims_, coord_type_tag<coord_type>());
  for (int iDim = 0; iDim < NumDims_; ++iDim) {
    ProcExtents.Begin(iDim) = GlobalExtents_.Begin(iDim) + coord_type(ProcLoc(iDim))*
      ProcSize_(iDim);
    ProcExtents.End(iDim) = Min(GlobalExtents_.Begin(iDim) + coord_type(ProcLoc(iDim)+1)*
      ProcSize_(iDim), GlobalExtents_.End(iDim));
  }

  tuple<coord_type> ProcBinSize = GetBinSize_(ProcExtents, ProcBinRange.Size());

  tuple<int> BinLoc = ClampToRange(ProcBinRange, MapToUniformGridCell(NumDims_, ProcExtents.Begin(),
    ProcBinSize, Point));

  int iBin = ProcBinIndexer.ToIndex(BinLoc);

  return {Rank, iBin};

}

template <typename RegionType> map<int,distributed_region_hash_retrieved_bins<RegionType>>
  distributed_region_hash<RegionType>::RetrieveBins(array_view<const elem<int,2>> BinIDs) const {

  set<int> BinRanks;
  for (auto &BinID : BinIDs) {
    int Rank = BinID(0);
    BinRanks.Insert(Rank);
  }

  array<int> RetrievingRanks = core::DynamicHandshake(Comm_, BinRanks);

  array<MPI_Request> Requests;

  Requests.Reserve(RetrievingRanks.Count()+BinRanks.Count());

  map<int,int> NumBinsSending;
  NumBinsSending.Reserve(RetrievingRanks.Count());

  for (int Rank : RetrievingRanks) {
    int &NumBins = NumBinsSending.Insert(Rank);
    MPI_Irecv(&NumBins, 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  map<int,int> NumBinsReceiving;
  NumBinsReceiving.Reserve(BinRanks.Count());

  for (auto &BinID : BinIDs) {
    int Rank = BinID(0);
    ++NumBinsReceiving.Fetch(Rank, 0);
  }

  for (int Rank : BinRanks) {
    int &NumBins = NumBinsReceiving(Rank);
    MPI_Isend(&NumBins, 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  Requests.Reserve(RetrievingRanks.Count()+BinRanks.Count());

  map<int,array<int>> BinsSending;
  BinsSending.Reserve(RetrievingRanks.Count());

  for (int Rank : RetrievingRanks) {
    int NumBins = NumBinsSending(Rank);
    array<int> &Bins = BinsSending.Insert(Rank);
    Bins.Resize({NumBins});
    MPI_Irecv(Bins.Data(), NumBins, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  map<int,array<int>> BinsReceiving;
  BinsReceiving.Reserve(BinRanks.Count());

  for (int Rank : BinRanks) {
    int NumBins = NumBinsReceiving(Rank);
    array<int> &Bins = BinsReceiving.Insert(Rank);
    Bins.Reserve(NumBins);
  }

  for (auto &BinID : BinIDs) {
    int Rank = BinID(0);
    int iBin = BinID(1);
    array<int> &Bins = BinsReceiving(Rank);
    Bins.Append(iBin);
  }

  for (int Rank : BinRanks) {
    const array<int> &Bins = BinsReceiving(Rank);
    MPI_Isend(Bins.Data(), Bins.Count(), MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  map<int,set<int>> SendingRegionIndices;
  SendingRegionIndices.Reserve(RetrievingRanks.Count());

  for (int Rank : RetrievingRanks) {
    set<int> &RegionIndices = SendingRegionIndices.Insert(Rank);
    const array<int> &Bins = BinsSending(Rank);
    for (int iBin : Bins) {
      int NumBinRegions = NumRegionsPerBin_[iBin];
      int BinRegionIndicesStart = BinRegionIndicesStarts_[iBin];
      for (int iBinRegion = 0; iBinRegion < NumBinRegions; ++iBinRegion) {
        int iRegion = BinRegionIndices_(BinRegionIndicesStart+iBinRegion);
        RegionIndices.Insert(iRegion);
      }
    }
  }

  Requests.Reserve(RetrievingRanks.Count()+BinRanks.Count());

  map<int,int> NumRecvRegions;
  NumRecvRegions.Reserve(BinRanks.Count());

  for (int Rank : BinRanks) {
    int &NumRegions = NumRecvRegions.Insert(Rank);
    MPI_Irecv(&NumRegions, 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  map<int,int> NumSendRegions;
  NumSendRegions.Reserve(RetrievingRanks.Count());

  for (int Rank : RetrievingRanks) {
    int &NumRegions = NumSendRegions.Insert(Rank);
    NumRegions = SendingRegionIndices(Rank).Count();
    MPI_Isend(&NumRegions, 1, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  Requests.Reserve(2*(RetrievingRanks.Count()+BinRanks.Count()));

  using mpi_region_type = typename mpi_traits::packed_type;
  auto RegionMPIType = mpi_traits::CreateMPIType();
  MPI_Type_commit(&RegionMPIType.Get());

  struct region_send_recv {
    array<mpi_region_type> Regions;
    array<int> Ranks;
  };

  map<int,region_send_recv> RecvRegionData;
  RecvRegionData.Reserve(BinRanks.Count());

  for (int Rank : BinRanks) {
    region_send_recv &Recv = RecvRegionData.Insert(Rank);
    int NumRegions = NumRecvRegions(Rank);
    Recv.Regions.Resize({NumRegions});
    Recv.Ranks.Resize({NumRegions});
    MPI_Irecv(Recv.Regions.Data(), NumRegions, RegionMPIType, Rank, 0, Comm_,
      &Requests.Append());
    MPI_Irecv(Recv.Ranks.Data(), NumRegions, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  map<int,region_send_recv> SendRegionData;
  SendRegionData.Reserve(RetrievingRanks.Count());

  for (int Rank : RetrievingRanks) {
    region_send_recv &Send = SendRegionData.Insert(Rank);
    const set<int> &RegionIndices = SendingRegionIndices(Rank);
    int NumRegions = RegionIndices.Count();
    Send.Regions.Resize({NumRegions});
    Send.Ranks.Resize({NumRegions});
    for (int iSendRegion = 0; iSendRegion < NumRegions; ++iSendRegion) {
      int iRegion = RegionIndices[iSendRegion];
      const region_data &Data = RegionData_(iRegion);
      Send.Regions(iSendRegion) = mpi_traits::Pack(Data.Region_);
      Send.Ranks(iSendRegion) = Data.Rank_;
    }
    MPI_Isend(Send.Regions.Data(), NumRegions, RegionMPIType, Rank, 0, Comm_,
      &Requests.Append());
    MPI_Isend(Send.Ranks.Data(), NumRegions, MPI_INT, Rank, 0, Comm_, &Requests.Append());
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  Requests.Reserve(RetrievingRanks.Count()+BinRanks.Count());

  struct bin_index_data {
    array<int> NumRegionsPerBin;
    array<int> BinRegionIndices;
  };

  map<int,bin_index_data> RecvBinIndexData;
  RecvBinIndexData.Reserve(BinRanks.Count());

  for (int Rank : BinRanks) {
    bin_index_data &BinIndexData = RecvBinIndexData.Insert(Rank);
    int NumBins = BinsReceiving(Rank).Count();
    BinIndexData.NumRegionsPerBin.Resize({NumBins});
    MPI_Irecv(BinIndexData.NumRegionsPerBin.Data(), NumBins, MPI_INT, Rank, 0, Comm_,
      &Requests.Append());
  }

  map<int,bin_index_data> SendBinIndexData;
  SendBinIndexData.Reserve(RetrievingRanks.Count());

  for (int Rank : RetrievingRanks) {
    bin_index_data &BinIndexData = SendBinIndexData.Insert(Rank);
    const array<int> &Bins = BinsSending(Rank);
    int NumBins = Bins.Count();
    BinIndexData.NumRegionsPerBin.Resize({NumBins});
    for (int iSendBin = 0; iSendBin < NumBins; ++iSendBin) {
      int iBin = Bins[iSendBin];
      BinIndexData.NumRegionsPerBin(iSendBin) = NumRegionsPerBin_[iBin];
    }
    MPI_Isend(BinIndexData.NumRegionsPerBin.Data(), NumBins, MPI_INT, Rank, 0, Comm_,
      &Requests.Append());
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  Requests.Reserve(2*(RetrievingRanks.Count()+BinRanks.Count()));

  for (int Rank : BinRanks) {
    bin_index_data &BinIndexData = RecvBinIndexData(Rank);
    long long TotalBinRegions = 0;
    for (int NumRegions : BinIndexData.NumRegionsPerBin) {
      TotalBinRegions += (long long)(NumRegions);
    }
    BinIndexData.BinRegionIndices.Resize({TotalBinRegions});
    MPI_Irecv(BinIndexData.BinRegionIndices.Data(), TotalBinRegions, MPI_INT, Rank, 0, Comm_,
      &Requests.Append());
  }

  for (int Rank : RetrievingRanks) {
    bin_index_data &BinIndexData = SendBinIndexData(Rank);
    const set<int> &RegionIndices = SendingRegionIndices(Rank);
    const array<int> &Bins = BinsSending(Rank);
    int NumBins = Bins.Count();
    array<int> BinRegionIndexToSendRegionIndex({RegionData_.Count()}, -1);
    for (int iSendRegion = 0; iSendRegion < RegionIndices.Count(); ++iSendRegion) {
      int iRegion = RegionIndices[iSendRegion];
      BinRegionIndexToSendRegionIndex(iRegion) = iSendRegion;
    }
    long long TotalBinRegions = 0;
    for (int iSendBin = 0; iSendBin < NumBins; ++iSendBin) {
      TotalBinRegions += (long long)(BinIndexData.NumRegionsPerBin(iSendBin));
    }
    BinIndexData.BinRegionIndices.Resize({TotalBinRegions});
    long long iSendBinRegionIndex = 0;
    for (int iSendBin = 0; iSendBin < NumBins; ++iSendBin) {
      int iBin = Bins[iSendBin];
      for (int iRegion = 0; iRegion < NumRegionsPerBin_[iBin]; ++iRegion) {
        long long iBinRegionIndex = BinRegionIndicesStarts_[iBin]+iRegion;
        BinIndexData.BinRegionIndices(iSendBinRegionIndex) = BinRegionIndexToSendRegionIndex(
          BinRegionIndices_[iBinRegionIndex]);
        ++iSendBinRegionIndex;
      }
    }
    MPI_Isend(BinIndexData.BinRegionIndices.Data(), TotalBinRegions, MPI_INT, Rank, 0, Comm_,
      &Requests.Append());
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  map<int,retrieved_bins> RetrievedBins;
  RetrievedBins.Reserve(BinRanks.Count());

  for (int Rank : BinRanks) {
    retrieved_bins &Bins = RetrievedBins.Insert(Rank);
    const region_send_recv &Recv = RecvRegionData(Rank);
    int NumRegions = Recv.Regions.Count();
    Bins.RegionData_.Resize({NumRegions});
    for (int iRegion = 0; iRegion < NumRegions; ++iRegion) {
      region_data &Data = Bins.RegionData_(iRegion);
      Data.Region_ = mpi_traits::Unpack(Recv.Regions(iRegion));
      Data.Rank_ = Recv.Ranks(iRegion);
    }
    bin_index_data &BinIndexData = RecvBinIndexData(Rank);
    const array<int> &BinIndices = BinsReceiving(Rank);
    int NumBins = BinIndices.Count();
    long long TotalBinRegions = BinIndexData.BinRegionIndices.Count();
    Bins.BinRegionIndicesIntervals_.Reserve(NumBins);
    Bins.BinRegionIndices_.Resize({TotalBinRegions});
    long long iBinRegionIndex = 0;
    for (int iRecvBin = 0; iRecvBin < NumBins; ++iRecvBin) {
      int iBin = BinIndices(iRecvBin);
      int NumRegions = BinIndexData.NumRegionsPerBin(iRecvBin);
      interval<long long> &Interval = Bins.BinRegionIndicesIntervals_.Insert(iBin);
      Interval = {iBinRegionIndex,iBinRegionIndex+NumRegions};
      iBinRegionIndex += NumRegions;
    }
    Bins.BinRegionIndices_ = std::move(BinIndexData.BinRegionIndices);
  }

  return RetrievedBins;

}

template <typename RegionType> template <typename IndexerType> set<typename IndexerType::index_type>
  distributed_region_hash<RegionType>::MapToBins_(const range &BinRange, const IndexerType
  &BinIndexer, const tuple<coord_type> &LowerCorner, const tuple<coord_type> &BinSize, const
  region_type &Region, maps_to_tag<hashable_region_maps_to::SET>) const {
  return region_traits::MapToBins(NumDims_, BinRange, BinIndexer, LowerCorner, BinSize,
    Region);
}

template <typename RegionType> template <typename IndexerType> set<typename IndexerType::index_type>
  distributed_region_hash<RegionType>::MapToBins_(const range &BinRange, const IndexerType
  &BinIndexer, const tuple<coord_type> &LowerCorner, const tuple<coord_type> &BinSize, const
  region_type &Region, maps_to_tag<hashable_region_maps_to::RANGE>) const {
  using index_type = typename IndexerType::index_type;
  range Bins = region_traits::MapToBins(NumDims_, BinRange, LowerCorner, BinSize, Region);
  set<index_type> BinsSet;
  BinsSet.Reserve(Bins.Count());
  for (int k = Bins.Begin(2); k < Bins.End(2); ++k) {
    for (int j = Bins.Begin(1); j < Bins.End(1); ++j) {
      for (int i = Bins.Begin(0); i < Bins.End(0); ++i) {
        index_type iBin = BinIndexer.ToIndex(i,j,k);
        BinsSet.Insert(iBin);
      }
    }
  }
  return BinsSet;
}

template <typename RegionType> interval<int,MAX_DIMS> distributed_region_hash<RegionType>::
  MakeEmptyExtents_(int NumDims, coord_type_tag<int>) {

  return MakeEmptyRange(NumDims);

}

template <typename RegionType> interval<double,MAX_DIMS> distributed_region_hash<RegionType>::
  MakeEmptyExtents_(int NumDims, coord_type_tag<double>) {

  return MakeEmptyBox(NumDims);

}

template <typename RegionType> interval<int,MAX_DIMS> distributed_region_hash<RegionType>::
  UnionExtents_(const interval<int,MAX_DIMS> &Left, const interval<int,MAX_DIMS> &Right) {

  return UnionRanges(Left, Right);

}

template <typename RegionType> interval<double,MAX_DIMS> distributed_region_hash<RegionType>::
  UnionExtents_(const interval<double,MAX_DIMS> &Left, const interval<double,MAX_DIMS> &Right) {

  return UnionBoxes(Left, Right);

}

template <typename RegionType> tuple<int> distributed_region_hash<RegionType>::BinDecomp_(int
  NumDims, const extents_type &Extents, int MaxBins) {

  tuple<int> NumBins = {1,1,1};

  if (NumDims == 1) {

    NumBins(0) = MaxBins;

  } else {

    tuple<double> Length;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      Length(iDim) = double(Extents.Size(iDim));
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

    extents_type ExtentsReduced = MakeEmptyExtents_(NumDims-1, coord_type_tag<coord_type>());

    int iReducedDim;

    iReducedDim = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        ExtentsReduced.Begin(iReducedDim) = Extents.Begin(iDim);
        ExtentsReduced.End(iReducedDim) = Extents.End(iDim);
        ++iReducedDim;
      }
    }

    int MaxBinsReduced = MaxBins/NumBins(iMinLengthDim);

    tuple<int> NumBinsReduced = BinDecomp_(NumDims-1, ExtentsReduced, MaxBinsReduced);

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

template <typename RegionType> tuple<int> distributed_region_hash<RegionType>::GetBinSize_(const
  interval<int,MAX_DIMS> &Extents, const tuple<int> &NumBins) {

  tuple<int> BinSize;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinSize(iDim) = (Extents.Size(iDim)+NumBins(iDim)-1)/NumBins(iDim);
  }

  return BinSize;

}

template <typename RegionType> tuple<double> distributed_region_hash<RegionType>::GetBinSize_(const
  interval<double,MAX_DIMS> &Extents, const tuple<int> &NumBins) {

  tuple<double> BinSize;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinSize(iDim) = Extents.Size(iDim)/double(NumBins(iDim));
  }

  return BinSize;

}

}}
