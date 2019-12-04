// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/DistributedRegionHash.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Box.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/RegionTraits.hpp"
#include "ovk/core/ScalarOps.hpp"
#include "ovk/core/Set.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <cmath>
#include <memory>
#include <numeric>
#include <type_traits>

namespace ovk {
namespace core {

template <typename CoordType> distributed_region_hash<CoordType>::distributed_region_hash(int
  NumDims, comm_view Comm):
  NumDims_(NumDims),
  Comm_(Comm),
  GlobalExtents_(traits::MakeEmptyRegion(NumDims)),
  ProcRange_(MakeEmptyRange(NumDims)),
  ProcIndexer_(ProcRange_),
  ProcSize_(MakeUniformTuple<coord_type>(NumDims, coord_type(0), coord_type(1))),
  AuxDataNumBytes_(0),
  AuxDataMPIType_(MPI_DATATYPE_NULL),
  BinRange_(MakeEmptyRange(NumDims))
{}

template <typename CoordType> distributed_region_hash<CoordType>::distributed_region_hash(int
  NumDims, comm_view Comm, int NumLocalRegions, array_view<const region_type> LocalRegions):
  distributed_region_hash(NumDims, Comm, NumLocalRegions, LocalRegions, array<const byte *>({
    NumLocalRegions}, nullptr), 0, MPI_DATATYPE_NULL)
{}

template <typename CoordType> distributed_region_hash<CoordType>::distributed_region_hash(int
  NumDims, comm_view Comm, int NumLocalRegions, array_view<const region_type> LocalRegions, const
  array<const byte *> &LocalRegionAuxData, long long AuxDataNumBytes, MPI_Datatype AuxDataMPIType):
  NumDims_(NumDims),
  Comm_(Comm),
  AuxDataNumBytes_(AuxDataNumBytes),
  AuxDataMPIType_(AuxDataMPIType)
{

  MPI_Datatype MPICoordType = GetMPIDataType<coord_type>();

  GlobalExtents_ = traits::MakeEmptyRegion(NumDims);
  for (auto &Region : LocalRegions) {
    GlobalExtents_ = traits::UnionRegions(GlobalExtents_, Region);
  }

  MPI_Allreduce(MPI_IN_PLACE, GlobalExtents_.Begin().Data(), NumDims_, MPICoordType, MPI_MIN, Comm_);
  MPI_Allreduce(MPI_IN_PLACE, GlobalExtents_.End().Data(), NumDims_, MPICoordType, MPI_MAX, Comm_);

  tuple<int> NumProcs = BinDecomp_(NumDims_, GlobalExtents_, Comm_.Size());

  ProcRange_ = range(NumProcs);
  ProcIndexer_ = range_indexer<int>(ProcRange_);

  ProcSize_ = GetBinSize_(GlobalExtents_, NumProcs);

  array<set<int>> LocalRegionOverlappedRanks({NumLocalRegions});

  for (int iRegion = 0; iRegion < NumLocalRegions; ++iRegion) {
    tuple<coord_type> LowerCorner = traits::GetRegionLowerCorner(LocalRegions(iRegion));
    tuple<coord_type> UpperCorner = traits::GetRegionUpperCorner(LocalRegions(iRegion));
    tuple<int> ProcLocLower = ClampToRange(ProcRange_, MapToUniformCell_(NumDims_,
      GlobalExtents_.Begin(), ProcSize_, LowerCorner));
    tuple<int> ProcLocUpper = ClampToRange(ProcRange_, MapToUniformCell_(NumDims_,
      GlobalExtents_.Begin(), ProcSize_, UpperCorner));
    range OverlappedProcRange = MakeEmptyRange(NumDims_);
    OverlappedProcRange = ExtendRange(OverlappedProcRange, ProcLocLower);
    OverlappedProcRange = ExtendRange(OverlappedProcRange, ProcLocUpper);
    for (int k = OverlappedProcRange.Begin(2); k < OverlappedProcRange.End(2); ++k) {
      for (int j = OverlappedProcRange.Begin(1); j < OverlappedProcRange.End(1); ++j) {
        for (int i = OverlappedProcRange.Begin(0); i < OverlappedProcRange.End(0); ++i) {
          int Rank = ProcIndexer_.ToIndex(i,j,k);
          LocalRegionOverlappedRanks(iRegion).Insert(Rank);
        }
      }
    }
  }

  set<int> SendToRanks;

  for (auto &Ranks : LocalRegionOverlappedRanks) {
    for (int Rank : Ranks) {
      SendToRanks.Insert(Rank);
    }
  }

  array<int> RecvFromRanks = core::DynamicHandshake(Comm_, SendToRanks);

  map<int,int> NumRegionsToRank;
  map<int,int> NumRegionsFromRank;

  for (auto &Ranks : LocalRegionOverlappedRanks) {
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

  Requests.Reserve(4*(NumSends + NumRecvs));

  RegionData_.Resize({NumRecvs});

  int iNextRegion = 0;
  for (int Rank : RecvFromRanks) {
    int NumRegions = NumRegionsFromRank(Rank);
    for (int iRegionFromRank = 0; iRegionFromRank < NumRegions; ++iRegionFromRank) {
      region_data &Data = RegionData_(iNextRegion);
      MPI_Irecv(Data.Region_.Begin().Data(), MAX_DIMS, MPICoordType, Rank, 0, Comm_,
        &Requests.Append());
      MPI_Irecv(Data.Region_.End().Data(), MAX_DIMS, MPICoordType, Rank, 0, Comm_,
        &Requests.Append());
      Data.Rank_ = Rank;
      Data.AuxData_.Resize({AuxDataNumBytes_});
      if (AuxDataNumBytes_ > 0) {
        MPI_Irecv(Data.AuxData_.Data(), 1, AuxDataMPIType_, Rank, 0, Comm_, &Requests.Append());
      }
      ++iNextRegion;
    }
  }

  for (int iRegion = 0; iRegion < NumLocalRegions; ++iRegion) {
    auto &Ranks = LocalRegionOverlappedRanks(iRegion);
    for (int Rank : Ranks) {
      MPI_Isend(LocalRegions(iRegion).Begin().Data(), MAX_DIMS, MPICoordType, Rank, 0, Comm_,
        &Requests.Append());
      MPI_Isend(LocalRegions(iRegion).End().Data(), MAX_DIMS, MPICoordType, Rank, 0, Comm_,
        &Requests.Append());
      if (AuxDataNumBytes_ > 0) {
        MPI_Isend(LocalRegionAuxData(iRegion), 1, AuxDataMPIType_, Rank, 0, Comm_,
          &Requests.Append());
      }
    }
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);
  Requests.Clear();

  int ProcToBinMultiplier = 1;

  if (RegionData_.Count() > 0) {

    double AvgBinRegionLength = 0.;
    for (auto &Data : RegionData_) {
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        AvgBinRegionLength += double(Data.Region_.Size(iDim));
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

  array<map<int,range>> LocalRegionOverlappedBinRanges({NumLocalRegions});

  for (int iRegion = 0; iRegion < NumLocalRegions; ++iRegion) {
    tuple<coord_type> LowerCorner = traits::GetRegionLowerCorner(LocalRegions(iRegion));
    tuple<coord_type> UpperCorner = traits::GetRegionUpperCorner(LocalRegions(iRegion));
    tuple<int> ProcLocLower = ClampToRange(ProcRange_, MapToUniformCell_(NumDims_,
      GlobalExtents_.Begin(), ProcSize_, LowerCorner));
    tuple<int> ProcLocUpper = ClampToRange(ProcRange_, MapToUniformCell_(NumDims_,
      GlobalExtents_.Begin(), ProcSize_, UpperCorner));
    range OverlappedProcRange = MakeEmptyRange(NumDims_);
    OverlappedProcRange = ExtendRange(OverlappedProcRange, ProcLocLower);
    OverlappedProcRange = ExtendRange(OverlappedProcRange, ProcLocUpper);
    for (int k = OverlappedProcRange.Begin(2); k < OverlappedProcRange.End(2); ++k) {
      for (int j = OverlappedProcRange.Begin(1); j < OverlappedProcRange.End(1); ++j) {
        for (int i = OverlappedProcRange.Begin(0); i < OverlappedProcRange.End(0); ++i) {
          tuple<int> ProcLoc = {i,j,k};
          int Rank = ProcIndexer_.ToIndex(ProcLoc);
          range ProcBinRange = MakeEmptyRange(NumDims_);
          for (int iDim = 0; iDim < NumDims_; ++iDim) {
            ProcBinRange.Begin(iDim) = 0;
            ProcBinRange.End(iDim) = ProcToBinMultipliers_(Rank);
          }
          region_type ProcExtents = traits::MakeEmptyRegion(NumDims_);
          for (int iDim = 0; iDim < NumDims_; ++iDim) {
            ProcExtents.Begin(iDim) = GlobalExtents_.Begin(iDim) + coord_type(ProcLoc(iDim))*
              ProcSize_(iDim);
            ProcExtents.End(iDim) = Min(GlobalExtents_.Begin(iDim) + coord_type(ProcLoc(iDim)+1)*
              ProcSize_(iDim), GlobalExtents_.End(iDim));
          }
          tuple<coord_type> ProcBinSize = GetBinSize_(ProcExtents, ProcBinRange.Size());
          tuple<int> BinLocLower = ClampToRange(ProcBinRange, MapToUniformCell_(NumDims_,
            ProcExtents.Begin(), ProcBinSize, LowerCorner));
          tuple<int> BinLocUpper = ClampToRange(ProcBinRange, MapToUniformCell_(NumDims_,
            ProcExtents.Begin(), ProcBinSize, UpperCorner));
          range &OverlappedBinRange = LocalRegionOverlappedBinRanges(iRegion).Insert(Rank);
          OverlappedBinRange = MakeEmptyRange(NumDims_);
          OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocLower);
          OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocUpper);
        }
      }
    }
  }

  Requests.Reserve(2*(NumSends + NumRecvs));

  array<range> BinRegionOverlappedBinRanges;
  BinRegionOverlappedBinRanges.Resize({RegionData_.Count()});

  iNextRegion = 0;
  for (int Rank : RecvFromRanks) {
    int NumRegions = NumRegionsFromRank(Rank);
    for (int iRegionFromRank = 0; iRegionFromRank < NumRegions; ++iRegionFromRank) {
      range &OverlappedBinRange = BinRegionOverlappedBinRanges(iNextRegion);
      MPI_Irecv(OverlappedBinRange.Begin().Data(), MAX_DIMS, MPI_INT, Rank, 0, Comm_,
        &Requests.Append());
      MPI_Irecv(OverlappedBinRange.End().Data(), MAX_DIMS, MPI_INT, Rank, 0, Comm_,
        &Requests.Append());
      ++iNextRegion;
    }
  }

  for (int iRegion = 0; iRegion < NumLocalRegions; ++iRegion) {
    auto &OverlappedBinRanges = LocalRegionOverlappedBinRanges(iRegion);
    for (auto &Entry : OverlappedBinRanges) {
      int Rank = Entry.Key();
      const range &OverlappedBinRange = Entry.Value();
      MPI_Isend(OverlappedBinRange.Begin().Data(), MAX_DIMS, MPI_INT, Rank, 0, Comm_,
        &Requests.Append());
      MPI_Isend(OverlappedBinRange.End().Data(), MAX_DIMS, MPI_INT, Rank, 0, Comm_,
        &Requests.Append());
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

    NumRegionsPerBin_.Resize(BinRange_, 0);

    for (auto &OverlappedBinRange : BinRegionOverlappedBinRanges) {
      for (int k = OverlappedBinRange.Begin(2); k < OverlappedBinRange.End(2); ++k) {
        for (int j = OverlappedBinRange.Begin(1); j < OverlappedBinRange.End(1); ++j) {
          for (int i = OverlappedBinRange.Begin(0); i < OverlappedBinRange.End(0); ++i) {
            tuple<int> Bin = {i,j,k};
            ++NumRegionsPerBin_(Bin);
          }
        }
      }
    }

    BinRegionIndicesStarts_.Resize(BinRange_);

    long long TotalBinRegionIndices = 0;
    for (int k = BinRange_.Begin(2); k < BinRange_.End(2); ++k) {
      for (int j = BinRange_.Begin(1); j < BinRange_.End(1); ++j) {
        for (int i = BinRange_.Begin(0); i < BinRange_.End(0); ++i) {
          tuple<int> Bin = {i,j,k};
          BinRegionIndicesStarts_(Bin) = TotalBinRegionIndices;
          TotalBinRegionIndices += NumRegionsPerBin_(Bin);
        }
      }
    }

    BinRegionIndices_.Resize({TotalBinRegionIndices});

    // Reset for filling in indices
    NumRegionsPerBin_.Fill(0);

    for (int iRegion = 0; iRegion < RegionData_.Count(); ++iRegion) {
      const range &OverlappedBinRange = BinRegionOverlappedBinRanges(iRegion);
      for (int k = OverlappedBinRange.Begin(2); k < OverlappedBinRange.End(2); ++k) {
        for (int j = OverlappedBinRange.Begin(1); j < OverlappedBinRange.End(1); ++j) {
          for (int i = OverlappedBinRange.Begin(0); i < OverlappedBinRange.End(0); ++i) {
            tuple<int> Bin = {i,j,k};
            long long iBinRegionIndex = BinRegionIndicesStarts_(Bin) + NumRegionsPerBin_(Bin);
            BinRegionIndices_(iBinRegionIndex) = iRegion;
            ++NumRegionsPerBin_(Bin);
          }
        }
      }
    }

  }

}

template <typename CoordType> elem<int,2> distributed_region_hash<CoordType>::MapToBin(const
  tuple<coord_type> &Point) const {

  tuple<int> ProcLoc = ClampToRange(ProcRange_, MapToUniformCell_(NumDims_, GlobalExtents_.Begin(),
    ProcSize_, Point));

  int Rank = ProcIndexer_.ToIndex(ProcLoc);

  range ProcBinRange = MakeEmptyRange(NumDims_);
  for (int iDim = 0; iDim < NumDims_; ++iDim) {
    ProcBinRange.Begin(iDim) = 0;
    ProcBinRange.End(iDim) = ProcToBinMultipliers_(Rank);
  }

  range_indexer_c<int> ProcBinIndexer(ProcBinRange);

  region_type ProcExtents = traits::MakeEmptyRegion(NumDims_);
  for (int iDim = 0; iDim < NumDims_; ++iDim) {
    ProcExtents.Begin(iDim) = GlobalExtents_.Begin(iDim) + coord_type(ProcLoc(iDim))*
      ProcSize_(iDim);
    ProcExtents.End(iDim) = Min(GlobalExtents_.Begin(iDim) + coord_type(ProcLoc(iDim)+1)*
      ProcSize_(iDim), GlobalExtents_.End(iDim));
  }

  tuple<coord_type> ProcBinSize = GetBinSize_(ProcExtents, ProcBinRange.Size());

  tuple<int> BinLoc = ClampToRange(ProcBinRange, MapToUniformCell_(NumDims_, ProcExtents.Begin(),
    ProcBinSize, Point));

  int iBin = ProcBinIndexer.ToIndex(BinLoc);

  return {Rank, iBin};

}

template <typename CoordType> map<int,distributed_region_hash_retrieved_bins<CoordType>>
  distributed_region_hash<CoordType>::RetrieveBins(array_view<const elem<int,2>> BinIDs) const {

  MPI_Datatype MPICoordType = GetMPIDataType<coord_type>();

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

  Requests.Reserve(3*(RetrievingRanks.Count()+BinRanks.Count()));

  struct region_send_recv {
    array<coord_type,3> Regions;
    array<int> Ranks;
    array<byte,2> AuxData;
  };

  map<int,region_send_recv> RecvRegionData;
  RecvRegionData.Reserve(BinRanks.Count());

  for (int Rank : BinRanks) {
    region_send_recv &Recv = RecvRegionData.Insert(Rank);
    int NumRegions = NumRecvRegions(Rank);
    Recv.Regions.Resize({{2,MAX_DIMS,NumRegions}});
    Recv.Ranks.Resize({NumRegions});
    Recv.AuxData.Resize({{NumRegions,AuxDataNumBytes_}});
    MPI_Irecv(Recv.Regions.Data(), 2*MAX_DIMS*NumRegions, MPICoordType, Rank, 0, Comm_,
      &Requests.Append());
    MPI_Irecv(Recv.Ranks.Data(), NumRegions, MPI_INT, Rank, 0, Comm_, &Requests.Append());
    if (AuxDataNumBytes_ > 0) {
      MPI_Irecv(Recv.AuxData.Data(), NumRegions, AuxDataMPIType_, Rank, 0, Comm_,
        &Requests.Append());
    }
  }

  map<int,region_send_recv> SendRegionData;
  SendRegionData.Reserve(RetrievingRanks.Count());

  for (int Rank : RetrievingRanks) {
    region_send_recv &Send = SendRegionData.Insert(Rank);
    const set<int> &RegionIndices = SendingRegionIndices(Rank);
    int NumRegions = RegionIndices.Count();
    Send.Regions.Resize({{2,MAX_DIMS,NumRegions}});
    Send.Ranks.Resize({NumRegions});
    Send.AuxData.Resize({{NumRegions,AuxDataNumBytes_}});
    for (int iSendRegion = 0; iSendRegion < NumRegions; ++iSendRegion) {
      int iRegion = RegionIndices[iSendRegion];
      const region_data &Data = RegionData_(iRegion);
      Send.Ranks(iSendRegion) = Data.Rank_;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Send.Regions(0,iDim,iSendRegion) = Data.Region_.Begin(iDim);
        Send.Regions(1,iDim,iSendRegion) = Data.Region_.End(iDim);
      }
      std::memcpy(Send.AuxData.Data(iSendRegion,0), Data.AuxData_.Data(), AuxDataNumBytes_);
    }
    MPI_Isend(Send.Regions.Data(), 2*MAX_DIMS*NumRegions, MPICoordType, Rank, 0, Comm_,
      &Requests.Append());
    MPI_Isend(Send.Ranks.Data(), NumRegions, MPI_INT, Rank, 0, Comm_, &Requests.Append());
    if (AuxDataNumBytes_ > 0) {
      MPI_Isend(Send.AuxData.Data(), NumRegions, AuxDataMPIType_, Rank, 0, Comm_,
        &Requests.Append());
    }
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
    int NumRegions = Recv.Ranks.Count();
    Bins.RegionData_.Resize({NumRegions});
    for (int iRegion = 0; iRegion < NumRegions; ++iRegion) {
      region_data &Data = Bins.RegionData_(iRegion);
      Data.Region_.Begin() = {
        Recv.Regions(0,0,iRegion),
        Recv.Regions(0,1,iRegion),
        Recv.Regions(0,2,iRegion)
      };
      Data.Region_.End() = {
        Recv.Regions(1,0,iRegion),
        Recv.Regions(1,1,iRegion),
        Recv.Regions(1,2,iRegion)
      };
      Data.Rank_ = Recv.Ranks(iRegion);
      Data.AuxData_.Resize({AuxDataNumBytes_});
      std::memcpy(Data.AuxData_.Data(), Recv.AuxData.Data(iRegion,0), AuxDataNumBytes_);
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

template <typename CoordType> tuple<int> distributed_region_hash<CoordType>::BinDecomp_(int NumDims,
  const region_type &GlobalExtents, int MaxBins) {

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

    region_type GlobalExtentsReduced = traits::MakeEmptyRegion(NumDims-1);

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

template class distributed_region_hash<int>;
template class distributed_region_hash<double>;

}}
