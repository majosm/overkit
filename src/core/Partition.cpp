// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Partition.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Halo.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/ScalarOps.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <map>
#include <set>
#include <string>
#include <utility>

namespace ovk {
namespace core {

namespace {

array<range> CreateSubregions(int NumDims, const range &Range, int NumSubregions);

}

partition::partition(std::shared_ptr<context> Context, const cart &Cart, comm_view Comm, const range
  &LocalRange, int ExtendAmount, int NumSubregions, const array<int> &NeighborRanks):
  Context_(std::move(Context)),
  Cart_(Cart),
  Comm_(Comm),
  LocalRange_(LocalRange),
  ExtendedRange_(ExtendLocalRange(Cart, LocalRange, ExtendAmount)),
  LocalSubregions_(CreateSubregions(Cart.Dimension(), LocalRange, NumSubregions)),
  ExtendedSubregions_(CreateSubregions(Cart.Dimension(), ExtendedRange_, NumSubregions)),
  Neighbors_(RetrievePartitionInfo(Comm_, NeighborRanks, LocalRange, ExtendedRange_)),
  Halo_(Context_, Cart, Comm_, LocalRange, ExtendedRange_, Neighbors_)
{}

range ExtendLocalRange(const cart &Cart, const range &LocalRange, int ExtendAmount) {

  const range &GlobalRange = Cart.Range();
  const tuple<bool> &Periodic = Cart.Periodic();

  range ExtendedRange = LocalRange;

  for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
    if (LocalRange.Begin(iDim) != GlobalRange.Begin(iDim) || Periodic[iDim]) {
      ExtendedRange.Begin(iDim) -= ExtendAmount;
    }
    if (LocalRange.End(iDim) != GlobalRange.End(iDim) || Periodic[iDim]) {
      ExtendedRange.End(iDim) += ExtendAmount;
    }
  }

  return ExtendedRange;

}

array<int> DetectNeighbors(const cart &Cart, comm_view Comm, const range &LocalRange, const
  partition_hash &Hash) {

  range ExtendedRange = ExtendLocalRange(Cart, LocalRange, 1);

  long long NumExtendedPoints = ExtendedRange.Count() - LocalRange.Count();

  array<int,2> ExtendedPoints({{MAX_DIMS,NumExtendedPoints}});

  long long iNextPoint = 0;
  for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
    for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
      for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        if (!LocalRange.Contains(Point)) {
          Point = Cart.PeriodicAdjust(Point);
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            ExtendedPoints(iDim,iNextPoint) = Point[iDim];
          }
          ++iNextPoint;
        }
      }
    }
  }

  array<int> ExtendedPointBinIndices({NumExtendedPoints});

  Hash.MapToBins(ExtendedPoints, ExtendedPointBinIndices);

  std::set<int> UniqueBinIndices;
  for (long long iPoint = 0; iPoint < NumExtendedPoints; ++iPoint) {
    UniqueBinIndices.insert(ExtendedPointBinIndices(iPoint));
  }

  std::map<int, core::partition_hash::bin> Bins;
  for (int iBin : UniqueBinIndices) {
    Bins.emplace(iBin, core::partition_hash::bin());
  }

  Hash.RetrieveBins(Bins);

  array<int> ExtendedPointRanks({NumExtendedPoints});

  Hash.FindPartitions(Bins, ExtendedPoints, ExtendedPointBinIndices, ExtendedPointRanks);

  Bins.clear();

  std::set<int> UniqueExtendedRanks;

  for (long long iPoint = 0; iPoint < NumExtendedPoints; ++iPoint) {
    UniqueExtendedRanks.insert(ExtendedPointRanks(iPoint));
  }
  UniqueExtendedRanks.erase(Comm.Rank());

  return {{int(UniqueExtendedRanks.size())}, UniqueExtendedRanks.begin()};

}

array<partition_info> RetrievePartitionInfo(comm_view Comm, const array<int> &Ranks, const range
  &LocalRange, const range &ExtendedRange) {

  const array<int> &RecvFromRanks = Ranks;
  array<int> SendToRanks = core::DynamicHandshake(Comm, RecvFromRanks);

  int NumSends = SendToRanks.Count();
  int NumRecvs = RecvFromRanks.Count();

  array<int,4> RetrievedRangeValues({{NumRecvs,2,2,MAX_DIMS}});

  array<MPI_Request> Requests;
  Requests.Reserve(NumSends+NumRecvs);

  for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
    int Rank = RecvFromRanks(iRecv);
    MPI_Request &Request = Requests.Append();
    MPI_Irecv(RetrievedRangeValues.Data(iRecv,0,0,0), 4*MAX_DIMS, MPI_INT, Rank, 0, Comm, &Request);
  }

  array<int,3> SelfRangeValues({{2,2,MAX_DIMS}});
//   static_array<int,4*MAX_DIMS,3> LocalRangeValues({{2,2,MAX_DIMS}});
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    SelfRangeValues(0,0,iDim) = LocalRange.Begin(iDim);
    SelfRangeValues(0,1,iDim) = LocalRange.End(iDim);
    SelfRangeValues(1,0,iDim) = ExtendedRange.Begin(iDim);
    SelfRangeValues(1,1,iDim) = ExtendedRange.End(iDim);
  }

  for (int iSend = 0; iSend < NumSends; ++iSend) {
    int Rank = SendToRanks(iSend);
    MPI_Request &Request = Requests.Append();
    MPI_Isend(SelfRangeValues.Data(), 4*MAX_DIMS, MPI_INT, Rank, 0, Comm, &Request);
  }

  array<partition_info> RetrievedPartitionInfo({NumRecvs});

  while (true) {
    int iRequest;
    MPI_Waitany(Requests.Count(), Requests.Data(), &iRequest, MPI_STATUSES_IGNORE);
    if (iRequest == MPI_UNDEFINED) break;
    if (iRequest < NumRecvs) {
      int iRecv = iRequest;
      partition_info &Info = RetrievedPartitionInfo(iRecv);
      Info.Rank = RecvFromRanks(iRecv);
      tuple<int> LocalBegin(RetrievedRangeValues.Data(iRecv,0,0,0));
      tuple<int> LocalEnd(RetrievedRangeValues.Data(iRecv,0,1,0));
      tuple<int> ExtendedBegin(RetrievedRangeValues.Data(iRecv,1,0,0));
      tuple<int> ExtendedEnd(RetrievedRangeValues.Data(iRecv,1,1,0));
      Info.LocalRange = {LocalBegin, LocalEnd};
      Info.ExtendedRange = {ExtendedBegin, ExtendedEnd};
    }
  }

  return RetrievedPartitionInfo;

}

namespace {

array<range> CreateSubregions(int NumDims, const range &Range, int NumSubregions) {

  int iMaxStrideDim = NumDims-1;

  int Size = Range.Size(iMaxStrideDim);
  int SubregionSize = Size/NumSubregions;
  int Remainder = Size - NumSubregions*SubregionSize;

  array<range> Subregions({NumSubregions});
  for (int iSubregion = 0; iSubregion < NumSubregions; ++iSubregion) {
    Subregions(iSubregion) = Range;
    Subregions(iSubregion).Begin(iMaxStrideDim) = Range.Begin(iMaxStrideDim) +
      SubregionSize*iSubregion + Min(Remainder, iSubregion);
    Subregions(iSubregion).End(iMaxStrideDim) = Range.Begin(iMaxStrideDim) +
      SubregionSize*(iSubregion+1) + Min(Remainder, iSubregion+1);
  }

  return Subregions;

}

}

}}
