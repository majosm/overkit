// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Partition.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/DistributedRegionHash.hpp"
#include "ovk/core/Field.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/ScalarOps.hpp"
#include "ovk/core/ScopeGuard.hpp"
#include "ovk/core/TextProcessing.hpp"
#include "ovk/core/TypeTraits.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

namespace {

array<range> CreateSubregions(int NumDims, const range &Range, int NumSubregions);

}

partition::partition(std::shared_ptr<context> Context, const cart &Cart, comm_view Comm, const range
  &LocalRange, const range &ExtendedRange, int NumSubregions, array_view<const int> NeighborRanks):
  Context_(std::move(Context)),
  Cart_(Cart),
  Comm_(Comm),
  LocalRange_(LocalRange),
  ExtendedRange_(ExtendedRange),
  NumSubregions_(NumSubregions),
  LocalSubregions_(CreateSubregions(Cart.Dimension(), LocalRange, NumSubregions)),
  ExtendedSubregions_(CreateSubregions(Cart.Dimension(), ExtendedRange_, NumSubregions)),
  Neighbors_(partition_internal::RetrieveDecompInfo(Comm_, NeighborRanks, LocalRange,
    ExtendedRange_)),
  Halo_(Context_, Cart, Comm_, LocalRange, ExtendedRange_, Neighbors_)
{}

namespace core {

partition_hash CreatePartitionHash(int NumDims, comm_view Comm, const range &LocalRange) {

  return {NumDims, Comm, 1, array<range>({1}, {LocalRange}), array<int>({1}, {1})};

}

range ExtendLocalRange(const cart &Cart, const range &LocalRange, int ExtendAmount) {

  range ExtendedRange = LocalRange;

  for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
    if (LocalRange.Begin(iDim) != Cart.Range().Begin(iDim) || Cart.Periodic(iDim)) {
      ExtendedRange.Begin(iDim) -= ExtendAmount;
    }
    if (LocalRange.End(iDim) != Cart.Range().End(iDim) || Cart.Periodic(iDim)) {
      ExtendedRange.End(iDim) += ExtendAmount;
    }
  }

  if (Cart.PeriodicStorage() == periodic_storage::UNIQUE) {
    for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
      if (LocalRange.End(iDim) == Cart.Range().End(iDim) && Cart.Periodic(iDim)) {
        ExtendedRange.End(iDim) += 1;
      }
    }
  }

  return ExtendedRange;

}

cart CartPointToCell(const cart &Cart) {

  cart CellCart = Cart;

  for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
    CellCart.Range().Begin(iDim) = Cart.Range().Begin(iDim);
    if (Cart.Periodic(iDim) && Cart.PeriodicStorage() == periodic_storage::UNIQUE) {
      CellCart.Range().End(iDim) = Cart.Range().End(iDim);
    } else {
      CellCart.Range().End(iDim) = Cart.Range().End(iDim)-1;
    }
  }

  CellCart.PeriodicStorage() = periodic_storage::UNIQUE;

  return CellCart;

}

range LocalRangePointToCell(const cart &Cart, const range &LocalRange) {

  range CellLocalRange = MakeEmptyRange(Cart.Dimension());

  for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
    CellLocalRange.Begin(iDim) = LocalRange.Begin(iDim);
    if (Cart.Periodic(iDim) && Cart.PeriodicStorage() == periodic_storage::UNIQUE &&
      LocalRange.End(iDim) == Cart.Range().End(iDim)) {
      CellLocalRange.End(iDim) = LocalRange.End(iDim);
    } else {
      CellLocalRange.End(iDim) = LocalRange.End(iDim)-1;
    }
  }

  return CellLocalRange;

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
            ExtendedPoints(iDim,iNextPoint) = Point(iDim);
          }
          ++iNextPoint;
        }
      }
    }
  }

  array<int> ExtendedPointBinIndices({NumExtendedPoints});
  set<int> UniqueBinIndices;

  for (long long iPoint = 0; iPoint < NumExtendedPoints; ++iPoint) {
    tuple<int> Point = {
      ExtendedPoints(0,iPoint),
      ExtendedPoints(1,iPoint),
      ExtendedPoints(2,iPoint)
    };
    tuple<int> BinLoc = Hash.MapPointToBin(Point);
    ExtendedPointBinIndices(iPoint) = Hash.BinIndexer().ToIndex(BinLoc);
    UniqueBinIndices.Insert(ExtendedPointBinIndices(iPoint));
  }

  map<int,partition_hash_bin> Bins;
  for (int iBin : UniqueBinIndices) {
    Bins.Insert(iBin);
  }

  Hash.RetrieveBins(Bins);

  set<int> UniqueExtendedRanks;

  for (long long iPoint = 0; iPoint < NumExtendedPoints; ++iPoint) {
    tuple<int> Point = {
      ExtendedPoints(0,iPoint),
      ExtendedPoints(1,iPoint),
      ExtendedPoints(2,iPoint)
    };
    const partition_hash_bin &Bin = Bins(ExtendedPointBinIndices(iPoint));
    for (int iRegion = 0; iRegion < Bin.Regions().Count(); ++iRegion) {
      const partition_hash_region_data &Region = Bin.Region(iRegion);
      if (Region.Extents.Contains(Point)) {
        UniqueExtendedRanks.Insert(Region.Rank);
        break;
      }
    }
  }
  UniqueExtendedRanks.Erase(Comm.Rank());

  Bins.Clear();

  return {UniqueExtendedRanks};

}

}

namespace partition_internal {

map<int,decomp_info> RetrieveDecompInfo(comm_view Comm, array_view<const int> Ranks, const range
  &LocalRange, const range &ExtendedRange) {

  array_view<const int> &RecvFromRanks = Ranks;
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

  map<int,decomp_info> RetrievedDecompInfo;

  while (true) {
    int iRequest;
    MPI_Waitany(Requests.Count(), Requests.Data(), &iRequest, MPI_STATUSES_IGNORE);
    if (iRequest == MPI_UNDEFINED) break;
    if (iRequest < NumRecvs) {
      int iRecv = iRequest;
      decomp_info &DecompInfo = RetrievedDecompInfo.Fetch(RecvFromRanks(iRecv));
      tuple<int> LocalBegin(RetrievedRangeValues.Data(iRecv,0,0,0));
      tuple<int> LocalEnd(RetrievedRangeValues.Data(iRecv,0,1,0));
      tuple<int> ExtendedBegin(RetrievedRangeValues.Data(iRecv,1,0,0));
      tuple<int> ExtendedEnd(RetrievedRangeValues.Data(iRecv,1,1,0));
      DecompInfo.LocalRange = {LocalBegin, LocalEnd};
      DecompInfo.ExtendedRange = {ExtendedBegin, ExtendedEnd};
    }
  }

  return RetrievedDecompInfo;

}

namespace halo_internal {

halo_map::halo_map(const cart &Cart, const range &LocalRange, const range &ExtendedRange,
  const map<int,decomp_info> &Neighbors) {

  const range &GlobalRange = Cart.Range();
  int NumNeighbors = Neighbors.Count();

  field_indexer ExtendedIndexer(ExtendedRange);

  NeighborRanks_.Resize({NumNeighbors});
  NeighborSendIndices_.Resize({NumNeighbors});
  NeighborRecvIndices_.Resize({NumNeighbors});

  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    NeighborRanks_(iNeighbor) = Neighbors[iNeighbor].Key();
  }

  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    const range &NeighborLocalRange = Neighbors[iNeighbor].Value().LocalRange;
    const range &NeighborExtendedRange = Neighbors[iNeighbor].Value().ExtendedRange;
    array<long long> &SendIndices = NeighborSendIndices_(iNeighbor);
    if (Cart.Range().Includes(NeighborExtendedRange)) {
      range SendRange = IntersectRanges(NeighborExtendedRange, LocalRange);
      SendIndices.Reserve(SendRange.Count());
      for (int k = SendRange.Begin(2); k < SendRange.End(2); ++k) {
        for (int j = SendRange.Begin(1); j < SendRange.End(1); ++j) {
          for (int i = SendRange.Begin(0); i < SendRange.End(0); ++i) {
            long long iPoint = ExtendedIndexer.ToIndex(i,j,k);
            SendIndices.Append(iPoint);
          }
        }
      }
    } else {
      field<bool> SendMask(NeighborExtendedRange, false);
      long long NumSendPoints = 0;
      for (int k = NeighborExtendedRange.Begin(2); k < NeighborExtendedRange.End(2); ++k) {
        for (int j = NeighborExtendedRange.Begin(1); j < NeighborExtendedRange.End(1); ++j) {
          for (int i = NeighborExtendedRange.Begin(0); i < NeighborExtendedRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (!NeighborLocalRange.Contains(Point)) {
              if (Cart.Range().Contains(Point)) {
                SendMask(Point) = LocalRange.Contains(Point);
              } else {
                tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
                SendMask(Point) = LocalRange.Contains(AdjustedPoint);
              }
              if (SendMask(Point)) {
                ++NumSendPoints;
              }
            }
          }
        }
      }
      SendIndices.Reserve(NumSendPoints);
      for (int k = NeighborExtendedRange.Begin(2); k < NeighborExtendedRange.End(2); ++k) {
        for (int j = NeighborExtendedRange.Begin(1); j < NeighborExtendedRange.End(1); ++j) {
          for (int i = NeighborExtendedRange.Begin(0); i < NeighborExtendedRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (SendMask(Point)) {
              long long iPoint;
              if (Cart.Range().Contains(Point)) {
                iPoint = ExtendedIndexer.ToIndex(Point);
              } else {
                tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
                iPoint = ExtendedIndexer.ToIndex(AdjustedPoint);
              }
              SendIndices.Append(iPoint);
            }
          }
        }
      }
    }
  }

  if (GlobalRange.Includes(ExtendedRange)) {
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      const range &NeighborLocalRange = Neighbors[iNeighbor].Value().LocalRange;
      array<long long> &RecvIndices = NeighborRecvIndices_(iNeighbor);
      range RecvRange = IntersectRanges(ExtendedRange, NeighborLocalRange);
      RecvIndices.Reserve(RecvRange.Count());
      for (int k = RecvRange.Begin(2); k < RecvRange.End(2); ++k) {
        for (int j = RecvRange.Begin(1); j < RecvRange.End(1); ++j) {
          for (int i = RecvRange.Begin(0); i < RecvRange.End(0); ++i) {
            long long iPoint = ExtendedIndexer.ToIndex(i,j,k);
            RecvIndices.Append(iPoint);
          }
        }
      }
    }
  } else {
    field<bool> RecvMask(ExtendedRange);
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      RecvMask.Fill(false);
      const range &NeighborLocalRange = Neighbors[iNeighbor].Value().LocalRange;
      array<long long> &RecvIndices = NeighborRecvIndices_(iNeighbor);
      long long NumRecvPoints = 0;
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (!LocalRange.Contains(Point)) {
              if (Cart.Range().Contains(Point)) {
                RecvMask(Point) = NeighborLocalRange.Contains(Point);
              } else {
                tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
                RecvMask(Point) = NeighborLocalRange.Contains(AdjustedPoint);
              }
              if (RecvMask(Point)) {
                ++NumRecvPoints;
              }
            }
          }
        }
      }
      RecvIndices.Reserve(NumRecvPoints);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (RecvMask(Point)) {
              long long iPoint = ExtendedIndexer.ToIndex(Point);
              RecvIndices.Append(iPoint);
            }
          }
        }
      }
    }
  }

  // Source and destination may sometimes be on the same rank (e.g. if for a given
  // direction periodic==true and #procs==1)
  if (!GlobalRange.Includes(ExtendedRange)) {
    field<bool> LocalToLocalMask(ExtendedRange, false);
    long long NumLocalToLocal = 0;
    for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
      for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
        for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          if (!LocalRange.Contains(Point)) {
            tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
            if (LocalRange.Contains(AdjustedPoint)) {
              ++NumLocalToLocal;
            }
          }
        }
      }
    }
    LocalToLocalSourceIndices_.Reserve(NumLocalToLocal);
    LocalToLocalDestIndices_.Reserve(NumLocalToLocal);
    for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
      for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
        for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          if (!LocalRange.Contains(Point)) {
            tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
            if (LocalRange.Contains(AdjustedPoint)) {
              long long iSourcePoint = ExtendedIndexer.ToIndex(AdjustedPoint);
              long long iDestPoint = ExtendedIndexer.ToIndex(Point);
              LocalToLocalSourceIndices_.Append(iSourcePoint);
              LocalToLocalDestIndices_.Append(iDestPoint);
            }
          }
        }
      }
    }
  }

}

}

halo::halo(std::shared_ptr<context> Context, const cart &Cart, comm_view Comm, const range
  &LocalRange, const range &ExtendedRange, const map<int,decomp_info> &Neighbors):
  Context_(std::move(Context)),
  Comm_(Comm)
{

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.StartSync(TOTAL_TIME, Comm);
  Profiler.Start(SETUP_TIME);

  HaloMap_ = halo_map(Cart, LocalRange, ExtendedRange, Neighbors);

  Profiler.Stop(SETUP_TIME);
  Profiler.Stop(TOTAL_TIME);

}

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

}
