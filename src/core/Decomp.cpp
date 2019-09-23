// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Decomp.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/DistributedRegionHash.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/ElemSet.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Interval.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

namespace ovk {
namespace core {

decomp_hash CreateDecompHash(int NumDims, comm_view Comm, const range &LocalRange) {
  return {NumDims, Comm, 1, array<range>({1}, {LocalRange}), array<int>({1}, {1})};
}

array<int> DetectNeighbors(const cart &Cart, comm_view Comm, const range &LocalRange, const
  decomp_hash &DecompHash) {

  range ExtendedRange = LocalRange;
  for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
    if (LocalRange.Begin(iDim) != Cart.Range().Begin(iDim) || Cart.Periodic(iDim)) {
      ExtendedRange.Begin(iDim) -= 1;
    }
    if (LocalRange.End(iDim) != Cart.Range().End(iDim) || Cart.Periodic(iDim)) {
      ExtendedRange.End(iDim) += 1;
    }
  }

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

  array<elem<int,2>> ExtendedPointBinIDs({NumExtendedPoints});
  elem_set<int,2> UniqueBinIDs;

  for (long long iPoint = 0; iPoint < NumExtendedPoints; ++iPoint) {
    tuple<int> Point = {
      ExtendedPoints(0,iPoint),
      ExtendedPoints(1,iPoint),
      ExtendedPoints(2,iPoint)
    };
    elem<int,2> BinID = DecompHash.MapToBin(Point);
    ExtendedPointBinIDs(iPoint) = BinID;
    UniqueBinIDs.Insert(BinID);
  }

  map<int,decomp_hash_retrieved_bins> RetrievedBins = DecompHash.RetrieveBins(UniqueBinIDs);

  set<int> UniqueExtendedRanks;

  for (long long iPoint = 0; iPoint < NumExtendedPoints; ++iPoint) {
    tuple<int> Point = {
      ExtendedPoints(0,iPoint),
      ExtendedPoints(1,iPoint),
      ExtendedPoints(2,iPoint)
    };
    const elem<int,2> &BinID = ExtendedPointBinIDs(iPoint);
    int BinRank = BinID(0);
    int iBin = BinID(1);
    const decomp_hash_retrieved_bins &Bins = RetrievedBins(BinRank);
    const interval<long long> &RegionIndicesInterval = Bins.BinRegionIndicesIntervals(iBin);
    for (long long iBinRegionIndex = RegionIndicesInterval.Begin(0); iBinRegionIndex <
      RegionIndicesInterval.End(0); ++iBinRegionIndex) {
      int iRegion = Bins.BinRegionIndices(iBinRegionIndex);
      const decomp_hash_region_data &RegionData = Bins.RegionData(iRegion);
      if (RegionData.Region.Contains(Point)) {
        UniqueExtendedRanks.Insert(RegionData.Rank);
        break;
      }
    }
  }
  UniqueExtendedRanks.Erase(Comm.Rank());

  return {UniqueExtendedRanks};

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

  return ExtendedRange;

}

cart CartPointToCell(const cart &Cart) {

  cart CellCart(Cart.Dimension());

  for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
    CellCart.Range().Begin(iDim) = Cart.Range().Begin(iDim);
    if (Cart.Periodic(iDim) && Cart.PeriodicStorage() == periodic_storage::UNIQUE) {
      CellCart.Range().End(iDim) = Cart.Range().End(iDim);
    } else {
      CellCart.Range().End(iDim) = Cart.Range().End(iDim)-1;
    }
  }

  CellCart.Periodic() = Cart.Periodic();
  CellCart.PeriodicStorage() = periodic_storage::UNIQUE;

  return CellCart;

}

range RangePointToCell(const cart &Cart, const range &Range) {

  range CellRange = MakeEmptyRange(Cart.Dimension());

  for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
    CellRange.Begin(iDim) = Range.Begin(iDim);
    if (Range.End(iDim) != Cart.Range().End(iDim) || (Cart.Periodic(iDim) && Cart.PeriodicStorage()
      == periodic_storage::UNIQUE)) {
      CellRange.End(iDim) = Range.End(iDim);
    } else {
      CellRange.End(iDim) = Range.End(iDim)-1;
    }
  }

  return CellRange;

}

cart CartIncludeExteriorPoint(const cart &Cart) {

  cart CartIncludingExteriorPoint(Cart.Dimension());

  for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
    if (!Cart.Periodic(iDim)) {
      CartIncludingExteriorPoint.Range().Begin(iDim) = Cart.Range().Begin(iDim)-1;
      CartIncludingExteriorPoint.Range().End(iDim) = Cart.Range().End(iDim)+1;
    } else {
      CartIncludingExteriorPoint.Range().Begin(iDim) = Cart.Range().Begin(iDim);
      CartIncludingExteriorPoint.Range().End(iDim) = Cart.Range().End(iDim);
    }
  }

  CartIncludingExteriorPoint.Periodic() = Cart.Periodic();
  CartIncludingExteriorPoint.PeriodicStorage() = Cart.PeriodicStorage();

  return CartIncludingExteriorPoint;

}

range RangeIncludeExteriorPoint(const cart &Cart, const range &Range) {

  range RangeIncludingExteriorPoint = MakeEmptyRange(Cart.Dimension());

  for (int iDim = 0; iDim < Cart.Dimension(); ++iDim) {
    if (!Cart.Periodic(iDim) && Range.Begin(iDim) == Cart.Range().Begin(iDim)) {
      RangeIncludingExteriorPoint.Begin(iDim) = Range.Begin(iDim)-1;
    } else {
      RangeIncludingExteriorPoint.Begin(iDim) = Range.Begin(iDim);
    }
    if (!Cart.Periodic(iDim) && Range.End(iDim) == Cart.Range().End(iDim)) {
      RangeIncludingExteriorPoint.End(iDim) = Range.End(iDim)+1;
    } else {
      RangeIncludingExteriorPoint.End(iDim) = Range.End(iDim);
    }
  }

  return RangeIncludingExteriorPoint;

}

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

}}
