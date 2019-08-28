// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Halo.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "support/Decomp.hpp"

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

using testing::ElementsAreArray;

class HaloTests : public tests::mpi_test {};

using support::CartesianDecomp;

TEST_F(HaloTests, Exchange) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize1 = CreateSubsetComm(TestComm(), TestComm().Rank() < 1);
  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  auto CreateCart = [](int NumDims, bool IsPeriodic, bool Duplicated) -> ovk::cart {
    ovk::range GlobalRange = ovk::MakeEmptyRange(NumDims);
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      GlobalRange.End(iDim) = 20;
    }
    ovk::tuple<bool> Periodic = {false,false,false};
    if (IsPeriodic) Periodic[NumDims-1] = true;
    ovk::periodic_storage PeriodicStorage = ovk::periodic_storage::UNIQUE;
    if (Duplicated) PeriodicStorage = ovk::periodic_storage::DUPLICATED;
    return {NumDims, GlobalRange, Periodic, PeriodicStorage};
  };

  auto CreateNeighbors = [](const ovk::cart &Cart, ovk::comm_view Comm,
    const ovk::range &LocalRange, const ovk::range &ExtendedRange) -> ovk::map<int,
    ovk::core::decomp_info> {
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);
    return ovk::core::RetrieveDecompInfo(Comm, NeighborRanks, LocalRange, ExtendedRange);
  };

  auto CreateBeforeDataInt = [](const ovk::cart &Cart, const ovk::range &LocalRange, const
    ovk::range &ExtendedRange) -> ovk::field<int> {
    ovk::field<int> BeforeData(ExtendedRange, -1);
    ovk::range_indexer_c<int> GlobalIndexer(Cart.Range());
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ovk::tuple<int> Point = {i,j,k};
          BeforeData(Point) = GlobalIndexer.ToIndex(Point);
        }
      }
    }
    return BeforeData;
  };

  auto CreateAfterDataInt = [](const ovk::cart &Cart, const ovk::range &ExtendedRange) ->
    ovk::field<int> {
    ovk::field<int> AfterData(ExtendedRange);
    ovk::range_indexer_c<int> GlobalIndexer(Cart.Range());
    for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
      for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
        for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
          ovk::tuple<int> Point = {i,j,k};
          if (Cart.Range().Contains(Point)) {
            AfterData(Point) = GlobalIndexer.ToIndex(Point);
          } else {
            ovk::tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
            AfterData(Point) = GlobalIndexer.ToIndex(AdjustedPoint);
          }
        }
      }
    }
    return AfterData;
  };

  auto CreateBeforeDataDouble = [](const ovk::cart &Cart, const ovk::range &LocalRange, const
    ovk::range &ExtendedRange) -> ovk::field<double> {
    ovk::field<double> BeforeData(ExtendedRange, -1.);
    ovk::range_indexer_c<int> GlobalIndexer(Cart.Range());
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ovk::tuple<int> Point = {i,j,k};
          BeforeData(Point) = double(GlobalIndexer.ToIndex(Point))/double(Cart.Range().Count()-1);
        }
      }
    }
    return BeforeData;
  };

  auto CreateAfterDataDouble = [](const ovk::cart &Cart, const ovk::range &ExtendedRange) ->
    ovk::field<double> {
    ovk::field<double> AfterData(ExtendedRange);
    ovk::range_indexer_c<int> GlobalIndexer(Cart.Range());
    for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
      for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
        for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
          ovk::tuple<int> Point = {i,j,k};
          if (Cart.Range().Contains(Point)) {
            AfterData(Point) = double(GlobalIndexer.ToIndex(Point))/double(Cart.Range().Count()-1);
          } else {
            ovk::tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
            AfterData(Point) = double(GlobalIndexer.ToIndex(AdjustedPoint))/double(Cart.Range()
              .Count()-1);
          }
        }
      }
    }
    return AfterData;
  };

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(TestComm())
  ));

  // Extended range same as local range
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = LocalRange;
    ovk::map<int,ovk::core::decomp_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Cart, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, ExtendedRange);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, non-periodic
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::map<int,ovk::core::decomp_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Cart, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, ExtendedRange);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, periodic, unique
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, true, false);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::map<int,ovk::core::decomp_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Cart, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, ExtendedRange);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, periodic, duplicated
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, true, true);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::map<int,ovk::core::decomp_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Cart, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, ExtendedRange);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Parallel, non-periodic
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::map<int,ovk::core::decomp_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Cart, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, ExtendedRange);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Parallel, periodic, unique
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, true, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::map<int,ovk::core::decomp_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Cart, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, ExtendedRange);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Parallel, periodic, duplicated
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, true, true);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::map<int,ovk::core::decomp_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Cart, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, ExtendedRange);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, multiple simultaneous exchanges
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::map<int,ovk::core::decomp_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data1 = CreateBeforeDataInt(Cart, LocalRange, ExtendedRange);
    ovk::field<double> Data2 = CreateBeforeDataDouble(Cart, LocalRange, ExtendedRange);
    ovk::array<ovk::request> Requests({2});
    Requests(0) = Halo.Exchange(Data1);
    Requests(1) = Halo.Exchange(Data2);
    ovk::WaitAll(Requests);
    ovk::field<int> ExpectedData1 = CreateAfterDataInt(Cart, ExtendedRange);
    ovk::field<double> ExpectedData2 = CreateAfterDataDouble(Cart, ExtendedRange);
    EXPECT_THAT(Data1, ElementsAreArray(ExpectedData1));
    EXPECT_THAT(Data2, ElementsAreArray(ExpectedData2));
  }

  // Parallel, multiple simultaneous exchanges
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::map<int,ovk::core::decomp_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data1 = CreateBeforeDataInt(Cart, LocalRange, ExtendedRange);
    ovk::field<double> Data2 = CreateBeforeDataDouble(Cart, LocalRange, ExtendedRange);
    ovk::array<ovk::request> Requests({2});
    Requests(0) = Halo.Exchange(Data1);
    Requests(1) = Halo.Exchange(Data2);
    ovk::WaitAll(Requests);
    ovk::field<int> ExpectedData1 = CreateAfterDataInt(Cart, ExtendedRange);
    ovk::field<double> ExpectedData2 = CreateAfterDataDouble(Cart, ExtendedRange);
    EXPECT_THAT(Data1, ElementsAreArray(ExpectedData1));
    EXPECT_THAT(Data2, ElementsAreArray(ExpectedData2));
  }

}
