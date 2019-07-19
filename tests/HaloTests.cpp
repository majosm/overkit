// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Halo.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "support/Decomp.hpp"

#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Partition.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <memory>

using testing::ElementsAre;
using testing::ElementsAreArray;

class HaloTests : public tests::mpi_test {};

using support::CartesianDecomp;

TEST_F(HaloTests, Exchange) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize1 = CreateSubsetComm(TestComm(), TestComm().Rank() < 1);
  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  auto CreateCart = [](int NumDims, bool IsPeriodic) -> ovk::cart {
    ovk::range GlobalRange = ovk::MakeEmptyRange(NumDims);
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      GlobalRange.End(iDim) = 20;
    }
    ovk::tuple<bool> Periodic = {false,false,false};
    if (IsPeriodic) Periodic[NumDims-1] = true;
    return {NumDims, GlobalRange, Periodic, ovk::periodic_storage::UNIQUE};
  };

  auto CreateNeighbors = [](const ovk::cart &Cart, ovk::comm_view Comm,
    const ovk::range &LocalRange, const ovk::range &ExtendedRange) -> ovk::array<
    ovk::core::partition_info> {
    ovk::core::partition_hash Hash(Cart.Dimension(), Comm, 1, ovk::array<ovk::range>({1},
      {LocalRange}), ovk::array<int>({1}, {1}));
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    return ovk::core::RetrievePartitionInfo(Comm, NeighborRanks, LocalRange, ExtendedRange);
  };

  auto CreateBeforeDataInt = [](ovk::comm_view Comm, const ovk::range &LocalRange, const
    ovk::range &ExtendedRange) -> ovk::field<int> {
    ovk::field<int> BeforeData(ExtendedRange, -1);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          BeforeData(i,j,k) = Comm.Rank();
        }
      }
    }
    return BeforeData;
  };

  auto CreateAfterDataInt = [](const ovk::cart &Cart, ovk::comm_view Comm,
    const ovk::range &LocalRange, const ovk::range &ExtendedRange, const ovk::array<
    ovk::core::partition_info> &Neighbors) -> ovk::field<int> {
    ovk::field<int> AfterData(ExtendedRange);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          AfterData(i,j,k) = Comm.Rank();
        }
      }
    }
    for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
      for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
        for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
          ovk::tuple<int> Point = {i,j,k};
          if (!LocalRange.Contains(Point)) {
            ovk::tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
            if (LocalRange.Contains(AdjustedPoint)) {
              AfterData(Point) = Comm.Rank();
            } else {
              for (auto &Neighbor : Neighbors) {
                if (Neighbor.LocalRange.Contains(AdjustedPoint)) {
                  AfterData(Point) = Neighbor.Rank;
                  break;
                }
              }
            }
          }
        }
      }
    }
    return AfterData;
  };

  auto CreateBeforeDataDouble = [](ovk::comm_view Comm, const ovk::range &LocalRange, const
    ovk::range &ExtendedRange) -> ovk::field<double> {
    ovk::field<double> BeforeData(ExtendedRange, 0.);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          BeforeData(i,j,k) = double(Comm.Rank())/double(ovk::Max(Comm.Size()-1,1));
        }
      }
    }
    return BeforeData;
  };

  auto CreateAfterDataDouble = [](const ovk::cart &Cart, ovk::comm_view Comm,
    const ovk::range &LocalRange, const ovk::range &ExtendedRange, const ovk::array<
    ovk::core::partition_info> &Neighbors) -> ovk::field<double> {
    ovk::field<double> AfterData(ExtendedRange);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          AfterData(i,j,k) = double(Comm.Rank())/double(ovk::Max(Comm.Size()-1,1));
        }
      }
    }
    for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
      for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
        for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
          ovk::tuple<int> Point = {i,j,k};
          if (!LocalRange.Contains(Point)) {
            ovk::tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
            if (LocalRange.Contains(AdjustedPoint)) {
              AfterData(Point) = double(Comm.Rank())/double(ovk::Max(Comm.Size()-1,1));
            } else {
              for (auto &Neighbor : Neighbors) {
                if (Neighbor.LocalRange.Contains(AdjustedPoint)) {
                  AfterData(Point) = double(Neighbor.Rank)/double(ovk::Max(Comm.Size()-1,1));
                  break;
                }
              }
            }
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
    ovk::cart Cart = CreateCart(2, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = LocalRange;
    ovk::array<ovk::core::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, non-periodic
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, false);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::core::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, periodic
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, true);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::core::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Parallel, non-periodic
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::core::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Parallel, periodic
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, true);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::core::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, multiple simultaneous exchanges
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, false);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::core::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data1 = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    ovk::field<double> Data2 = CreateBeforeDataDouble(Comm, LocalRange, ExtendedRange);
    ovk::array<ovk::request> Requests({2});
    Requests(0) = Halo.Exchange(Data1);
    Requests(1) = Halo.Exchange(Data2);
    ovk::WaitAll(Requests);
    ovk::field<int> ExpectedData1 = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    ovk::field<double> ExpectedData2 = CreateAfterDataDouble(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data1, ElementsAreArray(ExpectedData1));
    EXPECT_THAT(Data2, ElementsAreArray(ExpectedData2));
  }

  // Parallel, multiple simultaneous exchanges
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::core::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::core::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data1 = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    ovk::field<double> Data2 = CreateBeforeDataDouble(Comm, LocalRange, ExtendedRange);
    ovk::array<ovk::request> Requests({2});
    Requests(0) = Halo.Exchange(Data1);
    Requests(1) = Halo.Exchange(Data2);
    ovk::WaitAll(Requests);
    ovk::field<int> ExpectedData1 = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    ovk::field<double> ExpectedData2 = CreateAfterDataDouble(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data1, ElementsAreArray(ExpectedData1));
    EXPECT_THAT(Data2, ElementsAreArray(ExpectedData2));
  }

}
