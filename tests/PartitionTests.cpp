// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Partition.hpp>

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
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

using testing::ElementsAre;
using testing::ElementsAreArray;

class PartitionTests : public tests::mpi_test {};

using support::CartesianDecomp;
using support::TriangularDecomp;

TEST_F(PartitionTests, ExtendLocalRange) {

  if (TestComm().Rank() != 0) return;

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

  // Extend by 0, 2D
  {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range LocalRange({5,5,0}, {15,15,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 0);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(5,5,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(15,15,1));
  }

  // Extend by 0, 3D
  {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range LocalRange({5,5,5}, {15,15,15});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 0);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(5,5,5));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(15,15,15));
  }

  // Full range, non-periodic, 2D
  {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, Cart.Range(), 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,20,1));
  }

  // Full range, non-periodic, 3D
  {
    ovk::cart Cart = CreateCart(3, false, false);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, Cart.Range(), 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,20,20));
  }

  // Full range, periodic, 2D
  {
    ovk::cart Cart = CreateCart(2, true, false);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, Cart.Range(), 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,-2,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,23,1));
  }

  // Full range, periodic, 3D
  {
    ovk::cart Cart = CreateCart(3, true, false);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, Cart.Range(), 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,0,-2));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,20,23));
  }

  // Lower boundary, non-periodic, 2D
  {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range LocalRange({0,0,0}, {10,10,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(12,12,1));
  }

  // Lower boundary, non-periodic, 3D
  {
    ovk::cart Cart = CreateCart(3, false, false);
    ovk::range LocalRange({0,0,0}, {10,10,10});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(12,12,12));
  }

  // Lower boundary, periodic, unique, 2D
  {
    ovk::cart Cart = CreateCart(2, true, false);
    ovk::range LocalRange({0,0,0}, {10,10,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,-2,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(12,12,1));
  }

  // Lower boundary, periodic, unique, 3D
  {
    ovk::cart Cart = CreateCart(3, true, false);
    ovk::range LocalRange({0,0,0}, {10,10,10});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,0,-2));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(12,12,12));
  }

  // Lower boundary, periodic, duplicated, 2D
  {
    ovk::cart Cart = CreateCart(2, true, true);
    ovk::range LocalRange({0,0,0}, {10,10,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,-2,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(12,12,1));
  }

  // Lower boundary, periodic, duplicated, 3D
  {
    ovk::cart Cart = CreateCart(3, true, true);
    ovk::range LocalRange({0,0,0}, {10,10,10});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(0,0,-2));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(12,12,12));
  }

  // Upper boundary, non-periodic, 2D
  {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range LocalRange({10,10,0}, {20,20,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(8,8,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,20,1));
  }

  // Upper boundary, non-periodic, 3D
  {
    ovk::cart Cart = CreateCart(3, false, false);
    ovk::range LocalRange({10,10,10}, {20,20,20});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(8,8,8));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,20,20));
  }

  // Upper boundary, periodic, unique, 2D
  {
    ovk::cart Cart = CreateCart(2, true, false);
    ovk::range LocalRange({10,10,0}, {20,20,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(8,8,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,23,1));
  }

  // Upper boundary, periodic, unique, 3D
  {
    ovk::cart Cart = CreateCart(3, true, false);
    ovk::range LocalRange({10,10,10}, {20,20,20});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(8,8,8));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,20,23));
  }

  // Upper boundary, periodic, duplicated, 2D
  {
    ovk::cart Cart = CreateCart(2, true, true);
    ovk::range LocalRange({10,10,0}, {20,20,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(8,8,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,22,1));
  }

  // Upper boundary, periodic, duplicated, 3D
  {
    ovk::cart Cart = CreateCart(3, true, true);
    ovk::range LocalRange({10,10,10}, {20,20,20});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(8,8,8));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(20,20,22));
  }

  // Interior, non-periodic, 2D
  {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range LocalRange({5,5,0}, {15,15,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(3,3,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(17,17,1));
  }

  // Interior, non-periodic, 3D
  {
    ovk::cart Cart = CreateCart(3, false, false);
    ovk::range LocalRange({5,5,5}, {15,15,15});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(3,3,3));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(17,17,17));
  }

  // Interior, periodic, unique, 2D
  {
    ovk::cart Cart = CreateCart(2, true, false);
    ovk::range LocalRange({5,5,0}, {15,15,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(3,3,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(17,17,1));
  }

  // Interior, periodic, unique, 3D
  {
    ovk::cart Cart = CreateCart(3, true, false);
    ovk::range LocalRange({5,5,5}, {15,15,15});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(3,3,3));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(17,17,17));
  }

  // Interior, periodic, duplicated, 2D
  {
    ovk::cart Cart = CreateCart(2, true, true);
    ovk::range LocalRange({5,5,0}, {15,15,1});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(3,3,0));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(17,17,1));
  }

  // Interior, periodic, duplicated, 3D
  {
    ovk::cart Cart = CreateCart(3, true, true);
    ovk::range LocalRange({5,5,5}, {15,15,15});
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    EXPECT_THAT(ExtendedRange.Begin(), ElementsAre(3,3,3));
    EXPECT_THAT(ExtendedRange.End(), ElementsAre(17,17,17));
  }

}

TEST_F(PartitionTests, DetectNeighbors) {

  ASSERT_GE(TestComm().Size(), 14);

  ovk::comm CommOfSize1 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 1);
  ovk::comm CommOfSize6 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 6);
  ovk::comm CommOfSize12 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 12);
  ovk::comm CommOfSize14 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 14);

  auto CreateCart = [](int NumDims, bool IsPeriodic) -> ovk::cart {
    ovk::range GlobalRange = ovk::MakeEmptyRange(NumDims);
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      GlobalRange.End(iDim) = 20;
    }
    ovk::tuple<bool> Periodic = {false,false,false};
    if (IsPeriodic) Periodic[NumDims-1] = true;
    return {NumDims, GlobalRange, Periodic, ovk::periodic_storage::UNIQUE};
  };

  auto CreatePartitionHash = [](int NumDims, ovk::comm_view Comm, const ovk::range
    &LocalRange) -> ovk::core::partition_hash {
    ovk::core::partition_hash Hash(NumDims, Comm, 1, ovk::array<ovk::range>({1}, {LocalRange}),
      ovk::array<int>({1}, {1}));
    return Hash;
  };

  // 1 rank, non-periodic, 2D
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, false);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, Cart.Range());
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, Cart.Range(), Hash);
    EXPECT_EQ(NeighborRanks.Count(), 0);
  }

  // 1 rank, non-periodic, 3D
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(3, false);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, Cart.Range());
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, Cart.Range(), Hash);
    EXPECT_EQ(NeighborRanks.Count(), 0);
  }

  // 1 rank, periodic, 2D
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, true);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, Cart.Range());
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, Cart.Range(), Hash);
    EXPECT_EQ(NeighborRanks.Count(), 0);
  }

  // 1 rank, periodic, 3D
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(3, true);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, Cart.Range());
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, Cart.Range(), Hash);
    EXPECT_EQ(NeighborRanks.Count(), 0);
  }

  // Cartesian, non-periodic, 2D
  if (CommOfSize6) {
    ovk::cart Cart = CreateCart(2, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize6, 2, {2,3,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    switch (Comm.Rank()) {
    // Lower corner
    case 0:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,3,4}));
      break;
    // Middle
    case 1:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,2,3,4,5}));
      break;
    // Upper corner
    case 5:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2,4}));
      break;
    default:
      break;
    }
  }

  // Cartesian, non-periodic, 3D
  if (CommOfSize12) {
    ovk::cart Cart = CreateCart(3, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize12, 3, {2,2,3}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    switch (Comm.Rank()) {
    // Lower corner
    case 0:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,3,4,6,7,9,10}));
      break;
    // Middle
    case 7:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,2,3,4,5,6,8,9,10,11}));
      break;
    // Upper corner
    case 11:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2,4,5,7,8,10}));
      break;
    default:
      break;
    }
  }

  // Cartesian, periodic, 2D
  if (CommOfSize6) {
    ovk::cart Cart = CreateCart(2, true);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize6, 2, {2,3,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    switch (Comm.Rank()) {
    // Lower corner
    case 0:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2,3,4,5}));
      break;
    // Middle
    case 1:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,2,3,4,5}));
      break;
    // Upper corner
    case 5:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,2,3,4}));
      break;
    default:
      break;
    }
  }

  // Cartesian, periodic, 3D
  if (CommOfSize12) {
    ovk::cart Cart = CreateCart(3, true);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize12, 3, {2,2,3}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    switch (Comm.Rank()) {
    // Lower corner
    case 0:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2,3,4,5,6,7,8,9,10,11}));
      break;
    // Middle
    case 7:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,2,3,4,5,6,8,9,10,11}));
      break;
    // Upper corner
    case 11:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,2,3,4,5,6,7,8,9,10}));
      break;
    default:
      break;
    }
  }

  // Non-cartesian, non-periodic, 2D
  if (CommOfSize6) {
    ovk::comm_view Comm = CommOfSize6;
    ovk::cart Cart = CreateCart(2, false);
    ovk::range LocalRange = TriangularDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    switch (Comm.Rank()) {
    // Bottom
    case 0:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2}));
      break;
    // Middle
    case 2:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,4,5}));
      break;
    // Top
    case 4:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2,3,5}));
      break;
    default:
      break;
    }
  }

  // Non-cartesian, non-periodic, 3D
  if (CommOfSize14) {
    ovk::comm_view Comm = CommOfSize14;
    ovk::cart Cart = CreateCart(3, false);
    ovk::range LocalRange = TriangularDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    switch (Comm.Rank()) {
    // Bottom
    case 0:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2,3,4}));
      break;
    // Middle
    case 4:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,2,3,9,10,12,13}));
      break;
    // Top
    case 9:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2,3,4,5,6,7,8,10,11,12,13}));
      break;
    default:
      break;
    }
  }

  // Non-cartesian, periodic, 2D
  if (CommOfSize6) {
    ovk::comm_view Comm = CommOfSize6;
    ovk::cart Cart = CreateCart(2, true);
    ovk::range LocalRange = TriangularDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    switch (Comm.Rank()) {
    // Bottom
    case 0:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2,3,4,5}));
      break;
    // Middle
    case 2:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,4,5}));
      break;
    // Top
    case 4:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,2,3,5}));
      break;
    default:
      break;
    }
  }

  // Non-cartesian, periodic, 3D
  if (CommOfSize14) {
    ovk::comm_view Comm = CommOfSize14;
    ovk::cart Cart = CreateCart(3, true);
    ovk::range LocalRange = TriangularDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::core::partition_hash Hash = CreatePartitionHash(Cart.Dimension(), Comm, LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    switch (Comm.Rank()) {
    // Bottom
    case 0:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({1,2,3,4,5,6,7,8,9,10,11,12,13}));
      break;
    // Middle
    case 4:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,2,3,9,10,12,13}));
      break;
    // Top
    case 9:
      EXPECT_THAT(NeighborRanks, ElementsAreArray({0,1,2,3,4,5,6,7,8,10,11,12,13}));
      break;
    default:
      break;
    }
  }

}

TEST_F(PartitionTests, RetrievePartitionInfo) {

  ASSERT_GE(TestComm().Size(), 6);

  ovk::comm CommOfSize6 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 6);

  // Self
  if (CommOfSize6) {
    ovk::cart Cart(2, {{0,0,0}, {20,20,1}}, {false,false,false}, ovk::periodic_storage::UNIQUE);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize6, 2, {2,3,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 1);
    ovk::array<int> Ranks({1}, Comm.Rank());
    ovk::array<ovk::partition_info> PartitionInfo = ovk::partition_internal::
      RetrievePartitionInfo(Comm, Ranks, LocalRange, ExtendedRange);
    EXPECT_EQ(PartitionInfo.Count(), 1);
    EXPECT_EQ(PartitionInfo[0].Rank, Comm.Rank());
    EXPECT_THAT(PartitionInfo[0].LocalRange.Begin(), ElementsAreArray(LocalRange.Begin()));
    EXPECT_THAT(PartitionInfo[0].LocalRange.End(), ElementsAreArray(LocalRange.End()));
    EXPECT_THAT(PartitionInfo[0].ExtendedRange.Begin(), ElementsAreArray(ExtendedRange.Begin()));
    EXPECT_THAT(PartitionInfo[0].ExtendedRange.End(), ElementsAreArray(ExtendedRange.End()));
  }

  // Others
  if (CommOfSize6) {
    ovk::cart Cart(2, {{0,0,0}, {20,20,1}}, {false,false,false}, ovk::periodic_storage::UNIQUE);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize6, 2, {2,3,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 1);
    ovk::array<int> Ranks;
    if (Comm.Rank() == 3) {
      Ranks.Resize({2});
      Ranks[0] = 0;
      Ranks[1] = 2;
    }
    ovk::array<ovk::partition_info> PartitionInfo = ovk::partition_internal::
      RetrievePartitionInfo(Comm, Ranks, LocalRange, ExtendedRange);
    switch (Comm.Rank()) {
    case 3:
      EXPECT_EQ(PartitionInfo.Count(), 2);
      EXPECT_EQ(PartitionInfo[0].Rank, 0);
      EXPECT_THAT(PartitionInfo[0].LocalRange.Begin(), ElementsAre(0,0,0));
      EXPECT_THAT(PartitionInfo[0].LocalRange.End(), ElementsAre(10,7,1));
      EXPECT_THAT(PartitionInfo[0].ExtendedRange.Begin(), ElementsAre(0,0,0));
      EXPECT_THAT(PartitionInfo[0].ExtendedRange.End(), ElementsAre(11,8,1));
      EXPECT_EQ(PartitionInfo[1].Rank, 2);
      EXPECT_THAT(PartitionInfo[1].LocalRange.Begin(), ElementsAre(0,14,0));
      EXPECT_THAT(PartitionInfo[1].LocalRange.End(), ElementsAre(10,20,1));
      EXPECT_THAT(PartitionInfo[1].ExtendedRange.Begin(), ElementsAre(0,13,0));
      EXPECT_THAT(PartitionInfo[1].ExtendedRange.End(), ElementsAre(11,20,1));
      break;
    default:
      EXPECT_EQ(PartitionInfo.Count(), 0);
      break;
    }
  }

}

TEST_F(PartitionTests, HaloExchange) {

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
    const ovk::range &LocalRange, const ovk::range &ExtendedRange) -> ovk::array<
    ovk::partition_info> {
    ovk::core::partition_hash Hash(Cart.Dimension(), Comm, 1, ovk::array<ovk::range>({1},
      {LocalRange}), ovk::array<int>({1}, {1}));
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);
    return ovk::partition_internal::RetrievePartitionInfo(Comm, NeighborRanks, LocalRange,
      ExtendedRange);
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
    ovk::partition_info> &Neighbors) -> ovk::field<int> {
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
    ovk::partition_info> &Neighbors) -> ovk::field<double> {
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
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = LocalRange;
    ovk::array<ovk::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::partition_internal::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, non-periodic
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::partition_internal::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, periodic, unique
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, true, false);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::partition_internal::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, periodic, duplicated
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, true, true);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::partition_internal::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Parallel, non-periodic
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::partition_internal::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Parallel, periodic, unique
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, true, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::partition_internal::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Parallel, periodic, duplicated
  if (CommOfSize4) {
    ovk::cart Cart = CreateCart(2, true, true);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::partition_internal::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
    ovk::field<int> Data = CreateBeforeDataInt(Comm, LocalRange, ExtendedRange);
    Halo.Exchange(Data);
    ovk::field<int> ExpectedData = CreateAfterDataInt(Cart, Comm, LocalRange, ExtendedRange,
      Neighbors);
    EXPECT_THAT(Data, ElementsAreArray(ExpectedData));
  }

  // Serial, multiple simultaneous exchanges
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::range LocalRange = Cart.Range();
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::partition_internal::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
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
    ovk::cart Cart = CreateCart(2, false, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize4, 2, {2,2,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 2);
    ovk::array<ovk::partition_info> Neighbors = CreateNeighbors(Cart, Comm, LocalRange,
      ExtendedRange);
    ovk::partition_internal::halo Halo(Context, Cart, Comm, LocalRange, ExtendedRange, Neighbors);
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

TEST_F(PartitionTests, ConstructPartition) {

  ASSERT_GE(TestComm().Size(), 9);

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(TestComm())
  ));

  ovk::comm CommOfSize6 = CreateSubsetComm(TestComm(), TestComm().Rank() < 9);

  if (CommOfSize6) {

    ovk::cart Cart(2, {{-1,0,0},{21,20,1}}, {false,true,false}, ovk::periodic_storage::DUPLICATED);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize6, 2, {3,3,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);

    ovk::core::partition_hash Hash(Cart.Dimension(), Comm, 1, ovk::array<ovk::range>({1},
      {LocalRange}), ovk::array<int>({1}, {1}));
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, Hash);

    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 1);

    ovk::partition Partition(Context, Cart, Comm, LocalRange, ExtendedRange, 2, NeighborRanks);

    EXPECT_EQ(Partition.Cart().Dimension(), 2);
    EXPECT_THAT(Partition.Cart().Range().Begin(), ElementsAre(-1,0,0));
    EXPECT_THAT(Partition.Cart().Range().End(), ElementsAre(21,20,1));
    EXPECT_THAT(Partition.Cart().Periodic(), ElementsAre(false,true,false));
    EXPECT_EQ(Partition.Cart().PeriodicStorage(), ovk::periodic_storage::DUPLICATED);
    EXPECT_EQ(Partition.Comm().Get(), Comm.Get());
    EXPECT_THAT(Partition.GlobalRange().Begin(), ElementsAre(-1,0,0));
    EXPECT_THAT(Partition.GlobalRange().End(), ElementsAre(21,20,1));
    switch (Comm.Rank()) {
    // Lower corner
    case 0:
      EXPECT_THAT(Partition.LocalRange().Begin(), ElementsAre(-1,0,0));
      EXPECT_THAT(Partition.LocalRange().End(), ElementsAre(7,7,1));
      EXPECT_THAT(Partition.ExtendedRange().Begin(), ElementsAre(-1,-1,0));
      EXPECT_THAT(Partition.ExtendedRange().End(), ElementsAre(8,8,1));
      EXPECT_EQ(Partition.LocalSubregions().Count(), 2);
      EXPECT_THAT(Partition.LocalSubregions()[0].Begin(), ElementsAre(-1,0,0));
      EXPECT_THAT(Partition.LocalSubregions()[0].End(), ElementsAre(7,4,1));
      EXPECT_THAT(Partition.LocalSubregions()[1].Begin(), ElementsAre(-1,4,0));
      EXPECT_THAT(Partition.LocalSubregions()[1].End(), ElementsAre(7,7,1));
      EXPECT_THAT(Partition.ExtendedSubregions()[0].Begin(), ElementsAre(-1,-1,0));
      EXPECT_THAT(Partition.ExtendedSubregions()[0].End(), ElementsAre(8,4,1));
      EXPECT_THAT(Partition.ExtendedSubregions()[1].Begin(), ElementsAre(-1,4,0));
      EXPECT_THAT(Partition.ExtendedSubregions()[1].End(), ElementsAre(8,8,1));
      EXPECT_EQ(Partition.Neighbors().Count(), 5);
      break;
    // Middle
    case 4:
      EXPECT_THAT(Partition.LocalRange().Begin(), ElementsAre(7,7,0));
      EXPECT_THAT(Partition.LocalRange().End(), ElementsAre(14,14,1));
      EXPECT_THAT(Partition.ExtendedRange().Begin(), ElementsAre(6,6,0));
      EXPECT_THAT(Partition.ExtendedRange().End(), ElementsAre(15,15,1));
      EXPECT_EQ(Partition.LocalSubregions().Count(), 2);
      EXPECT_THAT(Partition.LocalSubregions()[0].Begin(), ElementsAre(7,7,0));
      EXPECT_THAT(Partition.LocalSubregions()[0].End(), ElementsAre(14,11,1));
      EXPECT_THAT(Partition.LocalSubregions()[1].Begin(), ElementsAre(7,11,0));
      EXPECT_THAT(Partition.LocalSubregions()[1].End(), ElementsAre(14,14,1));
      EXPECT_THAT(Partition.ExtendedSubregions()[0].Begin(), ElementsAre(6,6,0));
      EXPECT_THAT(Partition.ExtendedSubregions()[0].End(), ElementsAre(15,11,1));
      EXPECT_THAT(Partition.ExtendedSubregions()[1].Begin(), ElementsAre(6,11,0));
      EXPECT_THAT(Partition.ExtendedSubregions()[1].End(), ElementsAre(15,15,1));
      EXPECT_EQ(Partition.Neighbors().Count(), 8);
      break;
    // Upper corner
    case 8:
      EXPECT_THAT(Partition.LocalRange().Begin(), ElementsAre(14,14,0));
      EXPECT_THAT(Partition.LocalRange().End(), ElementsAre(21,20,1));
      EXPECT_THAT(Partition.ExtendedRange().Begin(), ElementsAre(13,13,0));
      EXPECT_THAT(Partition.ExtendedRange().End(), ElementsAre(21,21,1));
      EXPECT_EQ(Partition.LocalSubregions().Count(), 2);
      EXPECT_THAT(Partition.LocalSubregions()[0].Begin(), ElementsAre(14,14,0));
      EXPECT_THAT(Partition.LocalSubregions()[0].End(), ElementsAre(21,17,1));
      EXPECT_THAT(Partition.LocalSubregions()[1].Begin(), ElementsAre(14,17,0));
      EXPECT_THAT(Partition.LocalSubregions()[1].End(), ElementsAre(21,20,1));
      EXPECT_THAT(Partition.ExtendedSubregions()[0].Begin(), ElementsAre(13,13,0));
      EXPECT_THAT(Partition.ExtendedSubregions()[0].End(), ElementsAre(21,17,1));
      EXPECT_THAT(Partition.ExtendedSubregions()[1].Begin(), ElementsAre(13,17,0));
      EXPECT_THAT(Partition.ExtendedSubregions()[1].End(), ElementsAre(21,21,1));
      EXPECT_EQ(Partition.Neighbors().Count(), 5);
      break;
    default:
      break;
    };

    ovk::field<int> Data(ExtendedRange, -1);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          Data(i,j,k) = Comm.Rank();
        }
      }
    }

    Partition.Exchange(Data);

    int MinValue = Comm.Size();
    for (auto &Value : Data) {
      MinValue = ovk::Min(MinValue, Value);
    }
    EXPECT_GE(MinValue, 0);

  }

}
