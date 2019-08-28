// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Decomp.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "support/Decomp.hpp"

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

using testing::ElementsAre;
using testing::ElementsAreArray;

class DecompTests : public tests::mpi_test {};

using support::CartesianDecomp;
using support::TriangularDecomp;

TEST_F(DecompTests, DetectNeighbors) {

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

  // 1 rank, non-periodic, 2D
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, false);
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      Cart.Range());
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, Cart.Range(),
      DecompHash);
    EXPECT_EQ(NeighborRanks.Count(), 0);
  }

  // 1 rank, non-periodic, 3D
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(3, false);
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      Cart.Range());
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, Cart.Range(),
      DecompHash);
    EXPECT_EQ(NeighborRanks.Count(), 0);
  }

  // 1 rank, periodic, 2D
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(2, true);
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      Cart.Range());
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, Cart.Range(),
      DecompHash);
    EXPECT_EQ(NeighborRanks.Count(), 0);
  }

  // 1 rank, periodic, 3D
  if (CommOfSize1) {
    ovk::comm_view Comm = CommOfSize1;
    ovk::cart Cart = CreateCart(3, true);
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      Cart.Range());
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, Cart.Range(),
      DecompHash);
    EXPECT_EQ(NeighborRanks.Count(), 0);
  }

  // Cartesian, non-periodic, 2D
  if (CommOfSize6) {
    ovk::cart Cart = CreateCart(2, false);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize6, 2, {2,3,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);
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
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);
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
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);
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
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);
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
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);
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
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);
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
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);
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
    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);
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

TEST_F(DecompTests, ExtendLocalRange) {

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

TEST_F(DecompTests, RetrieveDecompInfo) {

  ASSERT_GE(TestComm().Size(), 6);

  ovk::comm CommOfSize6 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 6);

  // Self
  if (CommOfSize6) {
    ovk::cart Cart(2, {{0,0,0}, {20,20,1}}, {false,false,false}, ovk::periodic_storage::UNIQUE);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize6, 2, {2,3,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);
    ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 1);
    ovk::array<int> Ranks({1}, Comm.Rank());
    ovk::map<int,ovk::core::decomp_info> DecompInfo = ovk::core::RetrieveDecompInfo(Comm, Ranks,
      LocalRange, ExtendedRange);
    EXPECT_EQ(DecompInfo.Count(), 1);
    EXPECT_EQ(DecompInfo[0].Key(), Comm.Rank());
    EXPECT_THAT(DecompInfo[0].Value().LocalRange.Begin(), ElementsAreArray(LocalRange.Begin()));
    EXPECT_THAT(DecompInfo[0].Value().LocalRange.End(), ElementsAreArray(LocalRange.End()));
    EXPECT_THAT(DecompInfo[0].Value().ExtendedRange.Begin(), ElementsAreArray(ExtendedRange.Begin()));
    EXPECT_THAT(DecompInfo[0].Value().ExtendedRange.End(), ElementsAreArray(ExtendedRange.End()));
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
    ovk::map<int,ovk::core::decomp_info> DecompInfo = ovk::core::RetrieveDecompInfo(Comm, Ranks,
      LocalRange, ExtendedRange);
    switch (Comm.Rank()) {
    case 3:
      EXPECT_EQ(DecompInfo.Count(), 2);
      EXPECT_EQ(DecompInfo[0].Key(), 0);
      EXPECT_THAT(DecompInfo[0].Value().LocalRange.Begin(), ElementsAre(0,0,0));
      EXPECT_THAT(DecompInfo[0].Value().LocalRange.End(), ElementsAre(10,7,1));
      EXPECT_THAT(DecompInfo[0].Value().ExtendedRange.Begin(), ElementsAre(0,0,0));
      EXPECT_THAT(DecompInfo[0].Value().ExtendedRange.End(), ElementsAre(11,8,1));
      EXPECT_EQ(DecompInfo[1].Key(), 2);
      EXPECT_THAT(DecompInfo[1].Value().LocalRange.Begin(), ElementsAre(0,14,0));
      EXPECT_THAT(DecompInfo[1].Value().LocalRange.End(), ElementsAre(10,20,1));
      EXPECT_THAT(DecompInfo[1].Value().ExtendedRange.Begin(), ElementsAre(0,13,0));
      EXPECT_THAT(DecompInfo[1].Value().ExtendedRange.End(), ElementsAre(11,20,1));
      break;
    default:
      EXPECT_EQ(DecompInfo.Count(), 0);
      break;
    }
  }

}
