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
#include <ovk/core/Decomp.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

using testing::ElementsAre;

class PartitionTests : public tests::mpi_test {};

using support::CartesianDecomp;

TEST_F(PartitionTests, Create) {

  ASSERT_GE(TestComm().Size(), 9);

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(TestComm())
    .SetStatusLoggingThreshold(0)
  ));

  ovk::comm CommOfSize9 = CreateSubsetComm(TestComm(), TestComm().Rank() < 9);

  if (CommOfSize9) {

    ovk::cart Cart(2, {{-1,0,0},{21,20,1}}, {false,true,false}, ovk::periodic_storage::DUPLICATED);
    ovk::comm Comm = ovk::CreateCartComm(CommOfSize9, 2, {3,3,1}, Cart.Periodic());
    ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), Comm);

    ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), Comm,
      LocalRange);
    ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, Comm, LocalRange, DecompHash);

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

    // Sanity check
    int MinValue = Comm.Size();
    for (auto &Value : Data) {
      MinValue = ovk::Min(MinValue, Value);
    }
    EXPECT_EQ(MinValue, -1);

    Partition.Exchange(Data);

    // Simple check since halo is tested elsewhere
    MinValue = Comm.Size();
    for (auto &Value : Data) {
      MinValue = ovk::Min(MinValue, Value);
    }
    EXPECT_GE(MinValue, 0);

  }

}
