// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Range.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

using testing::ElementsAre;

class RangeTests : public tests::mpi_test {};

TEST_F(RangeTests, MakeEmpty) {

  if (TestComm().Rank() != 0) return;

  // 1D
  {
    ovk::range Range = ovk::MakeEmptyRange(1);
    EXPECT_THAT(Range.Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(Range.End(), ElementsAre(0,1,1));
  }

  // 2D
  {
    ovk::range Range = ovk::MakeEmptyRange(2);
    EXPECT_THAT(Range.Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(Range.End(), ElementsAre(0,0,1));
  }

  // 3D
  {
    ovk::range Range = ovk::MakeEmptyRange(3);
    EXPECT_THAT(Range.Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(Range.End(), ElementsAre(0,0,0));
  }

}

// TODO: Add remaining tests
