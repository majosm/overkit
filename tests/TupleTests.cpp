// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Tuple.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>

#include <mpi.h>

using testing::ElementsAre;

class TupleTests : public tests::mpi_test {};

TEST_F(TupleTests, MakeUniform) {

  if (TestComm().Rank() != 0) return;

  // Default pad value
  {
    ovk::tuple<int> Tuple = ovk::MakeUniformTuple<int>(2, 5);
    EXPECT_THAT(Tuple, ElementsAre(5,5,0));
  }

  // Explicit pad value
  {
    ovk::tuple<int> Tuple = ovk::MakeUniformTuple<int>(2, 5, 1);
    EXPECT_THAT(Tuple, ElementsAre(5,5,1));
  }

}
