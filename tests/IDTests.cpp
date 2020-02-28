// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/ID.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Set.hpp>

#include <mpi.h>

class IDTests : public tests::mpi_test {};

TEST_F(IDTests, NextAvailableID) {

  if (TestComm().Rank() != 0) return;

  // First ID > 0
  {
    ovk::set<int> Set = {1, 2};
    EXPECT_EQ(ovk::NextAvailableID(Set), 0);
  }

  // Gap
  {
    ovk::set<int> Set = {0, 2};
    EXPECT_EQ(ovk::NextAvailableID(Set), 1);
  }

  // End
  {
    ovk::set<int> Set = {0, 1};
    EXPECT_EQ(ovk::NextAvailableID(Set), 2);
  }

}
