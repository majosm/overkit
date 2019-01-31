// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "MPITest.hpp"

#include <ovk/core/Comm.hpp>

#include <mpi.h>

using ovk::core::comm;
using ovk::core::CreateSubsetComm;

class ExampleTests : public tests::mpi_test {};

TEST_F(ExampleTests, OneEqualsOneParallel) {

  ASSERT_GE(Comm().Size(), 4);
  comm TestComm = CreateSubsetComm(Comm(), Comm().Rank() < 4);
  if (!TestComm) return;

  EXPECT_EQ(1, 1);

}

TEST_F(ExampleTests, OneEqualsOneSerial) {

  if (Comm().Rank() > 0) return;

  EXPECT_EQ(1, 1);

}
