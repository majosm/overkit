// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>

#include <mpi.h>

using ovk::core::comm;
using ovk::core::CreateSubsetComm;

class ExampleTests : public tests::mpi_test {};

TEST_F(ExampleTests, OneEqualsOneParallel) {

  ASSERT_GE(TestComm().Size(), 4);
  comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);
  if (!Comm) return;

  EXPECT_EQ(1, 1);

}

TEST_F(ExampleTests, OneEqualsOneSerial) {

  if (TestComm().Rank() > 0) return;

  EXPECT_EQ(1, 1);

}
