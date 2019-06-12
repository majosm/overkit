// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Context.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>

#include <mpi.h>

using testing::ElementsAre;
using testing::ElementsAreArray;

class ContextTests : public tests::mpi_test {};

TEST_F(ContextTests, Create) {

  // Implicit error handling
  {
    ovk::context Context = ovk::CreateContext(ovk::context::params()
      .SetComm(TestComm())
    );
  }

  // Explicit error handling
  {
    ovk::error Error;
    auto MaybeContext = ovk::CreateContext(ovk::context::params()
      .SetComm(TestComm())
    , Error);
    ASSERT_EQ(Error, ovk::error::NONE);
    ASSERT_TRUE(MaybeContext.Present());
  }

}
