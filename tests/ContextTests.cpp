// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Context.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>

#include <mpi.h>

using testing::ElementsAre;
using testing::ElementsAreArray;

class ContextTests : public tests::mpi_test {};

TEST_F(ContextTests, Create) {

  // Implicit error handling
  {
    ovk::context Context = ovk::CreateContext(ovk::context::params()
      .SetComm(TestComm())
      .SetStatusLoggingThreshold(0)
    );
  }

  // Explicit error handling
  {
    ovk::captured_error Error;
    auto MaybeContext = ovk::CreateContext(ovk::context::params()
      .SetComm(TestComm())
      .SetStatusLoggingThreshold(0)
    , Error);
    ASSERT_FALSE(Error);
    ASSERT_TRUE(MaybeContext.Present());
  }

}
