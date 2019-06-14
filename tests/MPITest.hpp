// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef TESTS_MPI_TEST_HPP_INCLUDED
#define TESTS_MPI_TEST_HPP_INCLUDED

#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>

#include <mpi.h>

namespace tests {

class mpi_test : public testing::Test {

public:

  mpi_test() = default;

  virtual void SetUp() override {
    Comm_ = MPI_COMM_WORLD;
  }

  virtual void TearDown() override {
    Comm_.Reset();
  }

  const ovk::comm_view &TestComm() { return Comm_; }

private:

  ovk::comm_view Comm_;

};

}

#endif
