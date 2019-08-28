// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/DistributedFieldOps.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "support/Decomp.hpp"

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Decomp.hpp>
#include <ovk/core/DistributedField.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Partition.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

class DistributedFieldOpsTests : public tests::mpi_test {};

using support::CartesianDecomp;

namespace {
// Have to also return cart comm at the moment because partition only stores comm_view
std::shared_ptr<const ovk::partition> CreatePartition(int NumDims, ovk::comm_view Comm, const
  ovk::range &GlobalRange, const ovk::tuple<int> &CartDims, bool IsPeriodic, bool Duplicated,
  ovk::comm &CartComm) {

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(Comm)
  ));

  ovk::tuple<bool> Periodic = {false,false,false};
  if (IsPeriodic) Periodic(NumDims-1) = true;

  ovk::periodic_storage PeriodicStorage = Duplicated ? ovk::periodic_storage::DUPLICATED :
    ovk::periodic_storage::UNIQUE;

  ovk::cart Cart(NumDims, GlobalRange, Periodic, PeriodicStorage);

  CartComm = ovk::CreateCartComm(Comm, NumDims, CartDims, Periodic);

  ovk::range LocalRange = CartesianDecomp(NumDims, GlobalRange, CartComm);
  ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 1);

  ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(NumDims, CartComm, LocalRange);

  ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, CartComm, LocalRange,
    DecompHash);

  return std::make_shared<ovk::partition>(std::move(Context), Cart, CartComm, LocalRange,
    ExtendedRange, 1, NeighborRanks);

}
}

TEST_F(DistributedFieldOpsTests, CountDistributedMask) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill({{1,1,0}, {5,5,1}}, true);
    long long NumTrue = ovk::core::CountDistributedMask(Mask);
    EXPECT_EQ(NumTrue, 16);
  }

}
