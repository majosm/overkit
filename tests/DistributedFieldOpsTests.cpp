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
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

using testing::ElementsAre;
using testing::ElementsAreArray;

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

TEST_F(DistributedFieldOpsTests, DetectEdge) {

  ASSERT_GE(TestComm().Size(), 8);

  ovk::comm CommOfSize2 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 2);
  ovk::comm CommOfSize4 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 4);
  ovk::comm CommOfSize8 = ovk::CreateSubsetComm(TestComm(), TestComm().Rank() < 8);

  auto ToRange1D = [](int LowerCorner, int UpperCorner) -> ovk::range {
    return {{LowerCorner,0,0}, {UpperCorner+1,1,1}};
  };

  auto ToRange2D = [](const ovk::elem<int,2> &LowerCorner, const ovk::elem<int,2> &UpperCorner) ->
    ovk::range {
    return {{LowerCorner(0),LowerCorner(1),0}, {UpperCorner(0)+1,UpperCorner(1)+1,1}};
  };

  auto ToRange3D = [](const ovk::elem<int,3> &LowerCorner, const ovk::elem<int,3> &UpperCorner) ->
    ovk::range {
    return {{LowerCorner(0),LowerCorner(1),LowerCorner(2)}, {UpperCorner(0)+1,UpperCorner(1)+1,
      UpperCorner(2)+1}};
  };

  auto CreatePartitionIncludingExteriorPoint = [](const std::shared_ptr<const ovk::partition>
    &Partition) -> std::shared_ptr<const ovk::partition> {
    ovk::cart Cart = ovk::core::CartIncludeExteriorPoint(Partition->Cart());
    ovk::range LocalRange = ovk::core::RangeIncludeExteriorPoint(Partition->Cart(),
      Partition->LocalRange());
    ovk::range ExtendedRange = ovk::core::RangeIncludeExteriorPoint(Partition->Cart(),
      Partition->ExtendedRange());
    return std::make_shared<ovk::partition>(Partition->SharedContext(), Cart, Partition->Comm(),
      LocalRange, ExtendedRange, Partition->SubregionCount(), Partition->NeighborRanks());
  };

  // 1D, interior, inner edge
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(1, 4), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange1D(1, 4), true);
    ExpectedValues.Fill(ToRange1D(2, 3), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, interior, outer edge
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(1, 4), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, true);
    ExpectedValues.Fill(ToRange1D(1, 4), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, periodic boundary, inner edge, unique
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, true, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(4, 5), true);
    Mask.Fill(ToRange1D(0, 1), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange1D(4, 5), true);
    ExpectedValues.Fill(ToRange1D(5, 5), false);
    ExpectedValues.Fill(ToRange1D(0, 1), true);
    ExpectedValues.Fill(ToRange1D(0, 0), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, periodic boundary, inner edge, duplicated
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, true, true, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(4, 5), true);
    Mask.Fill(ToRange1D(0, 1), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange1D(4, 5), true);
    ExpectedValues.Fill(ToRange1D(5, 5), false);
    ExpectedValues.Fill(ToRange1D(0, 1), true);
    ExpectedValues.Fill(ToRange1D(0, 0), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, periodic boundary, outer edge, unique
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, true, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(4, 5), true);
    Mask.Fill(ToRange1D(0, 1), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, true);
    ExpectedValues.Fill(ToRange1D(4, 5), false);
    ExpectedValues.Fill(ToRange1D(0, 1), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, periodic boundary, outer edge, duplicated
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, true, true, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(4, 5), true);
    Mask.Fill(ToRange1D(0, 1), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, true);
    ExpectedValues.Fill(ToRange1D(4, 5), false);
    ExpectedValues.Fill(ToRange1D(0, 1), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, false boundary, inner edge
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(0, 3), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange1D(0, 3), true);
    ExpectedValues.Fill(ToRange1D(1, 2), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, false boundary, outer edge
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(0, 3), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange1D(-1, 4), true);
    ExpectedValues.Fill(ToRange1D(0, 3), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, true boundary, inner edge
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(0, 3), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::TRUE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange1D(0, 3), true);
    ExpectedValues.Fill(ToRange1D(0, 2), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, true boundary, outer edge
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(0, 3), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::TRUE, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange1D(0, 5), true);
    ExpectedValues.Fill(ToRange1D(0, 3), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, mirror boundary, inner edge
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(0, 3), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::MIRROR, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange1D(-1, 3), true);
    ExpectedValues.Fill(ToRange1D(-1, 2), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 1D, mirror boundary, outer edge
  if (CommOfSize2) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(1, CommOfSize2, {{6,1,1}}, {2,1,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange1D(0, 3), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::MIRROR, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange1D(-1, 4), true);
    ExpectedValues.Fill(ToRange1D(-1, 3), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, interior, inner edge
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,1}, {4,4}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange2D({1,1}, {4,4}), true);
    ExpectedValues.Fill(ToRange2D({2,2}, {3,3}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, interior, outer edge
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,1}, {4,4}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, true);
    ExpectedValues.Fill(ToRange2D({1,1}, {4,4}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, periodic boundary, inner edge, unique
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, true, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,4}, {4,5}), true);
    Mask.Fill(ToRange2D({1,0}, {4,1}), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange2D({1,4}, {4,5}), true);
    ExpectedValues.Fill(ToRange2D({2,5}, {3,5}), false);
    ExpectedValues.Fill(ToRange2D({1,0}, {4,1}), true);
    ExpectedValues.Fill(ToRange2D({2,0}, {3,0}), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, periodic boundary, inner edge, duplicated
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, true, true, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,4}, {4,5}), true);
    Mask.Fill(ToRange2D({1,0}, {4,1}), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange2D({1,4}, {4,5}), true);
    ExpectedValues.Fill(ToRange2D({2,5}, {3,5}), false);
    ExpectedValues.Fill(ToRange2D({1,0}, {4,1}), true);
    ExpectedValues.Fill(ToRange2D({2,0}, {3,0}), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, periodic boundary, outer edge, unique
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, true, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,4}, {4,5}), true);
    Mask.Fill(ToRange2D({1,0}, {4,1}), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, true);
    ExpectedValues.Fill(ToRange2D({1,4}, {4,5}), false);
    ExpectedValues.Fill(ToRange2D({1,0}, {4,1}), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, periodic boundary, outer edge, duplicated
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, true, true, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,4}, {4,5}), true);
    Mask.Fill(ToRange2D({1,0}, {4,1}), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, true);
    ExpectedValues.Fill(ToRange2D({1,4}, {4,5}), false);
    ExpectedValues.Fill(ToRange2D({1,0}, {4,1}), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, false boundary, inner edge
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,0}, {4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange2D({1,0}, {4,3}), true);
    ExpectedValues.Fill(ToRange2D({2,1}, {3,2}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, false boundary, outer edge
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,0}, {4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange2D({0,-1}, {5,4}), true);
    ExpectedValues.Fill(ToRange2D({1,0}, {4,3}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, true boundary, inner edge
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,0}, {4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::TRUE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange2D({1,0}, {4,3}), true);
    ExpectedValues.Fill(ToRange2D({2,0}, {3,2}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, true boundary, outer edge
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,0}, {4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::TRUE, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange2D({0,0}, {5,5}), true);
    ExpectedValues.Fill(ToRange2D({1,0}, {4,3}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, mirror boundary, inner edge
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,0}, {4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::MIRROR, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange2D({1,-1}, {4,3}), true);
    ExpectedValues.Fill(ToRange2D({2,-1}, {3,2}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 2D, mirror boundary, outer edge
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,0}, {4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::MIRROR, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange2D({0,-1}, {5,4}), true);
    ExpectedValues.Fill(ToRange2D({1,-1}, {4,3}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, interior, inner edge
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,1}, {4,4,4}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange3D({1,1,1}, {4,4,4}), true);
    ExpectedValues.Fill(ToRange3D({2,2,2}, {3,3,3}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, interior, outer edge
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,1}, {4,4,4}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, true);
    ExpectedValues.Fill(ToRange3D({1,1,1}, {4,4,4}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, periodic boundary, inner edge, unique
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, true, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,4}, {4,4,5}), true);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,1}), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange3D({1,1,4}, {4,4,5}), true);
    ExpectedValues.Fill(ToRange3D({2,2,5}, {3,3,5}), false);
    ExpectedValues.Fill(ToRange3D({1,1,0}, {4,4,1}), true);
    ExpectedValues.Fill(ToRange3D({2,2,0}, {3,3,0}), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, periodic boundary, inner edge, duplicated
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, true, true, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,4}, {4,4,5}), true);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,1}), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange3D({1,1,4}, {4,4,5}), true);
    ExpectedValues.Fill(ToRange3D({2,2,5}, {3,3,5}), false);
    ExpectedValues.Fill(ToRange3D({1,1,0}, {4,4,1}), true);
    ExpectedValues.Fill(ToRange3D({2,2,0}, {3,3,0}), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, periodic boundary, outer edge, unique
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, true, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,4}, {4,4,5}), true);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,1}), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, true);
    ExpectedValues.Fill(ToRange3D({1,1,4}, {4,4,5}), false);
    ExpectedValues.Fill(ToRange3D({1,1,0}, {4,4,1}), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, periodic boundary, outer edge, duplicated
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, true, true, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,4}, {4,4,5}), true);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,1}), true);
    Mask.Exchange();
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, true);
    ExpectedValues.Fill(ToRange3D({1,1,4}, {4,4,5}), false);
    ExpectedValues.Fill(ToRange3D({1,1,0}, {4,4,1}), false);
    ExpectedValues.Exchange();
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, false boundary, inner edge
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::FALSE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange3D({1,1,0}, {4,4,3}), true);
    ExpectedValues.Fill(ToRange3D({2,2,1}, {3,3,2}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, false boundary, outer edge
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange3D({0,0,-1}, {5,5,4}), true);
    ExpectedValues.Fill(ToRange3D({1,1,0}, {4,4,3}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, true boundary, inner edge
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::TRUE, false,
      EdgeMask);
    ovk::distributed_field<bool> ExpectedValues(Partition, false);
    ExpectedValues.Fill(ToRange3D({1,1,0}, {4,4,3}), true);
    ExpectedValues.Fill(ToRange3D({2,2,0}, {3,3,2}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, true boundary, outer edge
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::TRUE, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange3D({0,0,0}, {5,5,5}), true);
    ExpectedValues.Fill(ToRange3D({1,1,0}, {4,4,3}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, mirror boundary, inner edge
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::INNER, ovk::core::mask_bc::MIRROR, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange3D({1,1,-1}, {4,4,3}), true);
    ExpectedValues.Fill(ToRange3D({2,2,-1}, {3,3,2}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // 3D, mirror boundary, outer edge
  if (CommOfSize8) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(3, CommOfSize8, {{6,6,6}}, {2,2,2}, false, false, CartComm);
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange3D({1,1,0}, {4,4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::MIRROR, true,
      EdgeMask);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange3D({0,0,-1}, {5,5,4}), true);
    ExpectedValues.Fill(ToRange3D({1,1,-1}, {4,4,3}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

  // With partition pool
  if (CommOfSize4) {
    ovk::comm CartComm;
    auto Partition = CreatePartition(2, CommOfSize4, {{6,6,1}}, {2,2,1}, false, false, CartComm);
    ovk::core::partition_pool PartitionPool(Partition->SharedContext(), CartComm,
      Partition->NeighborRanks());
    ovk::distributed_field<bool> Mask(Partition, false);
    Mask.Fill(ToRange2D({1,0}, {4,3}), true);
    ovk::distributed_field<bool> EdgeMask;
    ovk::core::DetectEdge(Mask, ovk::core::edge_type::OUTER, ovk::core::mask_bc::FALSE, true,
      EdgeMask, &PartitionPool);
    auto EdgePartition = CreatePartitionIncludingExteriorPoint(Partition);
    ovk::distributed_field<bool> ExpectedValues(EdgePartition, false);
    ExpectedValues.Fill(ToRange2D({0,-1}, {5,4}), true);
    ExpectedValues.Fill(ToRange2D({1,0}, {4,3}), false);
    EXPECT_THAT(EdgeMask, ElementsAreArray(ExpectedValues));
  }

}
