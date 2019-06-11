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

  ovk::range Range = ovk::MakeEmptyRange(2);
  EXPECT_THAT(Range.Begin(), ElementsAre(0,0,0));
  EXPECT_THAT(Range.End(), ElementsAre(0,0,1));

}

TEST_F(RangeTests, Extend) {

  if (TestComm().Rank() != 0) return;

  // Empty
  {
    ovk::range Range = ovk::MakeEmptyRange(2);
    ovk::range ExtendRange = ovk::ExtendRange(Range, {2,3,0});
    EXPECT_THAT(ExtendRange.Begin(), ElementsAre(2,3,0));
    EXPECT_THAT(ExtendRange.End(), ElementsAre(3,4,1));
  }

  // Inside
  {
    ovk::range Range({1,2,0}, {4,5,1});
    ovk::range ExtendRange = ovk::ExtendRange(Range, {2,3,0});
    EXPECT_THAT(ExtendRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(ExtendRange.End(), ElementsAre(4,5,1));
  }

  // Outside, below
  {
    ovk::range Range({1,2,0}, {4,5,1});
    ovk::range ExtendRange = ovk::ExtendRange(Range, {0,0,0});
    EXPECT_THAT(ExtendRange.Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(ExtendRange.End(), ElementsAre(4,5,1));
  }

  // Outside, above
  {
    ovk::range Range({1,2,0}, {4,5,1});
    ovk::range ExtendRange = ovk::ExtendRange(Range, {5,6,0});
    EXPECT_THAT(ExtendRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(ExtendRange.End(), ElementsAre(6,7,1));
  }

}

TEST_F(RangeTests, Overlaps) {

  if (TestComm().Rank() != 0) return;

  // Self
  {
    ovk::range Range({1,2,0}, {4,5,1});
    EXPECT_TRUE(ovk::RangesOverlap(Range, Range));
  }

  // Left empty
  {
    ovk::range Range1 = ovk::MakeEmptyRange(2);
    ovk::range Range2({1,2,0}, {4,5,1});
    EXPECT_FALSE(ovk::RangesOverlap(Range1, Range2));
  }

  // Right empty
  {
    ovk::range Range1({1,2,0}, {4,5,1});
    ovk::range Range2 = ovk::MakeEmptyRange(2);
    EXPECT_FALSE(ovk::RangesOverlap(Range1, Range2));
  }

  // Non-overlapping
  {
    ovk::range Range1({1,2,0}, {4,5,1});
    ovk::range Range2({4,5,0}, {7,8,1});
    EXPECT_FALSE(ovk::RangesOverlap(Range1, Range2));
  }

  // Right subset
  {
    ovk::range Range1({1,2,0}, {4,5,1});
    ovk::range Range2({2,3,0}, {3,4,1});
    EXPECT_TRUE(ovk::RangesOverlap(Range1, Range2));
  }

  // Left subset
  {
    ovk::range Range1({2,3,0}, {3,4,1});
    ovk::range Range2({1,2,0}, {4,5,1});
    EXPECT_TRUE(ovk::RangesOverlap(Range1, Range2));
  }

  // Partial overlap
  {
    ovk::range Range2({1,2,0}, {4,5,1});
    ovk::range Range1({3,4,0}, {6,7,1});
    EXPECT_TRUE(ovk::RangesOverlap(Range1, Range2));
  }

}

TEST_F(RangeTests, Union) {

  if (TestComm().Rank() != 0) return;

  // Self
  {
    ovk::range Range({1,2,0}, {4,5,1});
    ovk::range UnionRange = ovk::UnionRanges(Range, Range);
    EXPECT_THAT(UnionRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(UnionRange.End(), ElementsAre(4,5,1));
  }

  // Left empty
  {
    ovk::range Range1 = ovk::MakeEmptyRange(2);
    ovk::range Range2({1,2,0}, {4,5,1});
    ovk::range UnionRange = ovk::UnionRanges(Range1, Range2);
    EXPECT_THAT(UnionRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(UnionRange.End(), ElementsAre(4,5,1));
  }

  // Right empty
  {
    ovk::range Range1({1,2,0}, {4,5,1});
    ovk::range Range2 = ovk::MakeEmptyRange(2);
    ovk::range UnionRange = ovk::UnionRanges(Range1, Range2);
    EXPECT_THAT(UnionRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(UnionRange.End(), ElementsAre(4,5,1));
  }

  // Non-overlapping
  {
    ovk::range Range1({1,2,0}, {4,5,1});
    ovk::range Range2({4,5,0}, {7,8,1});
    ovk::range UnionRange = ovk::UnionRanges(Range1, Range2);
    EXPECT_THAT(UnionRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(UnionRange.End(), ElementsAre(7,8,1));
  }

  // Right subset
  {
    ovk::range Range1({1,2,0}, {4,5,1});
    ovk::range Range2({2,3,0}, {3,4,1});
    ovk::range UnionRange = ovk::UnionRanges(Range1, Range2);
    EXPECT_THAT(UnionRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(UnionRange.End(), ElementsAre(4,5,1));
  }

  // Left subset
  {
    ovk::range Range1({2,3,0}, {3,4,1});
    ovk::range Range2({1,2,0}, {4,5,1});
    ovk::range UnionRange = ovk::UnionRanges(Range1, Range2);
    EXPECT_THAT(UnionRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(UnionRange.End(), ElementsAre(4,5,1));
  }

  // Partial overlap
  {
    ovk::range Range2({1,2,0}, {4,5,1});
    ovk::range Range1({3,4,0}, {6,7,1});
    ovk::range UnionRange = ovk::UnionRanges(Range1, Range2);
    EXPECT_THAT(UnionRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(UnionRange.End(), ElementsAre(6,7,1));
  }

}

TEST_F(RangeTests, Intersect) {

  if (TestComm().Rank() != 0) return;

  // Self
  {
    ovk::range Range({1,2,0}, {4,5,1});
    ovk::range IntersectRange = ovk::IntersectRanges(Range, Range);
    EXPECT_THAT(IntersectRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(IntersectRange.End(), ElementsAre(4,5,1));
  }

  // Left empty
  {
    ovk::range Range1 = ovk::MakeEmptyRange(2);
    ovk::range Range2({1,2,0}, {4,5,1});
    ovk::range IntersectRange = ovk::IntersectRanges(Range1, Range2);
    EXPECT_THAT(IntersectRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(IntersectRange.End(), ElementsAre(0,0,1));
  }

  // Right empty
  {
    ovk::range Range1({1,2,0}, {4,5,1});
    ovk::range Range2 = ovk::MakeEmptyRange(2);
    ovk::range IntersectRange = ovk::IntersectRanges(Range1, Range2);
    EXPECT_THAT(IntersectRange.Begin(), ElementsAre(1,2,0));
    EXPECT_THAT(IntersectRange.End(), ElementsAre(0,0,1));
  }

  // Non-overlapping
  {
    ovk::range Range1({1,2,0}, {4,5,1});
    ovk::range Range2({4,5,0}, {7,8,1});
    ovk::range IntersectRange = ovk::IntersectRanges(Range1, Range2);
    EXPECT_THAT(IntersectRange.Begin(), ElementsAre(4,5,0));
    EXPECT_THAT(IntersectRange.End(), ElementsAre(4,5,1));
  }

  // Right subset
  {
    ovk::range Range1({1,2,0}, {4,5,1});
    ovk::range Range2({2,3,0}, {3,4,1});
    ovk::range IntersectRange = ovk::IntersectRanges(Range1, Range2);
    EXPECT_THAT(IntersectRange.Begin(), ElementsAre(2,3,0));
    EXPECT_THAT(IntersectRange.End(), ElementsAre(3,4,1));
  }

  // Left subset
  {
    ovk::range Range1({2,3,0}, {3,4,1});
    ovk::range Range2({1,2,0}, {4,5,1});
    ovk::range IntersectRange = ovk::IntersectRanges(Range1, Range2);
    EXPECT_THAT(IntersectRange.Begin(), ElementsAre(2,3,0));
    EXPECT_THAT(IntersectRange.End(), ElementsAre(3,4,1));
  }

  // Partial overlap
  {
    ovk::range Range2({1,2,0}, {4,5,1});
    ovk::range Range1({3,4,0}, {6,7,1});
    ovk::range IntersectRange = ovk::IntersectRanges(Range1, Range2);
    EXPECT_THAT(IntersectRange.Begin(), ElementsAre(3,4,0));
    EXPECT_THAT(IntersectRange.End(), ElementsAre(4,5,1));
  }

}

TEST_F(RangeTests, Clamp) {

  if (TestComm().Rank() != 0) return;

  // Inside
  {
    ovk::range Range({1,2,0}, {4,5,1});
    ovk::tuple<int> ClampedPoint = ovk::ClampToRange(Range, {2,3,0});
    EXPECT_THAT(ClampedPoint, ElementsAre(2,3,0));
  }

  // Outside, below
  {
    ovk::range Range({1,2,0}, {4,5,1});
    ovk::tuple<int> ClampedPoint = ovk::ClampToRange(Range, {1,1,0});
    EXPECT_THAT(ClampedPoint, ElementsAre(1,2,0));
  }

  // Outside, above
  {
    ovk::range Range({1,2,0}, {4,5,1});
    ovk::tuple<int> ClampedPoint = ovk::ClampToRange(Range, {3,5,0});
    EXPECT_THAT(ClampedPoint, ElementsAre(3,4,0));
  }

}
