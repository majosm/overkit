// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Box.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

using testing::ElementsAre;

class BoxTests : public tests::mpi_test {};

TEST_F(BoxTests, MakeEmpty) {

  if (TestComm().Rank() != 0) return;

  ovk::box Box = ovk::MakeEmptyBox(2);
  EXPECT_THAT(Box.Begin(), ElementsAre(0.,0.,0.));
  EXPECT_THAT(Box.End(), ElementsAre(-1.,-1.,0.));

}

TEST_F(BoxTests, Extend) {

  if (TestComm().Rank() != 0) return;

  // Empty
  {
    ovk::box Box = ovk::MakeEmptyBox(2);
    ovk::box ExtendBox = ovk::ExtendBox(Box, {2.,3.,0.});
    EXPECT_THAT(ExtendBox.Begin(), ElementsAre(2.,3.,0.));
    EXPECT_THAT(ExtendBox.End(), ElementsAre(2.,3.,0.));
  }

  // Inside
  {
    ovk::box Box({1.,2.,0.}, {4.,5.,0.});
    ovk::box ExtendBox = ovk::ExtendBox(Box, {2.,3.,0.});
    EXPECT_THAT(ExtendBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(ExtendBox.End(), ElementsAre(4.,5.,0.));
  }

  // Outside, below
  {
    ovk::box Box({1.,2.,0.}, {4.,5.,0.});
    ovk::box ExtendBox = ovk::ExtendBox(Box, {0.,0.,0.});
    EXPECT_THAT(ExtendBox.Begin(), ElementsAre(0.,0.,0.));
    EXPECT_THAT(ExtendBox.End(), ElementsAre(4.,5.,0.));
  }

  // Outside, above
  {
    ovk::box Box({1.,2.,0.}, {4.,5.,0.});
    ovk::box ExtendBox = ovk::ExtendBox(Box, {5.,6.,0.});
    EXPECT_THAT(ExtendBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(ExtendBox.End(), ElementsAre(5.,6.,0.));
  }

}

TEST_F(BoxTests, Overlaps) {

  if (TestComm().Rank() != 0) return;

  // Self
  {
    ovk::box Box({1.,2.,0.}, {4.,5.,0.});
    EXPECT_TRUE(ovk::BoxesOverlap(Box, Box));
  }

  // Left empty
  {
    ovk::box Box1 = ovk::MakeEmptyBox(2);
    ovk::box Box2({1.,2.,0.}, {4.,5.,0.});
    EXPECT_FALSE(ovk::BoxesOverlap(Box1, Box2));
  }

  // Right empty
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2 = ovk::MakeEmptyBox(2);
    EXPECT_FALSE(ovk::BoxesOverlap(Box1, Box2));
  }

  // Non-overlapping
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2({5.,6.,0.}, {8.,9.,0.});
    EXPECT_FALSE(ovk::BoxesOverlap(Box1, Box2));
  }

  // Right subset
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2({2.,3.,0.}, {3.,4.,0.});
    EXPECT_TRUE(ovk::BoxesOverlap(Box1, Box2));
  }

  // Left subset
  {
    ovk::box Box1({2.,3.,0.}, {3.,4.,0.});
    ovk::box Box2({1.,2.,0.}, {4.,5.,0.});
    EXPECT_TRUE(ovk::BoxesOverlap(Box1, Box2));
  }

  // Partial overlap
  {
    ovk::box Box2({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box1({3.,4.,0.}, {6.,7.,0.});
    EXPECT_TRUE(ovk::BoxesOverlap(Box1, Box2));
  }

  // Single point overlap
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2({4.,5.,0.}, {7.,8.,0.});
    EXPECT_TRUE(ovk::BoxesOverlap(Box1, Box2));
  }

}

TEST_F(BoxTests, Union) {

  if (TestComm().Rank() != 0) return;

  // Self
  {
    ovk::box Box({1.,2.,0.}, {4.,5.,0.});
    ovk::box UnionBox = ovk::UnionBoxes(Box, Box);
    EXPECT_THAT(UnionBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(UnionBox.End(), ElementsAre(4.,5.,0.));
  }

  // Left empty
  {
    ovk::box Box1 = ovk::MakeEmptyBox(2);
    ovk::box Box2({1.,2.,0.}, {4.,5.,0.});
    ovk::box UnionBox = ovk::UnionBoxes(Box1, Box2);
    EXPECT_THAT(UnionBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(UnionBox.End(), ElementsAre(4.,5.,0.));
  }

  // Right empty
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2 = ovk::MakeEmptyBox(2);
    ovk::box UnionBox = ovk::UnionBoxes(Box1, Box2);
    EXPECT_THAT(UnionBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(UnionBox.End(), ElementsAre(4.,5.,0.));
  }

  // Non-overlapping
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2({5.,6.,0.}, {8.,9.,0.});
    ovk::box UnionBox = ovk::UnionBoxes(Box1, Box2);
    EXPECT_THAT(UnionBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(UnionBox.End(), ElementsAre(8.,9.,0.));
  }

  // Right subset
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2({2.,3.,0.}, {3.,4.,0.});
    ovk::box UnionBox = ovk::UnionBoxes(Box1, Box2);
    EXPECT_THAT(UnionBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(UnionBox.End(), ElementsAre(4.,5.,0.));
  }

  // Left subset
  {
    ovk::box Box1({2.,3.,0.}, {3.,4.,0.});
    ovk::box Box2({1.,2.,0.}, {4.,5.,0.});
    ovk::box UnionBox = ovk::UnionBoxes(Box1, Box2);
    EXPECT_THAT(UnionBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(UnionBox.End(), ElementsAre(4.,5.,0.));
  }

  // Partial overlap
  {
    ovk::box Box2({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box1({3.,4.,0.}, {6.,7.,0.});
    ovk::box UnionBox = ovk::UnionBoxes(Box1, Box2);
    EXPECT_THAT(UnionBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(UnionBox.End(), ElementsAre(6.,7.,0.));
  }

}

TEST_F(BoxTests, Intersect) {

  if (TestComm().Rank() != 0) return;

  // Self
  {
    ovk::box Box({1.,2.,0.}, {4.,5.,0.});
    ovk::box IntersectBox = ovk::IntersectBoxes(Box, Box);
    EXPECT_THAT(IntersectBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(IntersectBox.End(), ElementsAre(4.,5.,0.));
  }

  // Left empty
  {
    ovk::box Box1 = ovk::MakeEmptyBox(2);
    ovk::box Box2({1.,2.,0.}, {4.,5.,0.});
    ovk::box IntersectBox = ovk::IntersectBoxes(Box1, Box2);
    EXPECT_THAT(IntersectBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(IntersectBox.End(), ElementsAre(-1.,-1.,0.));
  }

  // Right empty
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2 = ovk::MakeEmptyBox(2);
    ovk::box IntersectBox = ovk::IntersectBoxes(Box1, Box2);
    EXPECT_THAT(IntersectBox.Begin(), ElementsAre(1.,2.,0.));
    EXPECT_THAT(IntersectBox.End(), ElementsAre(-1.,-1.,0.));
  }

  // Non-overlapping
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2({5.,6.,0.}, {8.,9.,0.});
    ovk::box IntersectBox = ovk::IntersectBoxes(Box1, Box2);
    EXPECT_THAT(IntersectBox.Begin(), ElementsAre(5.,6.,0.));
    EXPECT_THAT(IntersectBox.End(), ElementsAre(4.,5.,0.));
  }

  // Right subset
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2({2.,3.,0.}, {3.,4.,0.});
    ovk::box IntersectBox = ovk::IntersectBoxes(Box1, Box2);
    EXPECT_THAT(IntersectBox.Begin(), ElementsAre(2.,3.,0.));
    EXPECT_THAT(IntersectBox.End(), ElementsAre(3.,4.,0.));
  }

  // Left subset
  {
    ovk::box Box1({2.,3.,0.}, {3.,4.,0.});
    ovk::box Box2({1.,2.,0.}, {4.,5.,0.});
    ovk::box IntersectBox = ovk::IntersectBoxes(Box1, Box2);
    EXPECT_THAT(IntersectBox.Begin(), ElementsAre(2.,3.,0.));
    EXPECT_THAT(IntersectBox.End(), ElementsAre(3.,4.,0.));
  }

  // Partial overlap
  {
    ovk::box Box2({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box1({3.,4.,0.}, {6.,7.,0.});
    ovk::box IntersectBox = ovk::IntersectBoxes(Box1, Box2);
    EXPECT_THAT(IntersectBox.Begin(), ElementsAre(3.,4.,0.));
    EXPECT_THAT(IntersectBox.End(), ElementsAre(4.,5.,0.));
  }

  // Single point overlap
  {
    ovk::box Box1({1.,2.,0.}, {4.,5.,0.});
    ovk::box Box2({4.,5.,0.}, {7.,8.,0.});
    ovk::box IntersectBox = ovk::IntersectBoxes(Box1, Box2);
    EXPECT_THAT(IntersectBox.Begin(), ElementsAre(4.,5.,0.));
    EXPECT_THAT(IntersectBox.End(), ElementsAre(4.,5.,0.));
  }

}

TEST_F(BoxTests, Clamp) {

  if (TestComm().Rank() != 0) return;

  // Inside
  {
    ovk::box Box({1.,2.,0.}, {4.,5.,0.});
    ovk::tuple<int> ClampedPoint = ovk::ClampToBox(Box, {2.,3.,0.});
    EXPECT_THAT(ClampedPoint, ElementsAre(2.,3.,0.));
  }

  // Outside, below
  {
    ovk::box Box({1.,2.,0.}, {4.,5.,0.});
    ovk::tuple<int> ClampedPoint = ovk::ClampToBox(Box, {1.,1.,0.});
    EXPECT_THAT(ClampedPoint, ElementsAre(1.,2.,0.));
  }

  // Outside, above
  {
    ovk::box Box({1.,2.,0.}, {4.,5.,0.});
    ovk::tuple<int> ClampedPoint = ovk::ClampToBox(Box, {3.,6.,0.});
    EXPECT_THAT(ClampedPoint, ElementsAre(3.,5.,0.));
  }

}

TEST_F(BoxTests, Move) {

  if (TestComm().Rank() != 0) return;

  ovk::box Box({1.,2.,0.}, {4.,5.,0.});
  ovk::box MoveBox = ovk::MoveBox(Box, {1., 2., 0.});
  EXPECT_THAT(MoveBox.Begin(), ElementsAre(2.,4.,0.));
  EXPECT_THAT(MoveBox.End(), ElementsAre(5.,7.,0.));

}

TEST_F(BoxTests, Grow) {

  if (TestComm().Rank() != 0) return;

  ovk::box Box({1.,2.,0.}, {4.,5.,0.});
  ovk::box GrowBox = ovk::GrowBox(Box, {1., 2., 0.});
  EXPECT_THAT(GrowBox.Begin(), ElementsAre(0.,0.,0.));
  EXPECT_THAT(GrowBox.End(), ElementsAre(5.,7.,0.));

}

TEST_F(BoxTests, Scale) {

  if (TestComm().Rank() != 0) return;

  ovk::box Box({1.,2.,0.}, {4.,5.,0.});
  ovk::box ScaleBox = ovk::ScaleBox(Box, {1., 2., 1.});
  EXPECT_THAT(ScaleBox.Begin(), ElementsAre(1.,0.5,0.));
  EXPECT_THAT(ScaleBox.End(), ElementsAre(4.,6.5,0.));

}
