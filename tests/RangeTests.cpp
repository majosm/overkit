// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "MPITest.hpp"

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

using testing::ElementsAre;

class RangeTests : public tests::mpi_test {};

TEST_F(RangeTests, ConstructEmpty2D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(2);

  EXPECT_EQ(Range.Dimension(), 2);
  EXPECT_THAT(Range.Begin(), ElementsAre(0,0,0));
  EXPECT_THAT(Range.End(), ElementsAre(0,0,1));

}

TEST_F(RangeTests, ConstructEmpty3D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(3);

  EXPECT_EQ(Range.Dimension(), 3);
  EXPECT_THAT(Range.Begin(), ElementsAre(0,0,0));
  EXPECT_THAT(Range.End(), ElementsAre(0,0,0));

}

TEST_F(RangeTests, ConstructAssigned2D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{3,4});

  EXPECT_EQ(Range.Dimension(), 2);
  EXPECT_THAT(Range.Begin(), ElementsAre(1,2,0));
  EXPECT_EQ(Range.Begin(0), 1);
  EXPECT_EQ(Range.Begin(1), 2);
  EXPECT_EQ(Range.Begin(2), 0);
  const int *BeginData = Range.BeginData();
  EXPECT_EQ(BeginData[0], 1);
  EXPECT_EQ(BeginData[1], 2);
  EXPECT_EQ(BeginData[2], 0);
  EXPECT_THAT(Range.End(), ElementsAre(3,4,1));
  EXPECT_EQ(Range.End(0), 3);
  EXPECT_EQ(Range.End(1), 4);
  EXPECT_EQ(Range.End(2), 1);
  const int *EndData = Range.EndData();
  EXPECT_EQ(EndData[0], 3);
  EXPECT_EQ(EndData[1], 4);
  EXPECT_EQ(EndData[2], 1);

}

TEST_F(RangeTests, ConstructAssigned3D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(3, ovk::elem<int,3>{1,2,3}, ovk::elem<int,3>{4,5,6});

  EXPECT_EQ(Range.Dimension(), 3);
  EXPECT_THAT(Range.Begin(), ElementsAre(1,2,3));
  EXPECT_EQ(Range.Begin(0), 1);
  EXPECT_EQ(Range.Begin(1), 2);
  EXPECT_EQ(Range.Begin(2), 3);
  const int *BeginData = Range.BeginData();
  EXPECT_EQ(BeginData[0], 1);
  EXPECT_EQ(BeginData[1], 2);
  EXPECT_EQ(BeginData[2], 3);
  EXPECT_THAT(Range.End(), ElementsAre(4,5,6));
  EXPECT_EQ(Range.End(0), 4);
  EXPECT_EQ(Range.End(1), 5);
  EXPECT_EQ(Range.End(2), 6);
  const int *EndData = Range.EndData();
  EXPECT_EQ(EndData[0], 4);
  EXPECT_EQ(EndData[1], 5);
  EXPECT_EQ(EndData[2], 6);

}

TEST_F(RangeTests, CopyConstruct2D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range1(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{3,4});
  ovk::range Range2(Range1);

  EXPECT_EQ(Range2.Dimension(), 2);
  EXPECT_THAT(Range2.Begin(), ElementsAre(1,2,0));
  EXPECT_THAT(Range2.End(), ElementsAre(3,4,1));

}

TEST_F(RangeTests, CopyConstruct3D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range1(3, ovk::elem<int,3>{1,2,3}, ovk::elem<int,3>{4,5,6});
  ovk::range Range2(Range1);

  EXPECT_EQ(Range2.Dimension(), 3);
  EXPECT_THAT(Range2.Begin(), ElementsAre(1,2,3));
  EXPECT_THAT(Range2.End(), ElementsAre(4,5,6));

}

TEST_F(RangeTests, CopyAssign2D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range1(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{3,4});
  ovk::range Range2;

  Range2 = Range1;

  EXPECT_EQ(Range2.Dimension(), 2);
  EXPECT_THAT(Range2.Begin(), ElementsAre(1,2,0));
  EXPECT_THAT(Range2.End(), ElementsAre(3,4,1));

}

TEST_F(RangeTests, CopyAssign3D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range1(3, ovk::elem<int,3>{1,2,3}, ovk::elem<int,3>{4,5,6});
  ovk::range Range2;

  Range2 = Range1;

  EXPECT_EQ(Range2.Dimension(), 3);
  EXPECT_THAT(Range2.Begin(), ElementsAre(1,2,3));
  EXPECT_THAT(Range2.End(), ElementsAre(4,5,6));

}

TEST_F(RangeTests, Equality) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range1, Range2;

  Range1 = ovk::range(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{3,4});

  // Self
  EXPECT_TRUE(Range1 == Range1);

  // Other with same values
  Range2 = Range1;
  EXPECT_TRUE(Range1 == Range2);

  // Different dimension
  Range2 = ovk::range(3, ovk::elem<int,3>{1,2,0}, ovk::elem<int,3>{3,4,1});
  EXPECT_TRUE(Range1 != Range2);

  // Different Begin
  Range2 = ovk::range(2, ovk::elem<int,2>{1,3}, ovk::elem<int,2>{3,4});
  EXPECT_TRUE(Range1 != Range2);

  // Different End
  Range2 = ovk::range(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{2,4});
  EXPECT_TRUE(Range1 != Range2);

}

TEST_F(RangeTests, Size2D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range;

  // Empty
  Range = ovk::range(2);
  EXPECT_THAT(Range.Size(), ElementsAre(0,0,1));
  EXPECT_EQ(Range.Size(0), 0);
  EXPECT_EQ(Range.Size(1), 0);
  EXPECT_EQ(Range.Size(2), 1);

  // Non-empty
  Range = ovk::range(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{4,3});
  EXPECT_THAT(Range.Size(), ElementsAre(3,1,1));
  EXPECT_EQ(Range.Size(0), 3);
  EXPECT_EQ(Range.Size(1), 1);
  EXPECT_EQ(Range.Size(2), 1);

}

TEST_F(RangeTests, Size3D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range;

  // Empty
  Range = ovk::range(3);
  EXPECT_THAT(Range.Size(), ElementsAre(0,0,0));
  EXPECT_EQ(Range.Size(0), 0);
  EXPECT_EQ(Range.Size(1), 0);
  EXPECT_EQ(Range.Size(2), 0);

  // Non-empty
  Range = ovk::range(3, ovk::elem<int,3>{1,2,3}, ovk::elem<int,3>{6,5,4});
  EXPECT_THAT(Range.Size(), ElementsAre(5,3,1));
  EXPECT_EQ(Range.Size(0), 5);
  EXPECT_EQ(Range.Size(1), 3);
  EXPECT_EQ(Range.Size(2), 1);

}

TEST_F(RangeTests, Count2D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range;

  // Empty
  Range = ovk::range(2);
  EXPECT_EQ(Range.Count(), 0);

  // Non-empty
  Range = ovk::range(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{4,3});
  EXPECT_EQ(Range.Count(), 3);

}

TEST_F(RangeTests, Count3D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range;

  // Empty
  Range = ovk::range(3);
  EXPECT_EQ(Range.Count(), 0);

  // Non-empty
  Range = ovk::range(3, ovk::elem<int,3>{1,2,3}, ovk::elem<int,3>{6,5,4});
  EXPECT_EQ(Range.Count(), 15);

}

TEST_F(RangeTests, IsEmpty2D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range;

  // Empty
  Range = ovk::range(2);
  EXPECT_TRUE(Range.Empty());

  // Non-empty
  Range = ovk::range(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{4,3});
  EXPECT_FALSE(Range.Empty());

}

TEST_F(RangeTests, IsEmpty3D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range;

  // Empty
  Range = ovk::range(3);
  EXPECT_TRUE(Range.Empty());

  // Non-empty
  Range = ovk::range(3, ovk::elem<int,3>{1,2,3}, ovk::elem<int,3>{6,5,4});
  EXPECT_FALSE(Range.Empty());

}

TEST_F(RangeTests, Contains2D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range;
  ovk::elem<int, 2> Point;

  // Inside
  Range = ovk::range(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{4,5});
  Point = {2,3};
  EXPECT_TRUE(ovk::RangeContains(Range, Point));

  // Outside
  Range = ovk::range(2, ovk::elem<int,2>{1,2}, ovk::elem<int,2>{4,5});
  Point = {0,0};
  EXPECT_FALSE(ovk::RangeContains(Range, Point));

  // Empty
  Range = ovk::range(2);
  Point = {1,1};
  EXPECT_FALSE(ovk::RangeContains(Range, Point));

}

TEST_F(RangeTests, Contains3D) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range;
  ovk::elem<int, 3> Point;

  // Inside
  Range = ovk::range(3, ovk::elem<int,3>{1,2,3}, ovk::elem<int,3>{4,5,6});
  Point = {2,3,4};
  EXPECT_TRUE(ovk::RangeContains(Range, Point));

  // Outside
  Range = ovk::range(3, ovk::elem<int,3>{1,2,3}, ovk::elem<int,3>{4,5,6});
  Point = {0,0,0};
  EXPECT_FALSE(ovk::RangeContains(Range, Point));

  // Empty
  Range = ovk::range(3);
  Point = {1,1,1};
  EXPECT_FALSE(ovk::RangeContains(Range, Point));

}

// TODO: Add remaining tests
