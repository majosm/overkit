// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "MPITest.hpp"

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

TEST_F(RangeTests, TupleToIndex2DRowMajor) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(2, ovk::elem<int,2>({2,2}), ovk::elem<int,2>({4,5}));
  ovk::elem<int,2> Tuple;
  long long Index;

  // Lower i, lower j corner
  Tuple = {2,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 0);

  // Upper i, lower j corner
  Tuple = {3,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 3);

  // Lower i, upper j corner
  Tuple = {2,4};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 2);

  // Upper i, upper j corner
  Tuple = {3,4};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 5);

}

TEST_F(RangeTests, TupleToIndex2DColumnMajor) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(2, ovk::elem<int,2>({2,2}), ovk::elem<int,2>({4,5}));
  ovk::elem<int,2> Tuple;
  long long Index;

  // Lower i, lower j corner
  Tuple = {2,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 0);

  // Upper i, lower j corner
  Tuple = {3,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 1);

  // Lower i, upper j corner
  Tuple = {2,4};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 4);

  // Upper i, upper j corner
  Tuple = {3,4};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 5);

}

TEST_F(RangeTests, TupleToIndex3DRowMajor) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(3, ovk::elem<int,3>({2,2,2}), ovk::elem<int,3>({4,5,6}));
  ovk::elem<int,3> Tuple;
  long long Index;

  // Lower i, lower j, lower k corner
  Tuple = {2,2,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 0);

  // Upper i, lower j, lower k corner
  Tuple = {3,2,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 12);

  // Lower i, upper j, lower k corner
  Tuple = {2,4,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 8);

  // Upper i, upper j, lower k corner
  Tuple = {3,4,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 20);

  // Lower i, lower j, upper k corner
  Tuple = {2,2,5};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 3);

  // Upper i, lower j, upper k corner
  Tuple = {3,2,5};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 15);

  // Lower i, upper j, upper k corner
  Tuple = {2,4,5};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 11);

  // Upper i, upper j, upper k corner
  Tuple = {3,4,5};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::ROW_MAJOR, Tuple);
  EXPECT_EQ(Index, 23);

}

TEST_F(RangeTests, TupleToIndex3DColumnMajor) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(3, ovk::elem<int,3>({2,2,2}), ovk::elem<int,3>({4,5,6}));
  ovk::elem<int,3> Tuple;
  long long Index;

  // Lower i, lower j, lower k corner
  Tuple = {2,2,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 0);

  // Upper i, lower j, lower k corner
  Tuple = {3,2,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 1);

  // Lower i, upper j, lower k corner
  Tuple = {2,4,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 4);

  // Upper i, upper j, lower k corner
  Tuple = {3,4,2};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 5);

  // Lower i, lower j, upper k corner
  Tuple = {2,2,5};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 18);

  // Upper i, lower j, upper k corner
  Tuple = {3,2,5};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 19);

  // Lower i, upper j, upper k corner
  Tuple = {2,4,5};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 22);

  // Upper i, upper j, upper k corner
  Tuple = {3,4,5};
  Index = ovk::RangeTupleToIndex(Range, ovk::array_layout::COLUMN_MAJOR, Tuple);
  EXPECT_EQ(Index, 23);

}

TEST_F(RangeTests, IndexToTuple2DRowMajor) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(2, ovk::elem<int,2>({2,2}), ovk::elem<int,2>({4,5}));
  ovk::elem<int,OVK_MAX_DIMS> Tuple;
  ovk::elem<int,2> Tuple2D;

  // Lower i, lower j corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 0);
  EXPECT_THAT(Tuple, ElementsAre(2,2,0));
  Tuple2D = ovk::RangeIndexToTuple<2>(Range, ovk::array_layout::ROW_MAJOR, 0);
  EXPECT_THAT(Tuple2D, ElementsAre(2,2));

  // Upper i, lower j corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 3);
  EXPECT_THAT(Tuple, ElementsAre(3,2,0));
  Tuple2D = ovk::RangeIndexToTuple<2>(Range, ovk::array_layout::ROW_MAJOR, 3);
  EXPECT_THAT(Tuple2D, ElementsAre(3,2));

  // Lower i, upper j corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 2);
  EXPECT_THAT(Tuple, ElementsAre(2,4,0));
  Tuple2D = ovk::RangeIndexToTuple<2>(Range, ovk::array_layout::ROW_MAJOR, 2);
  EXPECT_THAT(Tuple2D, ElementsAre(2,4));

  // Upper i, upper j corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 5);
  EXPECT_THAT(Tuple, ElementsAre(3,4,0));
  Tuple2D = ovk::RangeIndexToTuple<2>(Range, ovk::array_layout::ROW_MAJOR, 5);
  EXPECT_THAT(Tuple2D, ElementsAre(3,4));

}

TEST_F(RangeTests, IndexToTuple2DColumnMajor) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(2, ovk::elem<int,2>({2,2}), ovk::elem<int,2>({4,5}));
  ovk::elem<int,OVK_MAX_DIMS> Tuple;
  ovk::elem<int,2> Tuple2D;

  // Lower i, lower j corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 0);
  EXPECT_THAT(Tuple, ElementsAre(2,2,0));
  Tuple2D = ovk::RangeIndexToTuple<2>(Range, ovk::array_layout::COLUMN_MAJOR, 0);
  EXPECT_THAT(Tuple2D, ElementsAre(2,2));

  // Upper i, lower j corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 1);
  EXPECT_THAT(Tuple, ElementsAre(3,2,0));
  Tuple2D = ovk::RangeIndexToTuple<2>(Range, ovk::array_layout::COLUMN_MAJOR, 1);
  EXPECT_THAT(Tuple2D, ElementsAre(3,2));

  // Lower i, upper j corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 4);
  EXPECT_THAT(Tuple, ElementsAre(2,4,0));
  Tuple2D = ovk::RangeIndexToTuple<2>(Range, ovk::array_layout::COLUMN_MAJOR, 4);
  EXPECT_THAT(Tuple2D, ElementsAre(2,4));

  // Upper i, upper j corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 5);
  EXPECT_THAT(Tuple, ElementsAre(3,4,0));
  Tuple2D = ovk::RangeIndexToTuple<2>(Range, ovk::array_layout::COLUMN_MAJOR, 5);
  EXPECT_THAT(Tuple2D, ElementsAre(3,4));

}

TEST_F(RangeTests, IndexToTuple3DRowMajor) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(3, ovk::elem<int,3>({2,2,2}), ovk::elem<int,3>({4,5,6}));
  ovk::elem<int,OVK_MAX_DIMS> Tuple;

  // Lower i, lower j, lower k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 0);
  EXPECT_THAT(Tuple, ElementsAre(2,2,2));

  // Upper i, lower j, lower k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 12);
  EXPECT_THAT(Tuple, ElementsAre(3,2,2));

  // Lower i, upper j, lower k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 8);
  EXPECT_THAT(Tuple, ElementsAre(2,4,2));

  // Upper i, upper j, lower k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 20);
  EXPECT_THAT(Tuple, ElementsAre(3,4,2));

  // Lower i, lower j, upper k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 3);
  EXPECT_THAT(Tuple, ElementsAre(2,2,5));

  // Upper i, lower j, upper k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 15);
  EXPECT_THAT(Tuple, ElementsAre(3,2,5));

  // Lower i, upper j, upper k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 11);
  EXPECT_THAT(Tuple, ElementsAre(2,4,5));

  // Upper i, upper j, upper k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::ROW_MAJOR, 23);
  EXPECT_THAT(Tuple, ElementsAre(3,4,5));

}

TEST_F(RangeTests, IndexToTuple3DColumnMajor) {

  if (TestComm().Rank() != 0) return;

  ovk::range Range(3, ovk::elem<int,3>({2,2,2}), ovk::elem<int,3>({4,5,6}));
  ovk::elem<int,OVK_MAX_DIMS> Tuple;

  // Lower i, lower j, lower k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 0);
  EXPECT_THAT(Tuple, ElementsAre(2,2,2));

  // Upper i, lower j, lower k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 1);
  EXPECT_THAT(Tuple, ElementsAre(3,2,2));

  // Lower i, upper j, lower k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 4);
  EXPECT_THAT(Tuple, ElementsAre(2,4,2));

  // Upper i, upper j, lower k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 5);
  EXPECT_THAT(Tuple, ElementsAre(3,4,2));

  // Lower i, lower j, upper k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 18);
  EXPECT_THAT(Tuple, ElementsAre(2,2,5));

  // Upper i, lower j, upper k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 19);
  EXPECT_THAT(Tuple, ElementsAre(3,2,5));

  // Lower i, upper j, upper k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 22);
  EXPECT_THAT(Tuple, ElementsAre(2,4,5));

  // Upper i, upper j, upper k corner
  Tuple = ovk::RangeIndexToTuple(Range, ovk::array_layout::COLUMN_MAJOR, 23);
  EXPECT_THAT(Tuple, ElementsAre(3,4,5));

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
