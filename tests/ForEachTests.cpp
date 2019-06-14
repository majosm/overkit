// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/ForEach.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/Vector.hpp>

#include <mpi.h>

using testing::ElementsAre;

class ForEachTests : public tests::mpi_test {};

TEST_F(ForEachTests, ForEachTupleInInterval) {

  if (TestComm().Rank() != 0) return;

  // One-dimensional, tuple
  {
    ovk::interval<int> Interval(1, 5);
    ovk::core::vector<ovk::elem<int>> Tuples;
    ovk::core::ForEach(Interval, [&](const ovk::elem<int> &Tuple) {
      Tuples.Append(Tuple);
    });
    EXPECT_EQ(Tuples.Count(), 4);
    EXPECT_THAT(Tuples[0], ElementsAre(1));
    EXPECT_THAT(Tuples[1], ElementsAre(2));
    EXPECT_THAT(Tuples[2], ElementsAre(3));
    EXPECT_THAT(Tuples[3], ElementsAre(4));
  }

  // One-dimensional, tuple element
  {
    ovk::interval<int> Interval(1, 5);
    ovk::core::vector<ovk::elem<int>> Tuples;
    ovk::core::ForEach(Interval, [&](int i) {
      Tuples.Append(ovk::elem<int>(i));
    });
    EXPECT_EQ(Tuples.Count(), 4);
    EXPECT_THAT(Tuples[0], ElementsAre(1));
    EXPECT_THAT(Tuples[1], ElementsAre(2));
    EXPECT_THAT(Tuples[2], ElementsAre(3));
    EXPECT_THAT(Tuples[3], ElementsAre(4));
  }

  // Multidimensional, row major, tuple
  {
    ovk::interval<int,3> Interval({1,2,3}, {2,4,6});
    ovk::core::vector<ovk::elem<int,3>> Tuples;
    ovk::core::ForEach<ovk::array_layout::ROW_MAJOR>(Interval, [&](const ovk::elem<int,3> &Tuple) {
      Tuples.Append(Tuple);
    });
    EXPECT_EQ(Tuples.Count(), 6);
    EXPECT_THAT(Tuples[0], ElementsAre(1,2,3));
    EXPECT_THAT(Tuples[1], ElementsAre(1,2,4));
    EXPECT_THAT(Tuples[2], ElementsAre(1,2,5));
    EXPECT_THAT(Tuples[3], ElementsAre(1,3,3));
    EXPECT_THAT(Tuples[4], ElementsAre(1,3,4));
    EXPECT_THAT(Tuples[5], ElementsAre(1,3,5));
  }

  // Multidimensional, row major, tuple elements
  {
    ovk::interval<int,3> Interval({1,2,3}, {2,4,6});
    ovk::core::vector<ovk::elem<int,3>> Tuples;
    ovk::core::ForEach<ovk::array_layout::ROW_MAJOR>(Interval, [&](int i, int j, int k) {
      Tuples.Append(ovk::elem<int,3>(i,j,k));
    });
    EXPECT_EQ(Tuples.Count(), 6);
    EXPECT_THAT(Tuples[0], ElementsAre(1,2,3));
    EXPECT_THAT(Tuples[1], ElementsAre(1,2,4));
    EXPECT_THAT(Tuples[2], ElementsAre(1,2,5));
    EXPECT_THAT(Tuples[3], ElementsAre(1,3,3));
    EXPECT_THAT(Tuples[4], ElementsAre(1,3,4));
    EXPECT_THAT(Tuples[5], ElementsAre(1,3,5));
  }

  // Multidimensional, column major, tuple
  {
    ovk::interval<int,3> Interval({1,2,3}, {2,4,6});
    ovk::core::vector<ovk::elem<int,3>> Tuples;
    ovk::core::ForEach<ovk::array_layout::COLUMN_MAJOR>(Interval, [&](const ovk::elem<int,3>
      &Tuple) {
      Tuples.Append(Tuple);
    });
    EXPECT_EQ(Tuples.Count(), 6);
    EXPECT_THAT(Tuples[0], ElementsAre(1,2,3));
    EXPECT_THAT(Tuples[1], ElementsAre(1,3,3));
    EXPECT_THAT(Tuples[2], ElementsAre(1,2,4));
    EXPECT_THAT(Tuples[3], ElementsAre(1,3,4));
    EXPECT_THAT(Tuples[4], ElementsAre(1,2,5));
    EXPECT_THAT(Tuples[5], ElementsAre(1,3,5));
  }

  // Multidimensional, column major, tuple elements
  {
    ovk::interval<int,3> Interval({1,2,3}, {2,4,6});
    ovk::core::vector<ovk::elem<int,3>> Tuples;
    ovk::core::ForEach<ovk::array_layout::COLUMN_MAJOR>(Interval, [&](int i, int j, int k) {
      Tuples.Append(ovk::elem<int,3>(i,j,k));
    });
    EXPECT_EQ(Tuples.Count(), 6);
    EXPECT_THAT(Tuples[0], ElementsAre(1,2,3));
    EXPECT_THAT(Tuples[1], ElementsAre(1,3,3));
    EXPECT_THAT(Tuples[2], ElementsAre(1,2,4));
    EXPECT_THAT(Tuples[3], ElementsAre(1,3,4));
    EXPECT_THAT(Tuples[4], ElementsAre(1,2,5));
    EXPECT_THAT(Tuples[5], ElementsAre(1,3,5));
  }

}
