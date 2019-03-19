// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Interval.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>
#include <ovk/core/Elem.hpp>

#include <mpi.h>

using testing::ElementsAre;

class IntervalTests : public tests::mpi_test {};

namespace ovk {
namespace core {
template <typename T, int N> class test_helper<interval<T,N>> {
public:
  static const elem<T,N> &GetBegin(const interval<T,N> &Interval) { return Interval.Begin_; }
  static const elem<T,N> &GetEnd(const interval<T,N> &Interval) { return Interval.End_; }
};
}}

TEST_F(IntervalTests, Meta) {

  if (TestComm().Rank() != 0) return;

  EXPECT_TRUE((std::is_same<typename ovk::interval<int,3>::value_type, int>::value));
  EXPECT_EQ(int(ovk::interval<int,3>::Rank), 3);
  EXPECT_TRUE((std::is_same<typename ovk::interval<int,3>::tuple_type, ovk::elem<int,3>>::value));

}

TEST_F(IntervalTests, Create) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::interval<int,3>>;

  // Default
  {
    // Just check that it compiles
    ovk::interval<int,3> Interval;
    // Suppress warnings about unused variable
    EXPECT_NE(&Interval, nullptr);
  }

  // Begin/end
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_THAT(helper::GetBegin(Interval), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetEnd(Interval), ElementsAre(4,5,6));
  }

}

TEST_F(IntervalTests, MakeEmpty) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::interval<int,3>>;

  ovk::interval<int,3> Interval = ovk::MakeEmptyInterval<int,3>();

  EXPECT_THAT(helper::GetBegin(Interval), ElementsAre(0,0,0));
  EXPECT_THAT(helper::GetEnd(Interval), ElementsAre(0,0,0));

}

TEST_F(IntervalTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::interval<int,3>>;

  // Copy construct
  {
    ovk::interval<int,3> Interval1({1,2,3}, {4,5,6});
    ovk::interval<int,3> Interval2(Interval1);
    EXPECT_THAT(helper::GetBegin(Interval2), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetEnd(Interval2), ElementsAre(4,5,6));
  }

  // Copy assign
  {
    ovk::interval<int,3> Interval1({1,2,3}, {4,5,6});
    ovk::interval<int,3> Interval2;
    Interval2 = Interval1;
    EXPECT_THAT(helper::GetBegin(Interval2), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetEnd(Interval2), ElementsAre(4,5,6));
  }

}

TEST_F(IntervalTests, Convert) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::interval<int,3>>;

  // Convert construct
  {
    ovk::interval<int,3> Interval1({1,2,3}, {4,5,6});
    ovk::interval<long long,3> Interval2(Interval1);
    EXPECT_THAT(helper::GetBegin(Interval2), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetEnd(Interval2), ElementsAre(4,5,6));
  }

  // Convert assign
  {
    ovk::interval<int,3> Interval1({1,2,3}, {4,5,6});
    ovk::interval<long long,3> Interval2;
    Interval2 = Interval1;
    EXPECT_THAT(helper::GetBegin(Interval2), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetEnd(Interval2), ElementsAre(4,5,6));
  }

}

TEST_F(IntervalTests, Equality) {

  if (TestComm().Rank() != 0) return;

  // Self
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_TRUE(Interval == Interval);
    EXPECT_FALSE(Interval != Interval);
  }

  // Other with same values
  {
    ovk::interval<int,3> Interval1({1,2,3}, {4,5,6});
    ovk::interval<int,3> Interval2({1,2,3}, {4,5,6});
    EXPECT_TRUE(Interval1 == Interval2);
    EXPECT_FALSE(Interval1 != Interval2);
  }

  // Different begin
  {
    ovk::interval<int,3> Interval1({1,2,3}, {4,5,6});
    ovk::interval<int,3> Interval2({1,3,2}, {4,5,6});
    EXPECT_FALSE(Interval1 == Interval2);
    EXPECT_TRUE(Interval1 != Interval2);
  }

  // Different end
  {
    ovk::interval<int,3> Interval1({1,2,3}, {4,5,6});
    ovk::interval<int,3> Interval2({1,2,3}, {5,4,6});
    EXPECT_FALSE(Interval1 == Interval2);
    EXPECT_TRUE(Interval1 != Interval2);
  }

}

TEST_F(IntervalTests, Begin) {

  if (TestComm().Rank() != 0) return;

  ovk::interval<int,3> Interval({1,2,3}, {4,5,6});

  EXPECT_THAT(Interval.Begin(), ElementsAre(1,2,3));
  EXPECT_EQ(Interval.Begin(0), 1);
  EXPECT_EQ(Interval.Begin(1), 2);
  EXPECT_EQ(Interval.Begin(2), 3);

}

TEST_F(IntervalTests, End) {

  if (TestComm().Rank() != 0) return;

  ovk::interval<int,3> Interval({1,2,3}, {4,5,6});

  EXPECT_THAT(Interval.End(), ElementsAre(4,5,6));
  EXPECT_EQ(Interval.End(0), 4);
  EXPECT_EQ(Interval.End(1), 5);
  EXPECT_EQ(Interval.End(2), 6);

}

TEST_F(IntervalTests, Size) {

  if (TestComm().Rank() != 0) return;

  // Empty
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,-5,6});
    EXPECT_THAT(Interval.Size(), ElementsAre(3,0,3));
    EXPECT_EQ(Interval.Size(0), 3);
    EXPECT_EQ(Interval.Size(1), 0);
    EXPECT_EQ(Interval.Size(2), 3);
  }

  // Non-empty
  {
    ovk::interval<int,3> Interval({1,2,3}, {6,5,4});
    EXPECT_THAT(Interval.Size(), ElementsAre(5,3,1));
    EXPECT_EQ(Interval.Size(0), 5);
    EXPECT_EQ(Interval.Size(1), 3);
    EXPECT_EQ(Interval.Size(2), 1);
  }

}

TEST_F(IntervalTests, Count) {

  if (TestComm().Rank() != 0) return;

  // Empty
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,-5,6});
    EXPECT_EQ(Interval.Count(), 0);
  }

  // Non-empty
  {
    ovk::interval<int,3> Interval({1,2,3}, {6,5,4});
    EXPECT_EQ(Interval.Count(), 15);
  }

}

TEST_F(IntervalTests, Volume) {

  if (TestComm().Rank() != 0) return;

  // Empty
  {
    ovk::interval<double,3> Interval({1.,2.,3.}, {4.,-5.,6.});
    EXPECT_EQ(Interval.Volume(), 0.);
  }

  // Non-empty
  {
    ovk::interval<double,3> Interval({1.,2.,3.}, {6.,5.,4.});
    EXPECT_EQ(Interval.Volume(), 15.);
  }

}

TEST_F(IntervalTests, Empty) {

  if (TestComm().Rank() != 0) return;

  // Empty
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,-5,6});
    EXPECT_TRUE(Interval.Empty());
  }

  // Non-empty
  {
    ovk::interval<int,3> Interval({1,2,3}, {6,5,4});
    EXPECT_FALSE(Interval.Empty());
  }

}

TEST_F(IntervalTests, Contains) {

  if (TestComm().Rank() != 0) return;

  // Inside
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_TRUE(Interval.Contains({2,3,4}));
  }

  // Outside
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_FALSE(Interval.Contains({0,0,0}));
  }

  // Empty
  {
    ovk::interval<int,3> Interval = ovk::MakeEmptyInterval<int,3>();
    EXPECT_FALSE(Interval.Contains({0,0,0}));
  }

}

TEST_F(IntervalTests, Includes) {

  if (TestComm().Rank() != 0) return;

  // Equal
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_TRUE(Interval.Includes({{1,2,3}, {4,5,6}}));
  }

  // Inside
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_TRUE(Interval.Includes({{2,3,4}, {3,4,5}}));
  }

  // Outside
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_FALSE(Interval.Includes({{6,7,8}, {9,10,11}}));
  }

  // One dimension partially outside, lower
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_FALSE(Interval.Includes({{1,2,1}, {4,5,4}}));
  }

  // One dimension partially outside, upper
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_FALSE(Interval.Includes({{1,2,5}, {4,5,8}}));
  }

  // Multiple dimensions partially outside, lower
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_FALSE(Interval.Includes({{-1,0,1}, {2,3,4}}));
  }

  // Multiple dimensions partially outside, upper
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_FALSE(Interval.Includes({{3,4,5}, {6,7,8}}));
  }

  // Enclosing interval empty
  {
    ovk::interval<int,3> Interval = ovk::MakeEmptyInterval<int,3>();
    EXPECT_FALSE(Interval.Includes({{1,2,3}, {4,5,6}}));
  }

  // Enclosed interval empty
  {
    ovk::interval<int,3> Interval({1,2,3}, {4,5,6});
    EXPECT_TRUE(Interval.Includes({{1,2,3}, {4,2,6}}));
  }

}

// TODO: Add remaining tests (are there any?)
