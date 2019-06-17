// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/ArrayOps.hpp>

#include "tests/MPITest.hpp"
#include "tests/mocks/MultidimArray.hpp"
#include "tests/mocks/Noncopyable.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>

#include <mpi.h>

#include <type_traits>
#include <utility>
#include <vector>

using testing::ElementsAre;

class ArrayOpsTests : public tests::mpi_test {};

using tests::multidim_array_row;
using tests::multidim_array_col;
using tests::noncopyable;

TEST_F(ArrayOpsTests, Fill) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;
  using array_view_noncopyable = ovk::array_view<noncopyable<int>>;

  // Constant value, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    ovk::ArrayFill(Array, 1);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 27);
  }

  // Constant value, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    ovk::ArrayFill(View, 1);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += View[i];
    EXPECT_EQ(Sum, 27);
  }

  // Initializer list, array
  {
    multidim_array Array({{1,2,3}}, 0);
    ovk::ArrayFill(Array, {0,1,2,3,4,5});
    int Sum = 0;
    for (int i = 0; i < 6; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 15);
  }

  // Initializer list, array view
  {
    multidim_array Array({{1,2,3}}, 0);
    array_view View(Array);
    ovk::ArrayFill(View, {0,1,2,3,4,5});
    int Sum = 0;
    for (int i = 0; i < 6; ++i) Sum += View[i];
    EXPECT_EQ(Sum, 15);
  }

  // Iterator, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    std::array<int,27> SourceArray;
    for (int i = 0; i < 27; ++i) SourceArray[i] = i;
    ovk::ArrayFill(Array, SourceArray.begin());
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 351);
  }

  // Iterator, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    std::array<int,27> SourceArray;
    for (int i = 0; i < 27; ++i) SourceArray[i] = i;
    ovk::ArrayFill(View, SourceArray.begin());
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += View[i];
    EXPECT_EQ(Sum, 351);
  }

  // Array view, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    std::array<int,27> SourceArray;
    for (int i = 0; i < 27; ++i) SourceArray[i] = i;
    array_view SourceView(SourceArray.data(), {{1,2,3}, {4,5,6}});
    ovk::ArrayFill(Array, SourceView);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 351);
  }

  // Array view, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    std::array<int,27> SourceArray;
    for (int i = 0; i < 27; ++i) SourceArray[i] = i;
    array_view SourceView(SourceArray.data(), {{1,2,3}, {4,5,6}});
    ovk::ArrayFill(View, SourceView);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += View[i];
    EXPECT_EQ(Sum, 351);
  }

  // Array, lvalue ref, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    multidim_array SourceArray({{1,2,3}, {4,5,6}});
    for (int i = 0; i < 27; ++i) SourceArray[i] = i;
    ovk::ArrayFill(Array, SourceArray);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 351);
  }

  // Array, lvalue ref, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    multidim_array SourceArray({{1,2,3}, {4,5,6}});
    for (int i = 0; i < 27; ++i) SourceArray[i] = i;
    ovk::ArrayFill(View, SourceArray);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += View[i];
    EXPECT_EQ(Sum, 351);
  }

  // Array, rvalue ref, array
  {
    std::array<noncopyable<int>,4> Array;
    std::array<noncopyable<int>,4> SourceArray;
    for (int i = 0; i < 4; ++i) SourceArray[i] = 1;
    ovk::ArrayFill(Array, std::move(SourceArray));
    int Sum = 0;
    for (int i = 0; i < 4; ++i) Sum += Array[i].Value();
    EXPECT_EQ(Sum, 4);
    Sum = 0;
    for (int i = 0; i < 4; ++i) Sum += SourceArray[i].Value();
    EXPECT_EQ(Sum, 0);
  }

  // Array, rvalue ref, array view
  {
    std::array<noncopyable<int>,4> Array;
    array_view_noncopyable View(Array);
    std::array<noncopyable<int>,4> SourceArray;
    for (int i = 0; i < 4; ++i) SourceArray[i] = 1;
    ovk::ArrayFill(View, std::move(SourceArray));
    int Sum = 0;
    for (int i = 0; i < 4; ++i) Sum += View[i].Value();
    EXPECT_EQ(Sum, 4);
    Sum = 0;
    for (int i = 0; i < 4; ++i) Sum += SourceArray[i].Value();
    EXPECT_EQ(Sum, 0);
  }

}

TEST_F(ArrayOpsTests, None) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Convertible to bool, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    EXPECT_TRUE(ovk::ArrayNone(Array));
    Array(2,3,4) = 1;
    EXPECT_FALSE(ovk::ArrayNone(Array));
  }

  // Convertible to bool, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayNone(View));
    View(2,3,4) = 1;
    EXPECT_FALSE(ovk::ArrayNone(View));
  }

  // With condition function, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_TRUE(ovk::ArrayNone(Array, [](int Value) -> bool { return Value == 5; }));
    Array(2,3,4) = 5;
    EXPECT_FALSE(ovk::ArrayNone(Array, [](int Value) -> bool { return Value == 5; }));
  }

  // With condition function, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayNone(View, [](int Value) -> bool { return Value == 5; }));
    View(2,3,4) = 5;
    EXPECT_FALSE(ovk::ArrayNone(View, [](int Value) -> bool { return Value == 5; }));
  }

}

TEST_F(ArrayOpsTests, Any) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Convertible to bool, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    EXPECT_FALSE(ovk::ArrayAny(Array));
    Array(2,3,4) = 1;
    EXPECT_TRUE(ovk::ArrayAny(Array));
  }

  // Convertible to bool, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayAny(View));
    View(2,3,4) = 1;
    EXPECT_TRUE(ovk::ArrayAny(View));
  }

  // With condition function, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_FALSE(ovk::ArrayAny(Array, [](int Value) -> bool { return Value == 5; }));
    Array(2,3,4) = 5;
    EXPECT_TRUE(ovk::ArrayAny(Array, [](int Value) -> bool { return Value == 5; }));
  }

  // With condition function, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayAny(View, [](int Value) -> bool { return Value == 5; }));
    View(2,3,4) = 5;
    EXPECT_TRUE(ovk::ArrayAny(View, [](int Value) -> bool { return Value == 5; }));
  }

}

TEST_F(ArrayOpsTests, NotAll) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Convertible to bool, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_FALSE(ovk::ArrayNotAll(Array));
    Array(2,3,4) = 0;
    EXPECT_TRUE(ovk::ArrayNotAll(Array));
  }

  // Convertible to bool, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayNotAll(View));
    View(2,3,4) = 0;
    EXPECT_TRUE(ovk::ArrayNotAll(View));
  }

  // With condition function, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_FALSE(ovk::ArrayNotAll(Array, [](int Value) -> bool { return Value == 1; }));
    Array(2,3,4) = 5;
    EXPECT_TRUE(ovk::ArrayNotAll(Array, [](int Value) -> bool { return Value == 1; }));
  }

  // With condition function, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayNotAll(View, [](int Value) -> bool { return Value == 1; }));
    View(2,3,4) = 5;
    EXPECT_TRUE(ovk::ArrayNotAll(View, [](int Value) -> bool { return Value == 1; }));
  }

}

TEST_F(ArrayOpsTests, All) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Convertible to bool, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_TRUE(ovk::ArrayAll(Array));
    Array(2,3,4) = 0;
    EXPECT_FALSE(ovk::ArrayAll(Array));
  }

  // Convertible to bool, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayAll(View));
    View(2,3,4) = 0;
    EXPECT_FALSE(ovk::ArrayAll(View));
  }

  // With condition function, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_TRUE(ovk::ArrayAll(Array, [](int Value) -> bool { return Value == 1; }));
    Array(2,3,4) = 5;
    EXPECT_FALSE(ovk::ArrayAll(Array, [](int Value) -> bool { return Value == 1; }));
  }

  // With condition function, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayAll(View, [](int Value) -> bool { return Value == 1; }));
    View(2,3,4) = 5;
    EXPECT_FALSE(ovk::ArrayAll(View, [](int Value) -> bool { return Value == 1; }));
  }

}

TEST_F(ArrayOpsTests, Min) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_EQ(ovk::ArrayMin(Array), 1);
    Array(2,3,4) = 5;
    EXPECT_EQ(ovk::ArrayMin(Array), 1);
  }

  // Array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_EQ(ovk::ArrayMin(View), 1);
    View(2,3,4) = 5;
    EXPECT_EQ(ovk::ArrayMin(View), 1);
  }

}

TEST_F(ArrayOpsTests, Max) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_EQ(ovk::ArrayMax(Array), 1);
    Array(2,3,4) = 5;
    EXPECT_EQ(ovk::ArrayMax(Array), 5);
  }

  // Array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_EQ(ovk::ArrayMax(View), 1);
    View(2,3,4) = 5;
    EXPECT_EQ(ovk::ArrayMax(View), 5);
  }

}

TEST_F(ArrayOpsTests, Sum) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Array
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    for (int i = 0; i < 27; ++i) Array[i] = i;
    EXPECT_EQ(ovk::ArraySum(Array), 351);
  }

  // Array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View(Array);
    for (int i = 0; i < 27; ++i) View[i] = i;
    EXPECT_EQ(ovk::ArraySum(View), 351);
  }

}
