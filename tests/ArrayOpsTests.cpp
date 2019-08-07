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
using testing::ElementsAreArray;

class ArrayOpsTests : public tests::mpi_test {};

using tests::multidim_array_row;
using tests::multidim_array_col;
using tests::noncopyable;

TEST_F(ArrayOpsTests, ForEach) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Array
  {
    multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    ovk::ArrayForEach(Array, [](int &Value) { Value += 1; });
    EXPECT_THAT(Array, ElementsAreArray({1,2,3,4,5,6,7,8}));
  }

  // Array view
  {
    multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    array_view View(Array);
    ovk::ArrayForEach(View, [](int &Value) { Value += 1; });
    EXPECT_THAT(View, ElementsAreArray({1,2,3,4,5,6,7,8}));
  }

//   // Function takes value and index, default index type, array
//   {
//     multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
//     ovk::ArrayForEach(Array, [](int &Value, long long iValue) { Value += (iValue % 2); });
//     EXPECT_THAT(Array, ElementsAreArray({0,2,2,4,4,6,6,8}));
//   }

//   // Function takes value and index, default index type, array view
//   {
//     multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
//     array_view View(Array);
//     ovk::ArrayForEach(View, [](int &Value, long long iValue) { Value += (iValue % 2); });
//     EXPECT_THAT(View, ElementsAreArray({0,2,2,4,4,6,6,8}));
//   }

//   // Function takes value and index, explicit index type, array
//   {
//     multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
//     ovk::ArrayForEach<int>(Array, [](int &Value, int iValue) { Value += (iValue % 2); });
//     EXPECT_THAT(Array, ElementsAreArray({0,2,2,4,4,6,6,8}));
//   }

//   // Function takes value and index, explicit index type, array view
//   {
//     multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
//     array_view View(Array);
//     ovk::ArrayForEach<int>(View, [](int &Value, int iValue) { Value += (iValue % 2); });
//     EXPECT_THAT(View, ElementsAreArray({0,2,2,4,4,6,6,8}));
//   }

//   // Function takes value and tuple, default tuple element type, array
//   {
//     multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
//     ovk::ArrayForEach(Array, [](int &Value, const ovk::elem<long long,3> &Tuple) {
//       Value += (Tuple(0)-1) + (Tuple(1)-2) + (Tuple(2)-3);
//     });
//     EXPECT_THAT(Array, ElementsAreArray({0,2,3,5,5,7,8,10}));
//   }

//   // Function takes value and tuple, default tuple element type, array view
//   {
//     multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
//     array_view View(Array);
//     ovk::ArrayForEach(View, [](int &Value, const ovk::elem<long long,3> &Tuple) {
//       Value += (Tuple(0)-1) + (Tuple(1)-2) + (Tuple(2)-3);
//     });
//     EXPECT_THAT(View, ElementsAreArray({0,2,3,5,5,7,8,10}));
//   }

//   // Function takes value and tuple, explicit tuple element type, array
//   {
//     multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
//     ovk::ArrayForEach<int>(Array, [](int &Value, const ovk::elem<int,3> &Tuple) {
//       Value += (Tuple(0)-1) + (Tuple(1)-2) + (Tuple(2)-3);
//     });
//     EXPECT_THAT(Array, ElementsAreArray({0,2,3,5,5,7,8,10}));
//   }

//   // Function takes value and tuple, explicit tuple element type, array view
//   {
//     multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
//     array_view View(Array);
//     ovk::ArrayForEach<int>(View, [](int &Value, const ovk::elem<int,3> &Tuple) {
//       Value += (Tuple(0)-1) + (Tuple(1)-2) + (Tuple(2)-3);
//     });
//     EXPECT_THAT(View, ElementsAreArray({0,2,3,5,5,7,8,10}));
//   }

}

TEST_F(ArrayOpsTests, Transform) {

  if (TestComm().Rank() != 0) return;

  using multidim_array_1 = multidim_array_row<int>;
  using multidim_array_2 = multidim_array_row<double>;
  using array_view_1 = ovk::array_view<int,3>;
  using array_view_2 = ovk::array_view<double,3>;

  // Array, array
  {
    multidim_array_1 Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    multidim_array_2 TransformedArray({{1,2,3}, {3,4,5}}, 0.);
    ovk::ArrayTransform(Array, TransformedArray, [](int Value) -> double {
      return double(Value+1);
    });
    EXPECT_THAT(TransformedArray, ElementsAreArray({1.,2.,3.,4.,5.,6.,7.,8.}));
  }

  // Array, array view
  {
    multidim_array_1 Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    multidim_array_2 TransformedArray({{1,2,3}, {3,4,5}}, 0.);
    array_view_2 TransformedView(TransformedArray);
    ovk::ArrayTransform(Array, TransformedView, [](int Value) -> double {
      return double(Value+1);
    });
    EXPECT_THAT(TransformedView, ElementsAreArray({1.,2.,3.,4.,5.,6.,7.,8.}));
  }

  // Array view, array
  {
    multidim_array_1 Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    multidim_array_2 TransformedArray({{1,2,3}, {3,4,5}}, 0.);
    array_view_1 View(Array);
    ovk::ArrayTransform(View, TransformedArray, [](int Value) -> double {
      return double(Value+1);
    });
    EXPECT_THAT(TransformedArray, ElementsAreArray({1.,2.,3.,4.,5.,6.,7.,8.}));
  }

  // Array view, array view
  {
    multidim_array_1 Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    multidim_array_2 TransformedArray({{1,2,3}, {3,4,5}}, 0.);
    array_view_1 View(Array);
    array_view_2 TransformedView(TransformedArray);
    ovk::ArrayTransform(View, TransformedView, [](int Value) -> double {
      return double(Value+1);
    });
    EXPECT_THAT(TransformedView, ElementsAreArray({1.,2.,3.,4.,5.,6.,7.,8.}));
  }

}

TEST_F(ArrayOpsTests, Reduce) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Array
  {
    multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    bool Result = ovk::ArrayReduce(Array, false, [](bool &Partial, int Value) {
      Partial = Partial ^ ((Value % 3) == 0);
    });
    EXPECT_TRUE(Result);
  }

  // Array view
  {
    multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    array_view View(Array);
    bool Result = ovk::ArrayReduce(View, false, [](bool &Partial, int Value) {
      Partial = Partial ^ ((Value % 3) == 0);
    });
    EXPECT_TRUE(Result);
  }

}

TEST_F(ArrayOpsTests, Collapse) {

  if (TestComm().Rank() != 0) return;

  using multidim_array = multidim_array_row<int>;
  using array_view = ovk::array_view<int,3>;

  // Array
  {
    multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    int Sum = ovk::ArrayCollapse(Array, [](int &Partial, int Value) { Partial += Value; });
    EXPECT_EQ(Sum, 28);
  }

  // Array view
  {
    multidim_array Array({{1,2,3}, {3,4,5}}, {0,1,2,3,4,5,6,7});
    array_view View(Array);
    int Sum = ovk::ArrayCollapse(View, [](int &Partial, int Value) { Partial += Value; });
    EXPECT_EQ(Sum, 28);
  }

}

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

  // Interval and constant value, array
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    ovk::ArrayFill(Array, {{2,3,4}, {4,5,6}}, 1);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 8);
  }

  // Interval and constant value, array view
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    ovk::ArrayFill(View, {{2,3,4}, {4,5,6}}, 1);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += View[i];
    EXPECT_EQ(Sum, 8);
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

  // Convertible to bool, array, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 1);
    EXPECT_TRUE(ovk::ArrayNone(Array));
  }

  // Convertible to bool, array, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    EXPECT_TRUE(ovk::ArrayNone(Array));
    Array(2,3,4) = 1;
    EXPECT_FALSE(ovk::ArrayNone(Array));
  }

  // Convertible to bool, array view, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 1);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayNone(View));
  }

  // Convertible to bool, array view, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayNone(View));
    View(2,3,4) = 1;
    EXPECT_FALSE(ovk::ArrayNone(View));
  }

  // With condition function, array, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 5);
    EXPECT_TRUE(ovk::ArrayNone(Array, [](int Value) -> bool { return Value == 5; }));
  }

  // With condition function, array, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_TRUE(ovk::ArrayNone(Array, [](int Value) -> bool { return Value == 5; }));
    Array(2,3,4) = 5;
    EXPECT_FALSE(ovk::ArrayNone(Array, [](int Value) -> bool { return Value == 5; }));
  }

  // With condition function, array view, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 5);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayNone(View, [](int Value) -> bool { return Value == 5; }));
  }

  // With condition function, array view, non-empty
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

  // Convertible to bool, array, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 1);
    EXPECT_FALSE(ovk::ArrayAny(Array));
  }

  // Convertible to bool, array, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    EXPECT_FALSE(ovk::ArrayAny(Array));
    Array(2,3,4) = 1;
    EXPECT_TRUE(ovk::ArrayAny(Array));
  }

  // Convertible to bool, array view, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 1);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayAny(View));
  }

  // Convertible to bool, array view, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayAny(View));
    View(2,3,4) = 1;
    EXPECT_TRUE(ovk::ArrayAny(View));
  }

  // With condition function, array, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 5);
    EXPECT_FALSE(ovk::ArrayAny(Array, [](int Value) -> bool { return Value == 5; }));
  }

  // With condition function, array, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_FALSE(ovk::ArrayAny(Array, [](int Value) -> bool { return Value == 5; }));
    Array(2,3,4) = 5;
    EXPECT_TRUE(ovk::ArrayAny(Array, [](int Value) -> bool { return Value == 5; }));
  }

  // With condition function, array view, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 5);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayAny(View, [](int Value) -> bool { return Value == 5; }));
  }

  // With condition function, array view, non-empty
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

  // Convertible to bool, array, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 0);
    EXPECT_FALSE(ovk::ArrayNotAll(Array));
  }

  // Convertible to bool, array, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_FALSE(ovk::ArrayNotAll(Array));
    Array(2,3,4) = 0;
    EXPECT_TRUE(ovk::ArrayNotAll(Array));
  }

  // Convertible to bool, array view, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 0);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayNotAll(View));
  }

  // Convertible to bool, array view, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayNotAll(View));
    View(2,3,4) = 0;
    EXPECT_TRUE(ovk::ArrayNotAll(View));
  }

  // With condition function, array, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 5);
    EXPECT_FALSE(ovk::ArrayNotAll(Array, [](int Value) -> bool { return Value == 1; }));
  }

  // With condition function, array, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_FALSE(ovk::ArrayNotAll(Array, [](int Value) -> bool { return Value == 1; }));
    Array(2,3,4) = 5;
    EXPECT_TRUE(ovk::ArrayNotAll(Array, [](int Value) -> bool { return Value == 1; }));
  }

  // With condition function, array view, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 5);
    array_view View(Array);
    EXPECT_FALSE(ovk::ArrayNotAll(View, [](int Value) -> bool { return Value == 1; }));
  }

  // With condition function, array view, non-empty
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

  // Convertible to bool, array, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 0);
    EXPECT_TRUE(ovk::ArrayAll(Array));
  }

  // Convertible to bool, array, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_TRUE(ovk::ArrayAll(Array));
    Array(2,3,4) = 0;
    EXPECT_FALSE(ovk::ArrayAll(Array));
  }

  // Convertible to bool, array view, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 0);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayAll(View));
  }

  // Convertible to bool, array view, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayAll(View));
    View(2,3,4) = 0;
    EXPECT_FALSE(ovk::ArrayAll(View));
  }

  // With condition function, array, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 5);
    EXPECT_TRUE(ovk::ArrayAll(Array, [](int Value) -> bool { return Value == 1; }));
  }

  // With condition function, array, non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    EXPECT_TRUE(ovk::ArrayAll(Array, [](int Value) -> bool { return Value == 1; }));
    Array(2,3,4) = 5;
    EXPECT_FALSE(ovk::ArrayAll(Array, [](int Value) -> bool { return Value == 1; }));
  }

  // With condition function, array view, empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}}, 5);
    array_view View(Array);
    EXPECT_TRUE(ovk::ArrayAll(View, [](int Value) -> bool { return Value == 1; }));
  }

  // With condition function, array view, non-empty
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
