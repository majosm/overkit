// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Elem.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Comm.hpp>

#include <mpi.h>

using testing::ElementsAre;

class ElemTests : public tests::mpi_test {};

namespace ovk {
namespace core {
template <typename T, int N> class test_helper<elem<T,N>> {
public:
  static const T *GetValues(const elem<T,N> &Elem) { return Elem.Values_; }
};
}}

TEST_F(ElemTests, Meta) {

  if (TestComm().Rank() != 0) return;

  EXPECT_TRUE((std::is_same<typename ovk::elem<int,3>::value_type, int>::value));
  EXPECT_EQ(int(ovk::elem<int,3>::Rank), 3);

}

TEST_F(ElemTests, Create) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::elem<int,3>>;

  // Default
  {
    // Just check that it compiles
    ovk::elem<int,3> Elem;
    // Suppress warnings about unused variable
    EXPECT_NE(&Elem, nullptr);
  }

  // Values
  {
    ovk::elem<int,3> Elem = {1,2,3};
    const int *Values = helper::GetValues(Elem);
    EXPECT_EQ(Values[0], 1);
    EXPECT_EQ(Values[1], 2);
    EXPECT_EQ(Values[2], 3);
  }

}

TEST_F(ElemTests, MakeUniform) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::elem<int,3>>;

  ovk::elem<int,3> Elem = ovk::MakeUniformElem<int,3>(1);

  const int *Values = helper::GetValues(Elem);
  EXPECT_EQ(Values[0], 1);
  EXPECT_EQ(Values[1], 1);
  EXPECT_EQ(Values[2], 1);

}

TEST_F(ElemTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::elem<int,3>>;

  // Copy construct
  {
    ovk::elem<int,3> Elem1 = {1,2,3};
    ovk::elem<int,3> Elem2 = Elem1;
    const int *Values = helper::GetValues(Elem2);
    EXPECT_EQ(Values[0], 1);
    EXPECT_EQ(Values[1], 2);
    EXPECT_EQ(Values[2], 3);
  }

  // Copy assign
  {
    ovk::elem<int,3> Elem1 = {1,2,3};
    ovk::elem<int,3> Elem2;
    Elem2 = Elem1;
    const int *Values = helper::GetValues(Elem2);
    EXPECT_EQ(Values[0], 1);
    EXPECT_EQ(Values[1], 2);
    EXPECT_EQ(Values[2], 3);
  }

}

// Can't put these in anonymous namespace because clang optimizes them out and emits warnings
namespace elem_tests_internal {

struct intlike {
  int Value;
  intlike() = default;
  explicit intlike(int Value_):
    Value(Value_)
  {}
};

std::true_type ConvertsToElemLongLong(const ovk::elem<long long,3> &) { return {}; }
std::false_type ConvertsToElemLongLong(...) { return {}; }
std::true_type ConvertsToElemIntlike(const ovk::elem<intlike,3> &) { return {}; }
std::false_type ConvertsToElemIntlike(...) { return {}; }
std::true_type ConvertsToElemRank2(const ovk::elem<int,2> &) { return {}; }
std::false_type ConvertsToElemRank2(...) { return {}; }

}

TEST_F(ElemTests, ConvertToOtherElem) {

  if (TestComm().Rank() != 0) return;

  using namespace elem_tests_internal;

  using helper_long_long = ovk::core::test_helper<ovk::elem<long long,3>>;
  using helper_intlike = ovk::core::test_helper<ovk::elem<intlike,3>>;

  // Implicit convert construct
  {
    ovk::elem<int,3> Elem1 = {1,2,3};
    ovk::elem<long long,3> Elem2 = Elem1;
    const long long *Values = helper_long_long::GetValues(Elem2);
    EXPECT_EQ(Values[0], 1);
    EXPECT_EQ(Values[1], 2);
    EXPECT_EQ(Values[2], 3);
  }

  // Explicit convert construct
  {
    ovk::elem<int,3> Elem1 = {1,2,3};
    ovk::elem<intlike,3> Elem2(Elem1);
    const intlike *Values = helper_intlike::GetValues(Elem2);
    EXPECT_EQ(Values[0].Value, 1);
    EXPECT_EQ(Values[1].Value, 2);
    EXPECT_EQ(Values[2].Value, 3);
  }

  // Implicit convert assign
  {
    ovk::elem<int,3> Elem1 = {1,2,3};
    ovk::elem<long long,3> Elem2;
    Elem2 = Elem1;
    const long long *Values = helper_long_long::GetValues(Elem2);
    EXPECT_EQ(Values[0], 1);
    EXPECT_EQ(Values[1], 2);
    EXPECT_EQ(Values[2], 3);
  }

  // Make sure implicit conversion is only allowed for value types that are implicitly convertible
  {
    ovk::elem<int,3> Elem;
    EXPECT_TRUE(decltype(ConvertsToElemLongLong(Elem))::value);
    EXPECT_FALSE(decltype(ConvertsToElemIntlike(Elem))::value);
    EXPECT_FALSE(decltype(ConvertsToElemRank2(Elem))::value);
    // Sanity check / suppress unused function warnings
    EXPECT_FALSE(decltype(ConvertsToElemLongLong(ovk::elem<intlike,3>()))::value);
    EXPECT_TRUE(decltype(ConvertsToElemIntlike(ovk::elem<intlike,3>()))::value);
    EXPECT_TRUE(decltype(ConvertsToElemRank2(ovk::elem<int,2>()))::value);
  }

}

TEST_F(ElemTests, Equality) {

  if (TestComm().Rank() != 0) return;

  // Self
  {
    ovk::elem<int,3> Elem = {1,2,3};
    EXPECT_TRUE(Elem == Elem);
    EXPECT_FALSE(Elem != Elem);
  }

  // Other with same values
  {
    ovk::elem<int,3> Elem1 = {1,2,3};
    ovk::elem<int,3> Elem2 = {1,2,3};
    EXPECT_TRUE(Elem1 == Elem2);
    EXPECT_FALSE(Elem1 != Elem2);
  }

  // Different values
  {
    ovk::elem<int,3> Elem1 = {1,2,3};
    ovk::elem<int,3> Elem2 = {1,3,2};
    EXPECT_FALSE(Elem1 == Elem2);
    EXPECT_TRUE(Elem1 != Elem2);
  }

}

TEST_F(ElemTests, ParenthesisOperator) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::elem<int,3>>;

  // Const
  {
    const ovk::elem<int,3> Elem = {1,2,3};
    EXPECT_EQ(Elem(0), 1);
    EXPECT_EQ(Elem(1), 2);
    EXPECT_EQ(Elem(2), 3);
  }

  // Non-const
  {
    ovk::elem<int,3> Elem;
    Elem(0) = 1;
    Elem(1) = 2;
    Elem(2) = 3;
    const int *Values = helper::GetValues(Elem);
    EXPECT_EQ(Values[0], 1);
    EXPECT_EQ(Values[1], 2);
    EXPECT_EQ(Values[2], 3);
  }

}

TEST_F(ElemTests, BracketOperator) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::elem<int,3>>;

  // Const
  {
    const ovk::elem<int,3> Elem = {1,2,3};
    EXPECT_EQ(Elem[0], 1);
    EXPECT_EQ(Elem[1], 2);
    EXPECT_EQ(Elem[2], 3);
  }

  // Non-const
  {
    ovk::elem<int,3> Elem;
    Elem[0] = 1;
    Elem[1] = 2;
    Elem[2] = 3;
    const int *Values = helper::GetValues(Elem);
    EXPECT_EQ(Values[0], 1);
    EXPECT_EQ(Values[1], 2);
    EXPECT_EQ(Values[2], 3);
  }

}

TEST_F(ElemTests, Data) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::elem<int,3>>;

  // Const
  {
    const ovk::elem<int,3> Elem = {1,2,3};
    const int *Values = helper::GetValues(Elem);
    const int *Data = Elem.Data();
    EXPECT_EQ(Data, Values);
  }

  // Non-const
  {
    ovk::elem<int,3> Elem = {1,2,3};
    const int *Values = helper::GetValues(Elem);
    int *Data = Elem.Data();
    EXPECT_EQ(Data, Values);
  }

}

TEST_F(ElemTests, BeginEnd) {

  if (TestComm().Rank() != 0) return;

  using helper = ovk::core::test_helper<ovk::elem<int,3>>;

  // Const
  {
    const ovk::elem<int,3> Elem = {1,2,3};
    const int *Values = helper::GetValues(Elem);
    EXPECT_EQ(Elem.Begin(), Values);
    EXPECT_EQ(Elem.End(), Values+3);
    int Sum = 0;
    for (auto &Value : Elem) Sum += Value;
    EXPECT_EQ(Sum, 6);
  }

  // Non-const
  {
    ovk::elem<int,3> Elem = {1,2,3};
    const int *Values = helper::GetValues(Elem);
    EXPECT_EQ(Elem.Begin(), Values);
    EXPECT_EQ(Elem.End(), Values+3);
    for (auto &Value : Elem) Value = 1;
    EXPECT_EQ(Values[0], 1);
    EXPECT_EQ(Values[1], 1);
    EXPECT_EQ(Values[2], 1);
  }

}

namespace elem_tests_internal {

std::true_type ConvertsToInt(const int &) { return {}; }
std::false_type ConvertsToInt(...) { return {}; }

}

TEST_F(ElemTests, ConvertToScalar) {

  if (TestComm().Rank() != 0) return;

  using namespace elem_tests_internal;

  using helper = ovk::core::test_helper<ovk::elem<int,1>>;

  // Const, rank 1
  {
    const ovk::elem<int,1> Elem = {1};
    const int &Value = Elem;
    const int *Values = helper::GetValues(Elem);
    EXPECT_EQ(&Value, Values);
  }

  // Non-const, rank 1
  {
    ovk::elem<int,1> Elem = {1};
    int &Value = Elem;
    const int *Values = helper::GetValues(Elem);
    EXPECT_EQ(&Value, Values);
  }

  // Rank > 1 should not convert
  {
    ovk::elem<int,3> Elem = {1,2,3};
    EXPECT_FALSE(decltype(ConvertsToInt(Elem))::value);
  }

}

TEST_F(ElemTests, Traits) {

  if (TestComm().Rank() != 0) return;

  using elem = ovk::elem<int,3>;

  EXPECT_TRUE(ovk::core::IsArray<elem>());
  EXPECT_TRUE((std::is_same<ovk::core::array_value_type<elem>, int>::value));
  EXPECT_EQ(ovk::core::ArrayRank<elem>(), 1);
  EXPECT_EQ(ovk::core::ArrayLayout<elem>(), ovk::array_layout::ROW_MAJOR);
  EXPECT_TRUE(ovk::core::ArrayHasStaticExtents<elem>());
  EXPECT_TRUE((ovk::core::StaticArrayHasExtentsBegin<elem, 0>()));
  EXPECT_TRUE((ovk::core::StaticArrayHasExtentsEnd<elem, 3>()));

  elem Elem = {1,2,3};
  EXPECT_THAT(ovk::core::ArrayExtents(Elem).Begin(), ElementsAre(0));
  EXPECT_THAT(ovk::core::ArrayExtents(Elem).End(), ElementsAre(3));
  EXPECT_EQ(ovk::core::ArrayData(Elem), Elem.Data());

}

TEST_F(ElemTests, Addition) {

  if (TestComm().Rank() != 0) return;

  ovk::elem<int,3> Elem1 = {1,2,3};
  ovk::elem<int,3> Elem2 = {3,3,3};
  ovk::elem<int,3> Sum = Elem1 + Elem2;

  EXPECT_THAT(Sum, ElementsAre(4,5,6));

}

TEST_F(ElemTests, Subtraction) {

  if (TestComm().Rank() != 0) return;

  ovk::elem<int,3> Elem1 = {4,5,6};
  ovk::elem<int,3> Elem2 = {3,3,3};
  ovk::elem<int,3> Sum = Elem1 - Elem2;

  EXPECT_THAT(Sum, ElementsAre(1,2,3));

}

TEST_F(ElemTests, Min) {

  if (TestComm().Rank() != 0) return;

  ovk::elem<int,3> Elem1 = {1,2,3};
  ovk::elem<int,3> Elem2 = {3,2,1};
  ovk::elem<int,3> Min = ovk::Min(Elem1, Elem2);

  EXPECT_THAT(Min, ElementsAre(1,2,1));

}

TEST_F(ElemTests, Max) {

  if (TestComm().Rank() != 0) return;

  ovk::elem<int,3> Elem1 = {1,2,3};
  ovk::elem<int,3> Elem2 = {3,2,1};
  ovk::elem<int,3> Max = ovk::Max(Elem1, Elem2);

  EXPECT_THAT(Max, ElementsAre(3,2,3));

}

TEST_F(ElemTests, Concat) {

  if (TestComm().Rank() != 0) return;

  ovk::elem<int,3> Elem1 = {1,2,3};
  ovk::elem<int,2> Elem2 = {4,5};
  ovk::elem<int,5> Concat = ovk::ConcatElems(Elem1, Elem2);

  EXPECT_THAT(Concat, ElementsAre(1,2,3,4,5));

}
