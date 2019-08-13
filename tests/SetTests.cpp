// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Set.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <mpi.h>

#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>

class SetTests : public tests::mpi_test {};

using testing::ElementsAre;

namespace ovk {
namespace core {
template <typename ValueType, typename CompareType> class test_helper<set<ValueType, CompareType>> {
public:
  using set_type = set<ValueType, CompareType>;
  using value_type = typename set_type::value_type;
  static const array<value_type> &GetValues(const set_type &Set) { return Set.Values_; }
  static array<value_type> &GetValues(set_type &Set) { return Set.Values_; }
};
}}

TEST_F(SetTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int,std::greater<int>>;

  EXPECT_TRUE((std::is_same<typename set::value_type, int>::value));
  EXPECT_TRUE((std::is_same<typename set::compare_type, std::greater<int>>::value));
  EXPECT_TRUE((std::is_same<typename set::index_type, long long>::value));
  EXPECT_TRUE((std::is_same<typename set::iterator::pointer, const int *>::value));
  EXPECT_TRUE((std::is_same<typename set::const_iterator::pointer, const int *>::value));

}

TEST_F(SetTests, Create) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  // Default
  {
    set Set;
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 0);
  }

  // Initializer list
  {
    set Set = {2, 3};
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 3);
  }

  // Iterators
  {
    std::array<int,2> SourceValues = {{2, 3}};
    set Set(SourceValues.begin(), SourceValues.end());
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 3);
  }

}

TEST_F(SetTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  // Copy construct
  {
    set Set1 = {2, 3};
    set Set2 = Set1;
    auto &Values1 = helper::GetValues(Set1);
    auto &Values2 = helper::GetValues(Set2);
    EXPECT_EQ(Values1.Count(), 2);
    EXPECT_EQ(Values1(0), 2);
    EXPECT_EQ(Values1(1), 3);
    EXPECT_EQ(Values2.Count(), 2);
    EXPECT_EQ(Values2(0), 2);
    EXPECT_EQ(Values2(1), 3);
  }

  // Copy assign
  {
    set Set1 = {2, 3};
    set Set2;
    Set2 = Set1;
    auto &Values1 = helper::GetValues(Set1);
    auto &Values2 = helper::GetValues(Set2);
    EXPECT_EQ(Values1.Count(), 2);
    EXPECT_EQ(Values1(0), 2);
    EXPECT_EQ(Values1(1), 3);
    EXPECT_EQ(Values2.Count(), 2);
    EXPECT_EQ(Values2(0), 2);
    EXPECT_EQ(Values2(1), 3);
  }

}

TEST_F(SetTests, Move) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  // Move construct
  {
    set Set1 = {2, 3};
    set Set2(std::move(Set1));
    auto &Values1 = helper::GetValues(Set1);
    auto &Values2 = helper::GetValues(Set2);
    EXPECT_EQ(Values1.Count(), 0);
    EXPECT_EQ(Values2.Count(), 2);
    EXPECT_EQ(Values2(0), 2);
    EXPECT_EQ(Values2(1), 3);
  }

  // Move assign
  {
    set Set1 = {2, 3};
    set Set2;
    Set2 = std::move(Set1);
    auto &Values1 = helper::GetValues(Set1);
    auto &Values2 = helper::GetValues(Set2);
    EXPECT_EQ(Values1.Count(), 0);
    EXPECT_EQ(Values2.Count(), 2);
    EXPECT_EQ(Values2(0), 2);
    EXPECT_EQ(Values2(1), 3);
  }

}

TEST_F(SetTests, Assign) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  // Initializer list, operator=
  {
    set Set = {1, 2};
    Set = {2, 3};
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 3);
  }

  // Initializer list, Assign
  {
    set Set = {1, 2};
    Set.Assign({2, 3});
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 3);
  }

  // Iterators, Assign
  {
    set Set = {1, 2};
    std::array<int,2> SourceValues = {{2, 3}};
    Set.Assign(SourceValues.begin(), SourceValues.end());
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 3);
  }

}

TEST_F(SetTests, Reserve) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  set Set;
  auto &Values = helper::GetValues(Set);
  long long Capacity = Values.Capacity();
  Set.Reserve(Capacity+1);
  EXPECT_GE(Values.Capacity(), Capacity+1);

}

TEST_F(SetTests, Insert) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  // Doesn't already exist
  {
    set Set = {2, 4};
    Set.Insert(3);
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 3);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 3);
    EXPECT_EQ(Values(2), 4);
  }

  // Already exists
  {
    set Set = {2, 4};
    Set.Insert(4);
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 4);
  }

  // With lower bound iterator, doesn't already exist
  {
    set Set = {2, 4};
    auto &Values = helper::GetValues(Set);
    auto Iter = Set.Insert(Set.Begin()+1, 3);
    EXPECT_EQ(Values.Count(), 3);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 3);
    EXPECT_EQ(Values(2), 4);
    EXPECT_EQ(Iter, Set.Begin()+1);
  }

  // With lower bound iterator, already exists
  {
    set Set = {2, 4};
    auto &Values = helper::GetValues(Set);
    auto Iter = Set.Insert(Set.Begin()+1, 4);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 4);
    EXPECT_EQ(Iter, Set.Begin()+1);
  }

}

TEST_F(SetTests, Erase) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  // Key, exists
  {
    set Set = {2, 3, 4};
    Set.Erase(3);
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 4);
  }

  // Key, doesn't exist
  {
    set Set = {2, 4};
    Set.Erase(3);
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 4);
  }

  // Iterator
  {
    set Set = {2, 3, 4};
    Set.Erase(Set.Begin()+1);
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 2);
    EXPECT_EQ(Values(1), 4);
  }

}

TEST_F(SetTests, EraseIf) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  // Match exists
  {
    set Set = {2, 3, 4, 5};
    Set.EraseIf([](int Value) -> bool { return Value <= 3; });
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 4);
    EXPECT_EQ(Values(1), 5);
  }

  // Match doesn't exist
  {
    set Set = {4, 5};
    Set.EraseIf([](int Value) -> bool { return Value <= 3; });
    auto &Values = helper::GetValues(Set);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_EQ(Values(0), 4);
    EXPECT_EQ(Values(1), 5);
  }

}

TEST_F(SetTests, Clear) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  // Empty
  {
    set Set;
    Set.Clear();
    auto &Values = helper::GetValues(Set);
    EXPECT_TRUE(Values.Empty());
  }

  // Non-empty
  {
    set Set = {2, 3};
    Set.Clear();
    auto &Values = helper::GetValues(Set);
    EXPECT_TRUE(Values.Empty());
  }

}

TEST_F(SetTests, Contains) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;

  // Exists
  {
    set Set = {2, 3, 4};
    EXPECT_TRUE(Set.Contains(3));
  }

  // Doesn't exist
  {
    set Set = {2, 4};
    EXPECT_FALSE(Set.Contains(3));
  }

}

TEST_F(SetTests, Find) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;

  // Exists
  {
    set Set = {2, 3, 4};
    auto Iter = Set.Find(3);
    EXPECT_EQ(Iter, Set.Begin()+1);
  }

  // Doesn't exist
  {
    set Set = {2, 4};
    auto Iter = Set.Find(3);
    EXPECT_EQ(Iter, Set.End());
  }

}

TEST_F(SetTests, LowerBound) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;

  // Exists
  {
    set Set = {2, 3, 4};
    auto Iter = Set.LowerBound(3);
    EXPECT_EQ(Iter, Set.Begin()+1);
  }

  // Doesn't exist
  {
    set Set = {2, 4};
    auto Iter = Set.LowerBound(3);
    EXPECT_EQ(Iter, Set.Begin()+1);
  }

}

TEST_F(SetTests, UpperBound) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;

  // Exists
  {
    set Set = {2, 3, 4};
    auto Iter = Set.UpperBound(3);
    EXPECT_EQ(Iter, Set.Begin()+2);
  }

  // Doesn't exist
  {
    set Set = {2, 4};
    auto Iter = Set.UpperBound(3);
    EXPECT_EQ(Iter, Set.Begin()+1);
  }

}

TEST_F(SetTests, Count) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;

  set Set = {2, 3};
  EXPECT_EQ(Set.Count(), 2);

}

TEST_F(SetTests, Empty) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;

  // Empty
  {
    set Set;
    EXPECT_TRUE(Set.Empty());
  }

  // Not empty
  {
    set Set = {2, 3};
    EXPECT_FALSE(Set.Empty());
  }

}

TEST_F(SetTests, Capacity) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  set Set = {2, 3};
  auto &Values = helper::GetValues(Set);
  EXPECT_EQ(Set.Capacity(), Values.Capacity());

}

TEST_F(SetTests, BracketOperator) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  set Set = {2, 3};
  auto &Values = helper::GetValues(Set);
  EXPECT_EQ(&Set[1], Values.Data(1));

}

TEST_F(SetTests, Data) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  set Set = {2, 3};
  auto &Values = helper::GetValues(Set);
  EXPECT_EQ(Set.Data(), Values.Data());

}

TEST_F(SetTests, BeginEnd) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;
  using helper = ovk::core::test_helper<set>;

  set Set = {2, 3};
  auto &Values = helper::GetValues(Set);
  EXPECT_EQ(Set.Begin().Pointer(), Values.Data());
  EXPECT_EQ(Set.End().Pointer(), Values.Data()+2);

  int Sum = 0;
  for (int Value : Set) Sum += Value;
  EXPECT_EQ(Sum, 5);

}

TEST_F(SetTests, Compare) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int,std::greater<int>>;

  // Get compare instance
  {
    set Set;
    auto &Compare = Set.Compare();
    EXPECT_TRUE((std::is_same<decltype(Compare), const std::greater<int> &>::value));
    EXPECT_TRUE(Compare(3, 2));
  }

  // Do comparison directly
  {
    set Set;
    EXPECT_TRUE(Set.Compare(3, 2));
  }

}

TEST_F(SetTests, ArrayTraits) {

  if (TestComm().Rank() != 0) return;

  using set = ovk::set<int>;

  EXPECT_TRUE(ovk::core::IsArray<set>());
  EXPECT_TRUE((std::is_same<ovk::core::array_value_type<set>, typename set::value_type>::value));
  EXPECT_EQ(ovk::core::ArrayRank<set>(), 1);
  EXPECT_TRUE(ovk::core::ArrayHasRuntimeExtents<set>());

  set Set = {2, 3};
  EXPECT_THAT(ovk::core::ArrayExtents(Set).Begin(), ElementsAre(0));
  EXPECT_THAT(ovk::core::ArrayExtents(Set).End(), ElementsAre(2));
  EXPECT_EQ(ovk::core::ArrayData(Set), Set.Data());

}
