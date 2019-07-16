// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/IDSet.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Elem.hpp>

#include <array>
#include <mpi.h>

#include <utility>

class IDSetTests : public tests::mpi_test {};

using testing::ElementsAre;

namespace ovk {
namespace core {
template <int Rank> class test_helper<id_set<Rank>> {
public:
  using id_set_type = id_set<Rank>;
  using value_type = typename id_set<Rank>::value_type;
  static const array<value_type> &GetValues(const id_set_type &IDSet) { return IDSet.Values_; }
  static array<value_type> &GetValues(id_set_type &IDSet) { return IDSet.Values_; }
};
}}

TEST_F(IDSetTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;

  EXPECT_EQ(int(id_set::Rank), 2);
  EXPECT_TRUE((std::is_same<typename id_set::value_type, ovk::elem<int,2>>::value));
  EXPECT_TRUE((std::is_same<typename id_set::iterator::pointer, const ovk::elem<int,2> *>::value));
  EXPECT_TRUE((std::is_same<typename id_set::const_iterator::pointer, const ovk::elem<int,2>
    *>::value));

}

TEST_F(IDSetTests, Create) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  // Default
  {
    id_set IDSet;
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 0);
  }

  // Initializer list
  {
    id_set IDSet = {{1,2}, {1,3}};
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(1,3));
  }

  // Iterators
  {
    std::array<ovk::elem<int,2>,2> SourceValues = {{{1,2}, {1,3}}};
    id_set IDSet(SourceValues.begin(), SourceValues.end());
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(1,3));
  }

}

TEST_F(IDSetTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  // Copy construct
  {
    id_set IDSet1 = {{1,2}, {1,3}};
    id_set IDSet2 = IDSet1;
    auto &Values1 = helper::GetValues(IDSet1);
    auto &Values2 = helper::GetValues(IDSet2);
    EXPECT_EQ(Values1.Count(), 2);
    EXPECT_THAT(Values1(0), ElementsAre(1,2));
    EXPECT_THAT(Values1(1), ElementsAre(1,3));
    EXPECT_EQ(Values2.Count(), 2);
    EXPECT_THAT(Values2(0), ElementsAre(1,2));
    EXPECT_THAT(Values2(1), ElementsAre(1,3));
  }

  // Copy assign
  {
    id_set IDSet1 = {{1,2}, {1,3}};
    id_set IDSet2;
    IDSet2 = IDSet1;
    auto &Values1 = helper::GetValues(IDSet1);
    auto &Values2 = helper::GetValues(IDSet2);
    EXPECT_EQ(Values1.Count(), 2);
    EXPECT_THAT(Values1(0), ElementsAre(1,2));
    EXPECT_THAT(Values1(1), ElementsAre(1,3));
    EXPECT_EQ(Values2.Count(), 2);
    EXPECT_THAT(Values2(0), ElementsAre(1,2));
    EXPECT_THAT(Values2(1), ElementsAre(1,3));
  }

}

TEST_F(IDSetTests, Move) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  // Move construct
  {
    id_set IDSet1 = {{1,2}, {1,3}};
    id_set IDSet2(std::move(IDSet1));
    auto &Values1 = helper::GetValues(IDSet1);
    auto &Values2 = helper::GetValues(IDSet2);
    EXPECT_EQ(Values1.Count(), 0);
    EXPECT_EQ(Values2.Count(), 2);
    EXPECT_THAT(Values2(0), ElementsAre(1,2));
    EXPECT_THAT(Values2(1), ElementsAre(1,3));
  }

  // Move assign
  {
    id_set IDSet1 = {{1,2}, {1,3}};
    id_set IDSet2;
    IDSet2 = std::move(IDSet1);
    auto &Values1 = helper::GetValues(IDSet1);
    auto &Values2 = helper::GetValues(IDSet2);
    EXPECT_EQ(Values1.Count(), 0);
    EXPECT_EQ(Values2.Count(), 2);
    EXPECT_THAT(Values2(0), ElementsAre(1,2));
    EXPECT_THAT(Values2(1), ElementsAre(1,3));
  }

}

TEST_F(IDSetTests, Assign) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  // Initializer list, operator=
  {
    id_set IDSet = {{0,1}, {0,2}};
    IDSet = {{1,2}, {1,3}};
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(1,3));
  }

  // Initializer list, Assign
  {
    id_set IDSet = {{0,1}, {0,2}};
    IDSet.Assign({{1,2}, {1,3}});
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(1,3));
  }

  // Iterators, Assign
  {
    id_set IDSet = {{0,1}, {0,2}};
    std::array<ovk::elem<int,2>,2> SourceValues = {{{1,2}, {1,3}}};
    IDSet.Assign(SourceValues.begin(), SourceValues.end());
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(1,3));
  }

}

TEST_F(IDSetTests, Reserve) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  id_set IDSet;
  auto &Values = helper::GetValues(IDSet);
  long long Capacity = Values.Capacity();
  IDSet.Reserve(Capacity+1);
  EXPECT_GE(Values.Capacity(), Capacity+1);

}

TEST_F(IDSetTests, Insert) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  // Value, doesn't already exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    IDSet.Insert({1,3});
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 3);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(1,3));
    EXPECT_THAT(Values(2), ElementsAre(2,1));
  }

  // Value, already exists
  {
    id_set IDSet = {{1,2}, {2,1}};
    IDSet.Insert({2,1});
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(2,1));
  }

  // Separate IDs, doesn't already exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    IDSet.Insert(1,3);
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 3);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(1,3));
    EXPECT_THAT(Values(2), ElementsAre(2,1));
  }

  // Separate IDs, already exists
  {
    id_set IDSet = {{1,2}, {2,1}};
    IDSet.Insert(2,1);
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(2,1));
  }

  // With lower bound iterator, doesn't already exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    auto &Values = helper::GetValues(IDSet);
    IDSet.Insert(IDSet.Begin()+1, {1,3});
    EXPECT_EQ(Values.Count(), 3);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(1,3));
    EXPECT_THAT(Values(2), ElementsAre(2,1));
  }

  // With lower bound iterator, already exists
  {
    id_set IDSet = {{1,2}, {2,1}};
    auto &Values = helper::GetValues(IDSet);
    IDSet.Insert(IDSet.Begin()+1, {2,1});
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(2,1));
  }

}

TEST_F(IDSetTests, Erase) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  // Value, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    IDSet.Erase({1,3});
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(2,1));
  }

  // Value, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    IDSet.Erase({1,3});
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(2,1));
  }

  // Separate IDs, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    IDSet.Erase(1,3);
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(2,1));
  }

  // Separate IDs, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    IDSet.Erase(1,3);
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(2,1));
  }

  // Iterator
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    IDSet.Erase(IDSet.Begin()+1);
    auto &Values = helper::GetValues(IDSet);
    EXPECT_EQ(Values.Count(), 2);
    EXPECT_THAT(Values(0), ElementsAre(1,2));
    EXPECT_THAT(Values(1), ElementsAre(2,1));
  }

}

TEST_F(IDSetTests, Clear) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  // Empty
  {
    id_set IDSet;
    IDSet.Clear();
    auto &Values = helper::GetValues(IDSet);
    EXPECT_TRUE(Values.Empty());
  }

  // Non-empty
  {
    id_set IDSet = {{1,2}, {1,3}};
    IDSet.Clear();
    auto &Values = helper::GetValues(IDSet);
    EXPECT_TRUE(Values.Empty());
  }

}

TEST_F(IDSetTests, Contains) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;

  // Value, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    EXPECT_TRUE(IDSet.Contains({1,3}));
  }

  // Value, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    EXPECT_FALSE(IDSet.Contains({1,3}));
  }

  // Separate IDs, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    EXPECT_TRUE(IDSet.Contains(1,3));
  }

  // Separate IDs, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    EXPECT_FALSE(IDSet.Contains(1,3));
  }

}

TEST_F(IDSetTests, Find) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;

  // Value, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    auto Iter = IDSet.Find({1,3});
    EXPECT_EQ(Iter, IDSet.Begin()+1);
  }

  // Value, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    auto Iter = IDSet.Find({1,3});
    EXPECT_EQ(Iter, IDSet.End());
  }

  // Separate IDs, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    auto Iter = IDSet.Find(1,3);
    EXPECT_EQ(Iter, IDSet.Begin()+1);
  }

  // Separate IDs, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    auto Iter = IDSet.Find(1,3);
    EXPECT_EQ(Iter, IDSet.End());
  }

}

TEST_F(IDSetTests, LowerBound) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;

  // Value, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    auto Iter = IDSet.LowerBound({1,3});
    EXPECT_EQ(Iter, IDSet.Begin()+1);
  }

  // Value, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    auto Iter = IDSet.LowerBound({1,3});
    EXPECT_EQ(Iter, IDSet.Begin()+1);
  }

  // Separate IDs, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    auto Iter = IDSet.LowerBound(1,3);
    EXPECT_EQ(Iter, IDSet.Begin()+1);
  }

  // Separate IDs, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    auto Iter = IDSet.LowerBound(1,3);
    EXPECT_EQ(Iter, IDSet.Begin()+1);
  }

}

TEST_F(IDSetTests, UpperBound) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;

  // Value, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    auto Iter = IDSet.UpperBound({1,3});
    EXPECT_EQ(Iter, IDSet.Begin()+2);
  }

  // Value, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    auto Iter = IDSet.UpperBound({1,3});
    EXPECT_EQ(Iter, IDSet.Begin()+1);
  }

  // Separate IDs, exists
  {
    id_set IDSet = {{1,2}, {1,3}, {2,1}};
    auto Iter = IDSet.UpperBound(1,3);
    EXPECT_EQ(Iter, IDSet.Begin()+2);
  }

  // Separate IDs, doesn't exist
  {
    id_set IDSet = {{1,2}, {2,1}};
    auto Iter = IDSet.UpperBound(1,3);
    EXPECT_EQ(Iter, IDSet.Begin()+1);
  }

}

TEST_F(IDSetTests, NextAvailableValue) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<1>;

  // First ID > 0
  {
    id_set IDSet = {1, 2};
    EXPECT_EQ(IDSet.NextAvailableValue(), 0);
  }

  // Gap
  {
    id_set IDSet = {0, 2};
    EXPECT_EQ(IDSet.NextAvailableValue(), 1);
  }

  // End
  {
    id_set IDSet = {0, 1};
    EXPECT_EQ(IDSet.NextAvailableValue(), 2);
  }

}

TEST_F(IDSetTests, Count) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;

  id_set IDSet = {{1,2}, {1,3}};
  EXPECT_EQ(IDSet.Count(), 2);

}

TEST_F(IDSetTests, Empty) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;

  // Empty
  {
    id_set IDSet;
    EXPECT_TRUE(IDSet.Empty());
  }

  // Not empty
  {
    id_set IDSet = {{1,2}, {1,3}};
    EXPECT_FALSE(IDSet.Empty());
  }

}

TEST_F(IDSetTests, Capacity) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  id_set IDSet = {{1,2}, {1,3}};
  auto &Values = helper::GetValues(IDSet);
  EXPECT_EQ(IDSet.Capacity(), Values.Capacity());

}

TEST_F(IDSetTests, BracketOperator) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  id_set IDSet = {{1,2}, {1,3}};
  auto &Values = helper::GetValues(IDSet);
  EXPECT_EQ(&IDSet[1], Values.Data(1));

}

TEST_F(IDSetTests, Data) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  id_set IDSet = {{1,2}, {1,3}};
  auto &Values = helper::GetValues(IDSet);
  EXPECT_EQ(IDSet.Data(), Values.Data());

}

TEST_F(IDSetTests, BeginEnd) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;
  using helper = ovk::core::test_helper<id_set>;

  id_set IDSet = {{1,2}, {1,3}};
  auto &Values = helper::GetValues(IDSet);
  EXPECT_EQ(IDSet.Begin().Pointer(), Values.Data());
  EXPECT_EQ(IDSet.End().Pointer(), Values.Data()+2);

  int Sum = 0;
  for (auto &Value : IDSet) Sum += Value(0)+Value(1);
  EXPECT_EQ(Sum, 7);

}

TEST_F(IDSetTests, Less) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;

  // Equal
  {
    EXPECT_FALSE(id_set::Less({1,2}, {1,2}));
  }

  // Same first index, second index less
  {
    EXPECT_TRUE(id_set::Less({1,1}, {1,2}));
  }

  // Same first index, second index greater
  {
    EXPECT_FALSE(id_set::Less({1,3}, {1,2}));
  }

  // First index less
  {
    EXPECT_TRUE(id_set::Less({0,2}, {1,2}));
  }

  // First index greater
  {
    EXPECT_FALSE(id_set::Less({2,2}, {1,2}));
  }

}

TEST_F(IDSetTests, ArrayTraits) {

  if (TestComm().Rank() != 0) return;

  using id_set = ovk::id_set<2>;

  EXPECT_TRUE(ovk::core::IsArray<id_set>());
  EXPECT_TRUE((std::is_same<ovk::core::array_value_type<id_set>, typename id_set::value_type>::
    value));
  EXPECT_EQ(ovk::core::ArrayRank<id_set>(), 1);
  EXPECT_TRUE(ovk::core::ArrayHasRuntimeExtents<id_set>());

  id_set IDSet = {{1,2}, {1,3}};
  EXPECT_THAT(ovk::core::ArrayExtents(IDSet).Begin(), ElementsAre(0));
  EXPECT_THAT(ovk::core::ArrayExtents(IDSet).End(), ElementsAre(2));
  EXPECT_EQ(ovk::core::ArrayData(IDSet), IDSet.Data());

}
