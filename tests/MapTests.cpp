// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Map.hpp>

#include "tests/MPITest.hpp"
#include "tests/mocks/Noncopyable.hpp"
#include "tests/mocks/Nondefaultconstructible.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <mpi.h>

#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>

// #include <ovk/core/Profiler.hpp>
// #include <cstdlib>
// #include <ctime>
// #include <map>

class MapEntryTests : public tests::mpi_test {};
class MapTests : public tests::mpi_test {};

using testing::ElementsAre;

using tests::noncopyable;
using tests::nondefaultconstructible;

namespace ovk {
namespace core {
template <typename KeyType, typename ValueType, bool Contiguous> class test_helper<map_entry<
  KeyType, ValueType, Contiguous>> {
public:
  using entry_type = map_entry<KeyType, ValueType, Contiguous>;
  using key_type = typename entry_type::key_type;
  using value_wrapper_type = typename entry_type::value_wrapper_type;
  static const key_type &GetKey(const entry_type &Entry) { return Entry.Key_; }
  static const value_wrapper_type &GetValue(const entry_type &Entry) { return Entry.Value_; }
  static value_wrapper_type &GetValue(entry_type &Entry) { return Entry.Value_; }
};
template <typename KeyType, typename ValueType, typename KeyCompareType, bool Contiguous> class
  test_helper<map<KeyType,ValueType,KeyCompareType,Contiguous>> {
public:
  using map_type = map<KeyType, ValueType, KeyCompareType, Contiguous>;
  using key_set_type = typename map_type::key_set_type;
  using entry = typename map_type::entry;
  static const key_set_type &GetKeys(const map_type &Map) { return Map.Keys_; }
  static key_set_type &GetKeys(map_type &Map) { return Map.Keys_; }
  static const array<entry> &GetEntries(const map_type &Map) { return Map.Entries_; }
  static array<entry> &GetEntries(map_type &Map) { return Map.Entries_; }
};
}}

namespace {

struct multiargument {
  int v1, v2;
  multiargument() = default;
  multiargument(int v1_, int v2_):
    v1(v1_),
    v2(v2_)
  {}
};

}

TEST_F(MapEntryTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using entry = ovk::map_entry<int,double>;

  EXPECT_TRUE((std::is_same<typename entry::key_type, int>::value));
  EXPECT_TRUE((std::is_same<typename entry::value_type, double>::value));
  EXPECT_TRUE(bool(entry::Contiguous));

  using entry_noncontig = ovk::map_entry_noncontig<int,double>;

  EXPECT_FALSE(bool(entry_noncontig::Contiguous));

}

TEST_F(MapEntryTests, Create) {

  if (TestComm().Rank() != 0) return;

  using entry_contig = ovk::map_entry_contig<int,int>;
  using entry_noncontig = ovk::map_entry_noncontig<int,int>;
  using entry_noncopyable_contig = ovk::map_entry_contig<int,noncopyable<int>>;
  using entry_noncopyable_noncontig = ovk::map_entry_noncontig<int,noncopyable<int>>;
  using entry_multiargument_contig = ovk::map_entry_contig<int,multiargument>;
  using entry_multiargument_noncontig = ovk::map_entry_noncontig<int,multiargument>;
  using entry_nondefaultconstructible_contig = ovk::map_entry_contig<int,nondefaultconstructible<int>>;
  using entry_nondefaultconstructible_noncontig = ovk::map_entry_noncontig<int,nondefaultconstructible<int>>;

  using helper_contig = ovk::core::test_helper<entry_contig>;
  using helper_noncontig = ovk::core::test_helper<entry_noncontig>;
  using helper_noncopyable_contig = ovk::core::test_helper<entry_noncopyable_contig>;
  using helper_noncopyable_noncontig = ovk::core::test_helper<entry_noncopyable_noncontig>;
  using helper_multiargument_contig = ovk::core::test_helper<entry_multiargument_contig>;
  using helper_multiargument_noncontig = ovk::core::test_helper<entry_multiargument_noncontig>;
  using helper_nondefaultconstructible_contig = ovk::core::test_helper<
    entry_nondefaultconstructible_contig>;
  using helper_nondefaultconstructible_noncontig = ovk::core::test_helper<
    entry_nondefaultconstructible_noncontig>;

  // lvalue ref, contiguous
  {
    int SourceValue = 1;
    entry_contig Entry(2, SourceValue);
    auto &Key = helper_contig::GetKey(Entry);
    auto &Value = helper_contig::GetValue(Entry);
    EXPECT_EQ(Key, 2);
    EXPECT_EQ(Value.Get(), 1);
  }

  // lvalue ref, non-contiguous
  {
    int SourceValue = 1;
    entry_noncontig Entry(2, SourceValue);
    auto &Key = helper_noncontig::GetKey(Entry);
    auto &Value = helper_noncontig::GetValue(Entry);
    EXPECT_EQ(Key, 2);
    EXPECT_EQ(Value.Get(), 1);
  }

  // rvalue ref, contiguous
  {
    noncopyable<int> SourceValue(1);
    entry_noncopyable_contig Entry(1, std::move(SourceValue));
    auto &Key = helper_noncopyable_contig::GetKey(Entry);
    auto &Value = helper_noncopyable_contig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // rvalue ref, non-contiguous
  {
    noncopyable<int> SourceValue(1);
    entry_noncopyable_noncontig Entry(1, std::move(SourceValue));
    auto &Key = helper_noncopyable_noncontig::GetKey(Entry);
    auto &Value = helper_noncopyable_noncontig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // In-place, contiguous
  {
    entry_multiargument_contig Entry(1, 1, 2);
    auto &Key = helper_multiargument_contig::GetKey(Entry);
    auto &Value = helper_multiargument_contig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().v1, 1);
    EXPECT_EQ(Value.Get().v2, 2);
  }

  // In-place, non-contiguous
  {
    entry_multiargument_noncontig Entry(1, 1, 2);
    auto &Key = helper_multiargument_noncontig::GetKey(Entry);
    auto &Value = helper_multiargument_noncontig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().v1, 1);
    EXPECT_EQ(Value.Get().v2, 2);
  }

  // Non-default-constructible, lvalue ref, contiguous
  {
    nondefaultconstructible<int> SourceValue(1);
    entry_nondefaultconstructible_contig Entry(1, SourceValue);
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, lvalue ref, non-contiguous
  {
    nondefaultconstructible<int> SourceValue(1);
    entry_nondefaultconstructible_noncontig Entry(1, SourceValue);
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, rvalue ref, contiguous
  {
    nondefaultconstructible<int> SourceValue(1);
    entry_nondefaultconstructible_contig Entry(1, std::move(SourceValue));
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, rvalue ref, non-contiguous
  {
    nondefaultconstructible<int> SourceValue(1);
    entry_nondefaultconstructible_noncontig Entry(1, std::move(SourceValue));
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, in-place, contiguous
  {
    entry_nondefaultconstructible_contig Entry(1, 1);
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, in-place, non-contiguous
  {
    entry_nondefaultconstructible_noncontig Entry(1, 1);
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

}

TEST_F(MapEntryTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using entry_contig = ovk::map_entry_contig<int,int>;
  using entry_noncontig = ovk::map_entry_noncontig<int,int>;
  using entry_nondefaultconstructible_contig = ovk::map_entry_contig<int,nondefaultconstructible<int>>;
  using entry_nondefaultconstructible_noncontig = ovk::map_entry_noncontig<int,nondefaultconstructible<int>>;

  using helper_contig = ovk::core::test_helper<entry_contig>;
  using helper_noncontig = ovk::core::test_helper<entry_noncontig>;
  using helper_nondefaultconstructible_contig = ovk::core::test_helper<
    entry_nondefaultconstructible_contig>;
  using helper_nondefaultconstructible_noncontig = ovk::core::test_helper<
    entry_nondefaultconstructible_noncontig>;

  // Contiguous
  {
    entry_contig Entry1(2, 1);
    entry_contig Entry2(Entry1);
    auto &Key = helper_contig::GetKey(Entry2);
    auto &Value = helper_contig::GetValue(Entry2);
    EXPECT_EQ(Key, 2);
    EXPECT_EQ(Value.Get(), 1);
  }

  // Non-contiguous
  {
    entry_noncontig Entry1(2, 1);
    entry_noncontig Entry2(Entry1);
    auto &Key = helper_noncontig::GetKey(Entry2);
    auto &Value = helper_noncontig::GetValue(Entry2);
    EXPECT_EQ(Key, 2);
    EXPECT_EQ(Value.Get(), 1);
  }

  // Non-default-constructible, contiguous
  {
    entry_nondefaultconstructible_contig Entry1(1, 1);
    entry_nondefaultconstructible_contig Entry2(Entry1);
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry2);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry2);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, non-contiguous
  {
    entry_nondefaultconstructible_noncontig Entry1(1, 1);
    entry_nondefaultconstructible_noncontig Entry2(Entry1);
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry2);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry2);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

}

TEST_F(MapEntryTests, Move) {

  if (TestComm().Rank() != 0) return;

  using entry_noncopyable_contig = ovk::map_entry_contig<int,noncopyable<int>>;
  using entry_noncopyable_noncontig = ovk::map_entry_noncontig<int,noncopyable<int>>;
  using entry_nondefaultconstructible_contig = ovk::map_entry_contig<int,nondefaultconstructible<int>>;
  using entry_nondefaultconstructible_noncontig = ovk::map_entry_noncontig<int,nondefaultconstructible<int>>;

  using helper_noncopyable_contig = ovk::core::test_helper<entry_noncopyable_contig>;
  using helper_noncopyable_noncontig = ovk::core::test_helper<entry_noncopyable_noncontig>;
  using helper_nondefaultconstructible_contig = ovk::core::test_helper<
    entry_nondefaultconstructible_contig>;
  using helper_nondefaultconstructible_noncontig = ovk::core::test_helper<
    entry_nondefaultconstructible_noncontig>;

  // Contiguous
  {
    entry_noncopyable_contig Entry1(1, 1);
    entry_noncopyable_contig Entry2(std::move(Entry1));
    auto &Key = helper_noncopyable_contig::GetKey(Entry2);
    auto &Value = helper_noncopyable_contig::GetValue(Entry2);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-contiguous
  {
    entry_noncopyable_noncontig Entry1(1, 1);
    entry_noncopyable_noncontig Entry2(std::move(Entry1));
    auto &Key = helper_noncopyable_noncontig::GetKey(Entry2);
    auto &Value = helper_noncopyable_noncontig::GetValue(Entry2);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, contiguous
  {
    entry_nondefaultconstructible_contig Entry1(1, 1);
    entry_nondefaultconstructible_contig Entry2(std::move(Entry1));
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry2);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry2);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, non-contiguous
  {
    entry_nondefaultconstructible_noncontig Entry1(1, 1);
    entry_nondefaultconstructible_noncontig Entry2(std::move(Entry1));
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry2);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry2);
    EXPECT_EQ(Key, 1);
    EXPECT_EQ(Value.Get().Value(), 1);
  }

}

TEST_F(MapEntryTests, Key) {

  if (TestComm().Rank() != 0) return;

  using entry_contig = ovk::map_entry_contig<int,int>;
  using entry_noncontig = ovk::map_entry_noncontig<int,int>;

  using helper_contig = ovk::core::test_helper<entry_contig>;
  using helper_noncontig = ovk::core::test_helper<entry_noncontig>;

  // Contiguous
  {
    entry_contig Entry(1, 1);
    auto &Key = helper_contig::GetKey(Entry);
    EXPECT_EQ(&Entry.Key(), &Key);
  }

  // Non-contiguous
  {
    entry_noncontig Entry(1, 1);
    auto &Key = helper_noncontig::GetKey(Entry);
    EXPECT_EQ(&Entry.Key(), &Key);
  }

}

TEST_F(MapEntryTests, Value) {

  if (TestComm().Rank() != 0) return;

  using entry_contig = ovk::map_entry_contig<int,int>;
  using entry_noncontig = ovk::map_entry_noncontig<int,int>;

  using helper_contig = ovk::core::test_helper<entry_contig>;
  using helper_noncontig = ovk::core::test_helper<entry_noncontig>;

  // Contiguous
  {
    entry_contig Entry(1, 1);
    auto &Value = helper_contig::GetValue(Entry);
    EXPECT_EQ(&Entry.Value(), &Value.Get());
  }

  // Non-contiguous
  {
    entry_noncontig Entry(1, 1);
    auto &Value = helper_noncontig::GetValue(Entry);
    EXPECT_EQ(&Entry.Value(), &Value.Get());
  }

}

TEST_F(MapTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map_contig<int,double,std::greater<int>>;

  EXPECT_TRUE((std::is_same<typename map::key_type, int>::value));
  EXPECT_TRUE((std::is_same<typename map::value_type, double>::value));
  EXPECT_TRUE((std::is_same<typename map::key_compare_type, std::greater<int>>::value));
  EXPECT_TRUE((std::is_same<typename map::key_set_type, ovk::set<int,std::greater<int>>>::value));
  EXPECT_TRUE(bool(map::Contiguous));
  EXPECT_TRUE((std::is_same<typename map::index_type, long long>::value));
  EXPECT_TRUE((std::is_same<typename map::entry, ovk::map_entry_contig<int,double>>::value));
  EXPECT_TRUE((std::is_same<typename map::iterator::pointer, typename map::entry *>::value));
  EXPECT_TRUE((std::is_same<typename map::const_iterator::pointer, const typename map::entry *
    >::value));

  using map_noncontig = ovk::map_noncontig<int,double,std::greater<int>>;

  EXPECT_FALSE(bool(map_noncontig::Contiguous));
  EXPECT_TRUE((std::is_same<typename map_noncontig::entry, ovk::map_entry_noncontig<int,
    double>>::value));

}

TEST_F(MapTests, Create) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using map_nondefaultconstructible = ovk::map<int,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<map>;
  using helper_nondefaultconstructible = ovk::core::test_helper<map_nondefaultconstructible>;

  // Default
  {
    map Map;
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 0);
    EXPECT_EQ(Entries.Count(), 0);
  }

  // Initializer list
  {
    map Map = {{2, 1}, {3, 2}};
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Iterators
  {
    std::array<map::entry,2> SourceEntries = {{{2, 1}, {3, 2}}};
    map Map(SourceEntries.begin(), SourceEntries.end());
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Non-default-constructible, default
  {
    map_nondefaultconstructible Map;
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 0);
    EXPECT_EQ(Entries.Count(), 0);
  }

  // Non-default-constructible, initializer list
  {
    map_nondefaultconstructible Map = {{1, {1}}, {2, {2}}};
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

  // Non-default-constructible, iterators
  {
    std::array<map_nondefaultconstructible::entry,2> SourceEntries = {{{1, {1}}, {2, {2}}}};
    map_nondefaultconstructible Map(SourceEntries.begin(), SourceEntries.end());
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

}

TEST_F(MapTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using map_nondefaultconstructible = ovk::map<int,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<map>;
  using helper_nondefaultconstructible = ovk::core::test_helper<map_nondefaultconstructible>;

  // Copy construct
  {
    map Map1 = {{2, 1}, {3, 2}};
    map Map2 = Map1;
    auto &Keys1 = helper::GetKeys(Map1);
    auto &Entries1 = helper::GetEntries(Map1);
    auto &Keys2 = helper::GetKeys(Map2);
    auto &Entries2 = helper::GetEntries(Map2);
    EXPECT_EQ(Keys1.Count(), 2);
    EXPECT_EQ(Keys1[0], 2);
    EXPECT_EQ(Keys1[1], 3);
    EXPECT_EQ(Entries1.Count(), 2);
    EXPECT_EQ(Entries1(0).Key(), 2);
    EXPECT_EQ(Entries1(0).Value(), 1);
    EXPECT_EQ(Entries1(1).Key(), 3);
    EXPECT_EQ(Entries1(1).Value(), 2);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_EQ(Keys2[0], 2);
    EXPECT_EQ(Keys2[1], 3);
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_EQ(Entries2(0).Key(), 2);
    EXPECT_EQ(Entries2(0).Value(), 1);
    EXPECT_EQ(Entries2(1).Key(), 3);
    EXPECT_EQ(Entries2(1).Value(), 2);
  }

  // Copy assign
  {
    map Map1 = {{2, 1}, {3, 2}};
    map Map2;
    Map2 = Map1;
    auto &Keys1 = helper::GetKeys(Map1);
    auto &Entries1 = helper::GetEntries(Map1);
    auto &Keys2 = helper::GetKeys(Map2);
    auto &Entries2 = helper::GetEntries(Map2);
    EXPECT_EQ(Keys1.Count(), 2);
    EXPECT_EQ(Keys1[0], 2);
    EXPECT_EQ(Keys1[1], 3);
    EXPECT_EQ(Entries1.Count(), 2);
    EXPECT_EQ(Entries1(0).Key(), 2);
    EXPECT_EQ(Entries1(0).Value(), 1);
    EXPECT_EQ(Entries1(1).Key(), 3);
    EXPECT_EQ(Entries1(1).Value(), 2);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_EQ(Keys2[0], 2);
    EXPECT_EQ(Keys2[1], 3);
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_EQ(Entries2(0).Key(), 2);
    EXPECT_EQ(Entries2(0).Value(), 1);
    EXPECT_EQ(Entries2(1).Key(), 3);
    EXPECT_EQ(Entries2(1).Value(), 2);
  }

  // Non-default-constructible, copy construct
  {
    map_nondefaultconstructible Map1 = {{2, {1}}, {3, {2}}};
    map_nondefaultconstructible Map2 = Map1;
    auto &Keys1 = helper_nondefaultconstructible::GetKeys(Map1);
    auto &Entries1 = helper_nondefaultconstructible::GetEntries(Map1);
    auto &Keys2 = helper_nondefaultconstructible::GetKeys(Map2);
    auto &Entries2 = helper_nondefaultconstructible::GetEntries(Map2);
    EXPECT_EQ(Keys1.Count(), 2);
    EXPECT_EQ(Keys1[0], 2);
    EXPECT_EQ(Keys1[1], 3);
    EXPECT_EQ(Entries1.Count(), 2);
    EXPECT_EQ(Entries1(0).Key(), 2);
    EXPECT_EQ(Entries1(0).Value().Value(), 1);
    EXPECT_EQ(Entries1(1).Key(), 3);
    EXPECT_EQ(Entries1(1).Value().Value(), 2);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_EQ(Keys2[0], 2);
    EXPECT_EQ(Keys2[1], 3);
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_EQ(Entries2(0).Key(), 2);
    EXPECT_EQ(Entries2(0).Value().Value(), 1);
    EXPECT_EQ(Entries2(1).Key(), 3);
    EXPECT_EQ(Entries2(1).Value().Value(), 2);
  }

  // Non-default-constructible, copy assign
  {
    map_nondefaultconstructible Map1 = {{2, {1}}, {3, {2}}};
    map_nondefaultconstructible Map2;
    Map2 = Map1;
    auto &Keys1 = helper_nondefaultconstructible::GetKeys(Map1);
    auto &Entries1 = helper_nondefaultconstructible::GetEntries(Map1);
    auto &Keys2 = helper_nondefaultconstructible::GetKeys(Map2);
    auto &Entries2 = helper_nondefaultconstructible::GetEntries(Map2);
    EXPECT_EQ(Keys1.Count(), 2);
    EXPECT_EQ(Keys1[0], 2);
    EXPECT_EQ(Keys1[1], 3);
    EXPECT_EQ(Entries1.Count(), 2);
    EXPECT_EQ(Entries1(0).Key(), 2);
    EXPECT_EQ(Entries1(0).Value().Value(), 1);
    EXPECT_EQ(Entries1(1).Key(), 3);
    EXPECT_EQ(Entries1(1).Value().Value(), 2);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_EQ(Keys2[0], 2);
    EXPECT_EQ(Keys2[1], 3);
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_EQ(Entries2(0).Key(), 2);
    EXPECT_EQ(Entries2(0).Value().Value(), 1);
    EXPECT_EQ(Entries2(1).Key(), 3);
    EXPECT_EQ(Entries2(1).Value().Value(), 2);
  }

}

TEST_F(MapTests, Move) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using map_nondefaultconstructible = ovk::map<int,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<map>;
  using helper_nondefaultconstructible = ovk::core::test_helper<map_nondefaultconstructible>;

  // Move construct
  {
    map Map1 = {{2, 1}, {3, 2}};
    map Map2 = std::move(Map1);
    auto &Keys1 = helper::GetKeys(Map1);
    auto &Entries1 = helper::GetEntries(Map1);
    auto &Keys2 = helper::GetKeys(Map2);
    auto &Entries2 = helper::GetEntries(Map2);
    EXPECT_EQ(Keys1.Count(), 0);
    EXPECT_EQ(Entries1.Count(), 0);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_EQ(Keys2[0], 2);
    EXPECT_EQ(Keys2[1], 3);
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_EQ(Entries2(0).Key(), 2);
    EXPECT_EQ(Entries2(0).Value(), 1);
    EXPECT_EQ(Entries2(1).Key(), 3);
    EXPECT_EQ(Entries2(1).Value(), 2);
  }

  // Move assign
  {
    map Map1 = {{2, 1}, {3, 2}};
    map Map2;
    Map2 = std::move(Map1);
    auto &Keys1 = helper::GetKeys(Map1);
    auto &Entries1 = helper::GetEntries(Map1);
    auto &Keys2 = helper::GetKeys(Map2);
    auto &Entries2 = helper::GetEntries(Map2);
    EXPECT_EQ(Keys1.Count(), 0);
    EXPECT_EQ(Entries1.Count(), 0);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_EQ(Keys2[0], 2);
    EXPECT_EQ(Keys2[1], 3);
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_EQ(Entries2(0).Key(), 2);
    EXPECT_EQ(Entries2(0).Value(), 1);
    EXPECT_EQ(Entries2(1).Key(), 3);
    EXPECT_EQ(Entries2(1).Value(), 2);
  }

  // Non-default-constructible, move construct
  {
    map_nondefaultconstructible Map1 = {{2, {1}}, {3, {2}}};
    map_nondefaultconstructible Map2 = std::move(Map1);
    auto &Keys1 = helper_nondefaultconstructible::GetKeys(Map1);
    auto &Entries1 = helper_nondefaultconstructible::GetEntries(Map1);
    auto &Keys2 = helper_nondefaultconstructible::GetKeys(Map2);
    auto &Entries2 = helper_nondefaultconstructible::GetEntries(Map2);
    EXPECT_EQ(Keys1.Count(), 0);
    EXPECT_EQ(Entries1.Count(), 0);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_EQ(Keys2[0], 2);
    EXPECT_EQ(Keys2[1], 3);
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_EQ(Entries2(0).Key(), 2);
    EXPECT_EQ(Entries2(0).Value().Value(), 1);
    EXPECT_EQ(Entries2(1).Key(), 3);
    EXPECT_EQ(Entries2(1).Value().Value(), 2);
  }

  // Non-default-constructible, move assign
  {
    map_nondefaultconstructible Map1 = {{2, {1}}, {3, {2}}};
    map_nondefaultconstructible Map2;
    Map2 = std::move(Map1);
    auto &Keys1 = helper_nondefaultconstructible::GetKeys(Map1);
    auto &Entries1 = helper_nondefaultconstructible::GetEntries(Map1);
    auto &Keys2 = helper_nondefaultconstructible::GetKeys(Map2);
    auto &Entries2 = helper_nondefaultconstructible::GetEntries(Map2);
    EXPECT_EQ(Keys1.Count(), 0);
    EXPECT_EQ(Entries1.Count(), 0);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_EQ(Keys2[0], 2);
    EXPECT_EQ(Keys2[1], 3);
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_EQ(Entries2(0).Key(), 2);
    EXPECT_EQ(Entries2(0).Value().Value(), 1);
    EXPECT_EQ(Entries2(1).Key(), 3);
    EXPECT_EQ(Entries2(1).Value().Value(), 2);
  }

}

TEST_F(MapTests, Assign) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using map_nondefaultconstructible = ovk::map<int,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<map>;
  using helper_nondefaultconstructible = ovk::core::test_helper<map_nondefaultconstructible>;

  // Initializer list, operator=
  {
    map Map = {{1, 0}, {2, 0}};
    Map = {{2, 1}, {3, 2}};
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Initializer list, Assign
  {
    map Map = {{1, 0}, {2, 0}};
    Map.Assign({{2, 1}, {3, 2}});
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Iterators, Assign
  {
    map Map = {{1, 0}, {2, 0}};
    std::array<map::entry,2> SourceEntries = {{{2, 1}, {3, 2}}};
    Map.Assign(SourceEntries.begin(), SourceEntries.end());
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Non-default-constructible, initializer list, operator=
  {
    map_nondefaultconstructible Map = {{0, {0}}, {1, {0}}};
    Map = {{1, {1}}, {2, {2}}};
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

  // Non-default-constructible, initializer list, Assign
  {
    map_nondefaultconstructible Map = {{0, {0}}, {1, {0}}};
    Map.Assign({{1, {1}}, {2, {2}}});
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

  // Non-default-constructible, iterators
  {
    map_nondefaultconstructible Map = {{0, {0}}, {1, {0}}};
    std::array<map_nondefaultconstructible::entry,2> SourceEntries = {{{1, {1}}, {2, {2}}}};
    Map.Assign(SourceEntries.begin(), SourceEntries.end());
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

}

TEST_F(MapTests, Reserve) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  map Map;
  auto &Keys = helper::GetKeys(Map);
  auto &Entries = helper::GetEntries(Map);
  long long Capacity = Keys.Capacity();
  Map.Reserve(Capacity+1);
  EXPECT_GE(Keys.Capacity(), Capacity+1);
  EXPECT_GE(Entries.Capacity(), Capacity+1);

}

TEST_F(MapTests, Insert) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using map_noncopyable = ovk::map<int,noncopyable<int>>;
  using map_multiargument = ovk::map<int,multiargument>;
  using map_nondefaultconstructible = ovk::map<int,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<map>;
  using helper_noncopyable = ovk::core::test_helper<map_noncopyable>;
  using helper_multiargument = ovk::core::test_helper<map_multiargument>;
  using helper_nondefaultconstructible = ovk::core::test_helper<map_nondefaultconstructible>;

  // lvalue ref, key doesn't already exist
  {
    map Map = {{2, 1}, {4, 3}};
    int SourceValue = 2;
    int &Value = Map.Insert(3, SourceValue);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Keys[2], 4);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 4);
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // lvalue ref, key already exists
  {
    map Map = {{2, 1}, {4, 3}};
    int SourceValue = 2;
    int &Value = Map.Insert(4, SourceValue);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 4);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 4);
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // rvalue ref, key doesn't already exist
  {
    map_noncopyable Map;
    auto &Keys = helper_noncopyable::GetKeys(Map);
    auto &Entries = helper_noncopyable::GetEntries(Map);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {3}});
    noncopyable<int> SourceValue(2);
    noncopyable<int> &Value = Map.Insert(2, std::move(SourceValue));
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // rvalue ref, key already exists
  {
    map_noncopyable Map;
    auto &Keys = helper_noncopyable::GetKeys(Map);
    auto &Entries = helper_noncopyable::GetEntries(Map);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {2}});
    noncopyable<int> SourceValue(3);
    noncopyable<int> &Value = Map.Insert(3, std::move(SourceValue));
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value().Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, key, key doesn't already exist
  {
    map Map = {{2, 1}, {4, 3}};
    int &Value = Map.Insert(3);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Keys[2], 4);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_EQ(Entries(2).Key(), 4);
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, key, key already exists
  {
    map Map = {{2, 1}, {4, 3}};
    int &Value = Map.Insert(4);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 4);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 4);
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, separate IDs, key doesn't already exist
  {
    map Map = {{2, 1}, {4, 3}};
    int &Value = Map.Insert(3);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Keys[2], 4);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_EQ(Entries(2).Key(), 4);
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, separate IDs, key already exists
  {
    map Map = {{2, 1}, {4, 3}};
    int &Value = Map.Insert(4);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 4);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 4);
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, multiple arguments, key doesn't already exist
  {
    map_multiargument Map = {{1, {1,2}}, {3, {3,4}}};
    multiargument &Value = Map.Insert(2, 2, 3);
    auto &Keys = helper_multiargument::GetKeys(Map);
    auto &Entries = helper_multiargument::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().v1, 1);
    EXPECT_EQ(Entries(0).Value().v2, 2);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().v1, 2);
    EXPECT_EQ(Entries(1).Value().v2, 3);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().v1, 3);
    EXPECT_EQ(Entries(2).Value().v2, 4);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, multiple arguments, key already exists
  {
    map_multiargument Map = {{1, {1,2}}, {3, {2,3}}};
    multiargument &Value = Map.Insert(3, 3, 4);
    auto &Keys = helper_multiargument::GetKeys(Map);
    auto &Entries = helper_multiargument::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().v1, 1);
    EXPECT_EQ(Entries(0).Value().v2, 2);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value().v1, 3);
    EXPECT_EQ(Entries(1).Value().v2, 4);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // Non-default-constructible, lvalue ref
  {
    map_nondefaultconstructible Map = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = Map.Insert(2, SourceValue);
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // Non-default-constructible, rvalue ref
  {
    map_nondefaultconstructible Map = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = Map.Insert(2, std::move(SourceValue));
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // Non-default-constructible, in-place
  {
    map_nondefaultconstructible Map = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> &Value = Map.Insert(2, 2);
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // With lower bound iterator, lvalue ref, key doesn't already exist
  {
    map Map = {{2, 1}, {4, 3}};
    int SourceValue = 2;
    auto Iter = Map.Insert(Map.Begin()+1, 3, SourceValue);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Keys[2], 4);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 4);
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // With lower bound iterator, lvalue ref, key already exists
  {
    map Map = {{2, 1}, {4, 3}};
    int SourceValue = 2;
    auto Iter = Map.Insert(Map.Begin()+1, 4, SourceValue);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 4);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 4);
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // With lower bound iterator, rvalue ref, key doesn't already exist
  {
    map_noncopyable Map;
    auto &Keys = helper_noncopyable::GetKeys(Map);
    auto &Entries = helper_noncopyable::GetEntries(Map);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {3}});
    noncopyable<int> SourceValue(2);
    auto Iter = Map.Insert(Map.Begin()+1, 2, std::move(SourceValue));
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // With lower bound iterator, rvalue ref, key already exists
  {
    map_noncopyable Map;
    auto &Keys = helper_noncopyable::GetKeys(Map);
    auto &Entries = helper_noncopyable::GetEntries(Map);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {2}});
    noncopyable<int> SourceValue(3);
    auto Iter = Map.Insert(Map.Begin()+1, 3, std::move(SourceValue));
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value().Value(), 3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // With lower bound iterator, in-place, key doesn't already exist
  {
    map_multiargument Map = {{1, {1,2}}, {3, {3,4}}};
    auto Iter = Map.Insert(Map.Begin()+1, 2, 2, 3);
    auto &Keys = helper_multiargument::GetKeys(Map);
    auto &Entries = helper_multiargument::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().v1, 1);
    EXPECT_EQ(Entries(0).Value().v2, 2);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().v1, 2);
    EXPECT_EQ(Entries(1).Value().v2, 3);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().v1, 3);
    EXPECT_EQ(Entries(2).Value().v2, 4);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // With lower bound iterator, in-place, key already exists
  {
    map_multiargument Map = {{1, {1,2}}, {3, {2,3}}};
    auto Iter = Map.Insert(Map.Begin()+1, 3, 3, 4);
    auto &Keys = helper_multiargument::GetKeys(Map);
    auto &Entries = helper_multiargument::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().v1, 1);
    EXPECT_EQ(Entries(0).Value().v2, 2);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value().v1, 3);
    EXPECT_EQ(Entries(1).Value().v2, 4);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // With lower bound iterator, non-default-constructible, lvalue ref
  {
    map_nondefaultconstructible Map = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    auto Iter = Map.Insert(Map.Begin()+1, 2, SourceValue);
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // With lower bound iterator, non-default-constructible, rvalue ref
  {
    map_nondefaultconstructible Map = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    auto Iter = Map.Insert(Map.Begin()+1, 2, std::move(SourceValue));
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // With lower bound iterator, non-default-constructible, in-place
  {
    map_nondefaultconstructible Map = {{1, {1}}, {3, {3}}};
    auto Iter = Map.Insert(Map.Begin()+1, 2, 2);
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

}

TEST_F(MapTests, Erase) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  // Key, exists
  {
    map Map = {{2, 1}, {3, 2}, {4, 3}};
    Map.Erase(3);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 4);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 4);
    EXPECT_EQ(Entries(1).Value(), 3);
  }

  // Key, doesn't exist
  {
    map Map = {{2, 1}, {4, 3}};
    Map.Erase(3);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 4);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 4);
    EXPECT_EQ(Entries(1).Value(), 3);
  }

  // Iterator
  {
    map Map = {{2, 1}, {3, 2}, {4, 3}};
    auto Iter = Map.Erase(Map.Begin()+1);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 4);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 4);
    EXPECT_EQ(Entries(1).Value(), 3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

}

TEST_F(MapTests, EraseIf) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  // Function takes entry, match exists
  {
    map Map = {{2, 1}, {3, 2}, {4, 3}, {5, 4}};
    Map.EraseIf([](const typename map::entry &Entry) -> bool {
      return Entry.Key() <= 2 || Entry.Value() == 3;
    });
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 3);
    EXPECT_EQ(Keys[1], 5);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 3);
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_EQ(Entries(1).Key(), 5);
    EXPECT_EQ(Entries(1).Value(), 4);
  }

  // Function takes entry, match doesn't exist
  {
    map Map = {{3, 2}, {5, 4}};
    Map.EraseIf([](const typename map::entry &Entry) -> bool {
      return Entry.Key() <= 2 || Entry.Value() == 3;
    });
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 3);
    EXPECT_EQ(Keys[1], 5);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 3);
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_EQ(Entries(1).Key(), 5);
    EXPECT_EQ(Entries(1).Value(), 4);
  }

  // Function takes key, match exists
  {
    map Map = {{2, 1}, {3, 2}, {4, 3}, {5, 4}};
    Map.EraseIf([](int Key) -> bool {
      return Key <= 2 || Key == 4;
    });
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 3);
    EXPECT_EQ(Keys[1], 5);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 3);
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_EQ(Entries(1).Key(), 5);
    EXPECT_EQ(Entries(1).Value(), 4);
  }

  // Function takes key, match doesn't exist
  {
    map Map = {{3, 2}, {5, 4}};
    Map.EraseIf([](int Key) -> bool {
      return Key <= 2 || Key == 4;
    });
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 3);
    EXPECT_EQ(Keys[1], 5);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 3);
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_EQ(Entries(1).Key(), 5);
    EXPECT_EQ(Entries(1).Value(), 4);
  }

}

TEST_F(MapTests, Clear) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  // Empty
  {
    map Map;
    Map.Clear();
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_TRUE(Keys.Empty());
    EXPECT_TRUE(Entries.Empty());
  }

  // Non-empty
  {
    map Map = {{2, 1}, {3, 2}};
    Map.Clear();
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_TRUE(Keys.Empty());
    EXPECT_TRUE(Entries.Empty());
  }

}

TEST_F(MapTests, Contains) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;

  // Exists
  {
    map Map = {{2, 1}, {3, 2}, {4, 3}};
    EXPECT_TRUE(Map.Contains(3));
  }

  // Doesn't exist
  {
    map Map = {{2, 1}, {4, 3}};
    EXPECT_FALSE(Map.Contains(3));
  }

}

TEST_F(MapTests, Find) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;

  // Exists
  {
    map Map = {{2, 1}, {3, 2}, {4, 3}};
    auto Iter = Map.Find(3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // Doesn't exist
  {
    map Map = {{2, 1}, {4, 3}};
    auto Iter = Map.Find(3);
    EXPECT_EQ(Iter, Map.End());
  }

}

TEST_F(MapTests, LowerBound) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;

  // Exists
  {
    map Map = {{2, 1}, {3, 2}, {4, 3}};
    auto Iter = Map.LowerBound(3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

  // Doesn't exist
  {
    map Map = {{2, 1}, {4, 3}};
    auto Iter = Map.LowerBound(3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

}

TEST_F(MapTests, UpperBound) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;

  // Exists
  {
    map Map = {{2, 1}, {3, 2}, {4, 3}};
    auto Iter = Map.UpperBound(3);
    EXPECT_EQ(Iter, Map.Begin()+2);
  }

  // Doesn't exist
  {
    map Map = {{2, 1}, {4, 3}};
    auto Iter = Map.UpperBound(3);
    EXPECT_EQ(Iter, Map.Begin()+1);
  }

}

TEST_F(MapTests, Count) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;

  map Map = {{2, 1}, {3, 2}};
  EXPECT_EQ(Map.Count(), 2);

}

TEST_F(MapTests, Empty) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;

  // Empty
  {
    map Map;
    EXPECT_TRUE(Map.Empty());
  }

  // Not empty
  {
    map Map = {{2, 1}, {3, 2}};
    EXPECT_FALSE(Map.Empty());
  }

}

TEST_F(MapTests, Capacity) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  map Map = {{2, 1}, {3, 2}};
  auto &Keys = helper::GetKeys(Map);
  EXPECT_EQ(Map.Capacity(), Keys.Capacity());

}

TEST_F(MapTests, Keys) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  map Map = {{2, 1}, {3, 2}};
  auto &Keys = helper::GetKeys(Map);
  EXPECT_EQ(&Map.Keys(), &Keys);

}

TEST_F(MapTests, ParenthesisOperator) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  // Key
  {
    map Map = {{2, 1}, {3, 2}};
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(&Map(3), &Entries(1).Value());
  }

  // Separate IDs
  {
    map Map = {{2, 1}, {3, 2}};
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(&Map(3), &Entries(1).Value());
  }

}

TEST_F(MapTests, Fetch) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using map_noncopyable = ovk::map<int,noncopyable<int>>;
  using map_multiargument = ovk::map<int,multiargument>;
  using map_nondefaultconstructible = ovk::map<int,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<map>;
  using helper_noncopyable = ovk::core::test_helper<map_noncopyable>;
  using helper_multiargument = ovk::core::test_helper<map_multiargument>;
  using helper_nondefaultconstructible = ovk::core::test_helper<map_nondefaultconstructible>;

  // lvalue ref, key doesn't already exist
  {
    map Map = {{2, 1}, {4, 3}};
    int SourceValue = 2;
    int &Value = Map.Fetch(3, SourceValue);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Keys[2], 4);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 4);
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // lvalue ref, key already exists
  {
    map Map = {{2, 1}, {4, 3}};
    int SourceValue = 2;
    int &Value = Map.Fetch(4, SourceValue);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 4);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 4);
    EXPECT_EQ(Entries(1).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // rvalue ref, key doesn't already exist
  {
    map_noncopyable Map;
    auto &Keys = helper_noncopyable::GetKeys(Map);
    auto &Entries = helper_noncopyable::GetEntries(Map);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {3}});
    noncopyable<int> SourceValue(2);
    noncopyable<int> &Value = Map.Fetch(2, std::move(SourceValue));
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // rvalue ref, key already exists
  {
    map_noncopyable Map;
    auto &Keys = helper_noncopyable::GetKeys(Map);
    auto &Entries = helper_noncopyable::GetEntries(Map);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {2}});
    noncopyable<int> SourceValue(3);
    noncopyable<int> &Value = Map.Fetch(3, std::move(SourceValue));
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, key doesn't already exist
  {
    map Map = {{2, 1}, {4, 3}};
    int &Value = Map.Fetch(3);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Keys[2], 4);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_EQ(Entries(2).Key(), 4);
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, key already exists
  {
    map Map = {{2, 1}, {4, 3}};
    int &Value = Map.Fetch(4);
    auto &Keys = helper::GetKeys(Map);
    auto &Entries = helper::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 2);
    EXPECT_EQ(Keys[1], 4);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 2);
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 4);
    EXPECT_EQ(Entries(1).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, multiple arguments, key doesn't already exist
  {
    map_multiargument Map = {{1, {1,2}}, {3, {3,4}}};
    multiargument &Value = Map.Fetch(2, 2, 3);
    auto &Keys = helper_multiargument::GetKeys(Map);
    auto &Entries = helper_multiargument::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().v1, 1);
    EXPECT_EQ(Entries(0).Value().v2, 2);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().v1, 2);
    EXPECT_EQ(Entries(1).Value().v2, 3);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().v1, 3);
    EXPECT_EQ(Entries(2).Value().v2, 4);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, multiple arguments, key already exists
  {
    map_multiargument Map = {{1, {1,2}}, {3, {2,3}}};
    multiargument &Value = Map.Fetch(3, 3, 4);
    auto &Keys = helper_multiargument::GetKeys(Map);
    auto &Entries = helper_multiargument::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 3);
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().v1, 1);
    EXPECT_EQ(Entries(0).Value().v2, 2);
    EXPECT_EQ(Entries(1).Key(), 3);
    EXPECT_EQ(Entries(1).Value().v1, 2);
    EXPECT_EQ(Entries(1).Value().v2, 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // Non-default-constructible, lvalue ref
  {
    map_nondefaultconstructible Map = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = Map.Fetch(2, SourceValue);
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // Non-default-constructible, rvalue ref
  {
    map_nondefaultconstructible Map = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = Map.Fetch(2, std::move(SourceValue));
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // Non-default-constructible, in-place
  {
    map_nondefaultconstructible Map = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> &Value = Map.Fetch(2, 2);
    auto &Keys = helper_nondefaultconstructible::GetKeys(Map);
    auto &Entries = helper_nondefaultconstructible::GetEntries(Map);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_EQ(Keys[0], 1);
    EXPECT_EQ(Keys[1], 2);
    EXPECT_EQ(Keys[2], 3);
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_EQ(Entries(0).Key(), 1);
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_EQ(Entries(1).Key(), 2);
    EXPECT_EQ(Entries(1).Value().Value(), 2);
    EXPECT_EQ(Entries(2).Key(), 3);
    EXPECT_EQ(Entries(2).Value().Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

}

TEST_F(MapTests, BracketOperator) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  map Map = {{2, 1}, {3, 2}};
  auto &Entries = helper::GetEntries(Map);
  EXPECT_EQ(&Map[1], Entries.Data(1));

}

TEST_F(MapTests, Data) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  map Map = {{2, 1}, {3, 2}};
  auto &Entries = helper::GetEntries(Map);
  EXPECT_EQ(Map.Data(), Entries.Data());

}

TEST_F(MapTests, BeginEnd) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;
  using helper = ovk::core::test_helper<map>;

  map Map = {{2, 1}, {3, 2}};
  auto &Entries = helper::GetEntries(Map);
  EXPECT_EQ(Map.Begin().Pointer(), Entries.Data());
  EXPECT_EQ(Map.End().Pointer(), Entries.Data()+2);

  int Sum = 0;
  for (auto &Entry : Map) Sum += Entry.Key()+Entry.Value();
  EXPECT_EQ(Sum, 8);

}

TEST_F(MapTests, KeyCompare) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int,std::greater<int>>;

  // Get compare instance
  {
    map Map;
    auto &Compare = Map.KeyCompare();
    EXPECT_TRUE((std::is_same<decltype(Compare), const std::greater<int> &>::value));
    EXPECT_TRUE(Compare(3, 2));
  }

  // Do comparison directly
  {
    map Map;
    EXPECT_TRUE(Map.KeyCompare(3, 2));
  }

}

TEST_F(MapTests, ArrayTraits) {

  if (TestComm().Rank() != 0) return;

  using map = ovk::map<int,int>;

  EXPECT_TRUE(ovk::core::IsArray<map>());
  EXPECT_TRUE((std::is_same<ovk::core::array_value_type<map>, typename map::entry>::value));
  EXPECT_EQ(ovk::core::ArrayRank<map>(), 1);
  EXPECT_TRUE(ovk::core::ArrayHasRuntimeExtents<map>());

  map Map = {{2, 1}, {3, 2}};
  EXPECT_THAT(ovk::core::ArrayExtents(Map).Begin(), ElementsAre(0));
  EXPECT_THAT(ovk::core::ArrayExtents(Map).End(), ElementsAre(2));
  EXPECT_EQ(ovk::core::ArrayData(Map), Map.Data());

}

// TEST_F(MapTests, Benchmark) {

//   if (TestComm().Rank() != 0) return;

//   using ovk::core::timer;

//   constexpr int N = 10000;

//   constexpr int NumRuns = 50;

//   ovk::array<int> Keys;
//   Keys.Reserve(N);

// //   std::srand(std::time(nullptr));
//   std::srand(0);

//   for (int i = 0; i < N; ++i) {
//     Keys.Append(std::rand());
//   }

//   double stdInsertTime, ovkInsertTime;
//   double stdLookupTime, ovkLookupTime;
//   double stdIterateTime, ovkIterateTime;

//   {
//     ovk::map<int,int> Map;
//     timer Timer;
//     Timer.Start();
// //     Map.Reserve(N);
//     for (int i = 0; i < N; ++i) {
//       Map.Insert(Keys(i), i);
//     }
//     Timer.Stop();
//     ovkInsertTime = Timer.Elapsed();
//     Timer.Start();
//     long long Sum = 0;
//     for (int Run = 0; Run < NumRuns; ++Run) {
//       for (int i = 0; i < N; ++i) {
//         int Value = Map(Keys(i));
//         Sum += (long long)(Value);
//       }
//     }
//     Timer.Stop();
//     printf("Sum = %llu\n", Sum);
//     ovkLookupTime = Timer.Elapsed();
//     Timer.Start();
//     Sum = 0;
//     for (int Run = 0; Run < NumRuns; ++Run) {
//       for (auto &Entry : Map) {
//         Sum += (long long)(Entry.Key()) + (long long)(Entry.Value());
//       }
//     }
//     Timer.Stop();
//     printf("Sum = %llu\n", Sum);
//     ovkIterateTime = Timer.Elapsed();
//   }

//   {
//     std::map<int,int> Map;
//     timer Timer;
//     Timer.Start();
//     for (int i = 0; i < N; ++i) {
//       Map.insert({Keys(i), i});
//     }
//     Timer.Stop();
//     stdInsertTime = Timer.Elapsed();
//     Timer.Start();
//     long long Sum = 0;
//     for (int Run = 0; Run < NumRuns; ++Run) {
//       for (int i = 0; i < N; ++i) {
//         int Value = Map[Keys(i)];
//         Sum += (long long)(Value);
//       }
//     }
//     Timer.Stop();
//     printf("Sum = %llu\n", Sum);
//     stdLookupTime = Timer.Elapsed();
//     Timer.Start();
//     Sum = 0;
//     for (int Run = 0; Run < NumRuns; ++Run) {
//       for (auto &Pair : Map) {
//         Sum += (long long)(Pair.first) + (long long)(Pair.second);
//       }
//     }
//     Timer.Stop();
//     printf("Sum = %llu\n", Sum);
//     stdIterateTime = Timer.Elapsed();
//   }

//   printf("Insert time: ovk = %12.6e, std = %12.6e, ovk/std = %12.6f\n", ovkInsertTime,
//     stdInsertTime, ovkInsertTime/stdInsertTime);
//   printf("Lookup time: ovk = %12.6e, std = %12.6e, ovk/std = %12.6f\n", ovkLookupTime,
//     stdLookupTime, ovkLookupTime/stdLookupTime);
//   printf("Iterate time: ovk = %12.6e, std = %12.6e, ovk/std = %12.6f\n", ovkIterateTime,
//     stdIterateTime, ovkIterateTime/stdIterateTime);

// }
