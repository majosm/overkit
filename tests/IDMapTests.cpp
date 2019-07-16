// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/IDMap.hpp>

#include "tests/MPITest.hpp"
#include "tests/mocks/Noncopyable.hpp"
#include "tests/mocks/Nondefaultconstructible.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Elem.hpp>

#include <mpi.h>

#include <iterator>
#include <utility>

class IDMapEntryTests : public tests::mpi_test {};
class IDMapTests : public tests::mpi_test {};

using testing::ElementsAre;

using tests::noncopyable;
using tests::nondefaultconstructible;

namespace {
template <int Rank, typename T, bool Contiguous> using entry_type = ovk::id_map_internal::
  id_map_entry<Rank,T,Contiguous>;
}

namespace ovk {
namespace core {
template <int Rank, typename T, bool Contiguous> class test_helper<entry_type<Rank,T,Contiguous>> {
public:
  using entry = entry_type<Rank,T,Contiguous>;
  using key_type = typename entry::key_type;
  using value_wrapper_type = typename entry::value_wrapper_type;
  static const key_type &GetKey(const entry &Entry) { return Entry.Key_; }
  static const value_wrapper_type &GetValue(const entry &Entry) { return Entry.Value_; }
  static value_wrapper_type &GetValue(entry &Entry) { return Entry.Value_; }
};
template <int Rank, typename T, bool Contiguous> class test_helper<id_map<Rank,T,Contiguous>> {
public:
  using id_map_type = id_map<Rank,T,Contiguous>;
  using id_set_type = typename id_map_type::id_set_type;
  using entry = typename id_map_type::entry;
  static const id_set_type &GetKeys(const id_map_type &IDMap) { return IDMap.Keys_; }
  static id_set_type &GetKeys(id_map_type &IDMap) { return IDMap.Keys_; }
  static const array<entry> &GetEntries(const id_map_type &IDMap) { return IDMap.Entries_; }
  static array<entry> &GetEntries(id_map_type &IDMap) { return IDMap.Entries_; }
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

TEST_F(IDMapEntryTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using entry = typename ovk::id_map<2,int,true>::entry;

  EXPECT_EQ(int(entry::Rank), 2);
  EXPECT_TRUE((std::is_same<typename entry::key_type, ovk::elem<int,2>>::value));
  EXPECT_TRUE((std::is_same<typename entry::value_type, int>::value));
  EXPECT_TRUE(bool(entry::Contiguous));

  using entry_noncontig = typename ovk::id_map<2,int,false>::entry;

  EXPECT_FALSE(bool(entry_noncontig::Contiguous));

}

TEST_F(IDMapEntryTests, Create) {

  if (TestComm().Rank() != 0) return;

  using entry_contig = entry_type<2,int,true>;
  using entry_noncontig = entry_type<2,int,false>;
  using entry_noncopyable_contig = entry_type<1,noncopyable<int>,true>;
  using entry_noncopyable_noncontig = entry_type<1,noncopyable<int>,false>;
  using entry_multiargument_contig = entry_type<1,multiargument,true>;
  using entry_multiargument_noncontig = entry_type<1,multiargument,false>;
  using entry_nondefaultconstructible_contig = entry_type<1,nondefaultconstructible<int>,true>;
  using entry_nondefaultconstructible_noncontig = entry_type<1,nondefaultconstructible<int>,false>;

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
    entry_contig Entry({1,2}, SourceValue);
    auto &Key = helper_contig::GetKey(Entry);
    auto &Value = helper_contig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1,2));
    EXPECT_EQ(Value.Get(), 1);
  }

  // lvalue ref, non-contiguous
  {
    int SourceValue = 1;
    entry_noncontig Entry({1,2}, SourceValue);
    auto &Key = helper_noncontig::GetKey(Entry);
    auto &Value = helper_noncontig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1,2));
    EXPECT_EQ(Value.Get(), 1);
  }

  // rvalue ref, contiguous
  {
    noncopyable<int> SourceValue(1);
    entry_noncopyable_contig Entry(1, std::move(SourceValue));
    auto &Key = helper_noncopyable_contig::GetKey(Entry);
    auto &Value = helper_noncopyable_contig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // rvalue ref, non-contiguous
  {
    noncopyable<int> SourceValue(1);
    entry_noncopyable_noncontig Entry(1, std::move(SourceValue));
    auto &Key = helper_noncopyable_noncontig::GetKey(Entry);
    auto &Value = helper_noncopyable_noncontig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // In-place, contiguous
  {
    entry_multiargument_contig Entry(1, 1, 2);
    auto &Key = helper_multiargument_contig::GetKey(Entry);
    auto &Value = helper_multiargument_contig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().v1, 1);
    EXPECT_EQ(Value.Get().v2, 2);
  }

  // In-place, non-contiguous
  {
    entry_multiargument_noncontig Entry(1, 1, 2);
    auto &Key = helper_multiargument_noncontig::GetKey(Entry);
    auto &Value = helper_multiargument_noncontig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().v1, 1);
    EXPECT_EQ(Value.Get().v2, 2);
  }

  // Non-default-constructible, lvalue ref, contiguous
  {
    nondefaultconstructible<int> SourceValue(1);
    entry_nondefaultconstructible_contig Entry(1, SourceValue);
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, lvalue ref, non-contiguous
  {
    nondefaultconstructible<int> SourceValue(1);
    entry_nondefaultconstructible_noncontig Entry(1, SourceValue);
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, rvalue ref, contiguous
  {
    nondefaultconstructible<int> SourceValue(1);
    entry_nondefaultconstructible_contig Entry(1, std::move(SourceValue));
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, rvalue ref, non-contiguous
  {
    nondefaultconstructible<int> SourceValue(1);
    entry_nondefaultconstructible_noncontig Entry(1, std::move(SourceValue));
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, in-place, contiguous
  {
    entry_nondefaultconstructible_contig Entry(1, 1);
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, in-place, non-contiguous
  {
    entry_nondefaultconstructible_noncontig Entry(1, 1);
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

}

TEST_F(IDMapEntryTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using entry_contig = entry_type<2,int,true>;
  using entry_noncontig = entry_type<2,int,false>;
  using entry_nondefaultconstructible_contig = entry_type<1,nondefaultconstructible<int>,true>;
  using entry_nondefaultconstructible_noncontig = entry_type<1,nondefaultconstructible<int>,false>;

  using helper_contig = ovk::core::test_helper<entry_contig>;
  using helper_noncontig = ovk::core::test_helper<entry_noncontig>;
  using helper_nondefaultconstructible_contig = ovk::core::test_helper<
    entry_nondefaultconstructible_contig>;
  using helper_nondefaultconstructible_noncontig = ovk::core::test_helper<
    entry_nondefaultconstructible_noncontig>;

  // Contiguous
  {
    entry_contig Entry1({1,2}, 1);
    entry_contig Entry2(Entry1);
    auto &Key = helper_contig::GetKey(Entry2);
    auto &Value = helper_contig::GetValue(Entry2);
    EXPECT_THAT(Key, ElementsAre(1,2));
    EXPECT_EQ(Value.Get(), 1);
  }

  // Non-contiguous
  {
    entry_noncontig Entry1({1,2}, 1);
    entry_noncontig Entry2(Entry1);
    auto &Key = helper_noncontig::GetKey(Entry2);
    auto &Value = helper_noncontig::GetValue(Entry2);
    EXPECT_THAT(Key, ElementsAre(1,2));
    EXPECT_EQ(Value.Get(), 1);
  }

  // Non-default-constructible, contiguous
  {
    entry_nondefaultconstructible_contig Entry1(1, 1);
    entry_nondefaultconstructible_contig Entry2(Entry1);
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry2);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry2);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, non-contiguous
  {
    entry_nondefaultconstructible_noncontig Entry1(1, 1);
    entry_nondefaultconstructible_noncontig Entry2(Entry1);
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry2);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry2);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

}

TEST_F(IDMapEntryTests, Move) {

  if (TestComm().Rank() != 0) return;

  using entry_noncopyable_contig = entry_type<1,noncopyable<int>,true>;
  using entry_noncopyable_noncontig = entry_type<1,noncopyable<int>,false>;
  using entry_nondefaultconstructible_contig = entry_type<1,nondefaultconstructible<int>,true>;
  using entry_nondefaultconstructible_noncontig = entry_type<1,nondefaultconstructible<int>,false>;

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
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-contiguous
  {
    entry_noncopyable_noncontig Entry1(1, 1);
    entry_noncopyable_noncontig Entry2(std::move(Entry1));
    auto &Key = helper_noncopyable_noncontig::GetKey(Entry2);
    auto &Value = helper_noncopyable_noncontig::GetValue(Entry2);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, contiguous
  {
    entry_nondefaultconstructible_contig Entry1(1, 1);
    entry_nondefaultconstructible_contig Entry2(std::move(Entry1));
    auto &Key = helper_nondefaultconstructible_contig::GetKey(Entry2);
    auto &Value = helper_nondefaultconstructible_contig::GetValue(Entry2);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

  // Non-default-constructible, non-contiguous
  {
    entry_nondefaultconstructible_noncontig Entry1(1, 1);
    entry_nondefaultconstructible_noncontig Entry2(std::move(Entry1));
    auto &Key = helper_nondefaultconstructible_noncontig::GetKey(Entry2);
    auto &Value = helper_nondefaultconstructible_noncontig::GetValue(Entry2);
    EXPECT_THAT(Key, ElementsAre(1));
    EXPECT_EQ(Value.Get().Value(), 1);
  }

}

TEST_F(IDMapEntryTests, Key) {

  if (TestComm().Rank() != 0) return;

  using entry_contig = entry_type<1,int,true>;
  using entry_noncontig = entry_type<1,int,false>;

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

TEST_F(IDMapEntryTests, Value) {

  if (TestComm().Rank() != 0) return;

  using entry_contig = entry_type<1,int,true>;
  using entry_noncontig = entry_type<1,int,false>;

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

TEST_F(IDMapTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int,true>;

  EXPECT_EQ(int(id_map::Rank), 2);
  EXPECT_TRUE((std::is_same<typename id_map::key_type, ovk::elem<int,2>>::value));
  EXPECT_TRUE((std::is_same<typename id_map::value_type, int>::value));
  EXPECT_TRUE(bool(id_map::Contiguous));
  EXPECT_TRUE((std::is_same<typename id_map::iterator::pointer, typename id_map::entry *>::value));
  EXPECT_TRUE((std::is_same<typename id_map::const_iterator::pointer, const typename id_map::entry *
    >::value));

}

TEST_F(IDMapTests, Create) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using id_map_nondefaultconstructible = ovk::id_map<1,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<id_map>;
  using helper_nondefaultconstructible = ovk::core::test_helper<id_map_nondefaultconstructible>;

  // Default
  {
    id_map IDMap;
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 0);
    EXPECT_EQ(Entries.Count(), 0);
  }

  // Initializer list
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Iterators
  {
    std::array<id_map::entry,2> SourceEntries = {{{{1,2}, 1}, {{1,3}, 2}}};
    id_map IDMap(SourceEntries.begin(), SourceEntries.end());
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Non-default-constructible, default
  {
    id_map_nondefaultconstructible IDMap;
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 0);
    EXPECT_EQ(Entries.Count(), 0);
  }

  // Non-default-constructible, initializer list
  {
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {2, {2}}};
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1));
    EXPECT_THAT(Keys[1], ElementsAre(2));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1));
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2));
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

  // Non-default-constructible, iterators
  {
    std::array<id_map_nondefaultconstructible::entry,2> SourceEntries = {{{1, {1}}, {2, {2}}}};
    id_map_nondefaultconstructible IDMap(SourceEntries.begin(), SourceEntries.end());
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1));
    EXPECT_THAT(Keys[1], ElementsAre(2));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1));
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2));
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

}

TEST_F(IDMapTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using id_map_nondefaultconstructible = ovk::id_map<2,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<id_map>;
  using helper_nondefaultconstructible = ovk::core::test_helper<id_map_nondefaultconstructible>;

  // Copy construct
  {
    id_map IDMap1 = {{{1,2}, 1}, {{1,3}, 2}};
    id_map IDMap2 = IDMap1;
    auto &Keys1 = helper::GetKeys(IDMap1);
    auto &Entries1 = helper::GetEntries(IDMap1);
    auto &Keys2 = helper::GetKeys(IDMap2);
    auto &Entries2 = helper::GetEntries(IDMap2);
    EXPECT_EQ(Keys1.Count(), 2);
    EXPECT_THAT(Keys1[0], ElementsAre(1,2));
    EXPECT_THAT(Keys1[1], ElementsAre(1,3));
    EXPECT_EQ(Entries1.Count(), 2);
    EXPECT_THAT(Entries1(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries1(0).Value(), 1);
    EXPECT_THAT(Entries1(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries1(1).Value(), 2);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_THAT(Keys2[0], ElementsAre(1,2));
    EXPECT_THAT(Keys2[1], ElementsAre(1,3));
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_THAT(Entries2(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries2(0).Value(), 1);
    EXPECT_THAT(Entries2(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries2(1).Value(), 2);
  }

  // Copy assign
  {
    id_map IDMap1 = {{{1,2}, 1}, {{1,3}, 2}};
    id_map IDMap2;
    IDMap2 = IDMap1;
    auto &Keys1 = helper::GetKeys(IDMap1);
    auto &Entries1 = helper::GetEntries(IDMap1);
    auto &Keys2 = helper::GetKeys(IDMap2);
    auto &Entries2 = helper::GetEntries(IDMap2);
    EXPECT_EQ(Keys1.Count(), 2);
    EXPECT_THAT(Keys1[0], ElementsAre(1,2));
    EXPECT_THAT(Keys1[1], ElementsAre(1,3));
    EXPECT_EQ(Entries1.Count(), 2);
    EXPECT_THAT(Entries1(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries1(0).Value(), 1);
    EXPECT_THAT(Entries1(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries1(1).Value(), 2);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_THAT(Keys2[0], ElementsAre(1,2));
    EXPECT_THAT(Keys2[1], ElementsAre(1,3));
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_THAT(Entries2(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries2(0).Value(), 1);
    EXPECT_THAT(Entries2(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries2(1).Value(), 2);
  }

  // Non-default-constructible, copy construct
  {
    id_map_nondefaultconstructible IDMap1 = {{{1,2}, {1}}, {{1,3}, {2}}};
    id_map_nondefaultconstructible IDMap2 = IDMap1;
    auto &Keys1 = helper_nondefaultconstructible::GetKeys(IDMap1);
    auto &Entries1 = helper_nondefaultconstructible::GetEntries(IDMap1);
    auto &Keys2 = helper_nondefaultconstructible::GetKeys(IDMap2);
    auto &Entries2 = helper_nondefaultconstructible::GetEntries(IDMap2);
    EXPECT_EQ(Keys1.Count(), 2);
    EXPECT_THAT(Keys1[0], ElementsAre(1,2));
    EXPECT_THAT(Keys1[1], ElementsAre(1,3));
    EXPECT_EQ(Entries1.Count(), 2);
    EXPECT_THAT(Entries1(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries1(0).Value().Value(), 1);
    EXPECT_THAT(Entries1(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries1(1).Value().Value(), 2);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_THAT(Keys2[0], ElementsAre(1,2));
    EXPECT_THAT(Keys2[1], ElementsAre(1,3));
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_THAT(Entries2(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries2(0).Value().Value(), 1);
    EXPECT_THAT(Entries2(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries2(1).Value().Value(), 2);
  }

  // Non-default-constructible, copy assign
  {
    id_map_nondefaultconstructible IDMap1 = {{{1,2}, {1}}, {{1,3}, {2}}};
    id_map_nondefaultconstructible IDMap2;
    IDMap2 = IDMap1;
    auto &Keys1 = helper_nondefaultconstructible::GetKeys(IDMap1);
    auto &Entries1 = helper_nondefaultconstructible::GetEntries(IDMap1);
    auto &Keys2 = helper_nondefaultconstructible::GetKeys(IDMap2);
    auto &Entries2 = helper_nondefaultconstructible::GetEntries(IDMap2);
    EXPECT_EQ(Keys1.Count(), 2);
    EXPECT_THAT(Keys1[0], ElementsAre(1,2));
    EXPECT_THAT(Keys1[1], ElementsAre(1,3));
    EXPECT_EQ(Entries1.Count(), 2);
    EXPECT_THAT(Entries1(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries1(0).Value().Value(), 1);
    EXPECT_THAT(Entries1(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries1(1).Value().Value(), 2);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_THAT(Keys2[0], ElementsAre(1,2));
    EXPECT_THAT(Keys2[1], ElementsAre(1,3));
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_THAT(Entries2(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries2(0).Value().Value(), 1);
    EXPECT_THAT(Entries2(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries2(1).Value().Value(), 2);
  }

}

TEST_F(IDMapTests, Move) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using id_map_nondefaultconstructible = ovk::id_map<2,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<id_map>;
  using helper_nondefaultconstructible = ovk::core::test_helper<id_map_nondefaultconstructible>;

  // Move construct
  {
    id_map IDMap1 = {{{1,2}, 1}, {{1,3}, 2}};
    id_map IDMap2 = std::move(IDMap1);
    auto &Keys1 = helper::GetKeys(IDMap1);
    auto &Entries1 = helper::GetEntries(IDMap1);
    auto &Keys2 = helper::GetKeys(IDMap2);
    auto &Entries2 = helper::GetEntries(IDMap2);
    EXPECT_EQ(Keys1.Count(), 0);
    EXPECT_EQ(Entries1.Count(), 0);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_THAT(Keys2[0], ElementsAre(1,2));
    EXPECT_THAT(Keys2[1], ElementsAre(1,3));
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_THAT(Entries2(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries2(0).Value(), 1);
    EXPECT_THAT(Entries2(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries2(1).Value(), 2);
  }

  // Move assign
  {
    id_map IDMap1 = {{{1,2}, 1}, {{1,3}, 2}};
    id_map IDMap2;
    IDMap2 = std::move(IDMap1);
    auto &Keys1 = helper::GetKeys(IDMap1);
    auto &Entries1 = helper::GetEntries(IDMap1);
    auto &Keys2 = helper::GetKeys(IDMap2);
    auto &Entries2 = helper::GetEntries(IDMap2);
    EXPECT_EQ(Keys1.Count(), 0);
    EXPECT_EQ(Entries1.Count(), 0);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_THAT(Keys2[0], ElementsAre(1,2));
    EXPECT_THAT(Keys2[1], ElementsAre(1,3));
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_THAT(Entries2(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries2(0).Value(), 1);
    EXPECT_THAT(Entries2(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries2(1).Value(), 2);
  }

  // Non-default-constructible, move construct
  {
    id_map_nondefaultconstructible IDMap1 = {{{1,2}, {1}}, {{1,3}, {2}}};
    id_map_nondefaultconstructible IDMap2 = std::move(IDMap1);
    auto &Keys1 = helper_nondefaultconstructible::GetKeys(IDMap1);
    auto &Entries1 = helper_nondefaultconstructible::GetEntries(IDMap1);
    auto &Keys2 = helper_nondefaultconstructible::GetKeys(IDMap2);
    auto &Entries2 = helper_nondefaultconstructible::GetEntries(IDMap2);
    EXPECT_EQ(Keys1.Count(), 0);
    EXPECT_EQ(Entries1.Count(), 0);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_THAT(Keys2[0], ElementsAre(1,2));
    EXPECT_THAT(Keys2[1], ElementsAre(1,3));
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_THAT(Entries2(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries2(0).Value().Value(), 1);
    EXPECT_THAT(Entries2(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries2(1).Value().Value(), 2);
  }

  // Non-default-constructible, move assign
  {
    id_map_nondefaultconstructible IDMap1 = {{{1,2}, {1}}, {{1,3}, {2}}};
    id_map_nondefaultconstructible IDMap2;
    IDMap2 = std::move(IDMap1);
    auto &Keys1 = helper_nondefaultconstructible::GetKeys(IDMap1);
    auto &Entries1 = helper_nondefaultconstructible::GetEntries(IDMap1);
    auto &Keys2 = helper_nondefaultconstructible::GetKeys(IDMap2);
    auto &Entries2 = helper_nondefaultconstructible::GetEntries(IDMap2);
    EXPECT_EQ(Keys1.Count(), 0);
    EXPECT_EQ(Entries1.Count(), 0);
    EXPECT_EQ(Keys2.Count(), 2);
    EXPECT_THAT(Keys2[0], ElementsAre(1,2));
    EXPECT_THAT(Keys2[1], ElementsAre(1,3));
    EXPECT_EQ(Entries2.Count(), 2);
    EXPECT_THAT(Entries2(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries2(0).Value().Value(), 1);
    EXPECT_THAT(Entries2(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries2(1).Value().Value(), 2);
  }

}

TEST_F(IDMapTests, Assign) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using id_map_nondefaultconstructible = ovk::id_map<1,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<id_map>;
  using helper_nondefaultconstructible = ovk::core::test_helper<id_map_nondefaultconstructible>;

  // Initializer list, operator=
  {
    id_map IDMap = {{{0,1}, 0}, {{0,2}, 0}};
    IDMap = {{{1,2}, 1}, {{1,3}, 2}};
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Initializer list, Assign
  {
    id_map IDMap = {{{0,1}, 0}, {{0,2}, 0}};
    IDMap.Assign({{{1,2}, 1}, {{1,3}, 2}});
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Iterators, Assign
  {
    id_map IDMap = {{{0,1}, 0}, {{0,2}, 0}};
    std::array<id_map::entry,2> SourceEntries = {{{{1,2}, 1}, {{1,3}, 2}}};
    IDMap.Assign(SourceEntries.begin(), SourceEntries.end());
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 2);
  }

  // Non-default-constructible, initializer list, operator=
  {
    id_map_nondefaultconstructible IDMap = {{0, {0}}, {1, {0}}};
    IDMap = {{1, {1}}, {2, {2}}};
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1));
    EXPECT_THAT(Keys[1], ElementsAre(2));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1));
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2));
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

  // Non-default-constructible, initializer list, Assign
  {
    id_map_nondefaultconstructible IDMap = {{0, {0}}, {1, {0}}};
    IDMap.Assign({{1, {1}}, {2, {2}}});
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1));
    EXPECT_THAT(Keys[1], ElementsAre(2));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1));
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2));
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

  // Non-default-constructible, iterators
  {
    id_map_nondefaultconstructible IDMap = {{0, {0}}, {1, {0}}};
    std::array<id_map_nondefaultconstructible::entry,2> SourceEntries = {{{1, {1}}, {2, {2}}}};
    IDMap.Assign(SourceEntries.begin(), SourceEntries.end());
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1));
    EXPECT_THAT(Keys[1], ElementsAre(2));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1));
    EXPECT_EQ(Entries(0).Value().Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2));
    EXPECT_EQ(Entries(1).Value().Value(), 2);
  }

}

TEST_F(IDMapTests, Reserve) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  id_map IDMap;
  auto &Keys = helper::GetKeys(IDMap);
  auto &Entries = helper::GetEntries(IDMap);
  long long Capacity = Keys.Capacity();
  IDMap.Reserve(Capacity+1);
  EXPECT_GE(Keys.Capacity(), Capacity+1);
  EXPECT_GE(Entries.Capacity(), Capacity+1);

}

TEST_F(IDMapTests, Insert) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using id_map_noncopyable = ovk::id_map<1,noncopyable<int>>;
  using id_map_multiargument = ovk::id_map<1,multiargument>;
  using id_map_nondefaultconstructible = ovk::id_map<1,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<id_map>;
  using helper_noncopyable = ovk::core::test_helper<id_map_noncopyable>;
  using helper_multiargument = ovk::core::test_helper<id_map_multiargument>;
  using helper_nondefaultconstructible = ovk::core::test_helper<id_map_nondefaultconstructible>;

  // lvalue ref, key doesn't already exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int SourceValue = 2;
    int &Value = IDMap.Insert({1,3}, SourceValue);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_THAT(Keys[2], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_THAT(Entries(2).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // lvalue ref, key already exists
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int SourceValue = 2;
    int &Value = IDMap.Insert({2,1}, SourceValue);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // rvalue ref, key doesn't already exist
  {
    id_map_noncopyable IDMap;
    auto &Keys = helper_noncopyable::GetKeys(IDMap);
    auto &Entries = helper_noncopyable::GetEntries(IDMap);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {3}});
    noncopyable<int> SourceValue(2);
    noncopyable<int> &Value = IDMap.Insert(2, std::move(SourceValue));
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
    id_map_noncopyable IDMap;
    auto &Keys = helper_noncopyable::GetKeys(IDMap);
    auto &Entries = helper_noncopyable::GetEntries(IDMap);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {2}});
    noncopyable<int> SourceValue(3);
    noncopyable<int> &Value = IDMap.Insert(3, std::move(SourceValue));
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
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int &Value = IDMap.Insert({1,3});
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_THAT(Keys[2], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_THAT(Entries(2).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, key, key already exists
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int &Value = IDMap.Insert({2,1});
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, separate IDs, key doesn't already exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int &Value = IDMap.Insert(1,3);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_THAT(Keys[2], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_THAT(Entries(2).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, separate IDs, key already exists
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int &Value = IDMap.Insert(2,1);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, multiple arguments, key doesn't already exist
  {
    id_map_multiargument IDMap = {{1, {1,2}}, {3, {3,4}}};
    multiargument &Value = IDMap.Insert(2, 2, 3);
    auto &Keys = helper_multiargument::GetKeys(IDMap);
    auto &Entries = helper_multiargument::GetEntries(IDMap);
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
    id_map_multiargument IDMap = {{1, {1,2}}, {3, {2,3}}};
    multiargument &Value = IDMap.Insert(3, 3, 4);
    auto &Keys = helper_multiargument::GetKeys(IDMap);
    auto &Entries = helper_multiargument::GetEntries(IDMap);
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
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = IDMap.Insert(2, SourceValue);
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
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
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = IDMap.Insert(2, std::move(SourceValue));
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
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
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> &Value = IDMap.Insert(2, 2);
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
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
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int SourceValue = 2;
    int &Value = IDMap.Insert(IDMap.Begin()+1, {1,3}, SourceValue);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_THAT(Keys[2], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_THAT(Entries(2).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // With lower bound iterator, lvalue ref, key already exists
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int SourceValue = 2;
    int &Value = IDMap.Insert(IDMap.Begin()+1, {2,1}, SourceValue);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // With lower bound iterator, rvalue ref, key doesn't already exist
  {
    id_map_noncopyable IDMap;
    auto &Keys = helper_noncopyable::GetKeys(IDMap);
    auto &Entries = helper_noncopyable::GetEntries(IDMap);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {3}});
    noncopyable<int> SourceValue(2);
    noncopyable<int> &Value = IDMap.Insert(IDMap.Begin()+1, 2, std::move(SourceValue));
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

  // With lower bound iterator, rvalue ref, key already exists
  {
    id_map_noncopyable IDMap;
    auto &Keys = helper_noncopyable::GetKeys(IDMap);
    auto &Entries = helper_noncopyable::GetEntries(IDMap);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {2}});
    noncopyable<int> SourceValue(3);
    noncopyable<int> &Value = IDMap.Insert(IDMap.Begin()+1, 3, std::move(SourceValue));
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

  // With lower bound iterator, in-place, key doesn't already exist
  {
    id_map_multiargument IDMap = {{1, {1,2}}, {3, {3,4}}};
    multiargument &Value = IDMap.Insert(IDMap.Begin()+1, 2, 2, 3);
    auto &Keys = helper_multiargument::GetKeys(IDMap);
    auto &Entries = helper_multiargument::GetEntries(IDMap);
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

  // With lower bound iterator, in-place, key already exists
  {
    id_map_multiargument IDMap = {{1, {1,2}}, {3, {2,3}}};
    multiargument &Value = IDMap.Insert(IDMap.Begin()+1, 3, 3, 4);
    auto &Keys = helper_multiargument::GetKeys(IDMap);
    auto &Entries = helper_multiargument::GetEntries(IDMap);
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

  // With lower bound iterator, non-default-constructible, lvalue ref
  {
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = IDMap.Insert(IDMap.Begin()+1, 2, SourceValue);
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
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

  // With lower bound iterator, non-default-constructible, rvalue ref
  {
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = IDMap.Insert(IDMap.Begin()+1, 2, std::move(SourceValue));
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
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

  // With lower bound iterator, non-default-constructible, in-place
  {
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> &Value = IDMap.Insert(IDMap.Begin()+1, 2, 2);
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
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

TEST_F(IDMapTests, Erase) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  // Key, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    IDMap.Erase({1,3});
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 3);
  }

  // Key, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    IDMap.Erase({1,3});
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 3);
  }

  // Separate IDs, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    IDMap.Erase(1,3);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 3);
  }

  // Separate IDs, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    IDMap.Erase(1,3);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 3);
  }

  // Iterator
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    IDMap.Erase(IDMap.Begin()+1);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 3);
  }

}

TEST_F(IDMapTests, EraseIf) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  // Function takes entry, match exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}, {{1,4}, 4}};
    IDMap.EraseIf([](const typename id_map::entry &Entry) -> bool {
      return Entry.Key(1) == 2 || Entry.Value() == 3;
    });
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,3));
    EXPECT_THAT(Keys[1], ElementsAre(1,4));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,4));
    EXPECT_EQ(Entries(1).Value(), 4);
  }

  // Function takes entry, match doesn't exist
  {
    id_map IDMap = {{{1,3}, 2}, {{1,4}, 4}};
    IDMap.EraseIf([](const typename id_map::entry &Entry) -> bool {
      return Entry.Key(1) == 2 || Entry.Value() == 3;
    });
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,3));
    EXPECT_THAT(Keys[1], ElementsAre(1,4));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,4));
    EXPECT_EQ(Entries(1).Value(), 4);
  }

  // Function takes key, match exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}, {{1,4}, 4}};
    IDMap.EraseIf([](const ovk::elem<int,2> &Key) -> bool {
      return Key(0) == 2 || Key(1) == 2;
    });
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,3));
    EXPECT_THAT(Keys[1], ElementsAre(1,4));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,4));
    EXPECT_EQ(Entries(1).Value(), 4);
  }

  // Function takes key, match doesn't exist
  {
    id_map IDMap = {{{1,3}, 2}, {{1,4}, 4}};
    IDMap.EraseIf([](const ovk::elem<int,2> &Key) -> bool {
      return Key(0) == 2 || Key(1) == 2;
    });
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,3));
    EXPECT_THAT(Keys[1], ElementsAre(1,4));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,4));
    EXPECT_EQ(Entries(1).Value(), 4);
  }

  // Function takes separate IDs, match exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}, {{1,4},4}};
    IDMap.EraseIf([](int M, int N) -> bool {
      return M == 2 || N == 2;
    });
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,3));
    EXPECT_THAT(Keys[1], ElementsAre(1,4));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,4));
    EXPECT_EQ(Entries(1).Value(), 4);
  }

  // Function takes separate IDs, match doesn't exist
  {
    id_map IDMap = {{{1,3}, 2}, {{1,4}, 4}};
    IDMap.EraseIf([](int M, int N) -> bool {
      return M == 2 || N == 2;
    });
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,3));
    EXPECT_THAT(Keys[1], ElementsAre(1,4));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(0).Value(), 2);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,4));
    EXPECT_EQ(Entries(1).Value(), 4);
  }

}

TEST_F(IDMapTests, Clear) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  // Empty
  {
    id_map IDMap;
    IDMap.Clear();
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_TRUE(Keys.Empty());
    EXPECT_TRUE(Entries.Empty());
  }

  // Non-empty
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
    IDMap.Clear();
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_TRUE(Keys.Empty());
    EXPECT_TRUE(Entries.Empty());
  }

}

TEST_F(IDMapTests, Contains) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2, int>;

  // Key, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    EXPECT_TRUE(IDMap.Contains({1,3}));
  }

  // Key, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    EXPECT_FALSE(IDMap.Contains({1,3}));
  }

  // Separate IDs, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    EXPECT_TRUE(IDMap.Contains(1,3));
  }

  // Separate IDs, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    EXPECT_FALSE(IDMap.Contains(1,3));
  }

}

TEST_F(IDMapTests, Find) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;

  // Key, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    auto Iter = IDMap.Find({1,3});
    EXPECT_EQ(Iter, IDMap.Begin()+1);
  }

  // Key, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    auto Iter = IDMap.Find({1,3});
    EXPECT_EQ(Iter, IDMap.End());
  }

  // Separate IDs, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    auto Iter = IDMap.Find(1,3);
    EXPECT_EQ(Iter, IDMap.Begin()+1);
  }

  // Separate IDs, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    auto Iter = IDMap.Find(1,3);
    EXPECT_EQ(Iter, IDMap.End());
  }

}

TEST_F(IDMapTests, LowerBound) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;

  // Key, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    auto Iter = IDMap.LowerBound({1,3});
    EXPECT_EQ(Iter, IDMap.Begin()+1);
  }

  // Key, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    auto Iter = IDMap.LowerBound({1,3});
    EXPECT_EQ(Iter, IDMap.Begin()+1);
  }

  // Separate IDs, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    auto Iter = IDMap.LowerBound(1,3);
    EXPECT_EQ(Iter, IDMap.Begin()+1);
  }

  // Separate IDs, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    auto Iter = IDMap.LowerBound(1,3);
    EXPECT_EQ(Iter, IDMap.Begin()+1);
  }

}

TEST_F(IDMapTests, UpperBound) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;

  // Key, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    auto Iter = IDMap.UpperBound({1,3});
    EXPECT_EQ(Iter, IDMap.Begin()+2);
  }

  // Key, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    auto Iter = IDMap.UpperBound({1,3});
    EXPECT_EQ(Iter, IDMap.Begin()+1);
  }

  // Separate IDs, exists
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}, {{2,1}, 3}};
    auto Iter = IDMap.UpperBound(1,3);
    EXPECT_EQ(Iter, IDMap.Begin()+2);
  }

  // Separate IDs, doesn't exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    auto Iter = IDMap.UpperBound(1,3);
    EXPECT_EQ(Iter, IDMap.Begin()+1);
  }

}

TEST_F(IDMapTests, NextAvailableKey) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<1,int>;

  // First ID > 0
  {
    id_map IDMap = {{1, 1}, {2, 1}};
    EXPECT_EQ(IDMap.NextAvailableKey(), 0);
  }

  // Gap
  {
    id_map IDMap = {{0, 0}, {2, 2}};
    EXPECT_EQ(IDMap.NextAvailableKey(), 1);
  }

  // End
  {
    id_map IDMap = {{0, 0}, {1, 1}};
    EXPECT_EQ(IDMap.NextAvailableKey(), 2);
  }

}

TEST_F(IDMapTests, Count) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;

  id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
  EXPECT_EQ(IDMap.Count(), 2);

}

TEST_F(IDMapTests, Empty) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2, int>;

  // Empty
  {
    id_map IDMap;
    EXPECT_TRUE(IDMap.Empty());
  }

  // Not empty
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
    EXPECT_FALSE(IDMap.Empty());
  }

}

TEST_F(IDMapTests, Capacity) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
  auto &Keys = helper::GetKeys(IDMap);
  EXPECT_EQ(IDMap.Capacity(), Keys.Capacity());

}

TEST_F(IDMapTests, Keys) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
  auto &Keys = helper::GetKeys(IDMap);
  EXPECT_EQ(&IDMap.Keys(), &Keys);

}

TEST_F(IDMapTests, ParenthesisOperator) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  // Key
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
    auto &Entries = helper::GetEntries(IDMap);
    ovk::elem<int,2> Key = {1,3};
    EXPECT_EQ(&IDMap(Key), &Entries(1).Value());
  }

  // Separate IDs
  {
    id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(&IDMap(1,3), &Entries(1).Value());
  }

}

TEST_F(IDMapTests, Get) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using id_map_noncopyable = ovk::id_map<1,noncopyable<int>>;
  using id_map_multiargument = ovk::id_map<1,multiargument>;
  using id_map_nondefaultconstructible = ovk::id_map<1,nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<id_map>;
  using helper_noncopyable = ovk::core::test_helper<id_map_noncopyable>;
  using helper_multiargument = ovk::core::test_helper<id_map_multiargument>;
  using helper_nondefaultconstructible = ovk::core::test_helper<id_map_nondefaultconstructible>;

  // lvalue ref, key doesn't already exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int SourceValue = 2;
    int &Value = IDMap.Get({1,3}, SourceValue);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_THAT(Keys[2], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 2);
    EXPECT_THAT(Entries(2).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // lvalue ref, key already exists
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int SourceValue = 2;
    int &Value = IDMap.Get({2,1}, SourceValue);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // rvalue ref, key doesn't already exist
  {
    id_map_noncopyable IDMap;
    auto &Keys = helper_noncopyable::GetKeys(IDMap);
    auto &Entries = helper_noncopyable::GetEntries(IDMap);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {3}});
    noncopyable<int> SourceValue(2);
    noncopyable<int> &Value = IDMap.Get(2, std::move(SourceValue));
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
    id_map_noncopyable IDMap;
    auto &Keys = helper_noncopyable::GetKeys(IDMap);
    auto &Entries = helper_noncopyable::GetEntries(IDMap);
    Keys.Insert(1);
    Keys.Insert(3);
    Entries.Append({1, {1}});
    Entries.Append({3, {2}});
    noncopyable<int> SourceValue(3);
    noncopyable<int> &Value = IDMap.Get(3, std::move(SourceValue));
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

  // In-place, default, key, key doesn't already exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int &Value = IDMap.Get({1,3});
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_THAT(Keys[2], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_THAT(Entries(2).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, key, key already exists
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int &Value = IDMap.Get({2,1});
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, separate IDs, key doesn't already exist
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int &Value = IDMap.Get(1,3);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 3);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(1,3));
    EXPECT_THAT(Keys[2], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 3);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(1,3));
    EXPECT_EQ(Entries(1).Value(), 0);
    EXPECT_THAT(Entries(2).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(2).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, default, separate IDs, key already exists
  {
    id_map IDMap = {{{1,2}, 1}, {{2,1}, 3}};
    int &Value = IDMap.Get(2,1);
    auto &Keys = helper::GetKeys(IDMap);
    auto &Entries = helper::GetEntries(IDMap);
    EXPECT_EQ(Keys.Count(), 2);
    EXPECT_THAT(Keys[0], ElementsAre(1,2));
    EXPECT_THAT(Keys[1], ElementsAre(2,1));
    EXPECT_EQ(Entries.Count(), 2);
    EXPECT_THAT(Entries(0).Key(), ElementsAre(1,2));
    EXPECT_EQ(Entries(0).Value(), 1);
    EXPECT_THAT(Entries(1).Key(), ElementsAre(2,1));
    EXPECT_EQ(Entries(1).Value(), 3);
    EXPECT_EQ(&Value, &Entries(1).Value());
  }

  // In-place, multiple arguments, key doesn't already exist
  {
    id_map_multiargument IDMap = {{1, {1,2}}, {3, {3,4}}};
    multiargument &Value = IDMap.Get(2, 2, 3);
    auto &Keys = helper_multiargument::GetKeys(IDMap);
    auto &Entries = helper_multiargument::GetEntries(IDMap);
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
    id_map_multiargument IDMap = {{1, {1,2}}, {3, {2,3}}};
    multiargument &Value = IDMap.Get(3, 3, 4);
    auto &Keys = helper_multiargument::GetKeys(IDMap);
    auto &Entries = helper_multiargument::GetEntries(IDMap);
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
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = IDMap.Get(2, SourceValue);
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
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
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> SourceValue(2);
    nondefaultconstructible<int> &Value = IDMap.Get(2, std::move(SourceValue));
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
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
    id_map_nondefaultconstructible IDMap = {{1, {1}}, {3, {3}}};
    nondefaultconstructible<int> &Value = IDMap.Get(2, 2);
    auto &Keys = helper_nondefaultconstructible::GetKeys(IDMap);
    auto &Entries = helper_nondefaultconstructible::GetEntries(IDMap);
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

TEST_F(IDMapTests, BracketOperator) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
  auto &Entries = helper::GetEntries(IDMap);
  EXPECT_EQ(&IDMap[1], Entries.Data(1));

}

TEST_F(IDMapTests, Data) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
  auto &Entries = helper::GetEntries(IDMap);
  EXPECT_EQ(IDMap.Data(), Entries.Data());

}

TEST_F(IDMapTests, BeginEnd) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;
  using helper = ovk::core::test_helper<id_map>;

  id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
  auto &Entries = helper::GetEntries(IDMap);
  EXPECT_EQ(IDMap.Begin().Pointer(), Entries.Data());
  EXPECT_EQ(IDMap.End().Pointer(), Entries.Data()+2);

  int Sum = 0;
  for (auto &Entry : IDMap) Sum += Entry.Key(0)+Entry.Key(1)+Entry.Value();
  EXPECT_EQ(Sum, 10);

}

TEST_F(IDMapTests, KeyLess) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;

  // Equal
  {
    EXPECT_FALSE(id_map::KeyLess({1,2}, {1,2}));
  }

  // Same first index, second index less
  {
    EXPECT_TRUE(id_map::KeyLess({1,1}, {1,2}));
  }

  // Same first index, second index greater
  {
    EXPECT_FALSE(id_map::KeyLess({1,3}, {1,2}));
  }

  // First index less
  {
    EXPECT_TRUE(id_map::KeyLess({0,2}, {1,2}));
  }

  // First index greater
  {
    EXPECT_FALSE(id_map::KeyLess({2,2}, {1,2}));
  }

}

TEST_F(IDMapTests, ArrayTraits) {

  if (TestComm().Rank() != 0) return;

  using id_map = ovk::id_map<2,int>;

  EXPECT_TRUE(ovk::core::IsArray<id_map>());
  EXPECT_TRUE((std::is_same<ovk::core::array_value_type<id_map>, typename id_map::entry>::value));
  EXPECT_EQ(ovk::core::ArrayRank<id_map>(), 1);
  EXPECT_TRUE(ovk::core::ArrayHasRuntimeExtents<id_map>());

  id_map IDMap = {{{1,2}, 1}, {{1,3}, 2}};
  EXPECT_THAT(ovk::core::ArrayExtents(IDMap).Begin(), ElementsAre(0));
  EXPECT_THAT(ovk::core::ArrayExtents(IDMap).End(), ElementsAre(2));
  EXPECT_EQ(ovk::core::ArrayData(IDMap), IDMap.Data());

}
