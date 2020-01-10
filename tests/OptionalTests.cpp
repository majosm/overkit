// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Optional.hpp>

#include "tests/MPITest.hpp"
#include "tests/mocks/Noncopyable.hpp"
#include "tests/mocks/Nondefaultconstructible.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>

#include <mpi.h>

#include <type_traits>
#include <utility>

class OptionalTests : public tests::mpi_test {};

using tests::noncopyable;
using tests::nondefaultconstructible;

namespace ovk {
namespace core {
template <typename T> class test_helper<optional<T>> {
public:
  using value_storage_type = typename optional<T>::value_storage;
  static const value_storage_type &GetValueStorage(const optional<T> &Optional) {
    return Optional.ValueStorage_;
  }
  static value_storage_type &GetValueStorage(optional<T> &Optional) {
    return Optional.ValueStorage_;
  }
  static const bool &GetPresent(const optional<T> &Optional) {
    return Optional.Present_;
  }
  static bool &GetPresent(optional<T> &Optional) {
    return Optional.Present_;
  }
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

TEST_F(OptionalTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using optional = ovk::optional<int>;

  EXPECT_TRUE((std::is_same<typename optional::value_type, int>::value));

}

TEST_F(OptionalTests, Create) {

  if (TestComm().Rank() != 0) return;

  using optional = ovk::optional<int>;
  using optional_noncopyable = ovk::optional<noncopyable<int>>;
  using optional_nondefaultconstructible = ovk::optional<nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<optional>;
  using helper_noncopyable = ovk::core::test_helper<optional_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<optional_nondefaultconstructible>;

  // Default
  {
    optional Optional;
    auto &Present = helper::GetPresent(Optional);
    EXPECT_FALSE(Present);
  }

  // lvalue ref
  {
    int SourceValue = 1;
    optional Optional(SourceValue);
    auto &Present = helper::GetPresent(Optional);
    auto &ValueStorage = helper::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get(), 1);
  }

  // rvalue ref
  {
    noncopyable<int> SourceValue(1);
    optional_noncopyable Optional(std::move(SourceValue));
    auto &Present = helper_noncopyable::GetPresent(Optional);
    auto &ValueStorage = helper_noncopyable::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

  // Non-default-constructible, default
  {
    optional_nondefaultconstructible Optional;
    auto &Present = helper_nondefaultconstructible::GetPresent(Optional);
    EXPECT_FALSE(Present);
  }

  // Non-default-constructible, lvalue ref
  {
    nondefaultconstructible<int> SourceValue(1);
    optional_nondefaultconstructible Optional(SourceValue);
    auto &Present = helper_nondefaultconstructible::GetPresent(Optional);
    auto &ValueStorage = helper_nondefaultconstructible::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

  // Non-default-constructible, rvalue ref
  {
    nondefaultconstructible<int> SourceValue(1);
    optional_nondefaultconstructible Optional(std::move(SourceValue));
    auto &Present = helper_nondefaultconstructible::GetPresent(Optional);
    auto &ValueStorage = helper_nondefaultconstructible::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

}

TEST_F(OptionalTests, Assign) {

  if (TestComm().Rank() != 0) return;

  using optional = ovk::optional<int>;
  using optional_noncopyable = ovk::optional<noncopyable<int>>;
  using optional_nondefaultconstructible = ovk::optional<nondefaultconstructible<int>>;
  using optional_multiargument = ovk::optional<multiargument>;
  using helper = ovk::core::test_helper<optional>;
  using helper_noncopyable = ovk::core::test_helper<optional_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<optional_nondefaultconstructible>;
  using helper_multiargument = ovk::core::test_helper<optional_multiargument>;

  // lvalue ref, operator=
  {
    optional Optional;
    int SourceValue = 1;
    Optional = SourceValue;
    auto &Present = helper::GetPresent(Optional);
    auto &ValueStorage = helper::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get(), 1);
  }

  // lvalue ref, Assign
  {
    optional Optional;
    int SourceValue = 1;
    Optional.Assign(SourceValue);
    auto &Present = helper::GetPresent(Optional);
    auto &ValueStorage = helper::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get(), 1);
  }

  // rvalue ref, operator=
  {
    optional_noncopyable Optional;
    noncopyable<int> SourceValue(1);
    Optional = std::move(SourceValue);
    auto &Present = helper_noncopyable::GetPresent(Optional);
    auto &ValueStorage = helper_noncopyable::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

  // rvalue ref, Assign
  {
    optional_noncopyable Optional;
    noncopyable<int> SourceValue(1);
    Optional.Assign(std::move(SourceValue));
    auto &Present = helper_noncopyable::GetPresent(Optional);
    auto &ValueStorage = helper_noncopyable::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

  // In-place, Assign
  {
    optional_multiargument Optional;
    Optional.Assign(1, 2);
    auto &Present = helper_multiargument::GetPresent(Optional);
    auto &ValueStorage = helper_multiargument::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().v1, 1);
    EXPECT_EQ(ValueStorage.Get().v2, 2);
  }

  // Non-default-constructible, lvalue ref, operator=
  {
    optional_nondefaultconstructible Optional;
    nondefaultconstructible<int> SourceValue(1);
    Optional = SourceValue;
    auto &Present = helper_nondefaultconstructible::GetPresent(Optional);
    auto &ValueStorage = helper_nondefaultconstructible::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

  // Non-default-constructible, lvalue ref, Assign
  {
    optional_nondefaultconstructible Optional;
    nondefaultconstructible<int> SourceValue(1);
    Optional.Assign(SourceValue);
    auto &Present = helper_nondefaultconstructible::GetPresent(Optional);
    auto &ValueStorage = helper_nondefaultconstructible::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

  // Non-default-constructible, rvalue ref, operator=
  {
    optional_nondefaultconstructible Optional;
    nondefaultconstructible<int> SourceValue(1);
    Optional = std::move(SourceValue);
    auto &Present = helper_nondefaultconstructible::GetPresent(Optional);
    auto &ValueStorage = helper_nondefaultconstructible::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

  // Non-default-constructible, rvalue ref, Assign
  {
    optional_nondefaultconstructible Optional;
    nondefaultconstructible<int> SourceValue(1);
    Optional.Assign(std::move(SourceValue));
    auto &Present = helper_nondefaultconstructible::GetPresent(Optional);
    auto &ValueStorage = helper_nondefaultconstructible::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

  // Non-default-constructible, in-place, Assign
  {
    optional_nondefaultconstructible Optional;
    Optional.Assign(1);
    auto &Present = helper_nondefaultconstructible::GetPresent(Optional);
    auto &ValueStorage = helper_nondefaultconstructible::GetValueStorage(Optional);
    EXPECT_TRUE(Present);
    EXPECT_EQ(ValueStorage.Get().Value(), 1);
  }

}

TEST_F(OptionalTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using optional = ovk::optional<int>;
  using optional_nondefaultconstructible = ovk::optional<nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<optional>;
  using helper_nondefaultconstructible = ovk::core::test_helper<optional_nondefaultconstructible>;

  // Copy construct
  {
    optional Optional1(1);
    optional Optional2(Optional1);
    auto &Present1 = helper::GetPresent(Optional1);
    auto &Present2 = helper::GetPresent(Optional2);
    auto &ValueStorage1 = helper::GetValueStorage(Optional1);
    auto &ValueStorage2 = helper::GetValueStorage(Optional2);
    EXPECT_TRUE(Present1);
    EXPECT_EQ(ValueStorage1.Get(), 1);
    EXPECT_TRUE(Present2);
    EXPECT_EQ(ValueStorage2.Get(), 1);
  }

  // Copy assign
  {
    optional Optional1(1);
    optional Optional2;
    Optional2 = Optional1;
    auto &Present1 = helper::GetPresent(Optional1);
    auto &Present2 = helper::GetPresent(Optional2);
    auto &ValueStorage1 = helper::GetValueStorage(Optional1);
    auto &ValueStorage2 = helper::GetValueStorage(Optional2);
    EXPECT_TRUE(Present1);
    EXPECT_EQ(ValueStorage1.Get(), 1);
    EXPECT_TRUE(Present2);
    EXPECT_EQ(ValueStorage2.Get(), 1);
  }

  // Non-default-contructible, copy construct
  {
    optional_nondefaultconstructible Optional1(1);
    optional_nondefaultconstructible Optional2(Optional1);
    auto &Present1 = helper_nondefaultconstructible::GetPresent(Optional1);
    auto &Present2 = helper_nondefaultconstructible::GetPresent(Optional2);
    auto &ValueStorage1 = helper_nondefaultconstructible::GetValueStorage(Optional1);
    auto &ValueStorage2 = helper_nondefaultconstructible::GetValueStorage(Optional2);
    EXPECT_TRUE(Present1);
    EXPECT_EQ(ValueStorage1.Get().Value(), 1);
    EXPECT_TRUE(Present2);
    EXPECT_EQ(ValueStorage2.Get().Value(), 1);
  }

  // Non-default-contructible, copy assign
  {
    optional_nondefaultconstructible Optional1(1);
    optional_nondefaultconstructible Optional2;
    Optional2 = Optional1;
    auto &Present1 = helper_nondefaultconstructible::GetPresent(Optional1);
    auto &Present2 = helper_nondefaultconstructible::GetPresent(Optional2);
    auto &ValueStorage1 = helper_nondefaultconstructible::GetValueStorage(Optional1);
    auto &ValueStorage2 = helper_nondefaultconstructible::GetValueStorage(Optional2);
    EXPECT_TRUE(Present1);
    EXPECT_EQ(ValueStorage1.Get().Value(), 1);
    EXPECT_TRUE(Present2);
    EXPECT_EQ(ValueStorage2.Get().Value(), 1);
  }

}

TEST_F(OptionalTests, Move) {

  if (TestComm().Rank() != 0) return;

  using optional_noncopyable = ovk::optional<noncopyable<int>>;
  using optional_nondefaultconstructible = ovk::optional<nondefaultconstructible<int>>;
  using helper_noncopyable = ovk::core::test_helper<optional_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<optional_nondefaultconstructible>;

  // Move construct
  {
    optional_noncopyable Optional1(1);
    optional_noncopyable Optional2(std::move(Optional1));
    auto &Present1 = helper_noncopyable::GetPresent(Optional1);
    auto &Present2 = helper_noncopyable::GetPresent(Optional2);
    auto &ValueStorage2 = helper_noncopyable::GetValueStorage(Optional2);
    EXPECT_FALSE(Present1);
    EXPECT_TRUE(Present2);
    EXPECT_EQ(ValueStorage2.Get().Value(), 1);
  }

  // Move assign
  {
    optional_noncopyable Optional1(1);
    optional_noncopyable Optional2;
    Optional2 = std::move(Optional1);
    auto &Present1 = helper_noncopyable::GetPresent(Optional1);
    auto &Present2 = helper_noncopyable::GetPresent(Optional2);
    auto &ValueStorage2 = helper_noncopyable::GetValueStorage(Optional2);
    EXPECT_FALSE(Present1);
    EXPECT_TRUE(Present2);
    EXPECT_EQ(ValueStorage2.Get().Value(), 1);
  }

  // Non-default-contructible, move construct
  {
    optional_nondefaultconstructible Optional1(1);
    optional_nondefaultconstructible Optional2(std::move(Optional1));
    auto &Present1 = helper_nondefaultconstructible::GetPresent(Optional1);
    auto &Present2 = helper_nondefaultconstructible::GetPresent(Optional2);
    auto &ValueStorage2 = helper_nondefaultconstructible::GetValueStorage(Optional2);
    EXPECT_FALSE(Present1);
    EXPECT_TRUE(Present2);
    EXPECT_EQ(ValueStorage2.Get().Value(), 1);
  }

  // Non-default-contructible, move assign
  {
    optional_nondefaultconstructible Optional1(1);
    optional_nondefaultconstructible Optional2;
    Optional2 = std::move(Optional1);
    auto &Present1 = helper_nondefaultconstructible::GetPresent(Optional1);
    auto &Present2 = helper_nondefaultconstructible::GetPresent(Optional2);
    auto &ValueStorage2 = helper_nondefaultconstructible::GetValueStorage(Optional2);
    EXPECT_FALSE(Present1);
    EXPECT_TRUE(Present2);
    EXPECT_EQ(ValueStorage2.Get().Value(), 1);
  }

}

TEST_F(OptionalTests, Get) {

  if (TestComm().Rank() != 0) return;

  using optional = ovk::optional<int>;
  using optional_noncopyable = ovk::optional<noncopyable<int>>;
  using helper = ovk::core::test_helper<optional>;
  using helper_noncopyable = ovk::core::test_helper<optional_noncopyable>;

  // Get function, const
  {
    const optional Optional(1);
    auto &ValueStorage = helper::GetValueStorage(Optional);
    EXPECT_EQ(&Optional.Get(), &ValueStorage.Get());
    EXPECT_TRUE((std::is_same<decltype(Optional.Get()), const int &>::value));
    EXPECT_EQ(Optional.Get(), 1);
  }

  // Get function, non-const
  {
    optional Optional(1);
    auto &ValueStorage = helper::GetValueStorage(Optional);
    EXPECT_EQ(&Optional.Get(), &ValueStorage.Get());
    EXPECT_TRUE((std::is_same<decltype(Optional.Get()), int &>::value));
    EXPECT_EQ(Optional.Get(), 1);
  }

  // Asterisk operator, const
  {
    const optional Optional(1);
    auto &ValueStorage = helper::GetValueStorage(Optional);
    EXPECT_EQ(&(*Optional), &ValueStorage.Get());
    EXPECT_TRUE((std::is_same<decltype(*Optional), const int &>::value));
    EXPECT_EQ(*Optional, 1);
  }

  // Asterisk operator, non-const
  {
    optional Optional(1);
    auto &ValueStorage = helper::GetValueStorage(Optional);
    EXPECT_EQ(&(*Optional), &ValueStorage.Get());
    EXPECT_TRUE((std::is_same<decltype(*Optional), int &>::value));
    EXPECT_EQ(*Optional, 1);
  }

  // Arrow operator, const
  {
    const optional_noncopyable Optional(1);
    auto &ValueStorage = helper_noncopyable::GetValueStorage(Optional);
    EXPECT_EQ(Optional.operator->(), &ValueStorage.Get());
    EXPECT_TRUE((std::is_same<decltype(Optional.operator->()), const noncopyable<int> *>::value));
    EXPECT_EQ(Optional->Value(), 1);
  }

  // Arrow operator, non-const
  {
    optional_noncopyable Optional(1);
    auto &ValueStorage = helper_noncopyable::GetValueStorage(Optional);
    EXPECT_EQ(Optional.operator->(), &ValueStorage.Get());
    EXPECT_TRUE((std::is_same<decltype(Optional.operator->()), noncopyable<int> *>::value));
    EXPECT_EQ(Optional->Value(), 1);
  }

}

TEST_F(OptionalTests, Release) {

  if (TestComm().Rank() != 0) return;

  using optional_noncopyable = ovk::optional<noncopyable<int>>;
  using helper_noncopyable = ovk::core::test_helper<optional_noncopyable>;

  optional_noncopyable Optional(1);
  noncopyable<int> Value = Optional.Release();
  auto &Present = helper_noncopyable::GetPresent(Optional);
  EXPECT_EQ(Value.Value(), 1);
  EXPECT_FALSE(Present);

}

TEST_F(OptionalTests, Reset) {

  if (TestComm().Rank() != 0) return;

  using optional = ovk::optional<int>;
  using helper = ovk::core::test_helper<optional>;

  optional Optional(1);
  Optional.Reset();
  auto &Present = helper::GetPresent(Optional);
  EXPECT_FALSE(Present);

}
