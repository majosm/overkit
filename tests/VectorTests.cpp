// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Vector.hpp>

#include "tests/MPITest.hpp"
#include "tests/mocks/Noncopyable.hpp"
#include "tests/mocks/Nondefaultconstructible.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>

#include <mpi.h>

#include <iterator>
#include <utility>
#include <vector>

using testing::ElementsAre;
using testing::ElementsAreArray;

class VectorTests : public tests::mpi_test {};

using tests::noncopyable;
using tests::nondefaultconstructible;

namespace ovk {
namespace core {
template <typename T, typename Allocator> class test_helper<vector<T, Allocator>> {
public:
  using vector_type = vector<T, Allocator>;
  using storage_value_type = typename vector_type::storage_value_type;
  using storage_allocator_type = typename vector_type::storage_allocator_type;
  static const std::vector<storage_value_type, storage_allocator_type> &GetValues(const vector_type
    &Vector) {
    return Vector.Values_;
  }
  static std::vector<storage_value_type, storage_allocator_type> &GetValues(vector_type &Vector) {
    return Vector.Values_;
  }
};
}}

namespace {

using not_bool = ovk::core::vector_internal::not_bool;

struct multiargument {
  int v1, v2;
  multiargument() = default;
  multiargument(int v1_, int v2_):
    v1(v1_),
    v2(v2_)
  {}
};

}

TEST_F(VectorTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using vector = ovk::core::vector<int>;
  using helper = ovk::core::test_helper<vector>;

  EXPECT_TRUE((std::is_same<typename vector::value_type, int>::value));
  EXPECT_TRUE((std::is_same<typename vector::index_type, long long>::value));
  EXPECT_TRUE((std::is_same<typename vector::iterator::pointer, int *>::value));
  EXPECT_TRUE((std::is_same<typename vector::const_iterator::pointer, const int *>::value));
  EXPECT_TRUE((std::is_same<typename helper::storage_value_type, int>::value));

  using vector_bool = ovk::core::vector<bool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;

  EXPECT_TRUE((std::is_same<typename vector_bool::value_type, bool>::value));
  EXPECT_TRUE((std::is_same<typename vector_bool::iterator::pointer, bool *>::value));
  EXPECT_TRUE((std::is_same<typename vector_bool::const_iterator::pointer, const bool *>::value));
  EXPECT_TRUE((std::is_same<typename helper_bool::storage_value_type, not_bool>::value));

}

TEST_F(VectorTests, Create) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using vector_noncopyable = ovk::core::vector<noncopyable<int>>;
  using vector_nondefaultconstructible = ovk::core::vector<nondefaultconstructible<int>>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;
  using helper_noncopyable = ovk::core::test_helper<vector_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<vector_nondefaultconstructible>;

  // Non-bool, default
  {
    vector_nonbool Vector;
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 0);
  }

  // Non-bool, size
  {
    vector_nonbool Vector(5);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 5);
  }

  // Non-bool, size and value
  {
    vector_nonbool Vector(5, 1);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1}));
  }

  // Non-bool, initializer list
  {
    vector_nonbool Vector = {0,1,2,3,4};
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4}));
  }

  // Bool, default
  {
    vector_bool Vector;
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 0);
  }

  // Bool, size
  {
    vector_bool Vector(5);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 5);
  }

  // Bool, size and value
  {
    vector_bool Vector(5, true);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::TRUE,not_bool::TRUE,not_bool::TRUE,
      not_bool::TRUE,not_bool::TRUE}));
  }

  // Bool, initializer list
  {
    vector_bool Vector = {false,true,false,false,true};
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::FALSE,not_bool::TRUE,not_bool::FALSE,
      not_bool::FALSE,not_bool::TRUE}));
  }

  // Non-copyable, default
  {
    vector_noncopyable Vector;
    auto &Values = helper_noncopyable::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 0);
  }

  // Non-copyable, size
  {
    vector_noncopyable Vector(5);
    auto &Values = helper_noncopyable::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 5);
  }

  // Non-default-constructible, default
  {
    vector_nondefaultconstructible Vector;
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 0);
  }

  // Non-default-constructible, size and value
  {
    vector_nondefaultconstructible Vector(5, nondefaultconstructible<int>(1));
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 5);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 5);
  }

  // Non-default-constructible, initializer list
  {
    vector_nondefaultconstructible Vector = {nondefaultconstructible<int>(1),
      nondefaultconstructible<int>(2)};
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
  }

}

TEST_F(VectorTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using vector = ovk::core::vector<int>;
  using vector_nondefaultconstructible = ovk::core::vector<nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<vector>;
  using helper_nondefaultconstructible = ovk::core::test_helper<vector_nondefaultconstructible>;

  // Copy construct
  {
    vector Vector1 = {0,1,2,3,4};
    vector Vector2(Vector1);
    auto &Values1 = helper::GetValues(Vector1);
    auto &Values2 = helper::GetValues(Vector2);
    EXPECT_THAT(Values1, ElementsAreArray({0,1,2,3,4}));
    EXPECT_THAT(Values2, ElementsAreArray({0,1,2,3,4}));
  }

  // Copy assign
  {
    vector Vector1 = {0,1,2,3,4};
    vector Vector2;
    Vector2 = Vector1;
    auto &Values1 = helper::GetValues(Vector1);
    auto &Values2 = helper::GetValues(Vector2);
    EXPECT_THAT(Values1, ElementsAreArray({0,1,2,3,4}));
    EXPECT_THAT(Values2, ElementsAreArray({0,1,2,3,4}));
  }

  // Non-default-contructible, copy construct
  {
    vector_nondefaultconstructible Vector1 = {{0},{1},{2},{3},{4}};
    vector_nondefaultconstructible Vector2(Vector1);
    auto &Values1 = helper_nondefaultconstructible::GetValues(Vector1);
    auto &Values2 = helper_nondefaultconstructible::GetValues(Vector2);
    int Sum = 0;
    for (auto &Value : Values1) Sum += Value.Value();
    EXPECT_EQ(Sum, 10);
    Sum = 0;
    for (auto &Value : Values2) Sum += Value.Value();
    EXPECT_EQ(Sum, 10);
  }

  // Non-default-contructible, copy assign
  {
    vector_nondefaultconstructible Vector1 = {{0},{1},{2},{3},{4}};
    vector_nondefaultconstructible Vector2;
    Vector2 = Vector1;
    auto &Values1 = helper_nondefaultconstructible::GetValues(Vector1);
    auto &Values2 = helper_nondefaultconstructible::GetValues(Vector2);
    int Sum = 0;
    for (auto &Value : Values1) Sum += Value.Value();
    EXPECT_EQ(Sum, 10);
    Sum = 0;
    for (auto &Value : Values2) Sum += Value.Value();
    EXPECT_EQ(Sum, 10);
  }

}

TEST_F(VectorTests, Move) {

  if (TestComm().Rank() != 0) return;

  using vector = ovk::core::vector<int>;
  using vector_noncopyable = ovk::core::vector<noncopyable<int>>;
  using vector_nondefaultconstructible = ovk::core::vector<nondefaultconstructible<int>>;
  using helper = ovk::core::test_helper<vector>;
  using helper_noncopyable = ovk::core::test_helper<vector_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<vector_nondefaultconstructible>;

  // Move construct
  {
    vector Vector1 = {0,1,2,3,4};
    vector Vector2(std::move(Vector1));
    auto &Values1 = helper::GetValues(Vector1);
    auto &Values2 = helper::GetValues(Vector2);
    EXPECT_TRUE(Values1.empty());
    EXPECT_THAT(Values2, ElementsAreArray({0,1,2,3,4}));
  }

  // Move assign
  {
    vector Vector1 = {0,1,2,3,4};
    vector Vector2;
    Vector2 = std::move(Vector1);
    auto &Values1 = helper::GetValues(Vector1);
    auto &Values2 = helper::GetValues(Vector2);
    EXPECT_TRUE(Values1.empty());
    EXPECT_THAT(Values2, ElementsAreArray({0,1,2,3,4}));
  }

  // Non-copyable, move construct
  {
    vector_noncopyable Vector1;
    auto &Values1 = helper_noncopyable::GetValues(Vector1);
    for (int i = 0; i < 5; ++i) Values1.emplace_back(i);
    vector_noncopyable Vector2(std::move(Vector1));
    auto &Values2 = helper_noncopyable::GetValues(Vector2);
    EXPECT_TRUE(Values1.empty());
    int Sum = 0;
    for (auto &Value : Values2) Sum += Value.Value();
    EXPECT_EQ(Sum, 10);
  }

  // Non-copyable, move assign
  {
    vector_noncopyable Vector1;
    auto &Values1 = helper_noncopyable::GetValues(Vector1);
    for (int i = 0; i < 5; ++i) Values1.emplace_back(i);
    vector_noncopyable Vector2;
    Vector2 = std::move(Vector1);
    auto &Values2 = helper_noncopyable::GetValues(Vector2);
    EXPECT_TRUE(Values1.empty());
    int Sum = 0;
    for (auto &Value : Values2) Sum += Value.Value();
    EXPECT_EQ(Sum, 10);
  }

  // Non-default-contructible, move construct
  {
    vector_nondefaultconstructible Vector1 = {{0},{1},{2},{3},{4}};
    vector_nondefaultconstructible Vector2(std::move(Vector1));
    auto &Values1 = helper_nondefaultconstructible::GetValues(Vector1);
    auto &Values2 = helper_nondefaultconstructible::GetValues(Vector2);
    EXPECT_TRUE(Values1.empty());
    int Sum = 0;
    for (auto &Value : Values2) Sum += Value.Value();
    EXPECT_EQ(Sum, 10);
  }

  // Non-default-contructible, move assign
  {
    vector_nondefaultconstructible Vector1 = {{0},{1},{2},{3},{4}};
    vector_nondefaultconstructible Vector2;
    Vector2 = std::move(Vector1);
    auto &Values1 = helper_nondefaultconstructible::GetValues(Vector1);
    auto &Values2 = helper_nondefaultconstructible::GetValues(Vector2);
    EXPECT_TRUE(Values1.empty());
    int Sum = 0;
    for (auto &Value : Values2) Sum += Value.Value();
    EXPECT_EQ(Sum, 10);
  }

}

TEST_F(VectorTests, Assign) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using vector_noncopyable = ovk::core::vector<noncopyable<int>>;
  using vector_nondefaultconstructible = ovk::core::vector<nondefaultconstructible<int>>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;
  using helper_noncopyable = ovk::core::test_helper<vector_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<vector_nondefaultconstructible>;

  // Non-bool, size and value
  {
    vector_nonbool Vector(1, 0);
    Vector.Assign(5, 1);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1}));
  }

  // Non-bool, initializer list, operator=
  {
    vector_nonbool Vector(1, 0);
    Vector = {0,1,2,3,4};
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4}));
  }

  // Non-bool, initializer list, Assign
  {
    vector_nonbool Vector(1, 0);
    Vector.Assign({0,1,2,3,4});
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4}));
  }

  // Bool, size and value
  {
    vector_bool Vector(1, false);
    Vector.Assign(5, true);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::TRUE,not_bool::TRUE,not_bool::TRUE,
      not_bool::TRUE,not_bool::TRUE}));
  }

  // Bool, initializer list, operator=
  {
    vector_bool Vector(1, false);
    Vector = {false,true,false,false,true};
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::FALSE,not_bool::TRUE,not_bool::FALSE,
      not_bool::FALSE,not_bool::TRUE}));
  }

  // Bool, initializer list, Assign
  {
    vector_bool Vector(1, false);
    Vector.Assign({false,true,false,false,true});
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::FALSE,not_bool::TRUE,not_bool::FALSE,
      not_bool::FALSE,not_bool::TRUE}));
  }

  // Non-copyable, iterators
  {
    vector_noncopyable Vector(1);
    std::array<noncopyable<int>,2> SourceValues = {{{1}, {2}}};
    Vector.Assign(std::make_move_iterator(SourceValues.begin()), std::make_move_iterator(
      SourceValues.end()));
    auto &Values = helper_noncopyable::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
  }

  // Non-default-constructible, size and value
  {
    vector_nondefaultconstructible Vector(1, {0});
    Vector.Assign(5, {1});
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 5);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 5);
  }

  // Non-default-constructible, initializer list, operator=
  {
    vector_nondefaultconstructible Vector(1, {0});
    Vector = {{1}, {2}};
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
  }

  // Non-default-constructible, initializer list, Assign
  {
    vector_nondefaultconstructible Vector(1, {0});
    Vector.Assign({{1}, {2}});
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
  }

  // Non-default-constructible, iterators
  {
    vector_nondefaultconstructible Vector(1, {0});
    std::array<nondefaultconstructible<int>,2> SourceArray = {{{1}, {2}}};
    Vector.Assign(SourceArray.begin(), SourceArray.end());
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
  }

}

TEST_F(VectorTests, Reserve) {

  if (TestComm().Rank() != 0) return;

  using vector = ovk::core::vector<int>;
  using helper = ovk::core::test_helper<vector>;

  vector Vector;
  auto &Values = helper::GetValues(Vector);
  long long Capacity = (long long)(Values.capacity());
  Vector.Reserve(Capacity+1);
  EXPECT_GE((long long)(Values.capacity()), Capacity+1);

}

TEST_F(VectorTests, Resize) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using vector_noncopyable = ovk::core::vector<noncopyable<int>>;
  using vector_nondefaultconstructible = ovk::core::vector<nondefaultconstructible<int>>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;
  using helper_noncopyable = ovk::core::test_helper<vector_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<vector_nondefaultconstructible>;

  // Non-bool, size
  {
    vector_nonbool Vector;
    Vector.Resize(5);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 5);
  }

  // Non-bool, size and value
  {
    vector_nonbool Vector;
    Vector.Resize(5, 1);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1}));
  }

  // Bool, size
  {
    vector_bool Vector;
    Vector.Resize(5);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 5);
  }

  // Bool, size and value
  {
    vector_bool Vector;
    Vector.Resize(5, true);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::TRUE,not_bool::TRUE,not_bool::TRUE,
      not_bool::TRUE,not_bool::TRUE}));
  }

  // Non-copyable, size
  {
    vector_noncopyable Vector;
    Vector.Resize(5);
    auto &Values = helper_noncopyable::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 5);
  }

  // Non-default-constructible, size and value
  {
    vector_nondefaultconstructible Vector;
    Vector.Resize(5, 1);
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 5);
  }

}

TEST_F(VectorTests, Clear) {

  if (TestComm().Rank() != 0) return;

  using vector = ovk::core::vector<int>;
  using helper = ovk::core::test_helper<vector>;

  vector Vector(5);
  Vector.Clear();
  auto &Values = helper::GetValues(Vector);
  EXPECT_EQ(int(Values.size()), 0);

}

TEST_F(VectorTests, Append) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using vector_noncopyable = ovk::core::vector<noncopyable<int>>;
  using vector_multiargument = ovk::core::vector<multiargument>;
  using vector_nondefaultconstructible = ovk::core::vector<nondefaultconstructible<int>>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;
  using helper_noncopyable = ovk::core::test_helper<vector_noncopyable>;
  using helper_multiargument = ovk::core::test_helper<vector_multiargument>;
  using helper_nondefaultconstructible = ovk::core::test_helper<vector_nondefaultconstructible>;

  // Non-bool, lvalue ref
  {
    vector_nonbool Vector = {0,1,2,3};
    int SourceValue = 4;
    int &NewValue = Vector.Append(SourceValue);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4}));
    EXPECT_EQ(&NewValue, Values.data()+4);
  }

  // Bool, lvalue ref
  {
    vector_bool Vector = {false,true,false,false};
    bool SourceValue = true;
    bool &NewValue = Vector.Append(SourceValue);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::FALSE,not_bool::TRUE,not_bool::FALSE,
      not_bool::FALSE,not_bool::TRUE}));
    EXPECT_EQ(&NewValue, reinterpret_cast<bool *>(Values.data()+4));
  }

  // rvalue ref
  {
    vector_noncopyable Vector;
    noncopyable<int> Value(1);
    noncopyable<int> &NewValue = Vector.Append(std::move(Value));
    auto &Values = helper_noncopyable::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 1);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 1);
    EXPECT_EQ(&NewValue, Values.data());
  }

  // In-place
  {
    vector_multiargument Vector;
    multiargument &NewValue = Vector.Append(1, 2);
    auto &Values = helper_multiargument::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 1);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.v1 + Value.v2;
    EXPECT_EQ(Sum, 3);
    EXPECT_EQ(&NewValue, Values.data());
  }

  // Non-default-constructible, lvalue ref
  {
    vector_nondefaultconstructible Vector;
    nondefaultconstructible<int> SourceValue(1);
    nondefaultconstructible<int> &NewValue = Vector.Append(SourceValue);
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 1);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 1);
    EXPECT_EQ(&NewValue, Values.data());
  }

  // Non-default-constructible, rvalue ref
  {
    vector_nondefaultconstructible Vector;
    nondefaultconstructible<int> SourceValue(1);
    nondefaultconstructible<int> &NewValue = Vector.Append(std::move(SourceValue));
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 1);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 1);
    EXPECT_EQ(&NewValue, Values.data());
  }

  // Non-default-constructible, in-place
  {
    vector_nondefaultconstructible Vector;
    nondefaultconstructible<int> &NewValue = Vector.Append(1);
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 1);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 1);
    EXPECT_EQ(&NewValue, Values.data());
  }

}

TEST_F(VectorTests, Insert) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using vector_noncopyable = ovk::core::vector<noncopyable<int>>;
  using vector_multiargument = ovk::core::vector<multiargument>;
  using vector_nondefaultconstructible = ovk::core::vector<nondefaultconstructible<int>>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;
  using helper_noncopyable = ovk::core::test_helper<vector_noncopyable>;
  using helper_multiargument = ovk::core::test_helper<vector_multiargument>;
  using helper_nondefaultconstructible = ovk::core::test_helper<vector_nondefaultconstructible>;

  // Index, lvalue ref, non-bool
  {
    vector_nonbool Vector = {0,1,3,4};
    int SourceValue = 2;
    auto &NewValue = Vector.Insert(2, SourceValue);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4}));
    EXPECT_EQ(&NewValue, Values.data()+2);
  }

  // Index, lvalue ref, bool
  {
    vector_bool Vector = {false,true,false,true};
    bool SourceValue = false;
    auto &NewValue = Vector.Insert(2, SourceValue);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::FALSE,not_bool::TRUE,not_bool::FALSE,
      not_bool::FALSE,not_bool::TRUE}));
    EXPECT_EQ(&NewValue, reinterpret_cast<bool *>(Values.data()+2));
  }

  // Index, rvalue ref
  {
    vector_noncopyable Vector;
    Vector.Append(2);
    noncopyable<int> Value(1);
    auto &NewValue = Vector.Insert(0, std::move(Value));
    auto &Values = helper_noncopyable::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
    EXPECT_EQ(&NewValue, Values.data());
  }

  // Index, in-place
  {
    vector_multiargument Vector = {{3,4}};
    auto &NewValue = Vector.Insert(0, 1, 2);
    auto &Values = helper_multiargument::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.v1 + Value.v2;
    EXPECT_EQ(Sum, 10);
    EXPECT_EQ(&NewValue, Values.data());
  }

  // Iterator, lvalue ref, non-bool
  {
    vector_nonbool Vector = {0,1,3,4};
    int SourceValue = 2;
    auto Iter = Vector.Insert(Vector.Begin()+2, SourceValue);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4}));
    EXPECT_EQ(Iter, Vector.Begin()+2);
  }

  // Iterator, lvalue ref, bool
  {
    vector_bool Vector = {false,true,false,true};
    bool SourceValue = false;
    auto Iter = Vector.Insert(Vector.Begin()+2, SourceValue);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::FALSE,not_bool::TRUE,not_bool::FALSE,
      not_bool::FALSE,not_bool::TRUE}));
    EXPECT_EQ(Iter, Vector.Begin()+2);
  }

  // Iterator, rvalue ref
  {
    vector_noncopyable Vector;
    Vector.Append(2);
    noncopyable<int> Value(1);
    auto Iter = Vector.Insert(Vector.Begin(), std::move(Value));
    auto &Values = helper_noncopyable::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
    EXPECT_EQ(Iter, Vector.Begin());
  }

  // Iterator, in-place
  {
    vector_multiargument Vector = {{3,4}};
    auto Iter = Vector.Insert(Vector.Begin(), 1, 2);
    auto &Values = helper_multiargument::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.v1 + Value.v2;
    EXPECT_EQ(Sum, 10);
    EXPECT_EQ(Iter, Vector.Begin());
  }

  // Non-default-constructible, index, lvalue ref
  {
    vector_nondefaultconstructible Vector = {{2}};
    nondefaultconstructible<int> SourceValue(1);
    auto &NewValue = Vector.Insert(0, SourceValue);
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
    EXPECT_EQ(&NewValue, Values.data());
  }

  // Non-default-constructible, index, rvalue ref
  {
    vector_nondefaultconstructible Vector = {{2}};
    nondefaultconstructible<int> SourceValue(1);
    auto &NewValue = Vector.Insert(0, std::move(SourceValue));
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
    EXPECT_EQ(&NewValue, Values.data());
  }

  // Non-default-constructible, index, in-place
  {
    vector_nondefaultconstructible Vector = {{2}};
    auto &NewValue = Vector.Insert(0, 1);
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
    EXPECT_EQ(&NewValue, Values.data());
  }

  // Non-default-constructible, iterator, lvalue ref
  {
    vector_nondefaultconstructible Vector = {{2}};
    nondefaultconstructible<int> SourceValue(1);
    auto Iter = Vector.Insert(Vector.Begin(), SourceValue);
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
    EXPECT_EQ(Iter, Vector.Begin());
  }

  // Non-default-constructible, iterator, rvalue ref
  {
    vector_nondefaultconstructible Vector = {{2}};
    nondefaultconstructible<int> SourceValue(1);
    auto Iter = Vector.Insert(Vector.Begin(), std::move(SourceValue));
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
    EXPECT_EQ(Iter, Vector.Begin());
  }

  // Non-default-constructible, iterator, in-place
  {
    vector_nondefaultconstructible Vector = {{2}};
    auto Iter = Vector.Insert(Vector.Begin(), 1);
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 2);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
    EXPECT_EQ(Iter, Vector.Begin());
  }

}

TEST_F(VectorTests, Erase) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using vector_noncopyable = ovk::core::vector<noncopyable<int>>;
  using vector_nondefaultconstructible = ovk::core::vector<nondefaultconstructible<int>>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;
  using helper_noncopyable = ovk::core::test_helper<vector_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<vector_nondefaultconstructible>;

  // Index, non-bool
  {
    vector_nonbool Vector = {0,1,2,3,4};
    Vector.Erase(2);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({0,1,3,4}));
  }

  // Index, bool
  {
    vector_bool Vector = {false,true,false,false,true};
    Vector.Erase(2);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::FALSE,not_bool::TRUE,not_bool::FALSE,
      not_bool::TRUE}));
  }

  // Iterator, non-bool
  {
    vector_nonbool Vector = {0,1,2,3,4};
    auto Iter = Vector.Erase(Vector.Begin()+2);
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({0,1,3,4}));
    EXPECT_EQ(Iter, Vector.Begin()+2);
  }

  // Iterator, bool
  {
    vector_bool Vector = {false,true,false,false,true};
    auto Iter = Vector.Erase(Vector.Begin()+2);
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_THAT(Values, ElementsAreArray({not_bool::FALSE,not_bool::TRUE,not_bool::FALSE,
      not_bool::TRUE}));
    EXPECT_EQ(Iter, Vector.Begin()+2);
  }

  // Non-copyable, index
  {
    vector_noncopyable Vector;
    Vector.Append(1);
    Vector.Append(2);
    Vector.Erase(0);
    auto &Values = helper_noncopyable::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 1);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 2);
  }

  // Non-copyable, iterator
  {
    vector_noncopyable Vector;
    Vector.Append(1);
    Vector.Append(2);
    auto Iter = Vector.Erase(Vector.Begin());
    auto &Values = helper_noncopyable::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 1);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 2);
    EXPECT_EQ(Iter, Vector.Begin());
  }

  // Non-default-constructible, index
  {
    vector_nondefaultconstructible Vector = {{1},{2}};
    Vector.Erase(0);
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 1);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 2);
  }

  // Non-default-constructible, iterator
  {
    vector_nondefaultconstructible Vector = {{1},{2}};
    auto Iter = Vector.Erase(Vector.Begin());
    auto &Values = helper_nondefaultconstructible::GetValues(Vector);
    EXPECT_EQ(int(Values.size()), 1);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 2);
    EXPECT_EQ(Iter, Vector.Begin());
  }

}

TEST_F(VectorTests, Count) {

  if (TestComm().Rank() != 0) return;

  using vector = ovk::core::vector<int>;

  vector Vector(5);
  EXPECT_EQ(Vector.Count(), 5);

}

TEST_F(VectorTests, Capacity) {

  if (TestComm().Rank() != 0) return;

  using vector = ovk::core::vector<int>;
  using helper = ovk::core::test_helper<vector>;

  vector Vector;
  auto &Values = helper::GetValues(Vector);
  Values.reserve(5);
  EXPECT_EQ(Vector.Capacity(), (long long)(Values.capacity()));

}

TEST_F(VectorTests, Empty) {

  if (TestComm().Rank() != 0) return;

  using vector = ovk::core::vector<int>;

  // Empty
  {
    vector Vector;
    EXPECT_TRUE(Vector.Empty());
  }

  // Non-empty
  {
    vector Vector(5);
    EXPECT_FALSE(Vector.Empty());
  }

}

TEST_F(VectorTests, BracketOperator) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;

  // Non-bool, const
  {
    const vector_nonbool Vector(5);
    const int &Value = Vector[2];
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(&Value, Values.data()+2);
  }

  // Non-bool, non-const
  {
    vector_nonbool Vector(5);
    int &Value = Vector[2];
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(&Value, Values.data()+2);
  }

  // Bool, const
  {
    const vector_bool Vector(5);
    const bool &Value = Vector[2];
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(&Value, reinterpret_cast<const bool *>(Values.data()+2));
  }

  // Bool, non-const
  {
    vector_bool Vector(5);
    bool &Value = Vector[2];
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(&Value, reinterpret_cast<bool *>(Values.data()+2));
  }

}

TEST_F(VectorTests, Front) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;

  // Non-bool, const
  {
    const vector_nonbool Vector(5);
    const int &Value = Vector.Front();
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(&Value, Values.data());
  }

  // Non-bool, non-const
  {
    vector_nonbool Vector(5);
    int &Value = Vector.Front();
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(&Value, Values.data());
  }

  // Bool, const
  {
    const vector_bool Vector(5);
    const bool &Value = Vector.Front();
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(&Value, reinterpret_cast<const bool *>(Values.data()));
  }

  // Bool, non-const
  {
    vector_bool Vector(5);
    bool &Value = Vector.Front();
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(&Value, reinterpret_cast<bool *>(Values.data()));
  }

}

TEST_F(VectorTests, Back) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;

  // Non-bool, const
  {
    const vector_nonbool Vector(5);
    const int &Value = Vector.Back();
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(&Value, Values.data()+4);
  }

  // Non-bool, non-const
  {
    vector_nonbool Vector(5);
    int &Value = Vector.Back();
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(&Value, Values.data()+4);
  }

  // Bool, const
  {
    const vector_bool Vector(5);
    const bool &Value = Vector.Back();
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(&Value, reinterpret_cast<const bool *>(Values.data()+4));
  }

  // Bool, non-const
  {
    vector_bool Vector(5);
    bool &Value = Vector.Back();
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(&Value, reinterpret_cast<bool *>(Values.data()+4));
  }

}

TEST_F(VectorTests, Data) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;

  // Non-bool, const
  {
    const vector_nonbool Vector(5);
    const int *Data = Vector.Data();
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(Data, Values.data());
  }

  // Non-bool, non-const
  {
    vector_nonbool Vector(5);
    int *Data = Vector.Data();
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(Data, Values.data());
  }

  // Bool, const
  {
    const vector_bool Vector(5);
    const bool *Data = Vector.Data();
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(Data, reinterpret_cast<const bool *>(Values.data()));
  }

  // Bool, non-const
  {
    vector_bool Vector(5);
    bool *Data = Vector.Data();
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(Data, reinterpret_cast<bool *>(Values.data()));
  }

}

TEST_F(VectorTests, BeginEnd) {

  if (TestComm().Rank() != 0) return;

  using vector_nonbool = ovk::core::vector<int>;
  using vector_bool = ovk::core::vector<bool>;
  using helper_nonbool = ovk::core::test_helper<vector_nonbool>;
  using helper_bool = ovk::core::test_helper<vector_bool>;

  // Non-bool, const
  {
    const vector_nonbool Vector = {0,1,2,3,4};
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(Vector.Begin().Pointer(), Values.data());
    EXPECT_EQ(Vector.End().Pointer(), Values.data()+5);
    int Sum = 0;
    for (auto &Value : Vector) Sum += Value;
    EXPECT_EQ(Sum, 10);
  }

  // Non-bool, non-const
  {
    vector_nonbool Vector = {0,1,2,3,4};
    auto &Values = helper_nonbool::GetValues(Vector);
    EXPECT_EQ(Vector.Begin().Pointer(), Values.data());
    EXPECT_EQ(Vector.End().Pointer(), Values.data()+5);
    for (auto &Value : Vector) Value = 1;
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1}));
  }

  // Bool, const
  {
    const vector_bool Vector = {false,true,false,false,true};
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(Vector.Begin().Pointer(), reinterpret_cast<const bool *>(Values.data()));
    EXPECT_EQ(Vector.End().Pointer(), reinterpret_cast<const bool *>(Values.data()+5));
    int Sum = 0;
    for (auto &Value : Vector) Sum += int(Value);
    EXPECT_EQ(Sum, 2);
  }

  // Bool, non-const
  {
    vector_bool Vector = {false,true,false,false,true};
    auto &Values = helper_bool::GetValues(Vector);
    EXPECT_EQ(Vector.Begin().Pointer(), reinterpret_cast<bool *>(Values.data()));
    EXPECT_EQ(Vector.End().Pointer(), reinterpret_cast<bool *>(Values.data()+5));
    for (auto &Value : Vector) Value = false;
    EXPECT_THAT(Values, ElementsAreArray({not_bool::FALSE,not_bool::FALSE,not_bool::FALSE,
      not_bool::FALSE,not_bool::FALSE}));
  }

}

TEST_F(VectorTests, Traits) {

  if (TestComm().Rank() != 0) return;

  using vector = ovk::core::vector<int>;

  EXPECT_TRUE(ovk::core::IsArray<vector>());
  EXPECT_TRUE((std::is_same<ovk::core::array_value_type<vector>, int>::value));
  EXPECT_EQ(ovk::core::ArrayRank<vector>(), 1);
  EXPECT_EQ(ovk::core::ArrayLayout<vector>(), ovk::array_layout::ROW_MAJOR);
  EXPECT_TRUE(ovk::core::ArrayHasRuntimeExtents<vector>());

  vector Vector = {0,1,2,3,4};
  EXPECT_THAT(ovk::core::ArrayExtents(Vector).Begin(), ElementsAre(0));
  EXPECT_THAT(ovk::core::ArrayExtents(Vector).End(), ElementsAre(5));
  EXPECT_EQ(ovk::core::ArrayData(Vector), Vector.Data());

}
