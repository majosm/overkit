// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Array.hpp>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "MPITest.hpp"
#include "mocks/MultidimArray.hpp"
#include "mocks/Noncopyable.hpp"
#include "mocks/Nondefaultconstructible.hpp"

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Interval.hpp>

#include <mpi.h>

#include <array>
#include <initializer_list>
#include <type_traits>
#include <utility>
#include <vector>

using testing::ElementsAre;
using testing::ElementsAreArray;

class ArrayTests : public tests::mpi_test {};

using tests::multidim_array_row;
using tests::multidim_array_col;
using tests::noncopyable;
using tests::nondefaultconstructible;

namespace ovk {
namespace core {
template <typename T, int Rank, array_layout Layout> class test_helper<array<T,Rank,Layout>> {
public:
  using array_type = array<T,Rank,Layout>;
  using view_type = array_view<T,Rank,Layout>;
  static const vector<T> &GetValues(const array_type &Array) { return Array.Values_; }
  static vector<T> &GetValues(array_type &Array) { return Array.Values_; }
  static const view_type &GetView(const array_type &Array) { return Array.View_; }
  static view_type &GetView(array_type &Array) { return Array.View_; }
};
}}

TEST_F(ArrayTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;

  EXPECT_TRUE((std::is_same<typename array::value_type, int>::value));
  EXPECT_EQ(int(array::Rank), 3);
  EXPECT_EQ(ovk::array_layout(array::Layout), ovk::array_layout::ROW_MAJOR);
  EXPECT_TRUE((std::is_same<typename array::index_type, long long>::value));
  EXPECT_TRUE((std::is_same<typename array::tuple_element_type, long long>::value));
  EXPECT_TRUE((std::is_same<typename array::tuple_type, ovk::elem<long long,3>>::value));
  EXPECT_TRUE((std::is_same<typename array::interval_type, ovk::interval<long long,3>>::value));
  EXPECT_TRUE((std::is_same<typename array::indexer_type, ovk::indexer<long long, long long,
    3>>::value));
  EXPECT_TRUE((std::is_same<typename array::view_type, ovk::array_view<int,3>>::value));
  EXPECT_TRUE((std::is_same<typename array::const_view_type, ovk::array_view<const int,3>>::value));
  EXPECT_TRUE((std::is_same<typename array::iterator, int *>::value));
  EXPECT_TRUE((std::is_same<typename array::const_iterator, const int *>::value));

  using array_col = ovk::array<int,3,ovk::array_layout::COLUMN_MAJOR>;

  EXPECT_TRUE((std::is_same<typename array_col::indexer_type, ovk::indexer<long long, long long,
    3, ovk::array_layout::COLUMN_MAJOR>>::value));
  EXPECT_TRUE((std::is_same<typename array_col::view_type, ovk::array_view<int,3,
    ovk::array_layout::COLUMN_MAJOR>>::value));
  EXPECT_TRUE((std::is_same<typename array_col::const_view_type, ovk::array_view<const int,3,
    ovk::array_layout::COLUMN_MAJOR>>::value));
  EXPECT_EQ(ovk::array_layout(array_col::Layout), ovk::array_layout::COLUMN_MAJOR);

}

TEST_F(ArrayTests, Create) {

  if (TestComm().Rank() != 0) return;

  using array_1d = ovk::array<int>;
  using array_row = ovk::array<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_col = ovk::array<int,3,ovk::array_layout::COLUMN_MAJOR>;
  using array_noncopyable = ovk::array<noncopyable<int>>;
  using array_nondefaultconstructible = ovk::array<nondefaultconstructible<int>>;
  using array_view_1d = ovk::array_view<int>;
  using array_view_row = ovk::array_view<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_view_col = ovk::array_view<int,3,ovk::array_layout::COLUMN_MAJOR>;
  using array_view_nondefaultconstructible = ovk::array_view<nondefaultconstructible<int>>;
  using const_array_view_1d = ovk::array_view<const int>;
  using const_array_view_row = ovk::array_view<const int,3,ovk::array_layout::ROW_MAJOR>;
  using const_array_view_col = ovk::array_view<const int,3,ovk::array_layout::COLUMN_MAJOR>;
  using helper_1d = ovk::core::test_helper<array_1d>;
  using helper_row = ovk::core::test_helper<array_row>;
  using helper_col = ovk::core::test_helper<array_col>;
  using helper_noncopyable = ovk::core::test_helper<array_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<array_nondefaultconstructible>;

  // One-dimensional, default
  {
    array_1d Array;
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(0));
    EXPECT_EQ(Values.Count(), 0);
  }

  // One-dimensional, size
  {
    array_1d Array({4});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_EQ(Values.Count(), 4);
  }

  // One-dimensional, interval
  {
    array_1d Array({-2, 2});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_EQ(Values.Count(), 4);
  }

  // One-dimensional, size and value
  {
    array_1d Array({4}, 1);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1}));
  }

  // One-dimensional, interval and value
  {
    array_1d Array({-2, 2}, 1);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1}));
  }

  // One-dimensional, size and initializer list
  {
    array_1d Array({4}, {0,1,2,3});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and initializer list
  {
    array_1d Array({-2, 2}, {0,1,2,3});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, size and iterator
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({4}, SourceArray.begin());
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and iterator
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({-2, 2}, SourceArray.begin());
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array = array_view_1d(SourceArray, {-2, 2});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, size and array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({4}, array_view_1d(SourceArray, {-2, 2}));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({-2, 2}, array_view_1d(SourceArray));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, const array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array = const_array_view_1d(SourceArray, {-2, 2});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, size and const array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({4}, const_array_view_1d(SourceArray, {-2, 2}));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and const array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({-2, 2}, const_array_view_1d(SourceArray));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, array
  {
    const std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array = SourceArray;
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, size and array
  {
    const std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({4}, SourceArray);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and array
  {
    const std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({-2, 2}, SourceArray);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // Multidimensional, row major, default
  {
    array_row Array;
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(0,0,0));
    EXPECT_EQ(Values.Count(), 0);
  }

  // Multidimensional, row major, size
  {
    array_row Array({{1,2,3}});
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_EQ(Values.Count(), 6);
  }

  // Multidimensional, row major, interval
  {
    array_row Array({{1,2,3}, {2,4,6}});
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_EQ(Values.Count(), 6);
  }

  // Multidimensional, row major, size and value
  {
    array_row Array({{1,2,3}}, 1);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  }

  // Multidimensional, row major, interval and value
  {
    array_row Array({{1,2,3}, {2,4,6}}, 1);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  }

  // Multidimensional, row major, size and initializer list
  {
    array_row Array({{1,2,3}}, {0,1,2,3,4,5});
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and initializer list
  {
    array_row Array({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, size and iterator
  {
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    array_row Array({{1,2,3}}, SourceArray.begin());
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and iterator
  {
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    array_row Array({{1,2,3}, {2,4,6}}, SourceArray.begin());
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array = array_view_row(SourceArray);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, size and array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,2,3}}, array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_row Array({{1,2,3}, {2,4,6}}, array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, const array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array = const_array_view_row(SourceArray);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, size and const array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,2,3}}, const_array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and const array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_row Array({{1,2,3}, {2,4,6}}, const_array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, array
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array = SourceArray;
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, size and array
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,2,3}}, SourceArray);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and array
  {
    multidim_array_row<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_row Array({{1,2,3}, {2,4,6}}, SourceArray);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, default
  {
    array_col Array;
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(0,0,0));
    EXPECT_EQ(Values.Count(), 0);
  }

  // Multidimensional, column major, size
  {
    array_col Array({{1,2,3}});
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_EQ(Values.Count(), 6);
  }

  // Multidimensional, column major, interval
  {
    array_col Array({{1,2,3}, {2,4,6}});
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_EQ(Values.Count(), 6);
  }

  // Multidimensional, column major, size and value
  {
    array_col Array({{1,2,3}}, 1);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  }

  // Multidimensional, column major, interval and value
  {
    array_col Array({{1,2,3}, {2,4,6}}, 1);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  }

  // Multidimensional, column major, size and initializer list
  {
    array_col Array({{1,2,3}}, {0,1,2,3,4,5});
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and initializer list
  {
    array_col Array({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, size and iterator
  {
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    array_col Array({{1,2,3}}, SourceArray.begin());
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and iterator
  {
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    array_col Array({{1,2,3}, {2,4,6}}, SourceArray.begin());
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array = array_view_col(SourceArray);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, size and array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,2,3}}, array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_col Array({{1,2,3}, {2,4,6}}, array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, const array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array = const_array_view_col(SourceArray);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, size and const array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,2,3}}, const_array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and const array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_col Array({{1,2,3}, {2,4,6}}, const_array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, array
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array = SourceArray;
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, size and array
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,2,3}}, SourceArray);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and array
  {
    multidim_array_col<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_col Array({{1,2,3}, {2,4,6}}, SourceArray);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Non-copyable, default
  {
    array_noncopyable Array;
    auto &View = helper_noncopyable::GetView(Array);
    auto &Values = helper_noncopyable::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(0));
    EXPECT_EQ(Values.Count(), 0);
  }

  // Non-copyable, rvalue array
  {
    std::array<noncopyable<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_noncopyable Array = std::move(SourceArray);
    auto &View = helper_noncopyable::GetView(Array);
    auto &Values = helper_noncopyable::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_EQ(Values.Count(), 4);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
    Sum = 0;
    for (auto &Value : SourceArray) Sum += Value.Value();
    EXPECT_EQ(Sum, 0);
  }


  // Non-copyable, interval and rvalue array
  {
    std::array<noncopyable<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_noncopyable Array({-2, 2}, std::move(SourceArray));
    auto &View = helper_noncopyable::GetView(Array);
    auto &Values = helper_noncopyable::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_EQ(Values.Count(), 4);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
    Sum = 0;
    for (auto &Value : SourceArray) Sum += Value.Value();
    EXPECT_EQ(Sum, 0);
  }

  // Non-default-constructible, interval and value
  {
    array_nondefaultconstructible Array({-2, 2}, 1);
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 4);
  }

  // Non-default-constructible, interval and iterator
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({-2, 2}, SourceArray.begin());
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, array view
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array = array_view_nondefaultconstructible(SourceArray, {-2, 2});
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, interval and array view
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({-2, 2}, array_view_nondefaultconstructible(SourceArray));
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, array
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array = SourceArray;
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, interval and array
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({-2, 2}, SourceArray);
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

}

TEST_F(ArrayTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;
  using helper = ovk::core::test_helper<array>;

  // Copy construct
  {
    array Array1({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array Array2 = Array1;
    auto &View = helper::GetView(Array2);
    auto &Values = helper::GetValues(Array2);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Copy construct with size
  {
    array Array1({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array Array2({{1,2,3}}, Array1);
    auto &View = helper::GetView(Array2);
    auto &Values = helper::GetValues(Array2);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Copy construct with interval
  {
    array Array1({{1,2,3}}, {0,1,2,3,4,5});
    array Array2({{1,2,3}, {2,4,6}}, Array1);
    auto &View = helper::GetView(Array2);
    auto &Values = helper::GetValues(Array2);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Copy assign
  {
    array Array1({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array Array2;
    Array2 = Array1;
    auto &View = helper::GetView(Array2);
    auto &Values = helper::GetValues(Array2);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

}

TEST_F(ArrayTests, Move) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;
  using helper = ovk::core::test_helper<array>;

  // Move construct
  {
    array Array1({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array Array2 = std::move(Array1);
    auto &View1 = helper::GetView(Array1);
    auto &Values1 = helper::GetValues(Array1);
    auto &View2 = helper::GetView(Array2);
    auto &Values2 = helper::GetValues(Array2);
    EXPECT_EQ(View2.Data(), Values2.Data());
    EXPECT_THAT(View2.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View2.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values2, ElementsAreArray({0,1,2,3,4,5}));
    EXPECT_FALSE(static_cast<bool>(View1));
    EXPECT_TRUE(Values1.Empty());
  }

  // Move construct with size
  {
    array Array1({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array Array2({{1,2,3}}, std::move(Array1));
    auto &View1 = helper::GetView(Array1);
    auto &Values1 = helper::GetValues(Array1);
    auto &View2 = helper::GetView(Array2);
    auto &Values2 = helper::GetValues(Array2);
    EXPECT_EQ(View2.Data(), Values2.Data());
    EXPECT_THAT(View2.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View2.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values2, ElementsAreArray({0,1,2,3,4,5}));
    EXPECT_FALSE(static_cast<bool>(View1));
    EXPECT_TRUE(Values1.Empty());
  }

  // Move construct with interval
  {
    array Array1({{1,2,3}}, {0,1,2,3,4,5});
    array Array2({{1,2,3}, {2,4,6}}, std::move(Array1));
    auto &View1 = helper::GetView(Array1);
    auto &Values1 = helper::GetValues(Array1);
    auto &View2 = helper::GetView(Array2);
    auto &Values2 = helper::GetValues(Array2);
    EXPECT_EQ(View2.Data(), Values2.Data());
    EXPECT_THAT(View2.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View2.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values2, ElementsAreArray({0,1,2,3,4,5}));
    EXPECT_FALSE(static_cast<bool>(View1));
    EXPECT_TRUE(Values1.Empty());
  }

  // Move assign
  {
    array Array1({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array Array2;
    Array2 = std::move(Array1);
    auto &View1 = helper::GetView(Array1);
    auto &Values1 = helper::GetValues(Array1);
    auto &View2 = helper::GetView(Array2);
    auto &Values2 = helper::GetValues(Array2);
    EXPECT_EQ(View2.Data(), Values2.Data());
    EXPECT_THAT(View2.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View2.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values2, ElementsAreArray({0,1,2,3,4,5}));
    EXPECT_FALSE(static_cast<bool>(View1));
    EXPECT_TRUE(Values1.Empty());
  }

}

TEST_F(ArrayTests, Assign) {

  if (TestComm().Rank() != 0) return;

  using array_1d = ovk::array<int>;
  using array_row = ovk::array<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_col = ovk::array<int,3,ovk::array_layout::COLUMN_MAJOR>;
  using array_noncopyable = ovk::array<noncopyable<int>>;
  using array_nondefaultconstructible = ovk::array<nondefaultconstructible<int>>;
  using array_view_1d = ovk::array_view<int>;
  using array_view_row = ovk::array_view<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_view_col = ovk::array_view<int,3,ovk::array_layout::COLUMN_MAJOR>;
  using array_view_nondefaultconstructible = ovk::array_view<nondefaultconstructible<int>,1>;
  using const_array_view_1d = ovk::array_view<const int>;
  using const_array_view_row = ovk::array_view<const int,3,ovk::array_layout::ROW_MAJOR>;
  using const_array_view_col = ovk::array_view<const int,3,ovk::array_layout::COLUMN_MAJOR>;
  using helper_1d = ovk::core::test_helper<array_1d>;
  using helper_row = ovk::core::test_helper<array_row>;
  using helper_col = ovk::core::test_helper<array_col>;
  using helper_noncopyable = ovk::core::test_helper<array_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<array_nondefaultconstructible>;

  // One-dimensional, size and value
  {
    array_1d Array({1}, 0);
    Array.Assign({4}, 1);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1}));
  }

  // One-dimensional, interval and value
  {
    array_1d Array({1}, 0);
    Array.Assign({-2, 2}, 1);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1}));
  }

  // One-dimensional, size and initializer list
  {
    array_1d Array({1}, 0);
    Array.Assign({4}, {0,1,2,3});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and initializer list
  {
    array_1d Array({1}, 0);
    Array.Assign({-2, 2}, {0,1,2,3});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, size and iterator
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign({4}, SourceArray.begin());
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and iterator
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign({-2, 2}, SourceArray.begin());
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, array view, operator=
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array = array_view_1d(SourceArray, {-2, 2});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, array view, Assign
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign(array_view_1d(SourceArray, {-2, 2}));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, size and array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign({4}, array_view_1d(SourceArray, {-2, 2}));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign({-2, 2}, array_view_1d(SourceArray));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, const array view, operator=
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array = const_array_view_1d(SourceArray, {-2, 2});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, const array view, Assign
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign(const_array_view_1d(SourceArray, {-2, 2}));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, size and const array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign({4}, const_array_view_1d(SourceArray, {-2, 2}));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and const array view
  {
    std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign({-2, 2}, const_array_view_1d(SourceArray));
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, array, operator=
  {
    const std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array = SourceArray;
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, array, Assign
  {
    const std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign(SourceArray);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, size and array
  {
    const std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign({4}, SourceArray);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // One-dimensional, interval and array
  {
    const std::array<int,4> SourceArray = {{0,1,2,3}};
    array_1d Array({1}, 0);
    Array.Assign({-2, 2}, SourceArray);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3}));
  }

  // Multidimensional, row major, size and value
  {
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, 1);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  }

  // Multidimensional, row major, interval and value
  {
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, 1);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  }

  // Multidimensional, row major, size and initializer list
  {
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, {0,1,2,3,4,5});
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and initializer list
  {
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, size and iterator
  {
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, SourceArray.begin());
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and iterator
  {
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, SourceArray.begin());
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, array view, operator=
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array = array_view_row(SourceArray);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, array view, Assign
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array.Assign(array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, size and array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, const array view, operator=
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array = const_array_view_row(SourceArray);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, const array view, Assign
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array.Assign(const_array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, size and const array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, const_array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and const array view
  {
    multidim_array_row<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, const_array_view_row(SourceArray));
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, array, operator=
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array = SourceArray;
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, array, Assign
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array.Assign(SourceArray);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, size and array
  {
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, SourceArray);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, row major, interval and array
  {
    multidim_array_row<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_row Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, SourceArray);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, size and value
  {
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, 1);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  }

  // Multidimensional, column major, interval and value
  {
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, 1);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  }

  // Multidimensional, column major, size and initializer list
  {
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, {0,1,2,3,4,5});
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and initializer list
  {
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, size and iterator
  {
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, SourceArray.begin());
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and iterator
  {
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, SourceArray.begin());
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, array view, operator=
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array = array_view_col(SourceArray);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, array view, Assign
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array.Assign(array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, size and array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, const array view, operator=
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array = const_array_view_col(SourceArray);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, const array view, Assign
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array.Assign(const_array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, size and const array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, const_array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and const array view
  {
    multidim_array_col<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, const_array_view_col(SourceArray));
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, array, operator=
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array = SourceArray;
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, array, Assign
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array.Assign(SourceArray);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, size and array
  {
    multidim_array_col<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}}, SourceArray);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Multidimensional, column major, interval and array
  {
    multidim_array_col<int> SourceArray({{1,2,3}}, {0,1,2,3,4,5});
    array_col Array({{1,1,1}}, 0);
    Array.Assign({{1,2,3}, {2,4,6}}, SourceArray);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Non-copyable, rvalue array, operator=
  {
    std::array<noncopyable<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_noncopyable Array({1});
    Array = std::move(SourceArray);
    auto &View = helper_noncopyable::GetView(Array);
    auto &Values = helper_noncopyable::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_EQ(Values.Count(), 4);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
    Sum = 0;
    for (auto &Value : SourceArray) Sum += Value.Value();
    EXPECT_EQ(Sum, 0);
  }

  // Non-copyable, rvalue array, Assign
  {
    std::array<noncopyable<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_noncopyable Array({1});
    Array.Assign(std::move(SourceArray));
    auto &View = helper_noncopyable::GetView(Array);
    auto &Values = helper_noncopyable::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    EXPECT_EQ(Values.Count(), 4);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
    Sum = 0;
    for (auto &Value : SourceArray) Sum += Value.Value();
    EXPECT_EQ(Sum, 0);
  }

  // Non-copyable, interval and rvalue array
  {
    std::array<noncopyable<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_noncopyable Array({1});
    Array.Assign({-2, 2}, std::move(SourceArray));
    auto &View = helper_noncopyable::GetView(Array);
    auto &Values = helper_noncopyable::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    EXPECT_EQ(Values.Count(), 4);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
    Sum = 0;
    for (auto &Value : SourceArray) Sum += Value.Value();
    EXPECT_EQ(Sum, 0);
  }

  // Non-default-constructible, interval and value
  {
    array_nondefaultconstructible Array({1}, 0);
    Array.Assign({-2, 2}, 1);
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 4);
  }

  // Non-default-constructible, interval and iterator
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({1}, 0);
    Array.Assign({-2, 2}, SourceArray.begin());
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, array view, operator=
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({1}, 0);
    Array = array_view_nondefaultconstructible(SourceArray, {-2, 2});
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, array view, Assign
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({1}, 0);
    Array.Assign(array_view_nondefaultconstructible(SourceArray, {-2, 2}));
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, interval and array view
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({1}, 0);
    Array.Assign({-2, 2}, array_view_nondefaultconstructible(SourceArray));
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, array, operator=
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({1}, 0);
    Array = SourceArray;
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, array, Assign
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({1}, 0);
    Array.Assign(SourceArray);
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(4));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-default-constructible, interval and array
  {
    std::array<nondefaultconstructible<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    array_nondefaultconstructible Array({1}, 0);
    Array.Assign({-2, 2}, SourceArray);
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_EQ(View.Data(), Values.Data());
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-2));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2));
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

}

TEST_F(ArrayTests, Reserve) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;
  using helper = ovk::core::test_helper<array>;

  array Array({{1,2,3}, {2,4,6}}, 1);
  auto &View = helper::GetView(Array);
  auto &Values = helper::GetValues(Array);
  long long Capacity = Values.Capacity();
  Array.Reserve(Capacity+1);
  EXPECT_GE(Values.Capacity(), Capacity+1);
  EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  EXPECT_EQ(View.Data(), Values.Data());
  EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
  EXPECT_THAT(View.Extents().End(), ElementsAre(2,4,6));

}

TEST_F(ArrayTests, Resize) {

  if (TestComm().Rank() != 0) return;

  using array_1d = ovk::array<int>;
  using array_row = ovk::array<int,3>;
  using array_col = ovk::array<int,3,ovk::array_layout::COLUMN_MAJOR>;
  using array_noncopyable = ovk::array<noncopyable<int>>;
  using array_nondefaultconstructible = ovk::array<nondefaultconstructible<int>>;
  using helper_1d = ovk::core::test_helper<array_1d>;
  using helper_row = ovk::core::test_helper<array_row>;
  using helper_col = ovk::core::test_helper<array_col>;
  using helper_noncopyable = ovk::core::test_helper<array_noncopyable>;
  using helper_nondefaultconstructible = ovk::core::test_helper<array_nondefaultconstructible>;

  // One-dimensional, interval, change size
  {
    array_1d Array({4}, {0,1,2,3});
    Array.Resize({5});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(5));
    EXPECT_EQ(Values.Count(), 5);
    EXPECT_EQ(Values[0], 0);
    EXPECT_EQ(Values[1], 1);
    EXPECT_EQ(Values[2], 2);
    EXPECT_EQ(Values[3], 3);
  }

  // One-dimensional, interval, shift
  {
    array_1d Array({4}, {0,1,2,3});
    Array.Resize({-1,3});
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-1));
    EXPECT_THAT(View.Extents().End(), ElementsAre(3));
    EXPECT_EQ(Values.Count(), 4);
    EXPECT_EQ(Values[1], 0);
    EXPECT_EQ(Values[2], 1);
    EXPECT_EQ(Values[3], 2);
  }

  // One-dimensional, interval and value, change size
  {
    array_1d Array({4}, {0,1,2,3});
    Array.Resize({5}, -1);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(5));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,-1}));
  }

  // One-dimensional, interval and value, shift
  {
    array_1d Array({4}, {0,1,2,3});
    Array.Resize({-1,3}, -1);
    auto &View = helper_1d::GetView(Array);
    auto &Values = helper_1d::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-1));
    EXPECT_THAT(View.Extents().End(), ElementsAre(3));
    EXPECT_THAT(Values, ElementsAreArray({-1,0,1,2}));
  }

  // Multidimensional, row major, interval, change size in max-stride dimension
  {
    array_row Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{2,2,3}});
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,2,3));
    EXPECT_EQ(Values.Count(), 12);
    EXPECT_EQ(Values[0], 0);
    EXPECT_EQ(Values[1], 1);
    EXPECT_EQ(Values[2], 2);
    EXPECT_EQ(Values[3], 3);
    EXPECT_EQ(Values[4], 4);
    EXPECT_EQ(Values[5], 5);
  }

  // Multidimensional, row major, interval, change size in non-max-stride dimension
  {
    array_row Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{1,2,4}});
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,4));
    EXPECT_EQ(Values.Count(), 8);
    EXPECT_EQ(Values[0], 0);
    EXPECT_EQ(Values[1], 1);
    EXPECT_EQ(Values[2], 2);
    EXPECT_EQ(Values[4], 3);
    EXPECT_EQ(Values[5], 4);
    EXPECT_EQ(Values[6], 5);
  }

  // Multidimensional, row major, interval, shift
  {
    array_row Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{0,0,1}, {1,2,4}});
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,1));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,4));
    EXPECT_EQ(Values.Count(), 6);
    EXPECT_EQ(Values[0], 1);
    EXPECT_EQ(Values[1], 2);
    EXPECT_EQ(Values[3], 4);
    EXPECT_EQ(Values[4], 5);
  }

  // Multidimensional, row major, interval and value, change size in max-stride dimension
  {
    array_row Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{2,2,3}}, -1);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5,-1,-1,-1,-1,-1,-1}));
  }

  // Multidimensional, row major, interval and value, change size in non-max-stride dimension
  {
    array_row Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{1,2,4}}, -1);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,-1,3,4,5,-1}));
  }

  // Multidimensional, row major, interval and value, shift
  {
    array_row Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{0,0,1}, {1,2,4}}, -1);
    auto &View = helper_row::GetView(Array);
    auto &Values = helper_row::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,1));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,4));
    EXPECT_THAT(Values, ElementsAreArray({1,2,-1,4,5,-1}));
  }

  // Multidimensional, column major, interval, change size in max-stride dimension
  {
    array_col Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{1,2,4}});
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,4));
    EXPECT_EQ(Values.Count(), 8);
    EXPECT_EQ(Values[0], 0);
    EXPECT_EQ(Values[1], 1);
    EXPECT_EQ(Values[2], 2);
    EXPECT_EQ(Values[3], 3);
    EXPECT_EQ(Values[4], 4);
    EXPECT_EQ(Values[5], 5);
  }

  // Multidimensional, column major, interval, change size in non-max-stride dimension
  {
    array_col Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{2,2,3}});
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,2,3));
    EXPECT_EQ(Values.Count(), 12);
    EXPECT_EQ(Values[0], 0);
    EXPECT_EQ(Values[2], 1);
    EXPECT_EQ(Values[4], 2);
    EXPECT_EQ(Values[6], 3);
    EXPECT_EQ(Values[8], 4);
    EXPECT_EQ(Values[10], 5);
  }

  // Multidimensional, column major, interval, shift
  {
    array_col Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{0,0,1}, {1,2,4}});
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,1));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,4));
    EXPECT_EQ(Values.Count(), 6);
    EXPECT_EQ(Values[0], 2);
    EXPECT_EQ(Values[1], 3);
    EXPECT_EQ(Values[2], 4);
    EXPECT_EQ(Values[3], 5);
  }

  // Multidimensional, column major, interval and value, change size in max-stride dimension
  {
    array_col Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{1,2,4}}, -1);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,4));
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5,-1,-1}));
  }

  // Multidimensional, column major, interval and value, change size in non-max-stride dimension
  {
    array_col Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{2,2,3}}, -1);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(2,2,3));
    EXPECT_THAT(Values, ElementsAreArray({0,-1,1,-1,2,-1,3,-1,4,-1,5,-1}));
  }

  // Multidimensional, column major, interval and value, shift
  {
    array_col Array({{1,2,3}}, {0,1,2,3,4,5});
    Array.Resize({{0,0,1}, {1,2,4}}, -1);
    auto &View = helper_col::GetView(Array);
    auto &Values = helper_col::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,1));
    EXPECT_THAT(View.Extents().End(), ElementsAre(1,2,4));
    EXPECT_THAT(Values, ElementsAreArray({2,3,4,5,-1,-1}));
  }

  // Non-copyable, interval, change size
  {
    array_noncopyable Array({4});
    auto &View = helper_noncopyable::GetView(Array);
    auto &Values = helper_noncopyable::GetValues(Array);
    for (int i = 0; i < 4; ++i) Values[i] = {i};
    Array.Resize({5});
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(5));
    EXPECT_EQ(Values.Count(), 5);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
  }

  // Non-copyable, interval, shift
  {
    array_noncopyable Array({4});
    auto &View = helper_noncopyable::GetView(Array);
    auto &Values = helper_noncopyable::GetValues(Array);
    for (int i = 0; i < 4; ++i) Values[i] = {i};
    Array.Resize({-1,3});
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-1));
    EXPECT_THAT(View.Extents().End(), ElementsAre(3));
    EXPECT_EQ(Values.Count(), 4);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 3);
  }

  // Non-defaultconstructible, interval, change size
  {
    array_nondefaultconstructible Array({4}, {{0},{1},{2},{3}});
    Array.Resize({5}, {-1});
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(5));
    EXPECT_EQ(Values.Count(), 5);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 5);
  }

  // Non-defaultconstructible, interval, shift
  {
    array_nondefaultconstructible Array({4}, {{0},{1},{2},{3}});
    Array.Resize({-1,3}, -1);
    auto &View = helper_nondefaultconstructible::GetView(Array);
    auto &Values = helper_nondefaultconstructible::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(-1));
    EXPECT_THAT(View.Extents().End(), ElementsAre(3));
    EXPECT_EQ(Values.Count(), 4);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 2);
  }

}

TEST_F(ArrayTests, Clear) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;
  using helper = ovk::core::test_helper<array>;

  // Empty
  {
    array Array;
    Array.Clear();
    auto &View = helper::GetView(Array);
    auto &Values = helper::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(0,0,0));
    EXPECT_EQ(Values.Count(), 0);
  }

  // Non-empty
  {
    array Array({{1,2,3}, {2,4,6}}, 1);
    Array.Clear();
    auto &View = helper::GetView(Array);
    auto &Values = helper::GetValues(Array);
    EXPECT_THAT(View.Extents().Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(View.Extents().End(), ElementsAre(0,0,0));
    EXPECT_EQ(Values.Count(), 0);
  }

}

TEST_F(ArrayTests, Fill) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;
  using array_view = ovk::array_view<int,3>;
  using array_noncopyable = ovk::array<noncopyable<int>>;
  using helper = ovk::core::test_helper<array>;
  using helper_noncopyable = ovk::core::test_helper<array_noncopyable>;

  // Constant value
  {
    array Array({{1,2,3}, {2,4,6}}, 0);
    Array.Fill(1);
    auto &Values = helper::GetValues(Array);
    EXPECT_THAT(Values, ElementsAreArray({1,1,1,1,1,1}));
  }

  // Initializer list
  {
    array Array({{1,2,3}}, 0);
    Array.Fill({0,1,2,3,4,5});
    auto &Values = helper::GetValues(Array);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Iterator
  {
    array Array({{1,2,3}, {2,4,6}}, 0);
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    Array.Fill(SourceArray.begin());
    auto &Values = helper::GetValues(Array);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // View
  {
    array Array({{1,2,3}, {2,4,6}}, 0);
    std::array<int,6> SourceArray = {{0,1,2,3,4,5}};
    Array.Fill(array_view(SourceArray.data(), {{1,2,3}, {2,4,6}}));
    auto &Values = helper::GetValues(Array);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Array, const ref
  {
    array Array({{1,2,3}, {2,4,6}}, 0);
    multidim_array_row<int> SourceArray({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    Array.Fill(static_cast<const multidim_array_row<int> &>(SourceArray));
    auto &Values = helper::GetValues(Array);
    EXPECT_THAT(Values, ElementsAreArray({0,1,2,3,4,5}));
  }

  // Array, rvalue ref
  {
    array_noncopyable Array({4});
    std::array<noncopyable<int>,4> SourceArray = {{{0},{1},{2},{3}}};
    Array.Fill(std::move(SourceArray));
    auto &Values = helper_noncopyable::GetValues(Array);
    int Sum = 0;
    for (auto &Value : Values) Sum += Value.Value();
    EXPECT_EQ(Sum, 6);
    Sum = 0;
    for (auto &Value : SourceArray) Sum += Value.Value();
    EXPECT_EQ(Sum, 0);
  }

}

TEST_F(ArrayTests, Extents) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;

  array Array({{1,2,3}, {2,4,6}});
  EXPECT_THAT(Array.Extents().Begin(), ElementsAre(1,2,3));
  EXPECT_THAT(Array.Extents().End(), ElementsAre(2,4,6));

}

TEST_F(ArrayTests, Begin) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;

  array Array({{1,2,3}, {2,4,6}});
  EXPECT_THAT(Array.Begin(), ElementsAre(1,2,3));
  EXPECT_EQ(Array.Begin(0), 1);
  EXPECT_EQ(Array.Begin(1), 2);
  EXPECT_EQ(Array.Begin(2), 3);

}

TEST_F(ArrayTests, End) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;

  array Array({{1,2,3}, {2,4,6}});
  EXPECT_THAT(Array.End(), ElementsAre(2,4,6));
  EXPECT_EQ(Array.End(0), 2);
  EXPECT_EQ(Array.End(1), 4);
  EXPECT_EQ(Array.End(2), 6);

}

TEST_F(ArrayTests, Size) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;

  array Array({{1,2,3}, {2,4,6}});
  EXPECT_THAT(Array.Size(), ElementsAre(1,2,3));
  EXPECT_EQ(Array.Size(0), 1);
  EXPECT_EQ(Array.Size(1), 2);
  EXPECT_EQ(Array.Size(2), 3);

}

TEST_F(ArrayTests, Count) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;

  array Array({{1,2,3}, {2,4,6}});
  EXPECT_EQ(Array.Count(), 6);

}

TEST_F(ArrayTests, Empty) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;

  // Empty
  {
    array Array({{1,2,3}, {2,2,6}});
    EXPECT_TRUE(Array.Empty());
  }

  // Non-empty
  {
    array Array({{1,2,3}, {2,4,6}});
    EXPECT_FALSE(Array.Empty());
  }

}

TEST_F(ArrayTests, Capacity) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;
  using helper = ovk::core::test_helper<array>;

  array Array({{1,2,3}, {2,4,6}});
  auto &Values = helper::GetValues(Array);
  EXPECT_EQ(Array.Capacity(), Values.Capacity());

}

TEST_F(ArrayTests, Indexer) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;

  array Array({{1,2,3}, {4,5,6}});
  EXPECT_THAT(Array.Indexer().Begin(), ElementsAre(1,2,3));
  EXPECT_THAT(Array.Indexer().Stride(), ElementsAre(9,3,1));

}

TEST_F(ArrayTests, ParenthesisOperator) {

  if (TestComm().Rank() != 0) return;

  using array_row = ovk::array<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_col = ovk::array<int,3,ovk::array_layout::COLUMN_MAJOR>;
  using helper_row = ovk::core::test_helper<array_row>;
  using helper_col = ovk::core::test_helper<array_col>;

  // Row major, iterator
  {
    array_row Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_row::GetValues(Array);
    ovk::elem<int,3> Tuple = {2,3,4};
    EXPECT_EQ(&Array(Tuple.Data()), &Values[33]);
  }

  // Row major, array
  {
    array_row Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(&Array(ovk::elem<int,3>(2,3,4)), &Values[33]);
  }

  // Row major, separate
  {
    array_row Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(&Array(2,3,4), &Values[33]);
  }

  // Column major, iterator
  {
    array_col Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_col::GetValues(Array);
    ovk::elem<int,3> Tuple = {2,3,4};
    EXPECT_EQ(&Array(Tuple.Data()), &Values[43]);
  }

  // Column major, array
  {
    array_col Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(&Array(ovk::elem<int,3>(2,3,4)), &Values[43]);
  }

  // Column major, separate
  {
    array_col Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(&Array(2,3,4), &Values[43]);
  }

}

TEST_F(ArrayTests, BracketOperator) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;
  using helper = ovk::core::test_helper<array>;

  array Array({{1,1,1}, {4,5,6}});
  auto &Values = helper::GetValues(Array);
  EXPECT_EQ(&Array[33], &Values[33]);

}

TEST_F(ArrayTests, Data) {

  if (TestComm().Rank() != 0) return;

  using array_row = ovk::array<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_col = ovk::array<int,3,ovk::array_layout::COLUMN_MAJOR>;
  using helper_row = ovk::core::test_helper<array_row>;
  using helper_col = ovk::core::test_helper<array_col>;

  // Row major, no argument
  {
    array_row Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(Array.Data(), &Values[0]);
  }

  // Row major, iterator
  {
    array_row Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_row::GetValues(Array);
    ovk::elem<int,3> Tuple = {2,3,4};
    EXPECT_EQ(Array.Data(Tuple.Data()), &Values[33]);
  }

  // Row major, array
  {
    array_row Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(Array.Data(ovk::elem<int,3>(2,3,4)), &Values[33]);
  }

  // Row major, separate
  {
    array_row Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_row::GetValues(Array);
    EXPECT_EQ(Array.Data(2,3,4), &Values[33]);
  }

  // Column major, no argument
  {
    array_col Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(Array.Data(), &Values[0]);
  }

  // Column major, iterator
  {
    array_col Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_col::GetValues(Array);
    ovk::elem<int,3> Tuple = {2,3,4};
    EXPECT_EQ(Array.Data(Tuple.Data()), &Values[43]);
  }

  // Column major, array
  {
    array_col Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(Array.Data(ovk::elem<int,3>(2,3,4)), &Values[43]);
  }

  // Column major, separate
  {
    array_col Array({{1,1,1}, {4,5,6}});
    auto &Values = helper_col::GetValues(Array);
    EXPECT_EQ(Array.Data(2,3,4), &Values[43]);
  }

}

TEST_F(ArrayTests, LinearBeginEnd) {

  if (TestComm().Rank() != 0) return;

  using array = ovk::array<int,3>;
  using helper = ovk::core::test_helper<array>;

  // Const
  {
    const array Array({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    auto &Values = helper::GetValues(Array);
    EXPECT_EQ(Array.LinearBegin(), &Values[0]);
    EXPECT_EQ(Array.LinearEnd(), &Values[0] + 6);
    int Sum = 0;
    for (auto &Value : Array) Sum += Value;
    EXPECT_EQ(Sum, 15);
  }

  // Non-const
  {
    array Array({{1,2,3}, {2,4,6}}, {0,1,2,3,4,5});
    auto &Values = helper::GetValues(Array);
    EXPECT_EQ(Array.LinearBegin(), &Values[0]);
    EXPECT_EQ(Array.LinearEnd(), &Values[0] + 6);
    for (auto &Value : Array) Value = 1;
    int Sum = 0;
    for (int i = 0; i < 6; ++i) Sum += Values[i];
    EXPECT_EQ(Sum, 6);
  }

}

TEST_F(ArrayTests, Traits) {

  if (TestComm().Rank() != 0) return;

  using array_row = ovk::array<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_col = ovk::array<int,3,ovk::array_layout::COLUMN_MAJOR>;

  // Row major
  EXPECT_TRUE(ovk::core::IsArray<array_row>());
  EXPECT_TRUE((std::is_same<ovk::core::array_value_type<array_row>, int>::value));
  EXPECT_EQ(ovk::core::ArrayRank<array_row>(), 3);
  EXPECT_EQ(ovk::core::ArrayLayout<array_row>(), ovk::array_layout::ROW_MAJOR);
  EXPECT_TRUE(ovk::core::ArrayHasRuntimeExtents<array_row>());
  {
    array_row Array({{1,2,3}, {2,4,6}});
    EXPECT_THAT(ovk::core::ArrayBegin(Array), ElementsAre(1,2,3));
    EXPECT_THAT(ovk::core::ArrayEnd(Array), ElementsAre(2,4,6));
    EXPECT_EQ(ovk::core::ArrayData(Array), Array.Data());
  }

  // Column major
  EXPECT_TRUE(ovk::core::IsArray<array_col>());
  EXPECT_TRUE((std::is_same<ovk::core::array_value_type<array_col>, int>::value));
  EXPECT_EQ(ovk::core::ArrayRank<array_col>(), 3);
  EXPECT_EQ(ovk::core::ArrayLayout<array_col>(), ovk::array_layout::COLUMN_MAJOR);
  EXPECT_TRUE(ovk::core::ArrayHasRuntimeExtents<array_col>());
  {
    array_col Array({{1,2,3}, {2,4,6}});
    EXPECT_THAT(ovk::core::ArrayBegin(Array), ElementsAre(1,2,3));
    EXPECT_THAT(ovk::core::ArrayEnd(Array), ElementsAre(2,4,6));
    EXPECT_EQ(ovk::core::ArrayData(Array), Array.Data());
  }

}
