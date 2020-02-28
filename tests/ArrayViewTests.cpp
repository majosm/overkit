// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/ArrayView.hpp>

#include "tests/MPITest.hpp"
#include "tests/mocks/MultidimArray.hpp"
#include "tests/mocks/Noncopyable.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Interval.hpp>

#include <mpi.h>

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

using testing::ElementsAre;

class ArrayViewTests : public tests::mpi_test {};

using tests::multidim_array_row;
using tests::multidim_array_col;
using tests::noncopyable;

namespace ovk {
namespace core {
template <typename T, int Rank, array_layout Layout> class test_helper<array_view<T,Rank,Layout>> {
public:
  using array_view_type = array_view<T,Rank,Layout>;
  using interval_type = interval<long long,Rank>;
  using indexer_type = indexer<long long, long long, Rank, Layout>;
  static const T *GetPtr(const array_view_type &View) { return View.Ptr_; }
  static const interval_type &GetExtents(const array_view_type &View) { return View.Extents_; }
  static const indexer_type &GetIndexer(const array_view_type &View) { return View.Indexer_; }
};
}}

TEST_F(ArrayViewTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;

  EXPECT_TRUE((std::is_same<typename array_view::value_type, int>::value));
  EXPECT_EQ(int(array_view::Rank), 3);
  EXPECT_EQ(ovk::array_layout(array_view::Layout), ovk::array_layout::ROW_MAJOR);
  EXPECT_TRUE((std::is_same<typename array_view::index_type, long long>::value));
  EXPECT_TRUE((std::is_same<typename array_view::tuple_element_type, long long>::value));
  EXPECT_TRUE((std::is_same<typename array_view::tuple_type, ovk::elem<long long,3>>::value));
  EXPECT_TRUE((std::is_same<typename array_view::interval_type, ovk::interval<long long,3>>::value));
  EXPECT_TRUE((std::is_same<typename array_view::indexer_type, ovk::indexer<long long, long long,
    3>>::value));
  EXPECT_TRUE((std::is_same<typename array_view::iterator::pointer, int *>::value));

  using array_view_col = ovk::array_view<int,3,ovk::array_layout::COLUMN_MAJOR>;

  EXPECT_TRUE((std::is_same<typename array_view_col::indexer_type, ovk::indexer<long long,
    long long, 3, ovk::array_layout::COLUMN_MAJOR>>::value));
  EXPECT_EQ(ovk::array_layout(array_view_col::Layout), ovk::array_layout::COLUMN_MAJOR);

  using array_view_const = ovk::array_view<const int,3>;

  EXPECT_TRUE((std::is_same<typename array_view_const::value_type, const int>::value));
  EXPECT_TRUE((std::is_same<typename array_view_const::iterator::pointer, const int *>::value));

}

TEST_F(ArrayViewTests, Create) {

  if (TestComm().Rank() != 0) return;

  using array_view_1d = ovk::array_view<int>;
  using array_view_1d_const = ovk::array_view<const int>;
  using array_view_row = ovk::array_view<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_view_row_const = ovk::array_view<const int,3,ovk::array_layout::ROW_MAJOR>;
  using array_view_col = ovk::array_view<int,3,ovk::array_layout::COLUMN_MAJOR>;
  using array_view_col_const = ovk::array_view<const int,3,ovk::array_layout::COLUMN_MAJOR>;
  using helper_1d = ovk::core::test_helper<array_view_1d>;
  using helper_1d_const = ovk::core::test_helper<array_view_1d_const>;
  using helper_row = ovk::core::test_helper<array_view_row>;
  using helper_row_const = ovk::core::test_helper<array_view_row_const>;
  using helper_col = ovk::core::test_helper<array_view_col>;
  using helper_col_const = ovk::core::test_helper<array_view_col_const>;

  // One-dimensional, non-const, default
  {
    array_view_1d View;
    EXPECT_EQ(helper_1d::GetPtr(View), nullptr);
    EXPECT_THAT(helper_1d::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetExtents(View).End(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, non-const, pointer and size
  {
    int Array[4] = {};
    array_view_1d View(&Array[0], {4});
    EXPECT_EQ(helper_1d::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, non-const, pointer and interval
  {
    int Array[4] = {};
    array_view_1d View(&Array[0], {-2, 2});
    EXPECT_EQ(helper_1d::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d::GetExtents(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d::GetExtents(View).End(), ElementsAre(2));
    EXPECT_THAT(helper_1d::GetIndexer(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, non-const, array
  {
    std::array<int,4> Array = {{}};
    array_view_1d View(Array);
    EXPECT_EQ(helper_1d::GetPtr(View), Array.data());
    EXPECT_THAT(helper_1d::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, non-const, array and size
  {
    std::array<int,4> Array = {{}};
    array_view_1d View(Array, {4});
    EXPECT_EQ(helper_1d::GetPtr(View), Array.data());
    EXPECT_THAT(helper_1d::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, non-const, array and interval
  {
    std::array<int,4> Array = {{}};
    array_view_1d View(Array, {-2, 2});
    EXPECT_EQ(helper_1d::GetPtr(View), Array.data());
    EXPECT_THAT(helper_1d::GetExtents(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d::GetExtents(View).End(), ElementsAre(2));
    EXPECT_THAT(helper_1d::GetIndexer(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, const, default
  {
    array_view_1d_const View;
    EXPECT_EQ(helper_1d_const::GetPtr(View), nullptr);
    EXPECT_THAT(helper_1d_const::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetExtents(View).End(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, const, pointer and size
  {
    int Array[4] = {};
    array_view_1d_const View(&Array[0], {4});
    EXPECT_EQ(helper_1d_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d_const::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, const, pointer and interval
  {
    int Array[4] = {};
    array_view_1d_const View(&Array[0], {-2, 2});
    EXPECT_EQ(helper_1d_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d_const::GetExtents(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d_const::GetExtents(View).End(), ElementsAre(2));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, const, array
  {
    std::array<int,4> Array = {{}};
    array_view_1d_const View(Array);
    EXPECT_EQ(helper_1d_const::GetPtr(View), Array.data());
    EXPECT_THAT(helper_1d_const::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, const, array and size
  {
    std::array<int,4> Array = {{}};
    array_view_1d_const View(Array, {4});
    EXPECT_EQ(helper_1d_const::GetPtr(View), Array.data());
    EXPECT_THAT(helper_1d_const::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, const, array and interval
  {
    std::array<int,4> Array = {{}};
    array_view_1d_const View(Array, {-2, 2});
    EXPECT_EQ(helper_1d_const::GetPtr(View), Array.data());
    EXPECT_THAT(helper_1d_const::GetExtents(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d_const::GetExtents(View).End(), ElementsAre(2));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // Multidimensional, row major, non-const, default
  {
    array_view_row View;
    EXPECT_EQ(helper_row::GetPtr(View), nullptr);
    EXPECT_THAT(helper_row::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetExtents(View).End(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetIndexer(View).Stride(), ElementsAre(0,0,1));
  }

  // Multidimensional, row major, non-const, pointer and size
  {
    int Array[6] = {};
    array_view_row View(&Array[0], {{1,2,3}});
    EXPECT_EQ(helper_row::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetExtents(View).End(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetIndexer(View).Stride(), ElementsAre(6,3,1));
  }

  // Multidimensional, row major, non-const, pointer and interval
  {
    int Array[27] = {};
    array_view_row View(&Array[0], {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper_row::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, non-const, array
  {
    multidim_array_row<int> Array({{1,2,3}, {4,5,6}});
    array_view_row View(Array);
    EXPECT_EQ(helper_row::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, non-const, array and size
  {
    multidim_array_row<int> Array({{1,2,3}, {4,5,6}});
    array_view_row View(Array, {{3,3,3}});
    EXPECT_EQ(helper_row::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetExtents(View).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper_row::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, non-const, array and interval
  {
    multidim_array_row<int> Array({{0,0,0}, {3,3,3}});
    array_view_row View(Array, {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper_row::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, const, default
  {
    array_view_row_const View;
    EXPECT_EQ(helper_row_const::GetPtr(View), nullptr);
    EXPECT_THAT(helper_row_const::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row_const::GetExtents(View).End(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Stride(), ElementsAre(0,0,1));
  }

  // Multidimensional, row major, const, pointer and size
  {
    int Array[6] = {};
    array_view_row_const View(&Array[0], {{1,2,3}});
    EXPECT_EQ(helper_row_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row_const::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row_const::GetExtents(View).End(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Stride(), ElementsAre(6,3,1));
  }

  // Multidimensional, row major, const, pointer and interval
  {
    int Array[27] = {};
    array_view_row_const View(&Array[0], {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper_row_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, const, array
  {
    multidim_array_row<int> Array({{1,2,3}, {4,5,6}});
    array_view_row_const View(Array);
    EXPECT_EQ(helper_row_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, const, array and size
  {
    multidim_array_row<int> Array({{1,2,3}, {4,5,6}});
    array_view_row_const View(Array, {{3,3,3}});
    EXPECT_EQ(helper_row_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row_const::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row_const::GetExtents(View).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, const, array and interval
  {
    multidim_array_row<int> Array({{0,0,0}, {3,3,3}});
    array_view_row_const View(Array, {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper_row_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, column major, non-const, default
  {
    array_view_col View;
    EXPECT_EQ(helper_col::GetPtr(View), nullptr);
    EXPECT_THAT(helper_col::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetExtents(View).End(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetIndexer(View).Stride(), ElementsAre(1,0,0));
  }

  // Multidimensional, column major, non-const, pointer and size
  {
    int Array[6] = {};
    array_view_col View(&Array[0], {{1,2,3}});
    EXPECT_EQ(helper_col::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetExtents(View).End(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetIndexer(View).Stride(), ElementsAre(1,1,2));
  }

  // Multidimensional, column major, non-const, pointer and interval
  {
    int Array[27] = {};
    array_view_col View(&Array[0], {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper_col::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, non-const, array
  {
    multidim_array_col<int> Array({{1,2,3}, {4,5,6}});
    array_view_col View(Array);
    EXPECT_EQ(helper_col::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, non-const, array and size
  {
    multidim_array_col<int> Array({{1,2,3}, {4,5,6}});
    array_view_col View(Array, {{3,3,3}});
    EXPECT_EQ(helper_col::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetExtents(View).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper_col::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, non-const, array and interval
  {
    multidim_array_col<int> Array({{0,0,0}, {3,3,3}});
    array_view_col View(Array, {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper_col::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, const, default
  {
    array_view_col_const View;
    EXPECT_EQ(helper_col_const::GetPtr(View), nullptr);
    EXPECT_THAT(helper_col_const::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col_const::GetExtents(View).End(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Stride(), ElementsAre(1,0,0));
  }

  // Multidimensional, column major, const, pointer and size
  {
    int Array[6] = {};
    array_view_col_const View(&Array[0], {{1,2,3}});
    EXPECT_EQ(helper_col_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col_const::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col_const::GetExtents(View).End(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Stride(), ElementsAre(1,1,2));
  }

  // Multidimensional, column major, const, pointer and interval
  {
    int Array[27] = {};
    array_view_col_const View(&Array[0], {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper_col_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, const, array
  {
    multidim_array_col<int> Array({{1,2,3}, {4,5,6}});
    array_view_col_const View(Array);
    EXPECT_EQ(helper_col_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, const, array and size
  {
    multidim_array_col<int> Array({{1,2,3}, {4,5,6}});
    array_view_col_const View(Array, {{3,3,3}});
    EXPECT_EQ(helper_col_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col_const::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col_const::GetExtents(View).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, const, array and interval
  {
    multidim_array_col<int> Array({{0,0,0}, {3,3,3}});
    array_view_col_const View(Array, {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper_col_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

}

TEST_F(ArrayViewTests, MakeArrayView) {

  if (TestComm().Rank() != 0) return;

  using array_view_1d = ovk::array_view<int>;
  using array_view_1d_const = ovk::array_view<const int>;
  using array_view_row = ovk::array_view<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_view_row_const = ovk::array_view<const int,3,ovk::array_layout::ROW_MAJOR>;
  using array_view_col = ovk::array_view<int,3,ovk::array_layout::COLUMN_MAJOR>;
  using array_view_col_const = ovk::array_view<const int,3,ovk::array_layout::COLUMN_MAJOR>;
  using helper_1d = ovk::core::test_helper<array_view_1d>;
  using helper_1d_const = ovk::core::test_helper<array_view_1d_const>;
  using helper_row = ovk::core::test_helper<array_view_row>;
  using helper_row_const = ovk::core::test_helper<array_view_row_const>;
  using helper_col = ovk::core::test_helper<array_view_col>;
  using helper_col_const = ovk::core::test_helper<array_view_col_const>;

  // One-dimensional, non-const, array
  {
    std::array<int,4> Array = {{}};
    auto View = ovk::MakeArrayView(Array);
    EXPECT_TRUE((std::is_same<decltype(View), array_view_1d>::value));
    EXPECT_EQ(helper_1d::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, non-const, array and size
  {
    std::array<int,4> Array = {{}};
    auto View = ovk::MakeArrayView(Array, {4});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_1d>::value));
    EXPECT_EQ(helper_1d::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, non-const, array and interval
  {
    std::array<int,4> Array = {{}};
    auto View = ovk::MakeArrayView(Array, {-2, 2});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_1d>::value));
    EXPECT_EQ(helper_1d::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d::GetExtents(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d::GetExtents(View).End(), ElementsAre(2));
    EXPECT_THAT(helper_1d::GetIndexer(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, const, array
  {
    const std::array<int,4> Array = {{}};
    auto View = ovk::MakeArrayView(Array);
    EXPECT_TRUE((std::is_same<decltype(View), array_view_1d_const>::value));
    EXPECT_EQ(helper_1d_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d_const::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, const, array and size
  {
    const std::array<int,4> Array = {{}};
    auto View = ovk::MakeArrayView(Array, {4});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_1d_const>::value));
    EXPECT_EQ(helper_1d_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d_const::GetExtents(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetExtents(View).End(), ElementsAre(4));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Begin(), ElementsAre(0));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // One-dimensional, const, array and interval
  {
    const std::array<int,4> Array = {{}};
    auto View = ovk::MakeArrayView(Array, {-2, 2});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_1d_const>::value));
    EXPECT_EQ(helper_1d_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_1d_const::GetExtents(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d_const::GetExtents(View).End(), ElementsAre(2));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Begin(), ElementsAre(-2));
    EXPECT_THAT(helper_1d_const::GetIndexer(View).Stride(), ElementsAre(1));
  }

  // Multidimensional, row major, non-const, array
  {
    multidim_array_row<int> Array({{1,2,3}, {4,5,6}});
    auto View = ovk::MakeArrayView(Array);
    EXPECT_TRUE((std::is_same<decltype(View), array_view_row>::value));
    EXPECT_EQ(helper_row::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, non-const, array and size
  {
    multidim_array_row<int> Array({{1,2,3}, {4,5,6}});
    auto View = ovk::MakeArrayView(Array, {{3,3,3}});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_row>::value));
    EXPECT_EQ(helper_row::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetExtents(View).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper_row::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, non-const, array and interval
  {
    multidim_array_row<int> Array({{0,0,0}, {3,3,3}});
    auto View = ovk::MakeArrayView(Array, {{1,2,3}, {4,5,6}});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_row>::value));
    EXPECT_EQ(helper_row::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, const, array
  {
    const multidim_array_row<int> Array({{1,2,3}, {4,5,6}});
    auto View = ovk::MakeArrayView(Array);
    EXPECT_TRUE((std::is_same<decltype(View), array_view_row_const>::value));
    EXPECT_EQ(helper_row_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, const, array and size
  {
    const multidim_array_row<int> Array({{1,2,3}, {4,5,6}});
    auto View = ovk::MakeArrayView(Array, {{3,3,3}});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_row_const>::value));
    EXPECT_EQ(helper_row_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row_const::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row_const::GetExtents(View).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, row major, const, array and interval
  {
    const multidim_array_row<int> Array({{0,0,0}, {3,3,3}});
    auto View = ovk::MakeArrayView(Array, {{1,2,3}, {4,5,6}});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_row_const>::value));
    EXPECT_EQ(helper_row_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_row_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row_const::GetIndexer(View).Stride(), ElementsAre(9,3,1));
  }

  // Multidimensional, column major, non-const, array
  {
    multidim_array_col<int> Array({{1,2,3}, {4,5,6}});
    auto View = ovk::MakeArrayView(Array);
    EXPECT_TRUE((std::is_same<decltype(View), array_view_col>::value));
    EXPECT_EQ(helper_col::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, non-const, array and size
  {
    multidim_array_col<int> Array({{1,2,3}, {4,5,6}});
    auto View = ovk::MakeArrayView(Array, {{3,3,3}});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_col>::value));
    EXPECT_EQ(helper_col::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetExtents(View).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper_col::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, non-const, array and interval
  {
    multidim_array_col<int> Array({{0,0,0}, {3,3,3}});
    auto View = ovk::MakeArrayView(Array, {{1,2,3}, {4,5,6}});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_col>::value));
    EXPECT_EQ(helper_col::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, const, array
  {
    const multidim_array_col<int> Array({{1,2,3}, {4,5,6}});
    auto View = ovk::MakeArrayView(Array);
    EXPECT_TRUE((std::is_same<decltype(View), array_view_col_const>::value));
    EXPECT_EQ(helper_col_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, const, array and size
  {
    const multidim_array_col<int> Array({{1,2,3}, {4,5,6}});
    auto View = ovk::MakeArrayView(Array, {{3,3,3}});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_col_const>::value));
    EXPECT_EQ(helper_col_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col_const::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col_const::GetExtents(View).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

  // Multidimensional, column major, const, array and interval
  {
    const multidim_array_col<int> Array({{0,0,0}, {3,3,3}});
    auto View = ovk::MakeArrayView(Array, {{1,2,3}, {4,5,6}});
    EXPECT_TRUE((std::is_same<decltype(View), array_view_col_const>::value));
    EXPECT_EQ(helper_col_const::GetPtr(View), &Array[0]);
    EXPECT_THAT(helper_col_const::GetExtents(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetExtents(View).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col_const::GetIndexer(View).Stride(), ElementsAre(1,3,9));
  }

}

TEST_F(ArrayViewTests, Copy) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using helper = ovk::core::test_helper<array_view>;
  using multidim_array = multidim_array_row<int>;

  // Copy construct
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View1(Array);
    array_view View2 = View1;
    EXPECT_EQ(helper::GetPtr(View2), &Array[0]);
    EXPECT_THAT(helper::GetExtents(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetExtents(View2).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper::GetIndexer(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetIndexer(View2).Stride(), ElementsAre(9,3,1));
  }

  // Copy construct with size
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View1(Array);
    array_view View2(View1, {{3,3,3}});
    EXPECT_EQ(helper::GetPtr(View2), &Array[0]);
    EXPECT_THAT(helper::GetExtents(View2).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper::GetExtents(View2).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper::GetIndexer(View2).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper::GetIndexer(View2).Stride(), ElementsAre(9,3,1));
  }

  // Copy construct with interval
  {
    multidim_array Array({{0,0,0}, {3,3,3}});
    array_view View1(Array);
    array_view View2(View1, {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper::GetPtr(View2), &Array[0]);
    EXPECT_THAT(helper::GetExtents(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetExtents(View2).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper::GetIndexer(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetIndexer(View2).Stride(), ElementsAre(9,3,1));
  }

  // Copy assign
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View1(Array);
    array_view View2;
    View2 = View1;
    EXPECT_EQ(helper::GetPtr(View2), &Array[0]);
    EXPECT_THAT(helper::GetExtents(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetExtents(View2).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper::GetIndexer(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper::GetIndexer(View2).Stride(), ElementsAre(9,3,1));
  }

}

// Can't put these in anonymous namespace because clang optimizes them out and emits warnings
namespace array_view_tests_internal {

std::true_type ConvertsToNonConstView(ovk::array_view<int,3> View) { return {}; }
std::false_type ConvertsToNonConstView(...) { return {}; }

}

TEST_F(ArrayViewTests, Convert) {

  using namespace array_view_tests_internal;

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using array_view_const = ovk::array_view<const int,3>;
  using helper_const = ovk::core::test_helper<array_view_const>;
  using multidim_array = multidim_array_row<int>;

  // Convert to const construct
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View1(Array);
    array_view_const View2 = View1;
    EXPECT_EQ(helper_const::GetPtr(View2), &Array[0]);
    EXPECT_THAT(helper_const::GetExtents(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_const::GetExtents(View2).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_const::GetIndexer(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_const::GetIndexer(View2).Stride(), ElementsAre(9,3,1));
  }

  // Convert to const construct with size
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View1(Array);
    array_view_const View2(View1, {{3,3,3}});
    EXPECT_EQ(helper_const::GetPtr(View2), &Array[0]);
    EXPECT_THAT(helper_const::GetExtents(View2).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_const::GetExtents(View2).End(), ElementsAre(3,3,3));
    EXPECT_THAT(helper_const::GetIndexer(View2).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper_const::GetIndexer(View2).Stride(), ElementsAre(9,3,1));
  }

  // Convert to const construct with interval
  {
    multidim_array Array({{0,0,0}, {3,3,3}});
    array_view View1(Array);
    array_view_const View2(View1, {{1,2,3}, {4,5,6}});
    EXPECT_EQ(helper_const::GetPtr(View2), &Array[0]);
    EXPECT_THAT(helper_const::GetExtents(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_const::GetExtents(View2).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_const::GetIndexer(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_const::GetIndexer(View2).Stride(), ElementsAre(9,3,1));
  }

  // Convert to const assign
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View1(Array);
    array_view_const View2;
    View2 = View1;
    EXPECT_EQ(helper_const::GetPtr(View2), &Array[0]);
    EXPECT_THAT(helper_const::GetExtents(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_const::GetExtents(View2).End(), ElementsAre(4,5,6));
    EXPECT_THAT(helper_const::GetIndexer(View2).Begin(), ElementsAre(1,2,3));
    EXPECT_THAT(helper_const::GetIndexer(View2).Stride(), ElementsAre(9,3,1));
  }

  // Make sure non-const view can't be created from const array
  {
    const multidim_array Array({{1,2,3}, {4,5,6}});
    EXPECT_FALSE(decltype(ConvertsToNonConstView(Array))::value);
    // Sanity check / suppress unused function warnings
    EXPECT_TRUE(decltype(ConvertsToNonConstView(multidim_array({{1,2,3}, {4,5,6}})))::value);
  }

}

TEST_F(ArrayViewTests, Reset) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using helper = ovk::core::test_helper<array_view>;
  using multidim_array = multidim_array_row<int>;

  // Initially null
  {
    array_view View;
    View.Reset();
    EXPECT_EQ(helper::GetPtr(View), nullptr);
    EXPECT_THAT(helper::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper::GetExtents(View).End(), ElementsAre(0,0,0));
    EXPECT_THAT(helper::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper::GetIndexer(View).Stride(), ElementsAre(0,0,1));
  }

  // Initially not null
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View(Array);
    View.Reset();
    EXPECT_EQ(helper::GetPtr(View), nullptr);
    EXPECT_THAT(helper::GetExtents(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper::GetExtents(View).End(), ElementsAre(0,0,0));
    EXPECT_THAT(helper::GetIndexer(View).Begin(), ElementsAre(0,0,0));
    EXPECT_THAT(helper::GetIndexer(View).Stride(), ElementsAre(0,0,1));
  }

}

TEST_F(ArrayViewTests, Equality) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using array_view_const = ovk::array_view<const int,3>;
  using multidim_array = multidim_array_row<int>;

  // Self
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View(Array);
    EXPECT_TRUE(View == View);
    EXPECT_FALSE(View != View);
  }

  // Other with same values
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View1(Array);
    array_view View2(Array);
    EXPECT_TRUE(View1 == View2);
    EXPECT_FALSE(View1 != View2);
  }

  // Other with same values, different const
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View1(Array);
    array_view_const View2(Array);
    EXPECT_TRUE(View1 == View2);
    EXPECT_FALSE(View1 != View2);
    EXPECT_TRUE(View2 == View1);
    EXPECT_FALSE(View2 != View1);
  }

  // Different pointer
  {
    multidim_array Array1({{1,2,3}, {4,5,6}});
    multidim_array Array2({{1,2,3}, {4,5,6}});
    array_view View1(Array1);
    array_view View2(Array2);
    EXPECT_FALSE(View1 == View2);
    EXPECT_TRUE(View1 != View2);
  }

  // Different begin
  {
    multidim_array Array1({{1,2,3}, {4,5,6}});
    multidim_array Array2({{1,3,2}, {4,5,6}});
    array_view View1(Array1);
    array_view View2(Array2);
    EXPECT_FALSE(View1 == View2);
    EXPECT_TRUE(View1 != View2);
  }

  // Different end
  {
    multidim_array Array1({{1,2,3}, {4,5,6}});
    multidim_array Array2({{1,2,3}, {4,6,5}});
    array_view View1(Array1);
    array_view View2(Array2);
    EXPECT_FALSE(View1 == View2);
    EXPECT_TRUE(View1 != View2);
  }

  // Nullptr
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View1;
    array_view View2(Array);
    EXPECT_TRUE(View1 == nullptr);
    EXPECT_FALSE(View1 != nullptr);
    EXPECT_TRUE(nullptr == View1);
    EXPECT_FALSE(nullptr != View1);
    EXPECT_FALSE(View2 == nullptr);
    EXPECT_TRUE(View2 != nullptr);
    EXPECT_FALSE(nullptr == View2);
    EXPECT_TRUE(nullptr != View2);
  }

}

TEST_F(ArrayViewTests, ConvertToBool) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using multidim_array = multidim_array_row<int>;

  // Null
  {
    array_view View;
    EXPECT_FALSE(static_cast<bool>(View));
  }

  // Not null
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View(Array);
    EXPECT_TRUE(static_cast<bool>(View));
  }

}

TEST_F(ArrayViewTests, Extents) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using multidim_array = multidim_array_row<int>;

  multidim_array Array({{1,2,3}, {4,5,6}});
  array_view View(Array);
  EXPECT_THAT(View.Extents().Begin(), ElementsAre(1,2,3));
  EXPECT_THAT(View.Extents().End(), ElementsAre(4,5,6));

}

TEST_F(ArrayViewTests, Size) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using multidim_array = multidim_array_row<int>;

  multidim_array Array({{1,2,3}, {4,5,6}});
  array_view View(Array);
  EXPECT_THAT(View.Size(), ElementsAre(3,3,3));
  EXPECT_EQ(View.Size(0), 3);
  EXPECT_EQ(View.Size(1), 3);
  EXPECT_EQ(View.Size(2), 3);

}

TEST_F(ArrayViewTests, Count) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using multidim_array = multidim_array_row<int>;

  multidim_array Array({{1,2,3}, {4,5,6}});
  array_view View(Array);
  EXPECT_EQ(View.Count(), 27);

}

TEST_F(ArrayViewTests, Empty) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using multidim_array = multidim_array_row<int>;

  // Empty
  {
    multidim_array Array({{1,2,3}, {4,2,6}});
    array_view View(Array);
    EXPECT_TRUE(View.Empty());
  }

  // Non-empty
  {
    multidim_array Array({{1,2,3}, {4,5,6}});
    array_view View(Array);
    EXPECT_FALSE(View.Empty());
  }

}

TEST_F(ArrayViewTests, Indexer) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using multidim_array = multidim_array_row<int>;

  multidim_array Array({{1,2,3}, {4,5,6}});
  array_view View(Array);
  EXPECT_THAT(View.Indexer().Begin(), ElementsAre(1,2,3));
  EXPECT_THAT(View.Indexer().Stride(), ElementsAre(9,3,1));

}

TEST_F(ArrayViewTests, ParenthesisOperator) {

  if (TestComm().Rank() != 0) return;

  using array_view_row = ovk::array_view<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_view_col = ovk::array_view<int,3,ovk::array_layout::COLUMN_MAJOR>;

  // Row major, iterator
  {
    multidim_array_row<int> Array({{1,1,1}, {4,5,6}});
    array_view_row View(Array);
    ovk::elem<int,3> Tuple = {2,3,4};
    EXPECT_EQ(&View(Tuple.Data()), &Array[33]);
  }

  // Row major, array
  {
    multidim_array_row<int> Array({{1,1,1}, {4,5,6}});
    array_view_row View(Array);
    EXPECT_EQ(&View(ovk::elem<int,3>(2,3,4)), &Array[33]);
  }

  // Row major, separate
  {
    multidim_array_row<int> Array({{1,1,1}, {4,5,6}});
    array_view_row View(Array);
    EXPECT_EQ(&View(2,3,4), &Array[33]);
  }

  // Column major, iterator
  {
    multidim_array_col<int> Array({{1,1,1}, {4,5,6}});
    array_view_col View(Array);
    ovk::elem<int,3> Tuple = {2,3,4};
    EXPECT_EQ(&View(Tuple.Data()), &Array[43]);
  }

  // Column major, array
  {
    multidim_array_col<int> Array({{1,1,1}, {4,5,6}});
    array_view_col View(Array);
    EXPECT_EQ(&View(ovk::elem<int,3>(2,3,4)), &Array[43]);
  }

  // Column major, separate
  {
    multidim_array_col<int> Array({{1,1,1}, {4,5,6}});
    array_view_col View(Array);
    EXPECT_EQ(&View(2,3,4), &Array[43]);
  }

}

TEST_F(ArrayViewTests, BracketOperator) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;

  multidim_array_row<int> Array({{1,1,1}, {4,5,6}});
  array_view View(Array);
  EXPECT_EQ(&View[33], &Array[33]);

}

TEST_F(ArrayViewTests, Data) {

  if (TestComm().Rank() != 0) return;

  using array_view_row = ovk::array_view<int,3,ovk::array_layout::ROW_MAJOR>;
  using array_view_col = ovk::array_view<int,3,ovk::array_layout::COLUMN_MAJOR>;

  // Row major, no argument
  {
    multidim_array_row<int> Array({{1,1,1}, {4,5,6}});
    array_view_row View(Array);
    EXPECT_EQ(View.Data(), &Array[0]);
  }

  // Row major, iterator
  {
    multidim_array_row<int> Array({{1,1,1}, {4,5,6}});
    array_view_row View(Array);
    ovk::elem<int,3> Tuple = {2,3,4};
    EXPECT_EQ(View.Data(Tuple.Data()), &Array[33]);
  }

  // Row major, array
  {
    multidim_array_row<int> Array({{1,1,1}, {4,5,6}});
    array_view_row View(Array);
    EXPECT_EQ(View.Data(ovk::elem<int,3>(2,3,4)), &Array[33]);
  }

  // Row major, separate
  {
    multidim_array_row<int> Array({{1,1,1}, {4,5,6}});
    array_view_row View(Array);
    EXPECT_EQ(View.Data(2,3,4), &Array[33]);
  }

  // Column major, no argument
  {
    multidim_array_col<int> Array({{1,1,1}, {4,5,6}});
    array_view_col View(Array);
    EXPECT_EQ(View.Data(), &Array[0]);
  }

  // Column major, iterator
  {
    multidim_array_col<int> Array({{1,1,1}, {4,5,6}});
    array_view_col View(Array);
    ovk::elem<int,3> Tuple = {2,3,4};
    EXPECT_EQ(View.Data(Tuple.Data()), &Array[43]);
  }

  // Column major, array
  {
    multidim_array_col<int> Array({{1,1,1}, {4,5,6}});
    array_view_col View(Array);
    EXPECT_EQ(View.Data(ovk::elem<int,3>(2,3,4)), &Array[43]);
  }

  // Column major, separate
  {
    multidim_array_col<int> Array({{1,1,1}, {4,5,6}});
    array_view_col View(Array);
    EXPECT_EQ(View.Data(2,3,4), &Array[43]);
  }

}

TEST_F(ArrayViewTests, BeginEnd) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using array_view_const = ovk::array_view<const int,3>;
  using multidim_array = multidim_array_row<int>;

  // Const
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 1);
    array_view_const View(Array);
    EXPECT_EQ(View.Begin().Pointer(), &Array[0]);
    EXPECT_EQ(View.End().Pointer(), &Array[0] + 27);
    int Sum = 0;
    for (auto &Value : View) Sum += Value;
    EXPECT_EQ(Sum, 27);
  }

  // Non-const
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    EXPECT_EQ(View.Begin().Pointer(), &Array[0]);
    EXPECT_EQ(View.End().Pointer(), &Array[0] + 27);
    for (auto &Value : View) Value = 1;
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 27);
  }

}

TEST_F(ArrayViewTests, Fill) {

  if (TestComm().Rank() != 0) return;

  using array_view = ovk::array_view<int,3>;
  using array_view_noncopyable = ovk::array_view<noncopyable<int>>;
  using multidim_array = multidim_array_row<int>;

  // Constant value
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    View.Fill(1);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 27);
  }

  // Interval and constant value
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    View.Fill({{2,3,4}, {4,5,6}}, 1);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 8);
  }

  // Initializer list
  {
    multidim_array Array({{1,2,3}}, 0);
    array_view View(Array);
    View.Fill({0,1,2,3,4,5});
    int Sum = 0;
    for (int i = 0; i < 6; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 15);
  }

  // Iterator
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    std::array<int,27> SourceArray;
    for (int i = 0; i < 27; ++i) SourceArray[i] = i;
    View.Fill(SourceArray.begin());
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 351);
  }

  // View
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    std::array<int,27> SourceArray;
    for (int i = 0; i < 27; ++i) SourceArray[i] = i;
    array_view SourceView(SourceArray.data(), {{1,2,3}, {4,5,6}});
    View.Fill(SourceView);
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 351);
  }

  // Array, const ref
  {
    multidim_array Array({{1,2,3}, {4,5,6}}, 0);
    array_view View(Array);
    multidim_array SourceArray({{1,2,3}, {4,5,6}});
    for (int i = 0; i < 27; ++i) SourceArray[i] = i;
    View.Fill(static_cast<const multidim_array &>(SourceArray));
    int Sum = 0;
    for (int i = 0; i < 27; ++i) Sum += Array[i];
    EXPECT_EQ(Sum, 351);
  }

  // Array, rvalue ref
  {
    std::array<noncopyable<int>,4> Array;
    array_view_noncopyable View(Array);
    std::array<noncopyable<int>,4> SourceArray;
    for (int i = 0; i < 4; ++i) SourceArray[i] = 1;
    View.Fill(std::move(SourceArray));
    int Sum = 0;
    for (int i = 0; i < 4; ++i) Sum += Array[i].Value();
    EXPECT_EQ(Sum, 4);
    Sum = 0;
    for (int i = 0; i < 4; ++i) Sum += SourceArray[i].Value();
    EXPECT_EQ(Sum, 0);
  }

}
