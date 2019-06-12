// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISPERSE_BASE_HPP_INCLUDED
#define OVK_CORE_DISPERSE_BASE_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

template <array_layout Layout> class disperse_base {

public:

  disperse_base(std::shared_ptr<context> &&Context, const array<int,2> &Points, int Count, const
    range &FieldValuesRange);

  // Can't define these here due to issues with GCC < 6.3 and Intel < 17
  // implementations of extern template
//   disperse_base(const disperse_base &Other) = delete;
//   disperse_base(disperse_base &&Other) noexcept = default;

//   disperse_base &operator=(const disperse_base &Other) = delete;
//   disperse_base &operator=(disperse_base &&Other) noexcept = default;

protected:

  using range_indexer = indexer<long long, int, MAX_DIMS, Layout>;

  std::shared_ptr<context> Context_;

  array_view<const int,2> Points_;
  int Count_;
  range FieldValuesRange_;
  range_indexer FieldValuesIndexer_;

};

extern template class disperse_base<array_layout::ROW_MAJOR>;
extern template class disperse_base<array_layout::COLUMN_MAJOR>;

template <typename T, array_layout Layout> class disperse_base_for_type : public disperse_base<
  Layout> {

private:

  using parent_type = disperse_base<Layout>;

public:

  using value_type = T;

  disperse_base_for_type(std::shared_ptr<context> &&Context, const array<int,2> &Points, int Count,
    const range &FieldValuesRange);

protected:

  using typename parent_type::range_indexer;

  using parent_type::Context_;
  using parent_type::Points_;
  using parent_type::Count_;
  using parent_type::FieldValuesRange_;
  using parent_type::FieldValuesIndexer_;

  array<array_view<const value_type>> PackedValues_;
  array<array_view<value_type>> FieldValues_;

  void SetBufferViews(const void * const *PackedValuesVoid, void **FieldValuesVoid);

};

extern template class disperse_base_for_type<bool, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<bool, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<unsigned char, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<unsigned char, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<int, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<int, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<long, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<long, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<long long, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<long long, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<unsigned int, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<unsigned int, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<unsigned long, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<unsigned long, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<unsigned long long, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<unsigned long long, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<float, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<float, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<double, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<double, array_layout::COLUMN_MAJOR>;

}}

#endif
