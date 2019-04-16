// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Collect.hpp"

#include "ovk/core/CollectAll.hpp"
#include "ovk/core/CollectAny.hpp"
#include "ovk/core/CollectInterp.hpp"
#include "ovk/core/CollectNone.hpp"
#include "ovk/core/CollectNotAll.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

namespace ovk {
namespace core {

namespace {

template <typename T> using collect_none_row = collect_internal::collect_none<T,
  array_layout::ROW_MAJOR>;
template <typename T> using collect_none_col = collect_internal::collect_none<T,
  array_layout::COLUMN_MAJOR>;

collect MakeCollectNone(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_none_row<bool>();
    case data_type::BYTE: return collect_none_row<unsigned char>();
    case data_type::INT: return collect_none_row<int>();
    case data_type::LONG: return collect_none_row<long>();
    case data_type::LONG_LONG: return collect_none_row<long long>();
    case data_type::UNSIGNED_INT: return collect_none_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_none_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_none_row<unsigned long long>();
    case data_type::FLOAT: return collect_none_row<float>();
    case data_type::DOUBLE: return collect_none_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_none_col<bool>();
    case data_type::BYTE: return collect_none_col<unsigned char>();
    case data_type::INT: return collect_none_col<int>();
    case data_type::LONG: return collect_none_col<long>();
    case data_type::LONG_LONG: return collect_none_col<long long>();
    case data_type::UNSIGNED_INT: return collect_none_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_none_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_none_col<unsigned long long>();
    case data_type::FLOAT: return collect_none_col<float>();
    case data_type::DOUBLE: return collect_none_col<double>();
    }
  }

  return {};

}

template <typename T> using collect_any_row = collect_internal::collect_any<T,
  array_layout::ROW_MAJOR>;
template <typename T> using collect_any_col = collect_internal::collect_any<T,
  array_layout::COLUMN_MAJOR>;

collect MakeCollectAny(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_any_row<bool>();
    case data_type::BYTE: return collect_any_row<unsigned char>();
    case data_type::INT: return collect_any_row<int>();
    case data_type::LONG: return collect_any_row<long>();
    case data_type::LONG_LONG: return collect_any_row<long long>();
    case data_type::UNSIGNED_INT: return collect_any_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_any_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_any_row<unsigned long long>();
    case data_type::FLOAT: return collect_any_row<float>();
    case data_type::DOUBLE: return collect_any_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_any_col<bool>();
    case data_type::BYTE: return collect_any_col<unsigned char>();
    case data_type::INT: return collect_any_col<int>();
    case data_type::LONG: return collect_any_col<long>();
    case data_type::LONG_LONG: return collect_any_col<long long>();
    case data_type::UNSIGNED_INT: return collect_any_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_any_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_any_col<unsigned long long>();
    case data_type::FLOAT: return collect_any_col<float>();
    case data_type::DOUBLE: return collect_any_col<double>();
    }
  }

  return {};

}

template <typename T> using collect_not_all_row = collect_internal::collect_not_all<T,
  array_layout::ROW_MAJOR>;
template <typename T> using collect_not_all_col = collect_internal::collect_not_all<T,
  array_layout::COLUMN_MAJOR>;

collect MakeCollectNotAll(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_not_all_row<bool>();
    case data_type::BYTE: return collect_not_all_row<unsigned char>();
    case data_type::INT: return collect_not_all_row<int>();
    case data_type::LONG: return collect_not_all_row<long>();
    case data_type::LONG_LONG: return collect_not_all_row<long long>();
    case data_type::UNSIGNED_INT: return collect_not_all_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_not_all_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_not_all_row<unsigned long long>();
    case data_type::FLOAT: return collect_not_all_row<float>();
    case data_type::DOUBLE: return collect_not_all_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_not_all_col<bool>();
    case data_type::BYTE: return collect_not_all_col<unsigned char>();
    case data_type::INT: return collect_not_all_col<int>();
    case data_type::LONG: return collect_not_all_col<long>();
    case data_type::LONG_LONG: return collect_not_all_col<long long>();
    case data_type::UNSIGNED_INT: return collect_not_all_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_not_all_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_not_all_col<unsigned long long>();
    case data_type::FLOAT: return collect_not_all_col<float>();
    case data_type::DOUBLE: return collect_not_all_col<double>();
    }
  }

  return {};

}

template <typename T> using collect_all_row = collect_internal::collect_all<T,
  array_layout::ROW_MAJOR>;
template <typename T> using collect_all_col = collect_internal::collect_all<T,
  array_layout::COLUMN_MAJOR>;

collect MakeCollectAll(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_all_row<bool>();
    case data_type::BYTE: return collect_all_row<unsigned char>();
    case data_type::INT: return collect_all_row<int>();
    case data_type::LONG: return collect_all_row<long>();
    case data_type::LONG_LONG: return collect_all_row<long long>();
    case data_type::UNSIGNED_INT: return collect_all_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_all_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_all_row<unsigned long long>();
    case data_type::FLOAT: return collect_all_row<float>();
    case data_type::DOUBLE: return collect_all_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_all_col<bool>();
    case data_type::BYTE: return collect_all_col<unsigned char>();
    case data_type::INT: return collect_all_col<int>();
    case data_type::LONG: return collect_all_col<long>();
    case data_type::LONG_LONG: return collect_all_col<long long>();
    case data_type::UNSIGNED_INT: return collect_all_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_all_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_all_col<unsigned long long>();
    case data_type::FLOAT: return collect_all_col<float>();
    case data_type::DOUBLE: return collect_all_col<double>();
    }
  }

  return {};

}

template <typename T> using collect_interp_row = collect_internal::collect_interp<T,
  array_layout::ROW_MAJOR>;
template <typename T> using collect_interp_col = collect_internal::collect_interp<T,
  array_layout::COLUMN_MAJOR>;

collect MakeCollectInterp(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::FLOAT: return collect_interp_row<float>();
    case data_type::DOUBLE: return collect_interp_row<double>();
    default:
      OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::FLOAT: return collect_interp_col<float>();
    case data_type::DOUBLE: return collect_interp_col<double>();
    default:
      OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
    }
  }

  return {};

}

}

collect MakeCollect(collect_op CollectOp, data_type ValueType, array_layout Layout) {

  switch (CollectOp) {
  case collect_op::NONE:
    return MakeCollectNone(ValueType, Layout);
  case collect_op::ANY:
    return MakeCollectAny(ValueType, Layout);
  case collect_op::NOT_ALL:
    return MakeCollectNotAll(ValueType, Layout);
  case collect_op::ALL:
    return MakeCollectAll(ValueType, Layout);
  case collect_op::INTERPOLATE:
    return MakeCollectInterp(ValueType, Layout);
  }

  return {};

}

}}
