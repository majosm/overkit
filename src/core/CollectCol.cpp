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
namespace collect_internal {

namespace {

template <typename T> using collect_none_col = collect_none<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectNoneCol(data_type ValueType) {

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

  return {};

}

template <typename T> using collect_any_col = collect_any<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectAnyCol(data_type ValueType) {

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

  return {};

}

template <typename T> using collect_not_all_col = collect_not_all<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectNotAllCol(data_type ValueType) {

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

  return {};

}

template <typename T> using collect_all_col = collect_all<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectAllCol(data_type ValueType) {

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

  return {};

}

template <typename T> using collect_interp_col = collect_interp<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectInterpCol(data_type ValueType) {

  switch (ValueType) {
  case data_type::FLOAT: return collect_interp_col<float>();
  case data_type::DOUBLE: return collect_interp_col<double>();
  default:
    OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
  }

  return {};

}

}

collect MakeCollectCol(collect_op CollectOp, data_type ValueType) {

  switch (CollectOp) {
  case collect_op::NONE:
    return MakeCollectNoneCol(ValueType);
  case collect_op::ANY:
    return MakeCollectAnyCol(ValueType);
  case collect_op::NOT_ALL:
    return MakeCollectNotAllCol(ValueType);
  case collect_op::ALL:
    return MakeCollectAllCol(ValueType);
  case collect_op::INTERPOLATE:
    return MakeCollectInterpCol(ValueType);
  }

  return {};

}

}}}
