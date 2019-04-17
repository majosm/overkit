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

template <typename T> using collect_none_row = collect_none<T, array_layout::ROW_MAJOR>;

collect MakeCollectNoneRow(data_type ValueType) {

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

  return {};

}

template <typename T> using collect_any_row = collect_any<T, array_layout::ROW_MAJOR>;

collect MakeCollectAnyRow(data_type ValueType) {

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

  return {};

}

template <typename T> using collect_not_all_row = collect_not_all<T, array_layout::ROW_MAJOR>;

collect MakeCollectNotAllRow(data_type ValueType) {

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

  return {};

}

template <typename T> using collect_all_row = collect_all<T, array_layout::ROW_MAJOR>;

collect MakeCollectAllRow(data_type ValueType) {

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

  return {};

}

template <typename T> using collect_interp_row = collect_interp<T, array_layout::ROW_MAJOR>;

collect MakeCollectInterpRow(data_type ValueType) {

  switch (ValueType) {
  case data_type::FLOAT: return collect_interp_row<float>();
  case data_type::DOUBLE: return collect_interp_row<double>();
  default:
    OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
  }

  return {};

}

}

collect MakeCollectRow(collect_op CollectOp, data_type ValueType) {

  switch (CollectOp) {
  case collect_op::NONE:
    return MakeCollectNoneRow(ValueType);
  case collect_op::ANY:
    return MakeCollectAnyRow(ValueType);
  case collect_op::NOT_ALL:
    return MakeCollectNotAllRow(ValueType);
  case collect_op::ALL:
    return MakeCollectAllRow(ValueType);
  case collect_op::INTERPOLATE:
    return MakeCollectInterpRow(ValueType);
  }

  return {};

}

}}}
