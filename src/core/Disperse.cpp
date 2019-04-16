// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Disperse.hpp"

#include "ovk/core/DisperseOverwrite.hpp"
#include "ovk/core/DisperseAppend.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

namespace ovk {
namespace core {

namespace {

template <typename T> using disperse_overwrite_row = disperse_internal::disperse_overwrite<T,
  array_layout::ROW_MAJOR>;
template <typename T> using disperse_overwrite_col = disperse_internal::disperse_overwrite<T,
  array_layout::COLUMN_MAJOR>;

disperse MakeDisperseOverwrite(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return disperse_overwrite_row<bool>();
    case data_type::BYTE: return disperse_overwrite_row<unsigned char>();
    case data_type::INT: return disperse_overwrite_row<int>();
    case data_type::LONG: return disperse_overwrite_row<long>();
    case data_type::LONG_LONG: return disperse_overwrite_row<long long>();
    case data_type::UNSIGNED_INT: return disperse_overwrite_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return disperse_overwrite_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return disperse_overwrite_row<unsigned long long>();
    case data_type::FLOAT: return disperse_overwrite_row<float>();
    case data_type::DOUBLE: return disperse_overwrite_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return disperse_overwrite_col<bool>();
    case data_type::BYTE: return disperse_overwrite_col<unsigned char>();
    case data_type::INT: return disperse_overwrite_col<int>();
    case data_type::LONG: return disperse_overwrite_col<long>();
    case data_type::LONG_LONG: return disperse_overwrite_col<long long>();
    case data_type::UNSIGNED_INT: return disperse_overwrite_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return disperse_overwrite_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return disperse_overwrite_col<unsigned long long>();
    case data_type::FLOAT: return disperse_overwrite_col<float>();
    case data_type::DOUBLE: return disperse_overwrite_col<double>();
    }
  }

  return {};

}

template <typename T> using disperse_append_row = disperse_internal::disperse_append<T,
  array_layout::ROW_MAJOR>;
template <typename T> using disperse_append_col = disperse_internal::disperse_append<T,
  array_layout::COLUMN_MAJOR>;

disperse MakeDisperseAppend(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return disperse_append_row<bool>();
    case data_type::BYTE: return disperse_append_row<unsigned char>();
    case data_type::INT: return disperse_append_row<int>();
    case data_type::LONG: return disperse_append_row<long>();
    case data_type::LONG_LONG: return disperse_append_row<long long>();
    case data_type::UNSIGNED_INT: return disperse_append_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return disperse_append_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return disperse_append_row<unsigned long long>();
    case data_type::FLOAT: return disperse_append_row<float>();
    case data_type::DOUBLE: return disperse_append_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return disperse_append_col<bool>();
    case data_type::BYTE: return disperse_append_col<unsigned char>();
    case data_type::INT: return disperse_append_col<int>();
    case data_type::LONG: return disperse_append_col<long>();
    case data_type::LONG_LONG: return disperse_append_col<long long>();
    case data_type::UNSIGNED_INT: return disperse_append_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return disperse_append_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return disperse_append_col<unsigned long long>();
    case data_type::FLOAT: return disperse_append_col<float>();
    case data_type::DOUBLE: return disperse_append_col<double>();
    }
  }

  return {};

}

}

disperse MakeDisperse(disperse_op DisperseOp, data_type ValueType, array_layout Layout) {

  switch (DisperseOp) {
  case disperse_op::OVERWRITE:
    return MakeDisperseOverwrite(ValueType, Layout);
  case disperse_op::APPEND:
    return MakeDisperseAppend(ValueType, Layout);
  }

  return {};

}

}}
