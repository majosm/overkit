// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Disperse.hpp"

#include "ovk/core/DisperseOverwrite.hpp"
#include "ovk/core/DisperseAppend.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"

namespace ovk {
namespace core {

namespace {
template <typename T> using disperse_overwrite_row = disperse_internal::disperse_overwrite<T,
  array_layout::ROW_MAJOR>;
template <typename T> using disperse_overwrite_col = disperse_internal::disperse_overwrite<T,
  array_layout::COLUMN_MAJOR>;
}

disperse MakeDisperseOverwrite(const array<int,2> &Points, data_type ValueType, int Count,
  const range &FieldValuesRange, array_layout FieldValuesLayout, profiler &Profiler) {

  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      return disperse_overwrite_row<bool>(Points, Count, FieldValuesRange, Profiler);
    case data_type::BYTE:
      return disperse_overwrite_row<unsigned char>(Points, Count, FieldValuesRange, Profiler);
    case data_type::INT:
      return disperse_overwrite_row<int>(Points, Count, FieldValuesRange, Profiler);
    case data_type::LONG:
      return disperse_overwrite_row<long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::LONG_LONG:
      return disperse_overwrite_row<long long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_INT:
      return disperse_overwrite_row<unsigned int>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_LONG:
      return disperse_overwrite_row<unsigned long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_LONG_LONG:
      return disperse_overwrite_row<unsigned long long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::FLOAT:
      return disperse_overwrite_row<float>(Points, Count, FieldValuesRange, Profiler);
    case data_type::DOUBLE:
      return disperse_overwrite_row<double>(Points, Count, FieldValuesRange, Profiler);
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      return disperse_overwrite_col<bool>(Points, Count, FieldValuesRange, Profiler);
    case data_type::BYTE:
      return disperse_overwrite_col<unsigned char>(Points, Count, FieldValuesRange, Profiler);
    case data_type::INT:
      return disperse_overwrite_col<int>(Points, Count, FieldValuesRange, Profiler);
    case data_type::LONG:
      return disperse_overwrite_col<long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::LONG_LONG:
      return disperse_overwrite_col<long long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_INT:
      return disperse_overwrite_col<unsigned int>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_LONG:
      return disperse_overwrite_col<unsigned long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_LONG_LONG:
      return disperse_overwrite_col<unsigned long long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::FLOAT:
      return disperse_overwrite_col<float>(Points, Count, FieldValuesRange, Profiler);
    case data_type::DOUBLE:
      return disperse_overwrite_col<double>(Points, Count, FieldValuesRange, Profiler);
    }
  }

  return {};

}

namespace {
template <typename T> using disperse_append_row = disperse_internal::disperse_append<T,
  array_layout::ROW_MAJOR>;
template <typename T> using disperse_append_col = disperse_internal::disperse_append<T,
  array_layout::COLUMN_MAJOR>;
}

disperse MakeDisperseAppend(const array<int,2> &Points, data_type ValueType, int Count,
  const range &FieldValuesRange, array_layout FieldValuesLayout, profiler &Profiler) {

  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      return disperse_append_row<bool>(Points, Count, FieldValuesRange, Profiler);
    case data_type::BYTE:
      return disperse_append_row<unsigned char>(Points, Count, FieldValuesRange, Profiler);
    case data_type::INT:
      return disperse_append_row<int>(Points, Count, FieldValuesRange, Profiler);
    case data_type::LONG:
      return disperse_append_row<long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::LONG_LONG:
      return disperse_append_row<long long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_INT:
      return disperse_append_row<unsigned int>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_LONG:
      return disperse_append_row<unsigned long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_LONG_LONG:
      return disperse_append_row<unsigned long long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::FLOAT:
      return disperse_append_row<float>(Points, Count, FieldValuesRange, Profiler);
    case data_type::DOUBLE:
      return disperse_append_row<double>(Points, Count, FieldValuesRange, Profiler);
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      return disperse_append_col<bool>(Points, Count, FieldValuesRange, Profiler);
    case data_type::BYTE:
      return disperse_append_col<unsigned char>(Points, Count, FieldValuesRange, Profiler);
    case data_type::INT:
      return disperse_append_col<int>(Points, Count, FieldValuesRange, Profiler);
    case data_type::LONG:
      return disperse_append_col<long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::LONG_LONG:
      return disperse_append_col<long long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_INT:
      return disperse_append_col<unsigned int>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_LONG:
      return disperse_append_col<unsigned long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::UNSIGNED_LONG_LONG:
      return disperse_append_col<unsigned long long>(Points, Count, FieldValuesRange, Profiler);
    case data_type::FLOAT:
      return disperse_append_col<float>(Points, Count, FieldValuesRange, Profiler);
    case data_type::DOUBLE:
      return disperse_append_col<double>(Points, Count, FieldValuesRange, Profiler);
    }
  }

  return {};

}

}}
