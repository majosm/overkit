// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Disperse.hpp"

#include "ovk/core/DisperseOverwrite.hpp"
#include "ovk/core/DisperseAppend.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

namespace {
template <typename T> using disperse_overwrite_row = disperse_internal::disperse_overwrite<T,
  array_layout::ROW_MAJOR>;
template <typename T> using disperse_overwrite_col = disperse_internal::disperse_overwrite<T,
  array_layout::COLUMN_MAJOR>;
}

disperse CreateDisperseOverwrite(std::shared_ptr<context> Context, const array<int,2> &Points,
  data_type ValueType, int Count, const range &FieldValuesRange, array_layout FieldValuesLayout) {

  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      return disperse_overwrite_row<bool>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::BYTE:
      return disperse_overwrite_row<unsigned char>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::INT:
      return disperse_overwrite_row<int>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::LONG:
      return disperse_overwrite_row<long>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::LONG_LONG:
      return disperse_overwrite_row<long long>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::UNSIGNED_INT:
      return disperse_overwrite_row<unsigned int>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::UNSIGNED_LONG:
      return disperse_overwrite_row<unsigned long>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::UNSIGNED_LONG_LONG:
      return disperse_overwrite_row<unsigned long long>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::FLOAT:
      return disperse_overwrite_row<float>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::DOUBLE:
      return disperse_overwrite_row<double>(std::move(Context), Points, Count, FieldValuesRange);
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      return disperse_overwrite_col<bool>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::BYTE:
      return disperse_overwrite_col<unsigned char>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::INT:
      return disperse_overwrite_col<int>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::LONG:
      return disperse_overwrite_col<long>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::LONG_LONG:
      return disperse_overwrite_col<long long>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::UNSIGNED_INT:
      return disperse_overwrite_col<unsigned int>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::UNSIGNED_LONG:
      return disperse_overwrite_col<unsigned long>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::UNSIGNED_LONG_LONG:
      return disperse_overwrite_col<unsigned long long>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::FLOAT:
      return disperse_overwrite_col<float>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::DOUBLE:
      return disperse_overwrite_col<double>(std::move(Context), Points, Count, FieldValuesRange);
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

disperse CreateDisperseAppend(std::shared_ptr<context> Context, const array<int,2> &Points,
  data_type ValueType, int Count, const range &FieldValuesRange, array_layout FieldValuesLayout) {

  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      return disperse_append_row<bool>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::BYTE:
      return disperse_append_row<unsigned char>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::INT:
      return disperse_append_row<int>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::LONG:
      return disperse_append_row<long>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::LONG_LONG:
      return disperse_append_row<long long>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::UNSIGNED_INT:
      return disperse_append_row<unsigned int>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::UNSIGNED_LONG:
      return disperse_append_row<unsigned long>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::UNSIGNED_LONG_LONG:
      return disperse_append_row<unsigned long long>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::FLOAT:
      return disperse_append_row<float>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::DOUBLE:
      return disperse_append_row<double>(std::move(Context), Points, Count, FieldValuesRange);
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      return disperse_append_col<bool>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::BYTE:
      return disperse_append_col<unsigned char>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::INT:
      return disperse_append_col<int>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::LONG:
      return disperse_append_col<long>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::LONG_LONG:
      return disperse_append_col<long long>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::UNSIGNED_INT:
      return disperse_append_col<unsigned int>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::UNSIGNED_LONG:
      return disperse_append_col<unsigned long>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::UNSIGNED_LONG_LONG:
      return disperse_append_col<unsigned long long>(std::move(Context), Points, Count,
        FieldValuesRange);
    case data_type::FLOAT:
      return disperse_append_col<float>(std::move(Context), Points, Count, FieldValuesRange);
    case data_type::DOUBLE:
      return disperse_append_col<double>(std::move(Context), Points, Count, FieldValuesRange);
    }
  }

  return {};

}

}}
