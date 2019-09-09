// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Disperse.hpp"

#include "ovk/core/DisperseOverwrite.hpp"
#include "ovk/core/DisperseAppend.hpp"

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

disperse CreateDisperseOverwrite(std::shared_ptr<context> Context, const disperse_map &DisperseMap,
  data_type ValueType, int Count, const range &FieldValuesRange, array_layout FieldValuesLayout) {

  disperse Disperse;

  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      Disperse = disperse_overwrite_row<bool>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::BYTE:
      Disperse = disperse_overwrite_row<unsigned char>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::INT:
      Disperse = disperse_overwrite_row<int>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::LONG:
      Disperse = disperse_overwrite_row<long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::LONG_LONG:
      Disperse = disperse_overwrite_row<long long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_INT:
      Disperse = disperse_overwrite_row<unsigned int>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_LONG:
      Disperse = disperse_overwrite_row<unsigned long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_LONG_LONG:
      Disperse = disperse_overwrite_row<unsigned long long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::FLOAT:
      Disperse = disperse_overwrite_row<float>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::DOUBLE:
      Disperse = disperse_overwrite_row<double>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      break;
    }
    break;
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      Disperse = disperse_overwrite_col<bool>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::BYTE:
      Disperse = disperse_overwrite_col<unsigned char>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::INT:
      Disperse = disperse_overwrite_col<int>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::LONG:
      Disperse = disperse_overwrite_col<long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::LONG_LONG:
      Disperse = disperse_overwrite_col<long long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_INT:
      Disperse = disperse_overwrite_col<unsigned int>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_LONG:
      Disperse = disperse_overwrite_col<unsigned long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_LONG_LONG:
      Disperse = disperse_overwrite_col<unsigned long long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::FLOAT:
      Disperse = disperse_overwrite_col<float>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::DOUBLE:
      Disperse = disperse_overwrite_col<double>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      break;
    }
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  return Disperse;

}

namespace {
template <typename T> using disperse_append_row = disperse_internal::disperse_append<T,
  array_layout::ROW_MAJOR>;
template <typename T> using disperse_append_col = disperse_internal::disperse_append<T,
  array_layout::COLUMN_MAJOR>;
}

disperse CreateDisperseAppend(std::shared_ptr<context> Context, const disperse_map &DisperseMap,
  data_type ValueType, int Count, const range &FieldValuesRange, array_layout FieldValuesLayout) {

  disperse Disperse;

  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      Disperse = disperse_append_row<bool>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::BYTE:
      Disperse = disperse_append_row<unsigned char>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::INT:
      Disperse = disperse_append_row<int>(std::move(Context), DisperseMap, Count, FieldValuesRange);
      break;
    case data_type::LONG:
      Disperse = disperse_append_row<long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::LONG_LONG:
      Disperse = disperse_append_row<long long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_INT:
      Disperse = disperse_append_row<unsigned int>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_LONG:
      Disperse = disperse_append_row<unsigned long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_LONG_LONG:
      Disperse = disperse_append_row<unsigned long long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::FLOAT:
      Disperse = disperse_append_row<float>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::DOUBLE:
      Disperse = disperse_append_row<double>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      break;
    }
    break;
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL:
      Disperse = disperse_append_col<bool>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::BYTE:
      Disperse = disperse_append_col<unsigned char>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::INT:
      Disperse = disperse_append_col<int>(std::move(Context), DisperseMap, Count, FieldValuesRange);
      break;
    case data_type::LONG:
      Disperse = disperse_append_col<long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::LONG_LONG:
      Disperse = disperse_append_col<long long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_INT:
      Disperse = disperse_append_col<unsigned int>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_LONG:
      Disperse = disperse_append_col<unsigned long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::UNSIGNED_LONG_LONG:
      Disperse = disperse_append_col<unsigned long long>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::FLOAT:
      Disperse = disperse_append_col<float>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    case data_type::DOUBLE:
      Disperse = disperse_append_col<double>(std::move(Context), DisperseMap, Count,
        FieldValuesRange);
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      break;
    }
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  return Disperse;

}

}}
