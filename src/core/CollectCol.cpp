// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Collect.hpp"

#include "ovk/core/CollectAll.hpp"
#include "ovk/core/CollectAny.hpp"
#include "ovk/core/CollectInterp.hpp"
#include "ovk/core/CollectNone.hpp"
#include "ovk/core/CollectNotAll.hpp"

#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Range.hpp"

namespace ovk {
namespace core {
namespace collect_internal {

template <typename T> using collect_none_col = collect_none<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectNoneCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_none_col<bool>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::BYTE:
    return collect_none_col<unsigned char>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::INT:
    return collect_none_col<int>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG:
    return collect_none_col<long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG_LONG:
    return collect_none_col<long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_INT:
    return collect_none_col<unsigned int>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG:
    return collect_none_col<unsigned long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_none_col<unsigned long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::FLOAT:
    return collect_none_col<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::DOUBLE:
    return collect_none_col<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  }

  return {};

}

template <typename T> using collect_any_col = collect_any<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectAnyCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_any_col<bool>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::BYTE:
    return collect_any_col<unsigned char>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::INT:
    return collect_any_col<int>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG:
    return collect_any_col<long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG_LONG:
    return collect_any_col<long long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::UNSIGNED_INT:
    return collect_any_col<unsigned int>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG:
    return collect_any_col<unsigned long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_any_col<unsigned long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::FLOAT:
    return collect_any_col<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::DOUBLE:
    return collect_any_col<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  }

  return {};

}

template <typename T> using collect_not_all_col = collect_not_all<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectNotAllCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_not_all_col<bool>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::BYTE:
    return collect_not_all_col<unsigned char>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::INT:
    return collect_not_all_col<int>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG:
    return collect_not_all_col<long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG_LONG:
    return collect_not_all_col<long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_INT:
    return collect_not_all_col<unsigned int>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG:
    return collect_not_all_col<unsigned long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_not_all_col<unsigned long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::FLOAT:
    return collect_not_all_col<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::DOUBLE:
    return collect_not_all_col<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  }

  return {};

}

template <typename T> using collect_all_col = collect_all<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectAllCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_all_col<bool>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::BYTE:
    return collect_all_col<unsigned char>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::INT:
    return collect_all_col<int>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG:
    return collect_all_col<long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG_LONG:
    return collect_all_col<long long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::UNSIGNED_INT:
    return collect_all_col<unsigned int>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG:
    return collect_all_col<unsigned long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_all_col<unsigned long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::FLOAT:
    return collect_all_col<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::DOUBLE:
    return collect_all_col<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  }

  return {};

}

template <typename T> using collect_interp_col = collect_interp<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectInterpCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler, array_view<const double,3> InterpCoefs) {

  switch (ValueType) {
  case data_type::FLOAT:
    return collect_interp_col<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler, InterpCoefs);
  case data_type::DOUBLE:
    return collect_interp_col<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler, InterpCoefs);
  default:
    OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
  }

  return {};

}

}}}
