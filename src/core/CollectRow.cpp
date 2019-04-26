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

template <typename T> using collect_none_row = collect_none<T, array_layout::ROW_MAJOR>;

collect MakeCollectNoneRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_none_row<bool>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::BYTE:
    return collect_none_row<unsigned char>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::INT:
    return collect_none_row<int>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG:
    return collect_none_row<long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG_LONG:
    return collect_none_row<long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_INT:
    return collect_none_row<unsigned int>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG:
    return collect_none_row<unsigned long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_none_row<unsigned long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::FLOAT:
    return collect_none_row<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::DOUBLE:
    return collect_none_row<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  }

  return {};

}

template <typename T> using collect_any_row = collect_any<T, array_layout::ROW_MAJOR>;

collect MakeCollectAnyRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_any_row<bool>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::BYTE:
    return collect_any_row<unsigned char>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::INT:
    return collect_any_row<int>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG:
    return collect_any_row<long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG_LONG:
    return collect_any_row<long long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::UNSIGNED_INT:
    return collect_any_row<unsigned int>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG:
    return collect_any_row<unsigned long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_any_row<unsigned long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::FLOAT:
    return collect_any_row<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::DOUBLE:
    return collect_any_row<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  }

  return {};

}

template <typename T> using collect_not_all_row = collect_not_all<T, array_layout::ROW_MAJOR>;

collect MakeCollectNotAllRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_not_all_row<bool>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::BYTE:
    return collect_not_all_row<unsigned char>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::INT:
    return collect_not_all_row<int>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG:
    return collect_not_all_row<long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG_LONG:
    return collect_not_all_row<long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_INT:
    return collect_not_all_row<unsigned int>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG:
    return collect_not_all_row<unsigned long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_not_all_row<unsigned long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::FLOAT:
    return collect_not_all_row<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::DOUBLE:
    return collect_not_all_row<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  }

  return {};

}

template <typename T> using collect_all_row = collect_all<T, array_layout::ROW_MAJOR>;

collect MakeCollectAllRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_all_row<bool>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::BYTE:
    return collect_all_row<unsigned char>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::INT:
    return collect_all_row<int>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG:
    return collect_all_row<long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::LONG_LONG:
    return collect_all_row<long long>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::UNSIGNED_INT:
    return collect_all_row<unsigned int>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG:
    return collect_all_row<unsigned long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_all_row<unsigned long long>(Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, Profiler);
  case data_type::FLOAT:
    return collect_all_row<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  case data_type::DOUBLE:
    return collect_all_row<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler);
  }

  return {};

}

template <typename T> using collect_interp_row = collect_interp<T, array_layout::ROW_MAJOR>;

collect MakeCollectInterpRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler, array_view<const double,3> InterpCoefs) {

  switch (ValueType) {
  case data_type::FLOAT:
    return collect_interp_row<float>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler, InterpCoefs);
  case data_type::DOUBLE:
    return collect_interp_row<double>(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      Profiler, InterpCoefs);
  default:
    OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
  }

  return {};

}

}}}
