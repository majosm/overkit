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
#include "ovk/core/Context.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {
namespace collect_internal {

template <typename T> using collect_none_col = collect_none<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectNoneCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_none_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::BYTE:
    return collect_none_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::INT:
    return collect_none_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::LONG:
    return collect_none_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::LONG_LONG:
    return collect_none_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::UNSIGNED_INT:
    return collect_none_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::UNSIGNED_LONG:
    return collect_none_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_none_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
  case data_type::FLOAT:
    return collect_none_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::DOUBLE:
    return collect_none_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  }

  return {};

}

template <typename T> using collect_any_col = collect_any<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectAnyCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_any_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::BYTE:
    return collect_any_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::INT:
    return collect_any_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::LONG:
    return collect_any_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::LONG_LONG:
    return collect_any_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::UNSIGNED_INT:
    return collect_any_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::UNSIGNED_LONG:
    return collect_any_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_any_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
  case data_type::FLOAT:
    return collect_any_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::DOUBLE:
    return collect_any_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    }

  return {};

}

template <typename T> using collect_not_all_col = collect_not_all<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectNotAllCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_not_all_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::BYTE:
    return collect_not_all_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
  case data_type::INT:
    return collect_not_all_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::LONG:
    return collect_not_all_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::LONG_LONG:
    return collect_not_all_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::UNSIGNED_INT:
    return collect_not_all_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::UNSIGNED_LONG:
    return collect_not_all_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_not_all_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
  case data_type::FLOAT:
    return collect_not_all_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::DOUBLE:
    return collect_not_all_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  }

  return {};

}

template <typename T> using collect_all_col = collect_all<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectAllCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  switch (ValueType) {
  case data_type::BOOL:
    return collect_all_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::BYTE:
    return collect_all_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::INT:
    return collect_all_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::LONG:
    return collect_all_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::LONG_LONG:
    return collect_all_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::UNSIGNED_INT:
    return collect_all_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::UNSIGNED_LONG:
    return collect_all_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
  case data_type::UNSIGNED_LONG_LONG:
    return collect_all_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
  case data_type::FLOAT:
    return collect_all_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  case data_type::DOUBLE:
    return collect_all_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
  }

  return {};

}

template <typename T> using collect_interp_col = collect_interp<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectInterpCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, array_view<const double,3> InterpCoefs) {

  switch (ValueType) {
  case data_type::FLOAT:
    return collect_interp_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, InterpCoefs);
  case data_type::DOUBLE:
    return collect_interp_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, InterpCoefs);
  default:
    OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
  }

  return {};

}

}}}
