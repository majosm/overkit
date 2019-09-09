// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Collect.hpp"

#include "ovk/core/CollectAll.hpp"
#include "ovk/core/CollectAny.hpp"
#include "ovk/core/CollectInterp.hpp"
#include "ovk/core/CollectMax.hpp"
#include "ovk/core/CollectMin.hpp"
#include "ovk/core/CollectNone.hpp"
#include "ovk/core/CollectNotAll.hpp"

#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {
namespace collect_internal {

template <typename T> using collect_none_row = collect_none<T, array_layout::ROW_MAJOR>;

collect CreateCollectNoneRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_none_row<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_none_row<unsigned char>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_none_row<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_none_row<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_none_row<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_none_row<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_none_row<unsigned long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_none_row<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_none_row<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_none_row<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  }

  return Collect;

}

template <typename T> using collect_any_row = collect_any<T, array_layout::ROW_MAJOR>;

collect CreateCollectAnyRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_any_row<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_any_row<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_any_row<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_any_row<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_any_row<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_any_row<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_any_row<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_any_row<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_any_row<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_any_row<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  }

  return Collect;

}

template <typename T> using collect_not_all_row = collect_not_all<T, array_layout::ROW_MAJOR>;

collect CreateCollectNotAllRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_not_all_row<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_not_all_row<unsigned char>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_not_all_row<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_not_all_row<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_not_all_row<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_not_all_row<unsigned int>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_not_all_row<unsigned long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_not_all_row<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_not_all_row<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_not_all_row<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  }

  return Collect;

}

template <typename T> using collect_all_row = collect_all<T, array_layout::ROW_MAJOR>;

collect CreateCollectAllRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_all_row<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_all_row<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_all_row<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_all_row<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_all_row<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_all_row<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_all_row<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_all_row<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_all_row<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_all_row<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  }

  return Collect;

}

template <typename T> using collect_min_row = collect_min<T, array_layout::ROW_MAJOR>;

collect CreateCollectMinRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_min_row<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_min_row<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_min_row<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_min_row<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_min_row<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_min_row<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_min_row<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_min_row<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_min_row<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_min_row<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  }

  return Collect;

}

template <typename T> using collect_max_row = collect_max<T, array_layout::ROW_MAJOR>;

collect CreateCollectMaxRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_max_row<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_max_row<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_max_row<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_max_row<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_max_row<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_max_row<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_max_row<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_max_row<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_max_row<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_max_row<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  }

  return Collect;

}

template <typename T> using collect_interp_row = collect_interp<T, array_layout::ROW_MAJOR>;

collect CreateCollectInterpRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, floating_ref<const array<double,3>> InterpCoefs) {

  collect Collect;

  switch (ValueType) {
  case data_type::FLOAT:
    Collect = collect_interp_row<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange, InterpCoefs);
    break;
  case data_type::DOUBLE:
    Collect = collect_interp_row<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange, InterpCoefs);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
    break;
  }

  return Collect;

}

}}}
