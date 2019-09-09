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

template <typename T> using collect_none_col = collect_none<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectNoneCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_none_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_none_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_none_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_none_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_none_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_none_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_none_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_none_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_none_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_none_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  return Collect;

}

template <typename T> using collect_any_col = collect_any<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectAnyCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_any_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_any_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_any_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_any_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_any_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_any_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_any_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_any_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_any_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_any_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  return Collect;

}

template <typename T> using collect_not_all_col = collect_not_all<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectNotAllCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_not_all_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_not_all_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_not_all_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_not_all_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_not_all_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_not_all_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_not_all_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_not_all_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_not_all_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_not_all_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  return Collect;

}

template <typename T> using collect_all_col = collect_all<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectAllCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_all_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_all_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_all_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_all_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_all_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_all_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_all_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_all_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_all_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_all_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  return Collect;

}

template <typename T> using collect_min_col = collect_min<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectMinCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_min_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_min_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_min_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_min_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_min_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_min_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_min_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_min_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_min_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_min_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  return Collect;

}

template <typename T> using collect_max_col = collect_max<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectMaxCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange) {

  collect Collect;

  switch (ValueType) {
  case data_type::BOOL:
    Collect = collect_max_col<bool>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::BYTE:
    Collect = collect_max_col<unsigned char>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::INT:
    Collect = collect_max_col<int>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG:
    Collect = collect_max_col<long>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::LONG_LONG:
    Collect = collect_max_col<long long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_INT:
    Collect = collect_max_col<unsigned int>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG:
    Collect = collect_max_col<unsigned long>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Collect = collect_max_col<unsigned long long>(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, Count, FieldValuesRange);
    break;
  case data_type::FLOAT:
    Collect = collect_max_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  case data_type::DOUBLE:
    Collect = collect_max_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  return Collect;

}

template <typename T> using collect_interp_col = collect_interp<T, array_layout::COLUMN_MAJOR>;

collect CreateCollectInterpCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, floating_ref<const array<double,3>> InterpCoefs) {

  collect Collect;

  switch (ValueType) {
  case data_type::FLOAT:
    Collect = collect_interp_col<float>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange, InterpCoefs);
    break;
  case data_type::DOUBLE:
    Collect = collect_interp_col<double>(std::move(Context), Comm, Cart, LocalRange, CollectMap,
      Count, FieldValuesRange, InterpCoefs);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
    break;
  }

  return Collect;

}

}}}
