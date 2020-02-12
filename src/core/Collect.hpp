// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_HPP_INCLUDED
#define OVK_CORE_COLLECT_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

class collect {

public:

  collect() = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, collect>::value)>
    collect(T &&Collect):
    Collect_(new model<remove_cvref<T>>(std::forward<T>(Collect)))
  {}

  collect(const collect &Other) = delete;
  collect(collect &&Other) noexcept = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, collect>::value)>
    collect &operator=(T &&Collect) {
    Collect_.reset(new model<remove_cvref<T>>(std::forward<T>(Collect)));
    return *this;
  }

  collect &operator=(const collect &Other) = delete;
  collect &operator=(collect &&Other) noexcept = default;

  void Collect(const void *FieldValues, void *PackedValues) {
    Collect_->Collect(FieldValues, PackedValues);
  }

private:

  class concept {
  public:
    virtual ~concept() noexcept {}
    virtual void Collect(const void *FieldValues, void *PackedValues) = 0;
  };

  template <typename T> class model final : public concept {
  public:
    model(T Collect):
      Collect_(std::move(Collect))
    {}
    virtual void Collect(const void *FieldValues, void *PackedValues) override {
      Collect_.Collect(FieldValues, PackedValues);
    }
  private:
    T Collect_;
  };

  std::unique_ptr<concept> Collect_;

};

namespace collect_internal {
collect CreateCollectNoneRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
collect CreateCollectNoneCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
}
inline collect CreateCollectNone(std::shared_ptr<context> Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, array_layout FieldValuesLayout) {
  collect Collect;
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    Collect = collect_internal::CreateCollectNoneRow(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  case array_layout::COLUMN_MAJOR:
    Collect = collect_internal::CreateCollectNoneCol(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }
  return Collect;
}

namespace collect_internal {
collect CreateCollectAnyRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
collect CreateCollectAnyCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
}
inline collect CreateCollectAny(std::shared_ptr<context> Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, array_layout FieldValuesLayout) {
  collect Collect;
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    Collect = collect_internal::CreateCollectAnyRow(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  case array_layout::COLUMN_MAJOR:
    Collect = collect_internal::CreateCollectAnyCol(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }
  return Collect;
}

namespace collect_internal {
collect CreateCollectNotAllRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
collect CreateCollectNotAllCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
}
inline collect CreateCollectNotAll(std::shared_ptr<context> Context, comm_view Comm, const cart
  &Cart, const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count,
  const range &FieldValuesRange, array_layout FieldValuesLayout) {
  collect Collect;
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    Collect = collect_internal::CreateCollectNotAllRow(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  case array_layout::COLUMN_MAJOR:
    Collect = collect_internal::CreateCollectNotAllCol(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }
  return Collect;
}

namespace collect_internal {
collect CreateCollectAllRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
collect CreateCollectAllCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
}
inline collect CreateCollectAll(std::shared_ptr<context> Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, array_layout FieldValuesLayout) {
  collect Collect;
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    Collect = collect_internal::CreateCollectAllRow(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  case array_layout::COLUMN_MAJOR:
    Collect = collect_internal::CreateCollectAllCol(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }
  return Collect;
}

namespace collect_internal {
collect CreateCollectMinRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
collect CreateCollectMinCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
}
inline collect CreateCollectMin(std::shared_ptr<context> Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, array_layout FieldValuesLayout) {
  collect Collect;
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    Collect = collect_internal::CreateCollectMinRow(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  case array_layout::COLUMN_MAJOR:
    Collect = collect_internal::CreateCollectMinCol(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }
  return Collect;
}

namespace collect_internal {
collect CreateCollectMaxRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
collect CreateCollectMaxCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange);
}
inline collect CreateCollectMax(std::shared_ptr<context> Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, array_layout FieldValuesLayout) {
  collect Collect;
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    Collect = collect_internal::CreateCollectMaxRow(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  case array_layout::COLUMN_MAJOR:
    Collect = collect_internal::CreateCollectMaxCol(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }
  return Collect;
}

namespace collect_internal {
collect CreateCollectInterpRow(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, floating_ref<const array<double,3>> InterpCoefs);
collect CreateCollectInterpCol(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
  const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count, const
  range &FieldValuesRange, floating_ref<const array<double,3>> InterpCoefs);
}
inline collect CreateCollectInterp(std::shared_ptr<context> Context, comm_view Comm, const cart
  &Cart, const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int Count,
  const range &FieldValuesRange, array_layout FieldValuesLayout, floating_ref<const array<double,3>>
  InterpCoefs) {
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    return collect_internal::CreateCollectInterpRow(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange, InterpCoefs);
    break;
  case array_layout::COLUMN_MAJOR:
    return collect_internal::CreateCollectInterpCol(std::move(Context), Comm, Cart, LocalRange,
      CollectMap, ValueType, Count, FieldValuesRange, InterpCoefs);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }
  return {};
}

#ifdef OVK_HAVE_OPENMP
namespace collect_internal {
collect CreateCollectInterpThreadedRow(std::shared_ptr<context> &&Context, comm_view Comm, const
  cart &Cart, const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int
  Count, const range &FieldValuesRange, floating_ref<const array<double,3>> InterpCoefs);
collect CreateCollectInterpThreadedCol(std::shared_ptr<context> &&Context, comm_view Comm, const
  cart &Cart, const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int
  Count, const range &FieldValuesRange, floating_ref<const array<double,3>> InterpCoefs);
}
inline collect CreateCollectInterpThreaded(std::shared_ptr<context> Context, comm_view Comm, const
  cart &Cart, const range &LocalRange, const collect_map &CollectMap, data_type ValueType, int
  Count, const range &FieldValuesRange, array_layout FieldValuesLayout, floating_ref<const
  array<double,3>> InterpCoefs) {
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    return collect_internal::CreateCollectInterpThreadedRow(std::move(Context), Comm, Cart,
      LocalRange, CollectMap, ValueType, Count, FieldValuesRange, InterpCoefs);
    break;
  case array_layout::COLUMN_MAJOR:
    return collect_internal::CreateCollectInterpThreadedCol(std::move(Context), Comm, Cart,
      LocalRange, CollectMap, ValueType, Count, FieldValuesRange, InterpCoefs);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }
  return {};
}
#endif

}}

#endif
