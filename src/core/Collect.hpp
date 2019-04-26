// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_HPP_INCLUDED
#define OVK_CORE_COLLECT_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Connectivity.hpp>
#include <ovk/core/ConnectivityD.hpp>
#include <ovk/core/DataType.hpp>
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
    Collect_(new model<T>(std::forward<T>(Collect)))
  {}

  collect(const collect &Other) = delete;
  collect(collect &&Other) noexcept = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, collect>::value)>
    collect &operator=(T &&Collect) {
    Collect_.reset(new model<T>(std::forward<T>(Collect)));
    return *this;
  }

  collect &operator=(const collect &Other) = delete;
  collect &operator=(collect &&Other) noexcept = default;

  void Collect(const void * const *FieldValues, void **PackedValues) {
    Collect_->Collect(FieldValues, PackedValues);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Collect(const void * const *FieldValues, void **PackedValues) = 0;
  };

  template <typename T> class model : public concept {
  public:
    model(T Collect):
      Collect_(std::move(Collect))
    {}
    virtual void Collect(const void * const *FieldValues, void **PackedValues) override {
      Collect_.Collect(FieldValues, PackedValues);
    }
  private:
    T Collect_;
  };

  std::unique_ptr<concept> Collect_;

};

namespace collect_internal {
collect MakeCollectNoneRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler);
collect MakeCollectNoneCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler);
}
inline collect MakeCollectNone(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  array_layout FieldValuesLayout, profiler &Profiler) {
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    return collect_internal::MakeCollectNoneRow(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler);
  case array_layout::COLUMN_MAJOR:
    return collect_internal::MakeCollectNoneCol(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler);
  }
  return {};
}

namespace collect_internal {
collect MakeCollectAnyRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler);
collect MakeCollectAnyCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler);
}
inline collect MakeCollectAny(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  array_layout FieldValuesLayout, profiler &Profiler) {
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    return collect_internal::MakeCollectAnyRow(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler);
  case array_layout::COLUMN_MAJOR:
    return collect_internal::MakeCollectAnyCol(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler);
  }
  return {};
}

namespace collect_internal {
collect MakeCollectNotAllRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler);
collect MakeCollectNotAllCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler);
}
inline collect MakeCollectNotAll(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  array_layout FieldValuesLayout, profiler &Profiler) {
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    return collect_internal::MakeCollectNotAllRow(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler);
  case array_layout::COLUMN_MAJOR:
    return collect_internal::MakeCollectNotAllCol(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler);
  }
  return {};
}

namespace collect_internal {
collect MakeCollectAllRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler);
collect MakeCollectAllCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler);
}
inline collect MakeCollectAll(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  array_layout FieldValuesLayout, profiler &Profiler) {
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    return collect_internal::MakeCollectAllRow(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler);
  case array_layout::COLUMN_MAJOR:
    return collect_internal::MakeCollectAllCol(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler);
  }
  return {};
}

namespace collect_internal {
collect MakeCollectInterpRow(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler, array_view<const double,3> InterpCoefs);
collect MakeCollectInterpCol(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  profiler &Profiler, array_view<const double,3> InterpCoefs);
}
inline collect MakeCollectInterp(comm_view Comm, const cart &Cart, const range &LocalRange, const
  collect_map &CollectMap, data_type ValueType, int Count, const range &FieldValuesRange,
  array_layout FieldValuesLayout, profiler &Profiler, array_view<const double,3> InterpCoefs) {
  switch (FieldValuesLayout) {
  case array_layout::ROW_MAJOR:
    return collect_internal::MakeCollectInterpRow(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler, InterpCoefs);
  case array_layout::COLUMN_MAJOR:
    return collect_internal::MakeCollectInterpCol(Comm, Cart, LocalRange, CollectMap, ValueType,
      Count, FieldValuesRange, Profiler, InterpCoefs);
  }
  return {};
}

}}

#endif
