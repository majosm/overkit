// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Collect.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/CollectBase.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

namespace {

template <typename T, array_layout Layout> class collect_none : public collect_base_for_type<T,
  Layout> {

protected:

  using parent_type = collect_base_for_type<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_none() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");

    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");

    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(true);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) = DonorValues(iCount)(iDonor) && !DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;

};

template <typename T, array_layout Layout> class collect_any : public collect_base_for_type<T,
  Layout> {

protected:

  using parent_type = collect_base_for_type<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_any() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");
    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(false);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) = DonorValues(iCount)(iDonor) || DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;

};

template <typename T, array_layout Layout> class collect_not_all : public collect_base_for_type<T,
  Layout> {

protected:

  using parent_type = collect_base_for_type<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_not_all() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");
    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(false);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) = DonorValues(iCount)(iDonor) || !DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;

};

template <typename T, array_layout Layout> class collect_all : public collect_base_for_type<T,
  Layout> {

protected:

  using parent_type = collect_base_for_type<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_all() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");
    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(true);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) = DonorValues(iCount)(iDonor) && DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;

};

template <typename T, array_layout Layout> class collect_interp : public collect_base_for_type<T,
  Layout> {

protected:

  using parent_type = collect_base_for_type<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Donors_;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_interp() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    InterpCoefs_ = Donors_->InterpCoefs_;

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});
    DonorPointCoefs_.Resize({MaxPointsInCell_});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");
    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      donor_indexer Indexer(DonorSize);

      for (int k = 0; k < DonorSize[2]; ++k) {
        for (int j = 0; j < DonorSize[1]; ++j) {
          for (int i = 0; i < DonorSize[0]; ++i) {
            int iPointInCell = Indexer.ToIndex(i,j,k);
            DonorPointCoefs_(iPointInCell) =
              InterpCoefs_(0,i,iDonor) *
              InterpCoefs_(1,j,iDonor) *
              InterpCoefs_(2,k,iDonor);
          }
        }
      }

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(0);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) += DonorPointCoefs_(iPointInCell)*DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array_view<const double,3> InterpCoefs_;
  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;
  array<double> DonorPointCoefs_;

};

template <typename T> using collect_none_row = collect_none<T, array_layout::ROW_MAJOR>;
template <typename T> using collect_none_col = collect_none<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectNone(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_none_row<bool>();
    case data_type::BYTE: return collect_none_row<unsigned char>();
    case data_type::INT: return collect_none_row<int>();
    case data_type::LONG: return collect_none_row<long>();
    case data_type::LONG_LONG: return collect_none_row<long long>();
    case data_type::UNSIGNED_INT: return collect_none_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_none_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_none_row<unsigned long long>();
    case data_type::FLOAT: return collect_none_row<float>();
    case data_type::DOUBLE: return collect_none_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_none_col<bool>();
    case data_type::BYTE: return collect_none_col<unsigned char>();
    case data_type::INT: return collect_none_col<int>();
    case data_type::LONG: return collect_none_col<long>();
    case data_type::LONG_LONG: return collect_none_col<long long>();
    case data_type::UNSIGNED_INT: return collect_none_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_none_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_none_col<unsigned long long>();
    case data_type::FLOAT: return collect_none_col<float>();
    case data_type::DOUBLE: return collect_none_col<double>();
    }
  }

  return {};

}

template <typename T> using collect_any_row = collect_any<T, array_layout::ROW_MAJOR>;
template <typename T> using collect_any_col = collect_any<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectAny(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_any_row<bool>();
    case data_type::BYTE: return collect_any_row<unsigned char>();
    case data_type::INT: return collect_any_row<int>();
    case data_type::LONG: return collect_any_row<long>();
    case data_type::LONG_LONG: return collect_any_row<long long>();
    case data_type::UNSIGNED_INT: return collect_any_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_any_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_any_row<unsigned long long>();
    case data_type::FLOAT: return collect_any_row<float>();
    case data_type::DOUBLE: return collect_any_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_any_col<bool>();
    case data_type::BYTE: return collect_any_col<unsigned char>();
    case data_type::INT: return collect_any_col<int>();
    case data_type::LONG: return collect_any_col<long>();
    case data_type::LONG_LONG: return collect_any_col<long long>();
    case data_type::UNSIGNED_INT: return collect_any_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_any_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_any_col<unsigned long long>();
    case data_type::FLOAT: return collect_any_col<float>();
    case data_type::DOUBLE: return collect_any_col<double>();
    }
  }

  return {};

}

template <typename T> using collect_not_all_row = collect_not_all<T, array_layout::ROW_MAJOR>;
template <typename T> using collect_not_all_col = collect_not_all<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectNotAll(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_not_all_row<bool>();
    case data_type::BYTE: return collect_not_all_row<unsigned char>();
    case data_type::INT: return collect_not_all_row<int>();
    case data_type::LONG: return collect_not_all_row<long>();
    case data_type::LONG_LONG: return collect_not_all_row<long long>();
    case data_type::UNSIGNED_INT: return collect_not_all_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_not_all_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_not_all_row<unsigned long long>();
    case data_type::FLOAT: return collect_not_all_row<float>();
    case data_type::DOUBLE: return collect_not_all_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_not_all_col<bool>();
    case data_type::BYTE: return collect_not_all_col<unsigned char>();
    case data_type::INT: return collect_not_all_col<int>();
    case data_type::LONG: return collect_not_all_col<long>();
    case data_type::LONG_LONG: return collect_not_all_col<long long>();
    case data_type::UNSIGNED_INT: return collect_not_all_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_not_all_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_not_all_col<unsigned long long>();
    case data_type::FLOAT: return collect_not_all_col<float>();
    case data_type::DOUBLE: return collect_not_all_col<double>();
    }
  }

  return {};

}

template <typename T> using collect_all_row = collect_all<T, array_layout::ROW_MAJOR>;
template <typename T> using collect_all_col = collect_all<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectAll(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_all_row<bool>();
    case data_type::BYTE: return collect_all_row<unsigned char>();
    case data_type::INT: return collect_all_row<int>();
    case data_type::LONG: return collect_all_row<long>();
    case data_type::LONG_LONG: return collect_all_row<long long>();
    case data_type::UNSIGNED_INT: return collect_all_row<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_all_row<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_all_row<unsigned long long>();
    case data_type::FLOAT: return collect_all_row<float>();
    case data_type::DOUBLE: return collect_all_row<double>();
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::BOOL: return collect_all_col<bool>();
    case data_type::BYTE: return collect_all_col<unsigned char>();
    case data_type::INT: return collect_all_col<int>();
    case data_type::LONG: return collect_all_col<long>();
    case data_type::LONG_LONG: return collect_all_col<long long>();
    case data_type::UNSIGNED_INT: return collect_all_col<unsigned int>();
    case data_type::UNSIGNED_LONG: return collect_all_col<unsigned long>();
    case data_type::UNSIGNED_LONG_LONG: return collect_all_col<unsigned long long>();
    case data_type::FLOAT: return collect_all_col<float>();
    case data_type::DOUBLE: return collect_all_col<double>();
    }
  }

  return {};

}

template <typename T> using collect_interp_row = collect_interp<T, array_layout::ROW_MAJOR>;
template <typename T> using collect_interp_col = collect_interp<T, array_layout::COLUMN_MAJOR>;

collect MakeCollectInterp(data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    switch (ValueType) {
    case data_type::FLOAT: return collect_interp_row<float>();
    case data_type::DOUBLE: return collect_interp_row<double>();
    default:
      OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
    }
  case array_layout::COLUMN_MAJOR:
    switch (ValueType) {
    case data_type::FLOAT: return collect_interp_col<float>();
    case data_type::DOUBLE: return collect_interp_col<double>();
    default:
      OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
    }
  }

  return {};

}

}

collect MakeCollect(collect_op CollectOp, data_type ValueType, array_layout Layout) {

  switch (CollectOp) {
  case collect_op::NONE:
    return MakeCollectNone(ValueType, Layout);
  case collect_op::ANY:
    return MakeCollectAny(ValueType, Layout);
  case collect_op::NOT_ALL:
    return MakeCollectNotAll(ValueType, Layout);
  case collect_op::ALL:
    return MakeCollectAll(ValueType, Layout);
  case collect_op::INTERPOLATE:
    return MakeCollectInterp(ValueType, Layout);
  }

  return {};

}

}}
