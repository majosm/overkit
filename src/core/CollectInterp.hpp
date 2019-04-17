// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_INTERP_HPP_INCLUDED
#define OVK_CORE_COLLECT_INTERP_HPP_INCLUDED

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/CollectBase.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/ConnectivityD.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

namespace ovk {
namespace core {
namespace collect_internal {

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
  using parent_type::GridValues_;
  using parent_type::DonorValues_;

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

  void Collect(const void * const *GridValuesVoid, void **DonorValuesVoid) {

    parent_type::SetBufferViews(GridValuesVoid, DonorValuesVoid);
    parent_type::RetrieveRemoteDonorValues(GridValues_, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");
    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues_, RemoteDonorValues_, iDonor, DonorSize,
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
        DonorValues_(iCount)(iDonor) = value_type(0);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues_(iCount)(iDonor) += DonorPointCoefs_(iPointInCell)*DonorPointValues_(iCount,
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

}}}

#endif
