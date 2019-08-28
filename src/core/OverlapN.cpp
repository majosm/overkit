// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/OverlapN.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Event.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Logger.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

namespace overlap_n_internal {

overlap_n_base::overlap_n_base(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
  &&SourceGridInfo):
  Context_(std::move(Context)),
  Grid_(&Grid),
  SourceGridInfo_(std::move(SourceGridInfo)),
  Comm_(Grid_->Comm())
{
  MPI_Barrier(Comm_);
}

overlap_n_base::~overlap_n_base() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Logger.LogDebug(Comm_.Rank() == 0, 0, "Destroyed N side of overlap (%s,%s).",
      SourceGridInfo_.Name(), Grid_->Name());
  }

}

}

overlap_n::overlap_n(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
  &&SourceGridInfo):
  overlap_n_base(std::move(Context), Grid, std::move(SourceGridInfo)),
  NumDims_(Grid_->Dimension()),
  Count_(0),
  Mask_(Grid_->ExtendedRange(), false),
  Points_({{MAX_DIMS,0}}),
  Sources_({{MAX_DIMS,0}}),
  SourceRanks_({0})
{

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();
  Logger.LogDebug(Comm_.Rank() == 0, 0, "Created N side of overlap (%s,%s).",
    SourceGridInfo_.Name(), Grid_->Name());

}

overlap_n::~overlap_n() noexcept {

  if (Context_) {
    // Barrier before cleaning up
    MPI_Barrier(Comm_);
  }

}

overlap_n overlap_n::internal_Create(std::shared_ptr<context> &&Context, const grid &Grid,
  grid_info &&SourceGridInfo) {

  return {std::move(Context), Grid, std::move(SourceGridInfo)};

}

namespace core {

overlap_n CreateOverlapN(std::shared_ptr<context> Context, const grid &Grid, grid_info
  SourceGridInfo) {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return overlap_n::internal_Create(std::move(Context), Grid, std::move(SourceGridInfo));

}

}

void overlap_n::Resize(long long Count) {

  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");

  MPI_Barrier(Comm_);

  OVK_DEBUG_ASSERT(!PointsEditor_.Active(), "Cannot resize while editing points.");
  OVK_DEBUG_ASSERT(!SourcesEditor_.Active(), "Cannot resize while editing sources.");
  OVK_DEBUG_ASSERT(!SourceRanksEditor_.Active(), "Cannot resize while editing source ranks.");

  Count_ = Count;

  Points_.Resize({{MAX_DIMS,Count}});
  Sources_.Resize({{MAX_DIMS,Count}});
  SourceRanks_.Resize({Count});

  for (long long iPoint = 0; iPoint < Count_; ++iPoint) {
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      Points_(iDim,iPoint) = Grid_->GlobalRange().Begin(iDim)-1;
    }
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      Points_(iDim,iPoint) = 0;
    }
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      Sources_(iDim,iPoint) = SourceGridInfo_.CellGlobalRange().Begin(iDim)-1;
    }
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      Sources_(iDim,iPoint) = 0;
    }
    SourceRanks_(iPoint) = -1;
  }

  MPI_Barrier(Comm_);

  ResizeEvent_.Trigger();
  PointsEvent_.Trigger();
  SourcesEvent_.Trigger();
  SourceRanksEvent_.Trigger();

  MPI_Barrier(Comm_);

}

bool overlap_n::EditingPoints() const {

  return PointsEditor_.Active();

}

edit_handle<array<int,2>> overlap_n::EditPoints() {

  if (!PointsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<overlap_n> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      overlap_n &OverlapN = *FloatingRef;
      OverlapN.UpdateMask_();
      MPI_Barrier(OverlapN.Comm_);
      OverlapN.PointsEvent_.Trigger();
      MPI_Barrier(OverlapN.Comm_);
    };
    PointsEditor_.Activate(std::move(DeactivateFunc));
  }

  return PointsEditor_.Edit(Points_);

}

void overlap_n::RestorePoints() {

  OVK_DEBUG_ASSERT(PointsEditor_.Active(), "Unable to restore points; not currently being edited.");

  PointsEditor_.Restore();

}

void overlap_n::UpdateMask_() {

  const partition &Partition = Grid_->Partition();

  Mask_.Fill(false);

  tuple<int> GlobalBeginMinusOne = Grid_->GlobalRange().Begin();
  for (int iDim = 0; iDim < NumDims_; ++iDim) {
    GlobalBeginMinusOne(iDim) -= 1;
  }

  for (long long iPoint = 0; iPoint < Count_; ++iPoint) {
    tuple<int> Point = {
      Points_(0,iPoint),
      Points_(1,iPoint),
      Points_(2,iPoint)
    };
    if (Point != GlobalBeginMinusOne) {
      Mask_(Point) = true;
    }
  }

  Partition.Exchange(Mask_);

}

bool overlap_n::EditingSources() const {

  return SourcesEditor_.Active();

}

edit_handle<array<int,2>> overlap_n::EditSources() {

  if (!SourcesEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<overlap_n> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      overlap_n &OverlapN = *FloatingRef;
      MPI_Barrier(OverlapN.Comm_);
      OverlapN.SourcesEvent_.Trigger();
      MPI_Barrier(OverlapN.Comm_);
    };
    SourcesEditor_.Activate(std::move(DeactivateFunc));
  }

  return SourcesEditor_.Edit(Sources_);

}

void overlap_n::RestoreSources() {

  OVK_DEBUG_ASSERT(SourcesEditor_.Active(), "Unable to restore sources; not currently being "
    "edited.");

  SourcesEditor_.Restore();

}

bool overlap_n::EditingSourceRanks() const {

  return SourceRanksEditor_.Active();

}

edit_handle<array<int>> overlap_n::EditSourceRanks() {

  if (!SourceRanksEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<overlap_n> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      overlap_n &OverlapN = *FloatingRef;
      MPI_Barrier(OverlapN.Comm_);
      OverlapN.SourceRanksEvent_.Trigger();
      MPI_Barrier(OverlapN.Comm_);
    };
    SourceRanksEditor_.Activate(std::move(DeactivateFunc));
  }

  return SourceRanksEditor_.Edit(SourceRanks_);

}

void overlap_n::RestoreSourceRanks() {

  OVK_DEBUG_ASSERT(SourceRanksEditor_.Active(), "Unable to restore source ranks; not currently "
    "being edited.");

  SourceRanksEditor_.Restore();

}

}
