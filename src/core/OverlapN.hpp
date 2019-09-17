// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_OVERLAP_N_HPP_INCLUDED
#define OVK_CORE_OVERLAP_N_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Editor.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

namespace overlap_n_internal {

// For doing stuff before creation and after destruction
class overlap_n_base {

protected:

  overlap_n_base(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
    &&SourceGridInfo);

  overlap_n_base(const overlap_n_base &Other) = delete;
  overlap_n_base(overlap_n_base &&Other) noexcept = default;

  overlap_n_base &operator=(const overlap_n_base &Other) = delete;
  overlap_n_base &operator=(overlap_n_base &&Other) noexcept = default;

  ~overlap_n_base() noexcept;

  std::shared_ptr<context> Context_;

  const grid *Grid_;

  grid_info SourceGridInfo_;

  comm_view Comm_;

};

}

class overlap_n : private overlap_n_internal::overlap_n_base {

public:

  overlap_n(const overlap_n &Other) = delete;
  overlap_n(overlap_n &&Other) noexcept = default;

  overlap_n &operator=(const overlap_n &Other) = delete;
  overlap_n &operator=(overlap_n &&Other) = default;

  ~overlap_n() noexcept;

  floating_ref<const overlap_n> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<overlap_n> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const grid &Grid() const { return *Grid_; }

  const grid_info &SourceGridInfo() const { return SourceGridInfo_; }
  
  int Dimension() const { return NumDims_; }

  comm_view Comm() const { return Comm_; }

  long long Size() const { return NumPoints_; }

  void Resize(long long NumPoints);
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddResizeEventListener(F Listener) const {
    return ResizeEvent_.AddListener(std::move(Listener));
  }

  const field<bool> &Mask() const { return Mask_; }

  const array<int,2> &Points() const { return Points_; }
  bool EditingPoints() const;
  edit_handle<array<int,2>> EditPoints();
  void RestorePoints();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddPointsEventListener(F Listener) const {
    return PointsEvent_.AddListener(std::move(Listener));
  }

  const array<int,2> &Sources() const { return Sources_; }
  bool EditingSources() const;
  edit_handle<array<int,2>> EditSources();
  void RestoreSources();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddSourcesEventListener(F Listener) const {
    return SourcesEvent_.AddListener(std::move(Listener));
  }

  const array<int> &SourceRanks() const { return SourceRanks_; }
  bool EditingSourceRanks() const;
  edit_handle<array<int>> EditSourceRanks();
  void RestoreSourceRanks();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddSourceRanksEventListener(F Listener) const {
    return SourceRanksEvent_.AddListener(std::move(Listener));
  }

  static overlap_n internal_Create(std::shared_ptr<context> &&Context, const grid &Grid,
    grid_info &&SourceGridInfo);

private:

  floating_ref_generator FloatingRefGenerator_;

  int NumDims_;

  long long NumPoints_;
  mutable event<void()> ResizeEvent_;

  field<bool> Mask_;

  array<int,2> Points_;
  editor PointsEditor_;
  mutable event<void()> PointsEvent_;

  array<int,2> Sources_;
  editor SourcesEditor_;
  mutable event<void()> SourcesEvent_;

  array<int> SourceRanks_;
  editor SourceRanksEditor_;
  mutable event<void()> SourceRanksEvent_;

  overlap_n(std::shared_ptr<context> &&Context, const grid &Grid, grid_info &&SourceGridInfo);

  void UpdateMask_();

};

namespace core {
overlap_n CreateOverlapN(std::shared_ptr<context> Context, const grid &Grid, grid_info
  SourceGridInfo);
}

}

#endif
