// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONNECTIVITY_N_HPP_INCLUDED
#define OVK_CORE_CONNECTIVITY_N_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Editor.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

namespace connectivity_n_internal {

// For doing stuff before creation and after destruction
class connectivity_n_base {

protected:

  connectivity_n_base(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
    &&SourceGridInfo);

  connectivity_n_base(const connectivity_n_base &Other) = delete;
  connectivity_n_base(connectivity_n_base &&Other) noexcept = default;

  connectivity_n_base &operator=(const connectivity_n_base &Other) = delete;
  connectivity_n_base &operator=(connectivity_n_base &&Other) noexcept = default;

  ~connectivity_n_base() noexcept;

  std::shared_ptr<context> Context_;

  const grid *Grid_;

  grid_info SourceGridInfo_;

  comm_view Comm_;

};

}

class connectivity_n : private connectivity_n_internal::connectivity_n_base {

public:

  connectivity_n(const connectivity_n &Other) = delete;
  connectivity_n(connectivity_n &&Other) noexcept = default;

  connectivity_n &operator=(const connectivity_n &Other) = delete;
  connectivity_n &operator=(connectivity_n &&Other) = default;

  ~connectivity_n() noexcept;

  floating_ref<const connectivity_n> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<connectivity_n> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const grid &Grid() const { return *Grid_; }

  const grid_info &SourceGridInfo() const { return SourceGridInfo_; }
  
  int Dimension() const { return NumDims_; }

  comm_view Comm() const { return Comm_; }

  long long Size() const { return NumReceivers_; }

  void Resize(long long NumReceivers);
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddResizeEventListener(F Listener) const {
    return ResizeEvent_.AddListener(std::move(Listener));
  }

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

  static connectivity_n internal_Create(std::shared_ptr<context> &&Context, const grid &Grid,
    grid_info &&SourceGridInfo);

private:

  floating_ref_generator FloatingRefGenerator_;

  int NumDims_;

  long long NumReceivers_;
  mutable event<void()> ResizeEvent_;

  array<int,2> Points_;
  editor PointsEditor_;
  mutable event<void()> PointsEvent_;

  array<int,2> Sources_;
  editor SourcesEditor_;
  mutable event<void()> SourcesEvent_;

  array<int> SourceRanks_;
  editor SourceRanksEditor_;
  mutable event<void()> SourceRanksEvent_;

  connectivity_n(std::shared_ptr<context> &&Context, const grid &Grid, grid_info &&SourceGridInfo);

};

namespace core {
connectivity_n CreateConnectivityN(std::shared_ptr<context> Context, const grid &Grid, grid_info
  SourceGridInfo);
}

}

#endif
