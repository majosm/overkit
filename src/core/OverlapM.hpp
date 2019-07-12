// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_OVERLAP_M_HPP_INCLUDED
#define OVK_CORE_OVERLAP_M_HPP_INCLUDED

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

namespace overlap_m_internal {

// For doing stuff before creation and after destruction
class overlap_m_base {

protected:

  overlap_m_base(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
    &&DestinationGridInfo);

  overlap_m_base(const overlap_m_base &Other) = delete;
  overlap_m_base(overlap_m_base &&Other) noexcept = default;

  overlap_m_base &operator=(const overlap_m_base &Other) = delete;
  overlap_m_base &operator=(overlap_m_base &&Other) noexcept = default;

  ~overlap_m_base() noexcept;

  std::shared_ptr<context> Context_;

  const grid *Grid_;

  grid_info DestinationGridInfo_;

  comm_view Comm_;

};

}

class overlap_m : private overlap_m_internal::overlap_m_base {

public:

  overlap_m(const overlap_m &Other) = delete;
  overlap_m(overlap_m &&Other) noexcept = default;

  overlap_m &operator=(const overlap_m &Other) = delete;
  overlap_m &operator=(overlap_m &&Other) noexcept = default;

  ~overlap_m() noexcept;

  floating_ref<const overlap_m> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<overlap_m> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const grid &Grid() const { return *Grid_; }

  const grid_info &DestinationGridInfo() const { return DestinationGridInfo_; }

  int Dimension() const { return NumDims_; }

  comm_view Comm() const { return Comm_; }

  long long Count() const { return Count_; }

  void Resize(long long Count);
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddResizeEventListener(F Listener) const {
    return ResizeEvent_.AddListener(std::move(Listener));
  }

  const array<int,2> &Cells() const { return Cells_; }
  bool EditingCells() const;
  edit_handle<array<int,2>> EditCells();
  void RestoreCells();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddCellsEventListener(F Listener) const {
    return CellsEvent_.AddListener(std::move(Listener));
  }

  const array<double,2> &Coords() const { return Coords_; }
  bool EditingCoords() const;
  edit_handle<array<double,2>> EditCoords();
  void RestoreCoords();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddCoordsEventListener(F Listener) const {
    return CoordsEvent_.AddListener(std::move(Listener));
  }

  const array<int,2> &Destinations() const { return Destinations_; }
  bool EditingDestinations() const;
  edit_handle<array<int,2>> EditDestinations();
  void RestoreDestinations();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddDestinationsEventListener(F Listener) const {
    return DestinationsEvent_.AddListener(std::move(Listener));
  }

  const array<int> &DestinationRanks() const { return DestinationRanks_; }
  bool EditingDestinationRanks() const;
  edit_handle<array<int>> EditDestinationRanks();
  void RestoreDestinationRanks();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddDestinationRanksEventListener(F Listener) const {
    return DestinationRanksEvent_.AddListener(std::move(Listener));
  }

  static overlap_m internal_Create(std::shared_ptr<context> &&Context, const grid &Grid,
    grid_info &&DestinationGridInfo);

private:

  floating_ref_generator FloatingRefGenerator_;

  int NumDims_;

  long long Count_;
  mutable event<void()> ResizeEvent_;

  array<int,2> Cells_;
  editor CellsEditor_;
  mutable event<void()> CellsEvent_;

  array<double,2> Coords_;
  editor CoordsEditor_;
  mutable event<void()> CoordsEvent_;

  array<int,2> Destinations_;
  editor DestinationsEditor_;
  mutable event<void()> DestinationsEvent_;

  array<int> DestinationRanks_;
  editor DestinationRanksEditor_;
  mutable event<void()> DestinationRanksEvent_;

  overlap_m(std::shared_ptr<context> &&Context, const grid &Grid, grid_info &&DestinationGridInfo);

};

namespace core {
overlap_m CreateOverlapM(std::shared_ptr<context> Context, const grid &Grid, grid_info
  DestinationGridInfo);
}

}

#endif
