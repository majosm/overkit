// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONNECTIVITY_M_HPP_INCLUDED
#define OVK_CORE_CONNECTIVITY_M_HPP_INCLUDED

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

namespace connectivity_m_internal {

// For doing stuff before creation and after destruction
class connectivity_m_base {

protected:

  connectivity_m_base(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
    &&DestinationGridInfo);

  connectivity_m_base(const connectivity_m_base &Other) = delete;
  connectivity_m_base(connectivity_m_base &&Other) noexcept = default;

  connectivity_m_base &operator=(const connectivity_m_base &Other) = delete;
  connectivity_m_base &operator=(connectivity_m_base &&Other) noexcept = default;

  ~connectivity_m_base() noexcept;

  std::shared_ptr<context> Context_;

  const grid *Grid_;

  grid_info DestinationGridInfo_;

  comm_view Comm_;

};

}

class connectivity_m : private connectivity_m_internal::connectivity_m_base {

public:

  connectivity_m(const connectivity_m &Other) = delete;
  connectivity_m(connectivity_m &&Other) noexcept = default;

  connectivity_m &operator=(const connectivity_m &Other) = delete;
  connectivity_m &operator=(connectivity_m &&Other) noexcept = default;

  ~connectivity_m() noexcept;

  floating_ref<const connectivity_m> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<connectivity_m> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const grid &Grid() const { return *Grid_; }

  const grid_info &DestinationGridInfo() const { return DestinationGridInfo_; }

  int Dimension() const { return NumDims_; }

  comm_view Comm() const { return Comm_; }

  long long Count() const { return Count_; }
  int MaxSize() const { return MaxSize_; }

  void Resize(long long Count, int MaxSize);
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddResizeEventListener(F Listener) const {
    return ResizeEvent_.AddListener(std::move(Listener));
  }

  const array<int,3> &Extents() const { return Extents_; }
  bool EditingExtents() const;
  edit_handle<array<int,3>> EditExtents();
  void RestoreExtents();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddExtentsEventListener(F Listener) const {
    return ExtentsEvent_.AddListener(std::move(Listener));
  }

  const array<double,2> &Coords() const { return Coords_; }
  bool EditingCoords() const;
  edit_handle<array<double,2>> EditCoords();
  void RestoreCoords();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddCoordsEventListener(F Listener) const {
    return CoordsEvent_.AddListener(std::move(Listener));
  }

  const array<double,3> &InterpCoefs() const { return InterpCoefs_; }
  bool EditingInterpCoefs() const;
  edit_handle<array<double,3>> EditInterpCoefs();
  void RestoreInterpCoefs();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddInterpCoefsEventListener(F Listener) const {
    return InterpCoefsEvent_.AddListener(std::move(Listener));
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

  static connectivity_m internal_Create(std::shared_ptr<context> &&Context, const grid &Grid,
    grid_info &&DestinationGridInfo);

private:

  floating_ref_generator FloatingRefGenerator_;

  int NumDims_;

  long long Count_;
  int MaxSize_;
  mutable event<void()> ResizeEvent_;

  array<int,3> Extents_;
  editor ExtentsEditor_;
  mutable event<void()> ExtentsEvent_;

  array<double,2> Coords_;
  editor CoordsEditor_;
  mutable event<void()> CoordsEvent_;

  array<double,3> InterpCoefs_;
  editor InterpCoefsEditor_;
  mutable event<void()> InterpCoefsEvent_;

  array<int,2> Destinations_;
  editor DestinationsEditor_;
  mutable event<void()> DestinationsEvent_;

  array<int> DestinationRanks_;
  editor DestinationRanksEditor_;
  mutable event<void()> DestinationRanksEvent_;

  connectivity_m(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
    &&DestinationGridInfo);

};

namespace core {
connectivity_m CreateConnectivityM(std::shared_ptr<context> Context, const grid &Grid, grid_info
  DestinationGridInfo);
}

}

#endif
