// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONNECTIVITY_COMPONENT_HPP_INCLUDED
#define OVK_CORE_CONNECTIVITY_COMPONENT_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/ConnectivityComponent.h>
#include <ovk/core/ConnectivityM.hpp>
#include <ovk/core/ConnectivityN.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DomainBase.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/IDMap.hpp>
#include <ovk/core/IDSet.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/StringWrapper.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <functional>
#include <memory>
#include <string>

namespace ovk {

enum class connectivity_event_flags : int {
  NONE = OVK_CONNECTIVITY_EVENT_FLAGS_NONE,
  CREATE = OVK_CONNECTIVITY_EVENT_FLAGS_CREATE,
  DESTROY = OVK_CONNECTIVITY_EVENT_FLAGS_DESTROY,
  RESIZE_M = OVK_CONNECTIVITY_EVENT_FLAGS_RESIZE_M,
  EDIT_M_EXTENTS = OVK_CONNECTIVITY_EVENT_FLAGS_EDIT_M_EXTENTS,
  EDIT_M_COORDS = OVK_CONNECTIVITY_EVENT_FLAGS_EDIT_M_COORDS,
  EDIT_M_INTERP_COEFS = OVK_CONNECTIVITY_EVENT_FLAGS_EDIT_M_INTERP_COEFS,
  EDIT_M_DESTINATIONS = OVK_CONNECTIVITY_EVENT_FLAGS_EDIT_M_DESTINATIONS,
  RESIZE_N = OVK_CONNECTIVITY_EVENT_FLAGS_RESIZE_N,
  EDIT_N_POINTS = OVK_CONNECTIVITY_EVENT_FLAGS_EDIT_N_POINTS,
  EDIT_N_SOURCES = OVK_CONNECTIVITY_EVENT_FLAGS_EDIT_N_SOURCES,
  ALL_EDITS = OVK_CONNECTIVITY_EVENT_FLAGS_ALL_EDITS,
  ALL = OVK_CONNECTIVITY_EVENT_FLAGS_ALL
};

constexpr inline connectivity_event_flags operator|(connectivity_event_flags Left,
  connectivity_event_flags Right) {
  return connectivity_event_flags(int(Left) | int(Right));
}
constexpr inline connectivity_event_flags operator&(connectivity_event_flags Left,
  connectivity_event_flags Right) {
  return connectivity_event_flags(int(Left) & int(Right));
}
constexpr inline connectivity_event_flags operator^(connectivity_event_flags Left,
  connectivity_event_flags Right) {
  return connectivity_event_flags(int(Left) ^ int(Right));
}
constexpr inline connectivity_event_flags operator~(connectivity_event_flags EventFlags) {
  return connectivity_event_flags(~int(EventFlags));
}
inline connectivity_event_flags operator|=(connectivity_event_flags &Left,
  connectivity_event_flags Right) {
  return Left = Left | Right;
}
inline connectivity_event_flags operator&=(connectivity_event_flags &Left,
  connectivity_event_flags Right) {
  return Left = Left & Right;
}
inline connectivity_event_flags operator^=(connectivity_event_flags &Left,
  connectivity_event_flags Right) {
  return Left = Left ^ Right;
}

namespace connectivity_component_internal {

// For doing stuff before creation and after destruction
class connectivity_component_base {

public:

  connectivity_component_base(const core::domain_base &Domain, std::string &&Name);

  connectivity_component_base(const connectivity_component_base &Other) = delete;
  connectivity_component_base(connectivity_component_base &&Other) noexcept = default;

  connectivity_component_base &operator=(const connectivity_component_base &Other) = delete;
  connectivity_component_base &operator=(connectivity_component_base &&Other) noexcept = default;

  ~connectivity_component_base() noexcept;

  floating_ref<const context> Context_;
  floating_ref<const core::domain_base> Domain_;

  core::string_wrapper Name_;

};

}

class connectivity_component : private connectivity_component_internal::connectivity_component_base
  {

public:

  class params {
  public:
    params() {}
    const std::string &Name() const { return Name_; }
    params &SetName(std::string Name);
  private:
    core::string_wrapper Name_ = "ConnectivityComponent";
    friend class connectivity_component;
  };

  connectivity_component(const core::domain_base &Domain, params Params = {});

  connectivity_component(const connectivity_component &Other) = delete;
  connectivity_component(connectivity_component &&Other) noexcept = default;

  connectivity_component &operator=(const connectivity_component &Other) = delete;
  connectivity_component &operator=(connectivity_component &&Other) noexcept = default;

  floating_ref<const connectivity_component> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<connectivity_component> GetFloatingRef() {
    return FloatingRefGenerator_.Generate(*this);
  }

  const std::string &Name() const { return *Name_; }

  int ConnectivityCount() const;

  const id_set<2> &ConnectivityIDs() const;

  bool ConnectivityExists(int MGridID, int NGridID) const;

  void CreateConnectivity(int MGridID, int NGridID);
  void CreateConnectivities(array_view<const int> MGridIDs, array_view<const int> NGridIDs);

  void DestroyConnectivity(int MGridID, int NGridID);
  void DestroyConnectivities(array_view<const int> MGridIDs, array_view<const int> NGridIDs);

  void ClearConnectivities();

  int LocalConnectivityMCount() const;
  int LocalConnectivityMCountForGrid(int MGridID) const;

  const id_set<2> &LocalConnectivityMIDs() const;

  const connectivity_m &ConnectivityM(int MGridID, int NGridID) const;
  bool EditingConnectivityM(int MGridID, int NGridID) const;
  edit_handle<connectivity_m> EditConnectivityM(int MGridID, int NGridID);
  void RestoreConnectivityM(int MGridID, int NGridID);

  int LocalConnectivityNCount() const;
  int LocalConnectivityNCountForGrid(int NGridID) const;

  const id_set<2> &LocalConnectivityNIDs() const;

  const connectivity_n &ConnectivityN(int MGridID, int NGridID) const;
  bool EditingConnectivityN(int MGridID, int NGridID) const;
  edit_handle<connectivity_n> EditConnectivityN(int MGridID, int NGridID);
  void RestoreConnectivityN(int MGridID, int NGridID);

  void StartEdit();
  void EndEdit();

  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F, int, int,
    connectivity_event_flags, bool>())> event_listener_handle AddConnectivityEventListener(F
    Listener) const {
    return ConnectivityEvent_.AddListener(std::move(Listener));
  }

private:

  struct connectivity_record {};

  struct local_m {
    connectivity_m Connectivity;
    connectivity_event_flags EventFlags;
    floating_ref_generator FloatingRefGenerator;
    event_listener_handle ResizeEventListener;
    event_listener_handle ExtentsEventListener;
    event_listener_handle CoordsEventListener;
    event_listener_handle InterpCoefsEventListener;
    event_listener_handle DestinationsEventListener;
    event_listener_handle DestinationRanksEventListener;
    editor Editor;
    explicit local_m(connectivity_m Connectivity);
  };

  struct local_n {
    connectivity_n Connectivity;
    connectivity_event_flags EventFlags;
    floating_ref_generator FloatingRefGenerator;
    event_listener_handle ResizeEventListener;
    event_listener_handle PointsEventListener;
    event_listener_handle SourcesEventListener;
    event_listener_handle SourceRanksEventListener;
    editor Editor;
    explicit local_n(connectivity_n Connectivity);
  };

  floating_ref_generator FloatingRefGenerator_;

  id_map<1,grid_event_flags> GridEventFlags_;
  event_listener_handle GridEventListener_;

  id_map<2,connectivity_record> ConnectivityRecords_;
  id_map<2,local_m,false> LocalMs_;
  id_map<2,local_n,false> LocalNs_;

  mutable event<void(int, int, connectivity_event_flags, bool)> ConnectivityEvent_;

  void OnGridEvent_();
  void DestroyConnectivitiesForDyingGrids_();

  void SyncEdits_();

};

}

#endif
