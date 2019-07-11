// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_STATE_COMPONENT_HPP_INCLUDED
#define OVK_CORE_STATE_COMPONENT_HPP_INCLUDED

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DomainBase.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/IDMap.hpp>
#include <ovk/core/IDSet.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/State.hpp>
#include <ovk/core/StateComponent.h>
#include <ovk/core/StringWrapper.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <functional>
#include <memory>
#include <string>

namespace ovk {

enum class state_event_flags : int {
  NONE = OVK_STATE_EVENT_FLAGS_NONE,
  CREATE = OVK_STATE_EVENT_FLAGS_CREATE,
  DESTROY = OVK_STATE_EVENT_FLAGS_DESTROY,
  EDIT_FLAGS = OVK_STATE_EVENT_FLAGS_EDIT_FLAGS,
  ALL = OVK_STATE_EVENT_FLAGS_ALL
};

constexpr inline state_event_flags operator|(state_event_flags Left, state_event_flags Right) {
  return state_event_flags(int(Left) | int(Right));
}
constexpr inline state_event_flags operator&(state_event_flags Left, state_event_flags Right) {
  return state_event_flags(int(Left) & int(Right));
}
constexpr inline state_event_flags operator^(state_event_flags Left, state_event_flags Right) {
  return state_event_flags(int(Left) ^ int(Right));
}
constexpr inline state_event_flags operator~(state_event_flags EventFlags) {
  return state_event_flags(~int(EventFlags));
}
inline state_event_flags operator|=(state_event_flags &Left, state_event_flags Right) {
  return Left = Left | Right;
}
inline state_event_flags operator&=(state_event_flags &Left, state_event_flags Right) {
  return Left = Left & Right;
}
inline state_event_flags operator^=(state_event_flags &Left, state_event_flags Right) {
  return Left = Left ^ Right;
}

namespace state_component_internal {

// For doing stuff before creation and after destruction
class state_component_base {

public:

  state_component_base(const core::domain_base &Domain, std::string &&Name);

  state_component_base(const state_component_base &Other) = delete;
  state_component_base(state_component_base &&Other) noexcept = default;

  state_component_base &operator=(const state_component_base &Other) = delete;
  state_component_base &operator=(state_component_base &&Other) noexcept = default;

  ~state_component_base() noexcept;

  floating_ref<const context> Context_;
  floating_ref<const core::domain_base> Domain_;

  core::string_wrapper Name_;

};

}

class state_component : private state_component_internal::state_component_base {

public:

  class params {
  public:
    params() {}
    const std::string &Name() const { return Name_; }
    params &SetName(std::string Name);
  private:
    core::string_wrapper Name_ = "StateComponent";
    friend class state_component;
  };

  state_component(const core::domain_base &Domain, params Params={});

  state_component(const state_component &Other) = delete;
  state_component(state_component &&Other) noexcept = default;

  state_component &operator=(const state_component &Other) = delete;
  state_component &operator=(state_component &&Other) noexcept = default;

  floating_ref<const state_component> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<state_component> GetFloatingRef() {
    return FloatingRefGenerator_.Generate(*this);
  }

  const std::string &Name() const { return *Name_; }

  int StateCount() const;

  const id_set<1> &StateIDs() const;

  bool StateExists(int GridID) const;

  void CreateState(int GridID, optional<state::params> MaybeParams={});
  void CreateStates(array_view<const int> GridIDs, array<optional<state::params>> MaybeParams={});

  void DestroyState(int GridID);
  void DestroyStates(array_view<const int> GridIDs);

  void ClearStates();

  int LocalStateCount() const;

  const id_set<1> &LocalStateIDs() const;

  const state &State(int GridID) const;
  bool EditingState(int GridID) const;
  edit_handle<state> EditState(int GridID);
  void RestoreState(int GridID);

  void StartEdit();
  void EndEdit();

  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F, int, state_event_flags,
    bool>())> event_listener_handle AddStateEventListener(F Listener) const {
    return StateEvent_.AddListener(std::move(Listener));
  }

private:

  struct state_record {};

  struct local {
    state State;
    state_event_flags EventFlags;
    floating_ref_generator FloatingRefGenerator;
    event_listener_handle CoordsEventListener;
    editor Editor;
    explicit local(state State);
  };

  floating_ref_generator FloatingRefGenerator_;

  id_map<1,grid_event_flags> GridEventFlags_;
  event_listener_handle GridEventListener_;

  id_map<1,state_record> StateRecords_;
  id_map<1,local,false> Locals_;

  mutable event<void(int, state_event_flags, bool)> StateEvent_;

  void OnGridEvent_();
  void DestroyStatesForDyingGrids_();

  void SyncEdits_();

};

}

#endif
