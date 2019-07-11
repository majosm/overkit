// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_STATE_HPP_INCLUDED
#define OVK_CORE_STATE_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Editor.hpp>
#include <ovk/core/Event.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/State.h>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

enum class state_flags : int {
  NONE = OVK_STATE_FLAGS_NONE,
  ACTIVE = OVK_STATE_FLAGS_ACTIVE,
  DOMAIN_BOUNDARY = OVK_STATE_FLAGS_DOMAIN_BOUNDARY,
  INTERNAL_BOUNDARY = OVK_STATE_FLAGS_INTERNAL_BOUNDARY,
  OVERLAPPED = OVK_STATE_FLAGS_OVERLAPPED,
  INFERRED_DOMAIN_BOUNDARY = OVK_STATE_FLAGS_INFERRED_DOMAIN_BOUNDARY,
  BOUNDARY_HOLE = OVK_STATE_FLAGS_BOUNDARY_HOLE,
  OCCLUDED = OVK_STATE_FLAGS_OCCLUDED,
  FRINGE = OVK_STATE_FLAGS_FRINGE,
  OUTER_FRINGE = OVK_STATE_FLAGS_OUTER_FRINGE,
  INNER_FRINGE = OVK_STATE_FLAGS_INNER_FRINGE,
  OVERLAP_MINIMIZED = OVK_STATE_FLAGS_OVERLAP_MINIMIZED,
  RECEIVER = OVK_STATE_FLAGS_RECEIVER,
  ORPHAN = OVK_STATE_FLAGS_ORPHAN,
  DEBUG1 = OVK_STATE_FLAGS_DEBUG1,
  DEBUG2 = OVK_STATE_FLAGS_DEBUG2,
  DEBUG3 = OVK_STATE_FLAGS_DEBUG3,
  DEBUG4 = OVK_STATE_FLAGS_DEBUG4,
  DEBUG5 = OVK_STATE_FLAGS_DEBUG5,
  ALL = OVK_STATE_FLAGS_ALL
};

inline bool ValidStateFlags(state_flags StateFlags) {
  return ovkValidStateFlags(ovk_state_flags(StateFlags));
}

constexpr inline state_flags operator|(state_flags Left, state_flags Right) {
  return state_flags(int(Left) | int(Right));
}
constexpr inline state_flags operator&(state_flags Left, state_flags Right) {
  return state_flags(int(Left) & int(Right));
}
constexpr inline state_flags operator^(state_flags Left, state_flags Right) {
  return state_flags(int(Left) ^ int(Right));
}
constexpr inline state_flags operator~(state_flags StateFlags) {
  return state_flags(~int(StateFlags));
}
inline state_flags operator|=(state_flags &Left, state_flags Right) {
  return Left = Left | Right;
}
inline state_flags operator&=(state_flags &Left, state_flags Right) {
  return Left = Left & Right;
}
inline state_flags operator^=(state_flags &Left, state_flags Right) {
  return Left = Left ^ Right;
}

namespace state_internal {

// For doing stuff before creation and after destruction
class state_base {

protected:

  state_base(std::shared_ptr<context> &&Context, const grid &Grid);

  state_base(const state_base &Other) = delete;
  state_base(state_base &&Other) noexcept = default;

  state_base &operator=(const state_base &Other) = delete;
  state_base &operator=(state_base &&Other) noexcept = default;

  ~state_base() noexcept;

  std::shared_ptr<context> Context_;

  const grid *Grid_;

  comm_view Comm_;

};

}

class state : private state_internal::state_base {

public:

  class params {
  public:
    params() = default;
  private:
    friend class state;
  };

  state(const state &Other) = delete;
  state(state &&Other) noexcept = default;

  state &operator=(const state &Other) = delete;
  state &operator=(state &&Other) = default;

  ~state() noexcept;

  floating_ref<const state> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<state> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const grid &Grid() const { return *Grid_; }

  comm_view Comm() const { return Comm_; }

  const field<state_flags> &Flags() const { return Flags_; }
  bool EditingFlags() const;
  edit_handle<field<state_flags>> EditFlags();
  void RestoreFlags();
  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableWith<F>())> event_listener_handle
    AddFlagsEventListener(F Listener) const {
    return FlagsEvent_.AddListener(std::move(Listener));
  }

  static state internal_Create(std::shared_ptr<context> &&Context, const grid &Grid, params
    &&Params);

private:

  floating_ref_generator FloatingRefGenerator_;

  field<state_flags> Flags_;
  editor FlagsEditor_;
  mutable event<void()> FlagsEvent_;

  state(std::shared_ptr<context> &&Context, const grid &Grid, params &&Params);

};

namespace core {
state CreateState(std::shared_ptr<context> Context, const grid &Grid, state::params Params={});
}

}

#endif
