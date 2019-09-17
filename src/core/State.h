// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_STATE_H_INCLUDED
#define OVK_CORE_STATE_H_INCLUDED

#include <ovk/core/Global.h>

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_STATE_FLAGS_NONE = 0,
  OVK_STATE_FLAGS_ACTIVE =                   1 <<  0,
  OVK_STATE_FLAGS_DOMAIN_BOUNDARY =          1 <<  1,
  OVK_STATE_FLAGS_INTERNAL_BOUNDARY =        1 <<  2,
  OVK_STATE_FLAGS_OVERLAPPED =               1 <<  3,
  OVK_STATE_FLAGS_INFERRED_DOMAIN_BOUNDARY = 1 <<  4,
  OVK_STATE_FLAGS_BOUNDARY_HOLE =            1 <<  5,
  OVK_STATE_FLAGS_OCCLUDED =                 1 <<  6,
  OVK_STATE_FLAGS_FRINGE =                   1 <<  7,
  OVK_STATE_FLAGS_OUTER_FRINGE =             1 <<  8,
  OVK_STATE_FLAGS_INNER_FRINGE =             1 <<  9,
  OVK_STATE_FLAGS_OVERLAP_MINIMIZED =        1 << 10,
  OVK_STATE_FLAGS_RECEIVER =                 1 << 11,
  OVK_STATE_FLAGS_ORPHAN =                   1 << 12,
  OVK_STATE_FLAGS_DEBUG1 =                   1 << 13,
  OVK_STATE_FLAGS_DEBUG2 =                   1 << 14,
  OVK_STATE_FLAGS_DEBUG3 =                   1 << 15,
  OVK_STATE_FLAGS_DEBUG4 =                   1 << 16,
  OVK_STATE_FLAGS_DEBUG5 =                   1 << 17,
  OVK_STATE_FLAGS_ALL = 
    OVK_STATE_FLAGS_ACTIVE |
    OVK_STATE_FLAGS_DOMAIN_BOUNDARY |
    OVK_STATE_FLAGS_INTERNAL_BOUNDARY |
    OVK_STATE_FLAGS_OVERLAPPED |
    OVK_STATE_FLAGS_INFERRED_DOMAIN_BOUNDARY |
    OVK_STATE_FLAGS_BOUNDARY_HOLE |
    OVK_STATE_FLAGS_OCCLUDED |
    OVK_STATE_FLAGS_FRINGE |
    OVK_STATE_FLAGS_OUTER_FRINGE |
    OVK_STATE_FLAGS_INNER_FRINGE |
    OVK_STATE_FLAGS_OVERLAP_MINIMIZED |
    OVK_STATE_FLAGS_RECEIVER |
    OVK_STATE_FLAGS_ORPHAN |
    OVK_STATE_FLAGS_DEBUG1 |
    OVK_STATE_FLAGS_DEBUG2 |
    OVK_STATE_FLAGS_DEBUG3 |
    OVK_STATE_FLAGS_DEBUG4 |
    OVK_STATE_FLAGS_DEBUG5
} ovk_state_flags;

static inline bool ovkValidStateFlags(ovk_state_flags StateFlags) {

  return StateFlags >= OVK_STATE_FLAGS_NONE && StateFlags <= OVK_STATE_FLAGS_ALL;

}

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
constexpr inline ovk_state_flags operator|(ovk_state_flags Left, ovk_state_flags Right) {
  return ovk_state_flags(int(Left) | int(Right));
}
constexpr inline ovk_state_flags operator&(ovk_state_flags Left, ovk_state_flags Right) {
  return ovk_state_flags(int(Left) & int(Right));
}
constexpr inline ovk_state_flags operator^(ovk_state_flags Left, ovk_state_flags Right) {
  return ovk_state_flags(int(Left) ^ int(Right));
}
constexpr inline ovk_state_flags operator~(ovk_state_flags StateFlags) {
  return ovk_state_flags(~int(StateFlags));
}
inline ovk_state_flags operator|=(ovk_state_flags &Left, ovk_state_flags Right) {
  return Left = Left | Right;
}
inline ovk_state_flags operator&=(ovk_state_flags &Left, ovk_state_flags Right) {
  return Left = Left & Right;
}
inline ovk_state_flags operator^=(ovk_state_flags &Left, ovk_state_flags Right) {
  return Left = Left ^ Right;
}
#endif

#endif
