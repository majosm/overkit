// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DOMAIN_BASE_H_INCLUDED
#define OVK_CORE_DOMAIN_BASE_H_INCLUDED

#include <ovk/core/GlobalBase.h>

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_GRID_EVENT_FLAGS_NONE = 0,
  OVK_GRID_EVENT_FLAGS_CREATE = 1 << 0,
  OVK_GRID_EVENT_FLAGS_DESTROY = 1 << 1,
  OVK_GRID_EVENT_FLAGS_ALL =
    OVK_GRID_EVENT_FLAGS_CREATE |
    OVK_GRID_EVENT_FLAGS_DESTROY
} ovk_grid_event_flags;

static inline bool ovkValidGridEventFlags(ovk_grid_event_flags EventFlags) {

  return EventFlags >= OVK_GRID_EVENT_FLAGS_NONE && EventFlags <= OVK_GRID_EVENT_FLAGS_ALL;

}

typedef enum {
  OVK_COMPONENT_EVENT_FLAGS_NONE = 0,
  OVK_COMPONENT_EVENT_FLAGS_CREATE = 1 << 0,
  OVK_COMPONENT_EVENT_FLAGS_DESTROY = 1 << 1,
  OVK_COMPONENT_EVENT_FLAGS_EDIT = 1 << 2,
  OVK_COMPONENT_EVENT_FLAGS_ALL =
    OVK_COMPONENT_EVENT_FLAGS_CREATE |
    OVK_COMPONENT_EVENT_FLAGS_DESTROY |
    OVK_COMPONENT_EVENT_FLAGS_EDIT
} ovk_component_event_flags;

static inline bool ovkValidComponentEventFlags(ovk_component_event_flags EventFlags) {

  return EventFlags >= OVK_COMPONENT_EVENT_FLAGS_NONE && EventFlags <=
    OVK_COMPONENT_EVENT_FLAGS_ALL;

}

#ifdef __cplusplus
}
#endif

#endif
