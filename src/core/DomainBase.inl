// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifdef __cplusplus
extern "C" {
#endif

static inline bool ovkValidGridEventFlags(ovk_grid_event_flags EventFlags) {

  return EventFlags >= OVK_GRID_EVENT_FLAGS_NONE && EventFlags <= OVK_GRID_EVENT_FLAGS_ALL;

}

static inline bool ovkValidComponentEventFlags(ovk_component_event_flags EventFlags) {

  return EventFlags >= OVK_COMPONENT_EVENT_FLAGS_NONE && EventFlags <=
    OVK_COMPONENT_EVENT_FLAGS_ALL;

}

#ifdef __cplusplus
}
#endif
