// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifdef __cplusplus
extern "C" {
#endif

static inline bool ovkValidConnectivityEventFlags(ovk_connectivity_event_flags EventFlags) {

  return EventFlags >= OVK_CONNECTIVITY_EVENT_FLAGS_NONE && EventFlags <=
    OVK_CONNECTIVITY_EVENT_FLAGS_ALL;

}

#ifdef __cplusplus
}
#endif
