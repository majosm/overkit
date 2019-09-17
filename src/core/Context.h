// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONTEXT_H_INCLUDED
#define OVK_CORE_CONTEXT_H_INCLUDED

#include <ovk/core/Global.h>

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_LOG_NONE = 0,
  OVK_LOG_ERRORS = 1 << 0,
  OVK_LOG_WARNINGS = 1 << 1,
  OVK_LOG_STATUS = 1 << 2,
  OVK_LOG_DEBUG = 1 << 3,
  OVK_LOG_ALL =
    OVK_LOG_ERRORS |
    OVK_LOG_WARNINGS |
    OVK_LOG_STATUS |
    OVK_LOG_DEBUG
} ovk_log_level;

static inline bool ovkValidLogLevel(ovk_log_level LogLevel) {

  return LogLevel >= OVK_LOG_NONE && LogLevel <= OVK_LOG_ALL;

}

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
constexpr inline ovk_log_level operator|(ovk_log_level Left, ovk_log_level Right) {
  return ovk_log_level(int(Left) | int(Right));
}
constexpr inline ovk_log_level operator&(ovk_log_level Left, ovk_log_level Right) {
  return ovk_log_level(int(Left) & int(Right));
}
constexpr inline ovk_log_level operator^(ovk_log_level Left, ovk_log_level Right) {
  return ovk_log_level(int(Left) ^ int(Right));
}
constexpr inline ovk_log_level operator~(ovk_log_level LogLevel) {
  return ovk_log_level(~int(LogLevel));
}
inline ovk_log_level operator|=(ovk_log_level &Left, ovk_log_level Right) {
  return Left = Left | Right;
}
inline ovk_log_level operator&=(ovk_log_level &Left, ovk_log_level Right) {
  return Left = Left & Right;
}
inline ovk_log_level operator^=(ovk_log_level &Left, ovk_log_level Right) {
  return Left = Left ^ Right;
}
#endif

#endif
