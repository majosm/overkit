// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_LOGGER_INCLUDED
#define OVK_CORE_LOGGER_INCLUDED

#include "ovk/core/Global.h"

#ifdef __cplusplus
extern "C" {
#endif

enum { LOG_BUFFER_LENGTH = 1024 };

typedef struct {
  ovk_log_level level;
  int write_rank;
} t_logger;

void PRIVATE(CreateLogger)(t_logger **Logger, ovk_log_level LogLevel, int WriteRank);
#define CreateLogger(...) PRIVATE(CreateLogger)(__VA_ARGS__)

void PRIVATE(DestroyLogger)(t_logger **Logger);
#define DestroyLogger(...) PRIVATE(DestroyLogger)(__VA_ARGS__)

static inline bool LoggingStatus(const t_logger *Logger) {
  return (Logger->level & OVK_LOG_STATUS) != 0;
}
static inline bool LoggingWarnings(const t_logger *Logger) {
  return (Logger->level & OVK_LOG_WARNINGS) != 0;
}
static inline bool LoggingErrors(const t_logger *Logger) {
  return (Logger->level & OVK_LOG_ERRORS) != 0;
}

void PRIVATE(GetLogLevel)(const t_logger *Logger, ovk_log_level *LogLevel);
#define GetLogLevel(...) PRIVATE(GetLogLevel)(__VA_ARGS__)
void PRIVATE(SetLogLevel)(t_logger *Logger, ovk_log_level LogLevel);
#define SetLogLevel(...) PRIVATE(SetLogLevel)(__VA_ARGS__)

void PRIVATE(LogStatus)(t_logger *Logger, bool WriteCondition, int IncrementLevel,
  const char *Format, ...);
#define LogStatus(...) PRIVATE(LogStatus)(__VA_ARGS__)
void PRIVATE(LogWarning)(t_logger *Logger, bool WriteCondition, const char *Format, ...);
#define LogWarning(...) PRIVATE(LogWarning)(__VA_ARGS__)
void PRIVATE(LogError)(t_logger *Logger, bool WriteCondition, const char *Format, ...);
#define LogError(...) PRIVATE(LogError)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
