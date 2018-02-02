// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_LOGGER_INCLUDED
#define OVK_CORE_LOGGER_INCLUDED

#include "Global.h"

#ifdef __cplusplus
extern "C" {
#endif

enum { LOG_BUFFER_SIZE = 1024 };

typedef struct {
  int rank;
  ovk_log_level level;
  char buffer[LOG_BUFFER_SIZE];
} t_logger;

static inline int LogRank(const t_logger *Logger) { return Logger->rank; }
static inline ovk_log_level LogLevel(const t_logger *Logger) { return Logger->level; }

void PRIVATE(CreateLogger)(t_logger **Logger, int Rank, ovk_log_level LogLevel);
#define CreateLogger(...) PRIVATE(CreateLogger)(__VA_ARGS__)

void PRIVATE(DestroyLogger)(t_logger **Logger);
#define DestroyLogger(...) PRIVATE(DestroyLogger)(__VA_ARGS__)

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
