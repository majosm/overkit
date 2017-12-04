// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_LOGGER_INCLUDED
#define OVK_CORE_LOGGER_INCLUDED

#include "Debug.h"
#include "Global.h"

#include <stdarg.h>
#include <stdio.h>

enum { LOG_BUFFER_SIZE = 1024 };

typedef struct {
  ovk_log_level level;
  char buffer[LOG_BUFFER_SIZE];
} t_logger;

static inline bool ValidLogLevel(ovk_log_level LogLevel) {

  return LogLevel >= OVK_LOG_NONE && LogLevel <= OVK_LOG_ALL;

}

static inline void CreateLogger(t_logger **Logger_, ovk_log_level LogLevel) {

  *Logger_ = malloc(sizeof(t_logger));
  t_logger *Logger = *Logger_;

  Logger->level = LogLevel;

}

static inline void DestroyLogger(t_logger **Logger_) {

  free(*Logger_);
  *Logger_ = NULL;

}

static inline int LogBufferWriteVAList(t_logger *Logger, int Offset, const char *Format,
  va_list ArgList) {

  OVK_DEBUG_ASSERT(Offset >= 0 && Offset < LOG_BUFFER_SIZE, "Invalid log buffer offset.");

  int FullWriteSize = vsnprintf(Logger->buffer + Offset, LOG_BUFFER_SIZE-Offset, Format,
    ArgList);

  // Must be strictly less to leave room for null terminator
  OVK_DEBUG_ASSERT(FullWriteSize < LOG_BUFFER_SIZE-Offset, "Insufficient buffer length for "
    "log message.");

  return min(FullWriteSize, LOG_BUFFER_SIZE-Offset);

}

static inline int LogBufferWrite(t_logger *Logger, int Offset, const char *Format, ...) {

  va_list ArgList;
  va_start(ArgList, Format);
  int WriteSize = LogBufferWriteVAList(Logger, Offset, Format, ArgList);
  va_end(ArgList);

  return WriteSize;

}

static inline void LogBufferFlush(t_logger *Logger, FILE *Stream) {

  fprintf(Stream, "%s", Logger->buffer);
  fflush(Stream);

}

static inline void LogStatus(t_logger *Logger, bool WriteCondition, int IncrementLevel,
  const char *Format, ...) {

  if ((Logger->level & OVK_LOG_STATUS) && WriteCondition) {

    int Offset = 0;

    Offset += LogBufferWrite(Logger, Offset, "ovk :: ");

    int PrefixOffset = Offset;
    while (Offset < PrefixOffset + 2*(IncrementLevel-1)) {
      Offset += LogBufferWrite(Logger, Offset, "  ");
    }

    if (IncrementLevel > 0) {
      if ((IncrementLevel-1) % 2 == 0) {
        Offset += LogBufferWrite(Logger, Offset, "* ");
      } else {
        Offset += LogBufferWrite(Logger, Offset, "- ");
      }
    }

    va_list ArgList;
    va_start(ArgList, Format);
    Offset += LogBufferWriteVAList(Logger, Offset, Format, ArgList);
    va_end(ArgList);

    Offset += LogBufferWrite(Logger, Offset, "\n");

    LogBufferFlush(Logger, stdout);

  }

}

static inline void LogWarning(t_logger *Logger, bool WriteCondition, const char *Format, ...) {

  if ((Logger->level & OVK_LOG_WARNINGS) && WriteCondition) {

    int Offset = 0;

    Offset += LogBufferWrite(Logger, Offset, "ovk :: ");

    Offset += LogBufferWrite(Logger, Offset, "WARNING: ");

    va_list ArgList;
    va_start(ArgList, Format);
    Offset += LogBufferWriteVAList(Logger, Offset, Format, ArgList);
    va_end(ArgList);

    Offset += LogBufferWrite(Logger, Offset, "\n");

    LogBufferFlush(Logger, stderr);

  }

}

static inline void LogError(t_logger *Logger, bool WriteCondition, const char *Format, ...) {

  if ((Logger->level & OVK_LOG_ERRORS) && WriteCondition) {

    int Offset = 0;

    Offset += LogBufferWrite(Logger, Offset, "ovk :: ");

    Offset += LogBufferWrite(Logger, Offset, "ERROR: ");

    va_list ArgList;
    va_start(ArgList, Format);
    Offset += LogBufferWriteVAList(Logger, Offset, Format, ArgList);
    va_end(ArgList);

    Offset += LogBufferWrite(Logger, Offset, "\n");

    LogBufferFlush(Logger, stderr);

  }

}

#endif
