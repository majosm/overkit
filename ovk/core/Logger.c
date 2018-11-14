// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Logger.h"

#include "ovk/core/Global.h"
#include "ovk/core/TextUtils.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

static int LogBufferWrite(char *Buffer, int Offset, const char *Format, ...);
static int LogBufferWriteVAList(char *Buffer, int Offset, const char *Format, va_list ArgList);
static void LogBufferFlush(const char *Buffer, FILE *Stream);
static void ReplaceRank(char *Buffer, int Rank);

void PRIVATE(CreateLogger)(t_logger **Logger_, ovk_log_level LogLevel, int WriteRank) {

  *Logger_ = malloc(sizeof(t_logger));
  t_logger *Logger = *Logger_;

  Logger->level = LogLevel;
  Logger->write_rank = WriteRank;

}

void PRIVATE(DestroyLogger)(t_logger **Logger_) {

  free(*Logger_);
  *Logger_ = NULL;

}

void PRIVATE(GetLogLevel)(const t_logger *Logger, ovk_log_level *LogLevel) {

  OVK_DEBUG_ASSERT(Logger, "Invalid logger pointer.");
  OVK_DEBUG_ASSERT(LogLevel, "Invalid log level pointer.");

  *LogLevel = Logger->level;

}

void PRIVATE(SetLogLevel)(t_logger *Logger, ovk_log_level LogLevel) {

  OVK_DEBUG_ASSERT(Logger, "Invalid logger pointer.");

  Logger->level = LogLevel;

}

void PRIVATE(LogStatus)(t_logger *Logger, bool WriteCondition, int IncrementLevel,
  const char *Format, ...) {

  OVK_DEBUG_ASSERT(Logger, "Invalid logger pointer.");
  OVK_DEBUG_ASSERT(IncrementLevel >= 0, "Invalid increment level.");

  if (LoggingStatus(Logger) && WriteCondition) {

    char Buffer[LOG_BUFFER_LENGTH];
    int Offset = 0;

    Offset += LogBufferWrite(Buffer, Offset, "ovk :: ");

    int PrefixOffset = Offset;
    while (Offset < PrefixOffset + 2*(IncrementLevel-1)) {
      Offset += LogBufferWrite(Buffer, Offset, "  ");
    }

    if (IncrementLevel > 0) {
      if ((IncrementLevel-1) % 2 == 0) {
        Offset += LogBufferWrite(Buffer, Offset, "* ");
      } else {
        Offset += LogBufferWrite(Buffer, Offset, "- ");
      }
    }

    va_list ArgList;
    va_start(ArgList, Format);
    Offset += LogBufferWriteVAList(Buffer, Offset, Format, ArgList);
    va_end(ArgList);

    Offset += LogBufferWrite(Buffer, Offset, "\n");

    // Must be strictly less to leave room for null terminator
    OVK_DEBUG_ASSERT(Offset < LOG_BUFFER_LENGTH, "Insufficient buffer length for log message.");

    ReplaceRank(Buffer, Logger->write_rank);

    LogBufferFlush(Buffer, stdout);

  }

}

void PRIVATE(LogWarning)(t_logger *Logger, bool WriteCondition, const char *Format, ...) {

  OVK_DEBUG_ASSERT(Logger, "Invalid logger pointer.");

  if (LoggingWarnings(Logger) && WriteCondition) {

    char Buffer[LOG_BUFFER_LENGTH];
    int Offset = 0;

    Offset += LogBufferWrite(Buffer, Offset, "ovk :: ");

    Offset += LogBufferWrite(Buffer, Offset, "WARNING: ");

    va_list ArgList;
    va_start(ArgList, Format);
    Offset += LogBufferWriteVAList(Buffer, Offset, Format, ArgList);
    va_end(ArgList);

    Offset += LogBufferWrite(Buffer, Offset, "\n");

    // Must be strictly less to leave room for null terminator
    OVK_DEBUG_ASSERT(Offset < LOG_BUFFER_LENGTH, "Insufficient buffer length for log message.");

    ReplaceRank(Buffer, Logger->write_rank);

    LogBufferFlush(Buffer, stderr);

  }

}

void PRIVATE(LogError)(t_logger *Logger, bool WriteCondition, const char *Format, ...) {

  OVK_DEBUG_ASSERT(Logger, "Invalid logger pointer.");

  if (LoggingErrors(Logger) && WriteCondition) {

    char Buffer[LOG_BUFFER_LENGTH];
    int Offset = 0;

    Offset += LogBufferWrite(Buffer, Offset, "ovk :: ");

    Offset += LogBufferWrite(Buffer, Offset, "ERROR: ");

    va_list ArgList;
    va_start(ArgList, Format);
    Offset += LogBufferWriteVAList(Buffer, Offset, Format, ArgList);
    va_end(ArgList);

    Offset += LogBufferWrite(Buffer, Offset, "\n");

    // Must be strictly less to leave room for null terminator
    OVK_DEBUG_ASSERT(Offset < LOG_BUFFER_LENGTH, "Insufficient buffer length for log message.");

    ReplaceRank(Buffer, Logger->write_rank);

    LogBufferFlush(Buffer, stderr);

  }

}

static int LogBufferWrite(char *Buffer, int Offset, const char *Format, ...) {

  va_list ArgList;
  va_start(ArgList, Format);
  int WriteSize = LogBufferWriteVAList(Buffer, Offset, Format, ArgList);
  va_end(ArgList);

  return WriteSize;

}

static int LogBufferWriteVAList(char *Buffer, int Offset, const char *Format, va_list ArgList) {

  int FullWriteSize = vsnprintf(Buffer + Offset, LOG_BUFFER_LENGTH-Offset, Format,
    ArgList);

  return min(FullWriteSize, LOG_BUFFER_LENGTH-Offset);

}

static void LogBufferFlush(const char *Buffer, FILE *Stream) {

  fprintf(Stream, "%s", Buffer);
  fflush(Stream);

}

static void ReplaceRank(char *Buffer, int Rank) {

  bool HasRank = strstr(Buffer, "@rank@") != NULL;
  if (!HasRank) return;

  int Length = strlen(Buffer);

  char RankString[NUMBER_STRING_LENGTH];
  IntToString(Rank, RankString);
  int RankStringLength = strlen(RankString);

  char NewBuffer[LOG_BUFFER_LENGTH];

  while (true) {
    const char *RankLocConst = strstr(Buffer, "@rank@");
    if (!RankLocConst) break;
    char *RankLoc = Buffer + (RankLocConst-Buffer);
    *RankLoc = '\0';
    const char *BufferBefore = Buffer;
    const char *BufferAfter = RankLoc+6;
    snprintf(NewBuffer, LOG_BUFFER_LENGTH, "%s%s%s", BufferBefore, RankString, BufferAfter);
    Length += RankStringLength-6;
    // Must be strictly less to leave room for null terminator
    OVK_DEBUG_ASSERT(Length < LOG_BUFFER_LENGTH, "Insufficient buffer length for log message.");
    strcpy(Buffer, NewBuffer);
  }

}
