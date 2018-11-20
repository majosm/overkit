// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Logger.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <string>

namespace ovk {
namespace core {

void CreateLogger(logger &Logger, log_level LogLevel, int WriteRank) {

  OVK_DEBUG_ASSERT(ValidLogLevel(LogLevel), "Invalid log level.");
  OVK_DEBUG_ASSERT(WriteRank >= 0, "Invalid write rank.");

  Logger.Level_ = LogLevel;
  Logger.WriteRank_ = WriteRank;

}

void DestroyLogger(logger &Logger) {

  Logger.Level_ = log_level::NONE;
  Logger.WriteRank_ = -1;

}

void GetLogLevel(const logger &Logger, log_level &LogLevel) {

  LogLevel = Logger.Level_;

}

void SetLogLevel(logger &Logger, log_level LogLevel) {

  OVK_DEBUG_ASSERT(ValidLogLevel(LogLevel), "Invalid log level.");

  Logger.Level_ = LogLevel;

}

namespace logger_internal {

void ReplaceRank(std::string &Message, int Rank) {

  size_t RankPos = Message.find("@rank@");

  if (RankPos == std::string::npos) return;

  std::string RankString = FormatNumber(Rank);

  std::string NewMessage;

  size_t Pos = 0;
  while (RankPos != std::string::npos) {
    NewMessage.append(Message, Pos, RankPos-Pos);
    NewMessage += RankString;
    Pos = RankPos + 6;
    RankPos = Message.find("@rank@", Pos);
  }
  NewMessage.append(Message, Pos, std::string::npos);

  Message = NewMessage;

}

}

}}
