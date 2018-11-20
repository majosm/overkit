// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_LOGGER_HPP_INCLUDED
#define OVK_CORE_LOGGER_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/TextProcessing.hpp>

#include <string>

namespace ovk {
namespace core {

struct logger {
  log_level Level_;
  int WriteRank_;
};

void CreateLogger(logger &Logger, log_level LogLevel, int WriteRank);
void DestroyLogger(logger &Logger);

inline bool LoggingStatus(const logger &Logger) {
  return (Logger.Level_ & log_level::STATUS) != log_level::NONE;
}
inline bool LoggingWarnings(const logger &Logger) {
  return (Logger.Level_ & log_level::WARNINGS) != log_level::NONE;
}
inline bool LoggingErrors(const logger &Logger) {
  return (Logger.Level_ & log_level::ERRORS) != log_level::NONE;
}

void GetLogLevel(const logger &Logger, log_level &LogLevel);
void SetLogLevel(logger &Logger, log_level LogLevel);

template <typename... Ts> void LogStatus(logger &Logger, bool WriteCondition, int IncrementLevel,
  const std::string &Format, const Ts &... Args);
template <typename... Ts> void LogWarning(logger &Logger, bool WriteCondition, const std::string
  &Format, const Ts &... Args);
template <typename... Ts> void LogError(logger &Logger, bool WriteCondition, const std::string
  &Format, const Ts &... Args);

}}

#include <ovk/core/Logger.inl>

#endif
