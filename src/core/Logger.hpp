// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_LOGGER_HPP_INCLUDED
#define OVK_CORE_LOGGER_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/TextProcessing.hpp>

#include <string>

namespace ovk {
namespace core {

class logger {

public:

  logger(int Rank);

  logger(const logger &Other) = delete;
  logger(logger &&Other) noexcept = default;

  logger &operator=(const logger &Other) = delete;
  logger &operator=(logger &&Other) noexcept = default;

  bool LoggingErrors() const { return LoggingErrors_; }
  void EnableErrorLogging();
  void DisableErrorLogging();
  template <typename... Ts> void LogError(bool WriteCondition, const std::string &Format, const
    Ts &... Args) const;

  bool LoggingWarnings() const { return LoggingWarnings_; }
  void EnableWarningLogging();
  void DisableWarningLogging();
  template <typename... Ts> void LogWarning(bool WriteCondition, const std::string &Format, const
    Ts &... Args) const;

  bool LoggingStatus() const { return LoggingStatus_; }
  void EnableStatusLogging();
  void DisableStatusLogging();
  template <typename... Ts> void LogStatus(bool WriteCondition, int IncrementLevel, const
    std::string &Format, const Ts &... Args) const;

  bool LoggingDebug() const { return LoggingDebug_; }
  void EnableDebugLogging();
  void DisableDebugLogging();
  template <typename... Ts> void LogDebug(bool WriteCondition, int IncrementLevel, const
    std::string &Format, const Ts &... Args) const;

private:

  bool LoggingErrors_ = false;
  bool LoggingWarnings_ = false;
  bool LoggingStatus_ = false;
  bool LoggingDebug_ = false;
  int Rank_;

};

}}

#include <ovk/core/Logger.inl>

#endif
