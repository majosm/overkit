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

class logger {

public:

  // Get rid of this when context no longer needs it
  logger();

  logger(log_level Level, int Rank);

  logger(const logger &Other) = delete;
  logger(logger &&Other) noexcept = default;

  logger &operator=(const logger &Other) = delete;
  logger &operator=(logger &&Other) noexcept = default;

  log_level Level() const { return Level_; }

  void SetLevel(log_level Level);

  bool LoggingStatus() const { return (Level_ & log_level::STATUS) != log_level::NONE; }
  bool LoggingWarnings() const { return (Level_ & log_level::WARNINGS) != log_level::NONE; }
  bool LoggingErrors() const { return (Level_ & log_level::ERRORS) != log_level::NONE; }

  template <typename... Ts> void LogStatus(bool WriteCondition, int IncrementLevel, const
    std::string &Format, const Ts &... Args) const;
  template <typename... Ts> void LogWarning(bool WriteCondition, const std::string &Format, const
    Ts &... Args) const;
  template <typename... Ts> void LogError(bool WriteCondition, const std::string &Format, const
    Ts &... Args) const;

private:

  log_level Level_;
  int Rank_;

  void ReplaceRank_(std::string &Message) const;

};

}}

#include <ovk/core/Logger.inl>

#endif
