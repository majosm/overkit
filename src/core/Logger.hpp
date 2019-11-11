// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_LOGGER_HPP_INCLUDED
#define OVK_CORE_LOGGER_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/TextProcessing.hpp>

#include <string>
#include <utility>

#include <mpi.h>

namespace ovk {
namespace core {

class logger {

  class increase_handle {
  public:
    increase_handle() = default;
    increase_handle(const increase_handle &Other) = delete;
    increase_handle(increase_handle &&Other) noexcept:
      ValuePtr_(Other.ValuePtr_),
      Amount_(Other.Amount_)
    {
      Other.ValuePtr_ = nullptr;
      Other.Amount_ = 0;
    }
    ~increase_handle() noexcept {
      if (ValuePtr_) {
        *ValuePtr_ -= Amount_;
      }
    }
    increase_handle &operator=(increase_handle Other) noexcept {
      std::swap(ValuePtr_, Other.ValuePtr_);
      std::swap(Amount_, Other.Amount_);
      return *this;
    }
    increase_handle &Reset() {
      *this = increase_handle();
      return *this;
    }
  protected:
    int *ValuePtr_ = nullptr;
    int Amount_ = 0;
    explicit increase_handle(int &Value, int Amount):
      ValuePtr_(&Value),
      Amount_(Amount)
    {}
    friend class logger;
  };

public:

  class status_level_handle : private increase_handle {
  private:
    using increase_handle::increase_handle;
    friend class logger;
  };

  class status_indent_handle : private increase_handle {
  private:
    using increase_handle::increase_handle;
    friend class logger;
  };

  class status_level_and_indent_handle {
  public:
    status_level_and_indent_handle() = default;
    status_level_and_indent_handle &Reset() {
      LevelHandle_.Reset();
      IndentHandle_.Reset();
      return *this;
    }
  private:
    status_level_handle LevelHandle_;
    status_indent_handle IndentHandle_;
    status_level_and_indent_handle(int &Level, int &Indent, int Amount):
      LevelHandle_(Level, Amount),
      IndentHandle_(Indent, Amount)
    {}
    friend class logger;
  };

  logger(int Rank);

  logger(const logger &Other) = delete;
  logger(logger &&Other) noexcept = default;

  logger &operator=(const logger &Other) = delete;
  logger &operator=(logger &&Other) noexcept = default;

  logger &SyncIndicator(comm_view Comm);

  bool LoggingErrors() const { return LoggingErrors_; }
  logger &EnableErrors();
  logger &DisableErrors();
  template <typename... Ts> void LogError(bool WriteCondition, const std::string &Format, const
    Ts &... Args);

  bool LoggingWarnings() const { return LoggingWarnings_; }
  logger &EnableWarnings();
  logger &DisableWarnings();
  template <typename... Ts> void LogWarning(bool WriteCondition, const std::string &Format, const
    Ts &... Args);

  bool LoggingStatus() const { return StatusLevel_ <= StatusThreshold_; }
  int StatusThreshold() const { return StatusThreshold_; }
  logger &SetStatusThreshold(int StatusThreshold);
  template <typename... Ts> void LogStatus(bool WriteCondition, const std::string &Format, const Ts
    &... Args);
  status_level_handle IncreaseStatusLevel(int IncreaseAmount=1);
  status_indent_handle IndentStatus(int IndentAmount=1);
  status_level_and_indent_handle IncreaseStatusLevelAndIndent(int IncreaseAmount=1);

  bool LoggingDebug() const { return LoggingDebug_; }
  logger &EnableDebug();
  logger &DisableDebug();
  template <typename... Ts> void LogDebug(bool WriteCondition, const std::string &Format, const Ts
    &... Args);

private:

  int Rank_;
  bool Indicator_ = false;

  bool LoggingErrors_ = false;
  bool LoggingWarnings_ = false;

  int StatusThreshold_ = 0;
  int StatusLevel_ = 1;
  int StatusIndent_ = 0;

  bool LoggingDebug_ = false;

};

}}

#include <ovk/core/Logger.inl>

#endif
