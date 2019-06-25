// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONTEXT_HPP_INCLUDED
#define OVK_CORE_CONTEXT_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.h>
#include <ovk/core/Error.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Logger.hpp>
#include <ovk/core/Moveabool.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Profiler.hpp>

#include <mpi.h>

#include <memory>

namespace ovk {

enum class log_level {
  NONE = OVK_LOG_NONE,
  ERRORS = OVK_LOG_ERRORS,
  WARNINGS = OVK_LOG_WARNINGS,
  STATUS = OVK_LOG_STATUS,
  DEBUG = OVK_LOG_DEBUG,
  ALL = OVK_LOG_ALL
};

inline bool ValidLogLevel(log_level LogLevel) {
  return ovkValidLogLevel(ovk_log_level(LogLevel));
}

constexpr inline log_level operator|(log_level Left, log_level Right) {
  return log_level(int(Left) | int(Right));
}
constexpr inline log_level operator&(log_level Left, log_level Right) {
  return log_level(int(Left) & int(Right));
}
constexpr inline log_level operator^(log_level Left, log_level Right) {
  return log_level(int(Left) ^ int(Right));
}
constexpr inline log_level operator~(log_level LogLevel) {
  return log_level(~int(LogLevel));
}
inline log_level operator|=(log_level &Left, log_level Right) {
  return Left = Left | Right;
}
inline log_level operator&=(log_level &Left, log_level Right) {
  return Left = Left & Right;
}
inline log_level operator^=(log_level &Left, log_level Right) {
  return Left = Left ^ Right;
}

namespace context_internal {

// For doing stuff before creation and after destruction
class context_base {

protected:

  context_base(MPI_Comm Comm, log_level LogLevel);

  context_base(const context_base &Other) = delete;
  context_base(context_base &&Other) noexcept = default;

  context_base &operator=(const context_base &Other) = delete;
  context_base &operator=(context_base &&Other) noexcept = default;

  ~context_base() noexcept;

  core::moveabool Exists_;

  comm Comm_;
  // TODO: Maybe mutability should be encapsulated inside?
  mutable core::logger Logger_;

};

}

class context : private context_internal::context_base {

public:

  class params {
  public:
    params() = default;
    MPI_Comm Comm() const { return Comm_; }
    params &SetComm(MPI_Comm Comm);
    log_level LogLevel() const { return LogLevel_; }
    params &SetLogLevel(log_level LogLevel);
    bool Profiling() const { return Profiling_; }
    params &SetProfiling(bool Profiling);
  private:
    MPI_Comm Comm_ = MPI_COMM_NULL;
    log_level LogLevel_ = log_level::ERRORS | log_level::WARNINGS;
    bool Profiling_ = false;
    friend class context;
  };

  context(const context &Other) = delete;
  context(context &&Other) noexcept = default;

  context &operator=(const context &Other) = delete;
  context &operator=(context &&Other) noexcept = default;

  ~context() noexcept;

  floating_ref<const context> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<context> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  const comm &Comm() const { return Comm_; }

  log_level LogLevel() const;
  void SetLogLevel(log_level LogLevel);

  bool Profiling() const { return Profiler_.Enabled(); }
  void EnableProfiling();
  void DisableProfiling();

  core::logger &core_Logger() const { return Logger_; }
  core::profiler &core_Profiler() const { return Profiler_; }

  static context internal_Create(params &&Params);

private:

  floating_ref_generator FloatingRefGenerator_;

  // TODO: Maybe mutability should be encapsulated inside?
  mutable core::profiler Profiler_;

  context(params &&Params);

};

context CreateContext(context::params Params);
optional<context> CreateContext(context::params Params, error &Error);

}

#endif
