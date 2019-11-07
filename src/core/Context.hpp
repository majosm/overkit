// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONTEXT_HPP_INCLUDED
#define OVK_CORE_CONTEXT_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/Error.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Logger.hpp>
#include <ovk/core/Moveabool.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Profiler.hpp>

#include <mpi.h>

#include <memory>
#include <type_traits>

namespace ovk {

namespace context_internal {

// For doing stuff before creation and after destruction
class context_base {

protected:

  context_base(MPI_Comm Comm, bool LoggingErrors, bool LoggingWarnings, int StatusLoggingThreshold);

  context_base(const context_base &Other) = delete;
  context_base(context_base &&Other) noexcept = default;

  context_base &operator=(const context_base &Other) = delete;
  context_base &operator=(context_base &&Other) noexcept = default;

  ~context_base() noexcept;

  core::moveabool Exists_;

  comm Comm_;
  mutable core::logger Logger_;
  core::logger::status_level_and_indent_handle Level1_;

};

}

class context : private context_internal::context_base {

public:

  class params {
  public:
    params() = default;
    MPI_Comm Comm() const { return Comm_; }
    params &SetComm(MPI_Comm Comm);
    bool ErrorLogging() const { return ErrorLogging_; }
    params &SetErrorLogging(bool ErrorLogging);
    bool WarningLogging() const { return WarningLogging_; }
    params &SetWarningLogging(bool WarningLogging);
    int StatusLoggingThreshold() const { return StatusLoggingThreshold_; }
    params &SetStatusLoggingThreshold(int StatusLoggingThreshold);
    bool Profiling() const { return Profiling_; }
    params &SetProfiling(bool Profiling);
  private:
    MPI_Comm Comm_ = MPI_COMM_NULL;
    bool ErrorLogging_ = true;
    bool WarningLogging_ = true;
    int StatusLoggingThreshold_ = 1;
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

  bool LoggingErrors() const { return Logger_.LoggingErrors(); }
  void EnableErrorLogging();
  void DisableErrorLogging();

  bool LoggingWarnings() const { return Logger_.LoggingWarnings(); }
  void EnableWarningLogging();
  void DisableWarningLogging();

  int StatusLoggingThreshold() const { return Logger_.StatusThreshold(); }
  void SetStatusLoggingThreshold(int StatusLoggingThreshold);

  bool Profiling() const { return Profiler_.Enabled(); }
  void EnableProfiling();
  void DisableProfiling();
  std::string WriteProfile() const;

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
