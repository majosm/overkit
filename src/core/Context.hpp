// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONTEXT_HPP_INCLUDED
#define OVK_CORE_CONTEXT_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
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

  core::comm Comm_;
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

  floating_ref<const context> GetFloatingRef() const { return FloatingRefGenerator_.Generate(); }
  floating_ref<context> GetFloatingRef() { return FloatingRefGenerator_.Generate(); }

  MPI_Comm Comm() const { return Comm_; }
  int CommSize() const { return Comm_.Size(); }
  int CommRank() const { return Comm_.Rank(); }

  log_level LogLevel() const { return Logger_.Level(); }
  void SetLogLevel(log_level LogLevel);

  bool Profiling() const { return Profiler_.Enabled(); }
  void EnableProfiling();
  void DisableProfiling();

  const core::comm &core_Comm() const { return Comm_; }
  core::logger &core_Logger() const { return Logger_; }
  core::profiler &core_Profiler() const { return Profiler_; }

  static context internal_Create(params &&Params);

private:

  floating_ref_generator<context> FloatingRefGenerator_;

  // TODO: Maybe mutability should be encapsulated inside?
  mutable core::profiler Profiler_;

  context(params &&Params);

};

context CreateContext(context::params Params);
optional<context> CreateContext(context::params Params, error &Error);

}

#endif
