// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Context.h"

#include "ovk/core-c/Global.h"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <cstdio>
#include <memory>
#include <utility>

void ovkCreateContext(ovk_context **Context, ovk_context_params **Params, ovk_error *Error) {

  if (OVK_DEBUG) {
    int MPIInitialized;
    MPI_Initialized(&MPIInitialized);
    // Can't use OVK_DEBUG_ASSERT here because it calls MPI_Abort
    if (!MPIInitialized) {
      std::fprintf(stderr, "ERROR: MPI not initialized.\n");
      std::fflush(stderr);
      exit(1);
    }
  }

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::context::params *>(*Params);

  ovk::error ErrorCPP;
  ovk::optional<ovk::context> MaybeContextCPP = ovk::CreateContext(std::move(*ParamsCPPPtr),
    ErrorCPP);

  delete ParamsCPPPtr;

  if (ErrorCPP == ovk::error::NONE) {
    auto ContextCPPPtr = new ovk::context(MaybeContextCPP.Release());
    *Context = reinterpret_cast<ovk_context *>(ContextCPPPtr);
  } else {
    *Context = nullptr;
  }

  *Params = nullptr;
  *Error = ovk_error(ErrorCPP);

}

void ovkDestroyContext(ovk_context **Context) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(*Context, "Invalid context pointer.");

  auto ContextCPPPtr = reinterpret_cast<ovk::context *>(*Context);

  delete ContextCPPPtr;

  *Context = nullptr;

}

void ovkShareContext(ovk_context **Context, ovk_shared_context **SharedContext) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(*Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(SharedContext, "Invalid shared context pointer.");

  auto ContextCPPPtr = reinterpret_cast<ovk::context *>(*Context);

  auto SharedContextCPPPtr = new std::shared_ptr<ovk::context>(std::make_shared<ovk::context>(
    std::move(*ContextCPPPtr)));

  delete ContextCPPPtr;

  *Context = reinterpret_cast<ovk_context *>(SharedContextCPPPtr->get());
  *SharedContext = reinterpret_cast<ovk_shared_context *>(SharedContextCPPPtr);

}

void ovkResetSharedContext(ovk_shared_context **SharedContext) {

  OVK_DEBUG_ASSERT(SharedContext, "Invalid shared context pointer.");
  OVK_DEBUG_ASSERT(*SharedContext, "Invalid shared context pointer.");

  auto SharedContextCPPPtr = reinterpret_cast<std::shared_ptr<ovk::context> *>(*SharedContext);

  delete SharedContextCPPPtr;

  *SharedContext = nullptr;

}

void ovkGetContextFromSharedC(const ovk_shared_context *SharedContext, const ovk_context **Context)
  {

  OVK_DEBUG_ASSERT(SharedContext, "Invalid shared context pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto SharedContextCPPPtr = reinterpret_cast<const std::shared_ptr<const ovk::context> *>(
    SharedContext);
  *Context = reinterpret_cast<const ovk_context *>(SharedContextCPPPtr->get());

}

void ovkGetContextFromShared(ovk_shared_context *SharedContext, ovk_context **Context) {

  OVK_DEBUG_ASSERT(SharedContext, "Invalid shared context pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto SharedContextCPPPtr = reinterpret_cast<std::shared_ptr<ovk::context> *>(SharedContext);
  *Context = reinterpret_cast<ovk_context *>(SharedContextCPPPtr->get());

}

void ovkGetContextComm(const ovk_context *Context, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  *Comm = ContextCPP.Comm();

}

void ovkGetContextCommSize(const ovk_context *Context, int *CommSize) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  *CommSize = ContextCPP.Comm().Size();

}

void ovkGetContextCommRank(const ovk_context *Context, int *CommRank) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  *CommRank = ContextCPP.Comm().Rank();

}

void ovkGetContextErrorLogging(const ovk_context *Context, bool *LoggingErrors) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(LoggingErrors, "Invalid logging errors pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  *LoggingErrors = ContextCPP.LoggingErrors();

}

void ovkSetContextErrorLogging(ovk_context *Context, bool LoggingErrors) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ContextCPP = *reinterpret_cast<ovk::context *>(Context);
  if (LoggingErrors) {
    ContextCPP.EnableErrorLogging();
  } else {
    ContextCPP.DisableErrorLogging();
  }

}

void ovkGetContextWarningLogging(const ovk_context *Context, bool *LoggingWarnings) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(LoggingWarnings, "Invalid logging warnings pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  *LoggingWarnings = ContextCPP.LoggingWarnings();

}

void ovkSetContextWarningLogging(ovk_context *Context, bool LoggingWarnings) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ContextCPP = *reinterpret_cast<ovk::context *>(Context);
  if (LoggingWarnings) {
    ContextCPP.EnableWarningLogging();
  } else {
    ContextCPP.DisableWarningLogging();
  }

}

void ovkGetContextStatusLoggingThreshold(const ovk_context *Context, int *StatusLoggingThreshold) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(StatusLoggingThreshold, "Invalid status logging threshold pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  *StatusLoggingThreshold = ContextCPP.StatusLoggingThreshold();

}

void ovkSetContextStatusLoggingThreshold(ovk_context *Context, int StatusLoggingThreshold) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ContextCPP = *reinterpret_cast<ovk::context *>(Context);
  ContextCPP.SetStatusLoggingThreshold(StatusLoggingThreshold);

}

void ovkGetContextProfiling(const ovk_context *Context, bool *Profiling) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Profiling, "Invalid profiling pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  *Profiling = ContextCPP.Profiling();

}

void ovkSetContextProfiling(ovk_context *Context, bool Profiling) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ContextCPP = *reinterpret_cast<ovk::context *>(Context);
  if (Profiling) {
    ContextCPP.EnableProfiling();
  } else {
    ContextCPP.DisableProfiling();
  }

}

void ovkWriteProfile(const ovk_context *Context, FILE *File) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);

  std::string ProfileString = ContextCPP.WriteProfile();

  if (ContextCPP.Comm().Rank() == 0) {
    std::fprintf(File, "%s", ProfileString.c_str());
    std::fflush(File);
  }

}

void ovkCreateContextParams(ovk_context_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::context::params();

  *Params = reinterpret_cast<ovk_context_params *>(ParamsCPPPtr);

}

void ovkDestroyContextParams(ovk_context_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::context::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetContextParamComm(const ovk_context_params *Params, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::context::params *>(Params);
  *Comm = ParamsCPP.Comm();

}

void ovkSetContextParamComm(ovk_context_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::context::params *>(Params);
  ParamsCPP.SetComm(Comm);

}

void ovkGetContextParamErrorLogging(const ovk_context_params *Params, bool *ErrorLogging) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(ErrorLogging, "Invalid error logging pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::context::params *>(Params);
  *ErrorLogging = ParamsCPP.ErrorLogging();

}

void ovkSetContextParamErrorLogging(ovk_context_params *Params, bool ErrorLogging) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::context::params *>(Params);
  ParamsCPP.SetErrorLogging(ErrorLogging);

}

void ovkGetContextParamWarningLogging(const ovk_context_params *Params, bool *WarningLogging) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(WarningLogging, "Invalid warning logging pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::context::params *>(Params);
  *WarningLogging = ParamsCPP.WarningLogging();

}

void ovkSetContextParamWarningLogging(ovk_context_params *Params, bool WarningLogging) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::context::params *>(Params);
  ParamsCPP.SetWarningLogging(WarningLogging);

}

void ovkGetContextParamStatusLoggingThreshold(const ovk_context_params *Params, int
  *StatusLoggingThreshold) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(StatusLoggingThreshold, "Invalid status logging threshold pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::context::params *>(Params);
  *StatusLoggingThreshold = ParamsCPP.StatusLoggingThreshold();

}

void ovkSetContextParamStatusLoggingThreshold(ovk_context_params *Params, int
  StatusLoggingThreshold) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::context::params *>(Params);
  ParamsCPP.SetStatusLoggingThreshold(StatusLoggingThreshold);

}

void ovkGetContextParamProfiling(const ovk_context_params *Params, bool *Profiling) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Profiling, "Invalid profiling pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::context::params *>(Params);
  *Profiling = ParamsCPP.Profiling();

}

void ovkSetContextParamProfiling(ovk_context_params *Params, bool Profiling) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::context::params *>(Params);
  ParamsCPP.SetProfiling(Profiling);

}
