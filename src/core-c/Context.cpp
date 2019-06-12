// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Context.h"

#include "ovk/core-c/Constants.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Context.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

void ovkCreateContext(ovk_context **Context, ovk_context_params **Params, ovk_error *Error) {

  if (OVK_DEBUG) {
    int MPIInitialized;
    MPI_Initialized(&MPIInitialized);
    // Can't use OVK_DEBUG_ASSERT here because it calls MPI_Abort
    if (!MPIInitialized) {
      fprintf(stderr, "ERROR: MPI not initialized.\n"); fflush(stderr);
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

void ovkGetContextFromSharedC(const ovk_shared_context **SharedContext, const ovk_context **Context)
  {

  OVK_DEBUG_ASSERT(SharedContext, "Invalid shared context pointer.");
  OVK_DEBUG_ASSERT(*SharedContext, "Invalid shared context pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto SharedContextCPPPtr = reinterpret_cast<const std::shared_ptr<const ovk::context> *>(
    *SharedContext);
  *Context = reinterpret_cast<const ovk_context *>(SharedContextCPPPtr->get());

}

void ovkGetContextFromShared(ovk_shared_context **SharedContext, ovk_context **Context) {

  OVK_DEBUG_ASSERT(SharedContext, "Invalid shared context pointer.");
  OVK_DEBUG_ASSERT(*SharedContext, "Invalid shared context pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto SharedContextCPPPtr = reinterpret_cast<std::shared_ptr<ovk::context> *>(*SharedContext);
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
  *CommSize = ContextCPP.CommSize();

}

void ovkGetContextCommRank(const ovk_context *Context, int *CommRank) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  *CommRank = ContextCPP.CommRank();

}

void ovkGetContextLogLevel(const ovk_context *Context, ovk_log_level *LogLevel) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(LogLevel, "Invalid log level pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  *LogLevel = ovk_log_level(ContextCPP.LogLevel());

}

void ovkSetContextLogLevel(ovk_context *Context, ovk_log_level LogLevel) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(ovkValidLogLevel(LogLevel), "Invalid log level.");

  auto &ContextCPP = *reinterpret_cast<ovk::context *>(Context);
  ContextCPP.SetLogLevel(ovk::log_level(LogLevel));

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

void ovkGetContextParamLogLevel(const ovk_context_params *Params, ovk_log_level *LogLevel) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LogLevel, "Invalid log level pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::context::params *>(Params);
  *LogLevel = ovk_log_level(ParamsCPP.LogLevel());

}

void ovkSetContextParamLogLevel(ovk_context_params *Params, ovk_log_level LogLevel) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::context::params *>(Params);
  ParamsCPP.SetLogLevel(ovk::log_level(LogLevel));

}
