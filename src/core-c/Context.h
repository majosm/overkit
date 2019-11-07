// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_CONTEXT_H_INCLUDED
#define OVK_CORE_C_CONTEXT_H_INCLUDED

#include <ovk/core-c/Error.h>
#include <ovk/core-c/Global.h>

#include <mpi.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_context;
typedef struct ovk_context ovk_context;

struct ovk_shared_context;
typedef struct ovk_shared_context ovk_shared_context;

struct ovk_context_params;
typedef struct ovk_context_params ovk_context_params;

void ovkCreateContext(ovk_context **Context, ovk_context_params **Params, ovk_error *Error);
void ovkDestroyContext(ovk_context **Context);

void ovkShareContext(ovk_context **Context, ovk_shared_context **SharedContext);
void ovkResetSharedContext(ovk_shared_context **SharedContext);
void ovkGetContextFromSharedC(const ovk_shared_context *SharedContext, const ovk_context
  **Context);
void ovkGetContextFromShared(ovk_shared_context *SharedContext, ovk_context **Context);

void ovkGetContextComm(const ovk_context *Context, MPI_Comm *Comm);
void ovkGetContextCommSize(const ovk_context *Context, int *CommSize);
void ovkGetContextCommRank(const ovk_context *Context, int *CommRank);
void ovkGetContextErrorLogging(const ovk_context *Context, bool *LoggingErrors);
void ovkSetContextErrorLogging(ovk_context *Context, bool LoggingErrors);
void ovkGetContextWarningLogging(const ovk_context *Context, bool *LoggingWarnings);
void ovkSetContextWarningLogging(ovk_context *Context, bool LoggingWarnings);
void ovkGetContextStatusLoggingThreshold(const ovk_context *Context, int *StatusLoggingThreshold);
void ovkSetContextStatusLoggingThreshold(ovk_context *Context, int StatusLoggingThreshold);
void ovkGetContextProfiling(const ovk_context *Context, bool *Profiling);
void ovkSetContextProfiling(ovk_context *Context, bool Profiling);
void ovkWriteProfile(const ovk_context *Context, FILE *File);

void ovkCreateContextParams(ovk_context_params **Params);
void ovkDestroyContextParams(ovk_context_params **Params);
void ovkGetContextParamComm(const ovk_context_params *Params, MPI_Comm *Comm);
void ovkSetContextParamComm(ovk_context_params *Params, MPI_Comm Comm);
void ovkGetContextParamErrorLogging(const ovk_context_params *Params, bool *ErrorLogging);
void ovkSetContextParamErrorLogging(ovk_context_params *Params, bool ErrorLogging);
void ovkGetContextParamWarningLogging(const ovk_context_params *Params, bool *WarningLogging);
void ovkSetContextParamWarningLogging(ovk_context_params *Params, bool WarningLogging);
void ovkGetContextParamStatusLoggingThreshold(const ovk_context_params *Params, int
  *StatusLoggingThreshold);
void ovkSetContextParamStatusLoggingThreshold(ovk_context_params *Params, int
  StatusLoggingThreshold);
void ovkGetContextParamProfiling(const ovk_context_params *Params, bool *Profiling);
void ovkSetContextParamProfiling(ovk_context_params *Params, bool Profiling);

#ifdef __cplusplus
}
#endif

#endif
