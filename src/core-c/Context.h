// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_CONTEXT_H_INCLUDED
#define OVK_CORE_C_CONTEXT_H_INCLUDED

#include <ovk/core-c/Error.h>
#include <ovk/core-c/Global.h>
#include <ovk/core/Context.h>

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
void ovkGetContextLogLevel(const ovk_context *Context, ovk_log_level *LogLevel);
void ovkSetContextLogLevel(ovk_context *Context, ovk_log_level LogLevel);
void ovkGetContextProfiling(const ovk_context *Context, bool *Profiling);
void ovkSetContextProfiling(ovk_context *Context, bool Profiling);
void ovkWriteProfile(const ovk_context *Context, FILE *File);

void ovkCreateContextParams(ovk_context_params **Params);
void ovkDestroyContextParams(ovk_context_params **Params);
void ovkGetContextParamComm(const ovk_context_params *Params, MPI_Comm *Comm);
void ovkSetContextParamComm(ovk_context_params *Params, MPI_Comm Comm);
void ovkGetContextParamLogLevel(const ovk_context_params *Params, ovk_log_level *LogLevel);
void ovkSetContextParamLogLevel(ovk_context_params *Params, ovk_log_level LogLevel);
void ovkGetContextParamProfiling(const ovk_context_params *Params, bool *Profiling);
void ovkSetContextParamProfiling(ovk_context_params *Params, bool Profiling);

#ifdef __cplusplus
}
#endif

#endif
