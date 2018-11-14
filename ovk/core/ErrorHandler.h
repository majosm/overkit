// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ERROR_HANDLER_INCLUDED
#define OVK_CORE_ERROR_HANDLER_INCLUDED

#include "ovk/core/Global.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct t_error_handler {
  ovk_error_handler_type type;
} t_error_handler;

#define OVK_EH_INIT(Handler) ovk_error EH_Error = OVK_NO_ERROR

#define OVK_EH_HANDLE(Handler, Error) \
  do { \
    EH_Error = (Error); \
    goto eh_finalize; \
  } while (false)

#define OVK_EH_HANDLE_GOTO(Handler, Error, Label) \
  do { \
    EH_Error = (Error); \
    goto Label; \
  } while (false)

#define OVK_EH_CHECK(Handler, Error) \
  do { \
    if (Error != OVK_NO_ERROR) { \
      OVK_EH_HANDLE(Handler, Error); \
    } \
  } while (0)

#define OVK_EH_CHECK_GOTO(Handler, Error, Label) \
  do { \
    if (Error != OVK_NO_ERROR) { \
      OVK_EH_HANDLE_GOTO(Handler, Error, Label); \
    } \
  } while (0)

#define OVK_EH_CHECK_ALL(Handler, Error, Comm) \
  do { \
    ovk_error CE_Error = OVK_MAX_ERROR; \
    if (Error != OVK_NO_ERROR) CE_Error = Error; \
    ovk_error CE_MinError; \
    MPI_Allreduce(&CE_Error, &CE_MinError, 1, MPI_INT, MPI_MIN, Comm); \
    if (CE_MinError != OVK_MAX_ERROR) { \
      OVK_EH_HANDLE(Handler, CE_MinError); \
    } \
  } while (0)

#define OVK_EH_CHECK_ALL_GOTO(Handler, Error, Comm, Label) \
  do { \
    ovk_error CE_Error = OVK_MAX_ERROR; \
    if (Error != OVK_NO_ERROR) CE_Error = Error; \
    ovk_error CE_MinError; \
    MPI_Allreduce(&CE_Error, &CE_MinError, 1, MPI_INT, MPI_MIN, Comm); \
    if (CE_MinError != OVK_MAX_ERROR) { \
      OVK_EH_HANDLE_GOTO(Handler, CE_MinError, Label); \
    } \
  } while (0)

// Note: the first line is to suppress compiler warnings about unused labels
#define OVK_EH_FINALIZE(Handler) \
  if (EH_Error != OVK_NO_ERROR) goto eh_finalize; \
  eh_finalize: \
    if (EH_Error != OVK_NO_ERROR) { \
      switch (Handler->type) { \
      case OVK_ERROR_HANDLER_ABORT: \
        MPI_Abort(MPI_COMM_WORLD, EH_Error); \
      case OVK_ERROR_HANDLER_RETURN: \
        return EH_Error; \
      } \
    }

void PRIVATE(CreateErrorHandler)(t_error_handler **Handler, ovk_error_handler_type Type);
#define CreateErrorHandler(...) PRIVATE(CreateErrorHandler)(__VA_ARGS__)

void PRIVATE(DestroyErrorHandler)(t_error_handler **Handler);
#define DestroyErrorHandler(...) PRIVATE(DestroyErrorHandler)(__VA_ARGS__)

void PRIVATE(GetErrorHandlerType)(const t_error_handler *Handler, ovk_error_handler_type *Type);
#define GetErrorHandlerType(...) PRIVATE(GetErrorHandlerType)(__VA_ARGS__)

void PRIVATE(SetErrorHandlerType)(t_error_handler *Handler, ovk_error_handler_type Type);
#define SetErrorHandlerType(...) PRIVATE(SetErrorHandlerType)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
