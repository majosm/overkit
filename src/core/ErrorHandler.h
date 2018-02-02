// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ERROR_HANDLER_INCLUDED
#define OVK_CORE_ERROR_HANDLER_INCLUDED

#include "Global.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct t_error_handler {
  ovk_error_handler_type type;
} t_error_handler;

#define OVK_HANDLE_ERROR(Handler, Error) \
  do { \
    switch (Handler->type) { \
    case OVK_ERROR_HANDLER_ABORT: \
      MPI_Abort(MPI_COMM_WORLD, Error); \
    case OVK_ERROR_HANDLER_RETURN: \
      return Error; \
    default: \
      return OVK_ERROR_NONE; \
    } \
  } while (0)

#define OVK_CHECK_ERROR(Handler, Error) \
  do { \
    if (Error != OVK_ERROR_NONE) { \
      OVK_HANDLE_ERROR(Handler, Error); \
    } \
  } while (0)

#define OVK_CHECK_ERROR_GOTO(Handler, Error, Label) \
  do { \
    if (Error != OVK_ERROR_NONE) { \
      goto Label; \
    } \
  } while (0)

#define OVK_CHECK_ERROR_COLLECTIVE(Handler, Error, Comm) \
  do { \
    ovk_error OVK_CE_Error = OVK_MAX_ERROR; \
    if (Error != OVK_ERROR_NONE) OVK_CE_Error = Error; \
    ovk_error OVK_CE_MinError; \
    MPI_Allreduce(&OVK_CE_Error, &OVK_CE_MinError, 1, MPI_INT, MPI_MIN, Comm); \
    if (OVK_CE_MinError != OVK_ERROR_NONE) { \
      OVK_HANDLE_ERROR(Handler, OVK_CE_MinError); \
    } \
  } while (0)

#define OVK_CHECK_ERROR_GOTO_COLLECTIVE(Handler, Error, Label, Comm) \
  do { \
    ovk_error OVK_CE_Error = OVK_MAX_ERROR; \
    if (Error != OVK_ERROR_NONE) OVK_CE_Error = Error; \
    ovk_error OVK_CE_MinError; \
    MPI_Allreduce(&OVK_CE_Error, &OVK_CE_MinError, 1, MPI_INT, MPI_MIN, Comm); \
    if (OVK_CE_MinError != OVK_ERROR_NONE) { \
      goto Label; \
    } \
  } while (0)

void PRIVATE(CreateErrorHandler)(t_error_handler **Handler, ovk_error_handler_type Type);
#define CreateErrorHandler(...) PRIVATE(CreateErrorHandler)(__VA_ARGS__)

void PRIVATE(DestroyErrorHandler)(t_error_handler **Handler);
#define DestroyErrorHandler(...) PRIVATE(DestroyErrorHandler)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
