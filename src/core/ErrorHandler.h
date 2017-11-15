// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ERROR_HANDLER_INCLUDED
#define OVK_CORE_ERROR_HANDLER_INCLUDED

#include "Global.h"

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

static inline bool ValidErrorHandlerType(ovk_error_handler_type Type) {

  switch (Type) {
  case OVK_ERROR_HANDLER_ABORT:
  case OVK_ERROR_HANDLER_RETURN:
    return true;
  default:
    return false;
  }

}

static inline void CreateErrorHandler(t_error_handler **Handler_, ovk_error_handler_type Type) {

  *Handler_ = malloc(sizeof(t_error_handler));
  t_error_handler *Handler = *Handler_;

  Handler->type = Type;

}

static inline void DestroyErrorHandler(t_error_handler **Handler_) {

  free(*Handler_);
  *Handler_ = NULL;

}

#endif
