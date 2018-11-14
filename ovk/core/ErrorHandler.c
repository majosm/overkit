// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ErrorHandler.h"

#include "ovk/core/Global.h"

void PRIVATE(CreateErrorHandler)(t_error_handler **Handler_, ovk_error_handler_type Type) {

  *Handler_ = malloc(sizeof(t_error_handler));
  t_error_handler *Handler = *Handler_;

  Handler->type = Type;

}

void PRIVATE(DestroyErrorHandler)(t_error_handler **Handler_) {

  free(*Handler_);
  *Handler_ = NULL;

}

void PRIVATE(GetErrorHandlerType)(const t_error_handler *Handler, ovk_error_handler_type *Type) {

  OVK_DEBUG_ASSERT(Handler, "Invalid error handler pointer.");
  OVK_DEBUG_ASSERT(Type, "Invalid error handler type pointer.");

  *Type = Handler->type;

}

void PRIVATE(SetErrorHandlerType)(t_error_handler *Handler, ovk_error_handler_type Type) {

  OVK_DEBUG_ASSERT(Handler, "Invalid error handler pointer.");

  Handler->type = Type;

}
