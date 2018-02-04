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