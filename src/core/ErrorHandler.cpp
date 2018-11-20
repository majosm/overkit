// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ErrorHandler.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

namespace ovk {
namespace core {

void CreateErrorHandler(error_handler &Handler, error_handler_type Type) {

  OVK_DEBUG_ASSERT(ValidErrorHandlerType(Type), "Invalid error handler type.");

  Handler.Type_ = Type;

}

void DestroyErrorHandler(error_handler &Handler) {

  Handler.Type_ = error_handler_type::ABORT;

}

void GetErrorHandlerType(const error_handler &Handler, error_handler_type &Type) {

  Type = Handler.Type_;

}

void SetErrorHandlerType(error_handler &Handler, error_handler_type Type) {

  OVK_DEBUG_ASSERT(ValidErrorHandlerType(Type), "Invalid error handler type.");

  Handler.Type_ = Type;

}

}}
