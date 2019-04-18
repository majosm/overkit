// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ErrorHandler.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

namespace ovk {
namespace core {

error_handler::error_handler():
  Type_(error_handler_type::ABORT)
{}

error_handler::error_handler(error_handler_type Type):
  Type_(Type)
{
  OVK_DEBUG_ASSERT(ValidErrorHandlerType(Type_), "Invalid error handler type.");
}

void error_handler::SetType(error_handler_type Type) {

  OVK_DEBUG_ASSERT(ValidErrorHandlerType(Type), "Invalid error handler type.");

  Type_ = Type;

}

}}
