// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Logger.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <string>

namespace ovk {
namespace core {

logger::logger(log_level Level, int Rank):
  Level_(Level),
  Rank_(Rank)
{

  OVK_DEBUG_ASSERT(ValidLogLevel(Level_), "Invalid log level.");
  OVK_DEBUG_ASSERT(Rank_ >= 0, "Invalid log rank.");

}

void logger::SetLevel(log_level Level) {

  OVK_DEBUG_ASSERT(ValidLogLevel(Level), "Invalid log level.");

  Level_ = Level;

}

}}
