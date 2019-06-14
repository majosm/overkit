// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Logger.hpp"

#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/TextProcessing.hpp"

namespace ovk {
namespace core {

logger::logger(int Rank):
  Rank_(Rank)
{
  OVK_DEBUG_ASSERT(Rank_ >= 0, "Invalid log rank.");
}

void logger::EnableErrorLogging() {

  LoggingErrors_ = true;

}

void logger::DisableErrorLogging() {

  LoggingErrors_ = false;

}

void logger::EnableWarningLogging() {

  LoggingWarnings_ = true;

}

void logger::DisableWarningLogging() {

  LoggingWarnings_ = false;

}

void logger::EnableStatusLogging() {

  LoggingStatus_ = true;

}

void logger::DisableStatusLogging() {

  LoggingStatus_ = false;

}

void logger::EnableDebugLogging() {

  LoggingDebug_ = true;

}

void logger::DisableDebugLogging() {

  LoggingDebug_ = false;

}

}}
