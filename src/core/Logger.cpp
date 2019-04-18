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

logger::logger():
  Level_(log_level::NONE),
  Rank_(-1)
{}

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

void logger::ReplaceRank_(std::string &Message) const {

  size_t RankPos = Message.find("@rank@");

  if (RankPos == std::string::npos) return;

  std::string RankString = FormatNumber(Rank_);

  std::string NewMessage;

  size_t Pos = 0;
  while (RankPos != std::string::npos) {
    NewMessage.append(Message, Pos, RankPos-Pos);
    NewMessage += RankString;
    Pos = RankPos + 6;
    RankPos = Message.find("@rank@", Pos);
  }
  NewMessage.append(Message, Pos, std::string::npos);

  Message = NewMessage;

}

}}
