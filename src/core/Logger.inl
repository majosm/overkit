// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename... Ts> void logger::LogStatus(bool WriteCondition, int IncrementLevel, const
  std::string &Format, const Ts &... Args) const {

  if (LoggingStatus() && WriteCondition) {
    std::string Message = "ovk :: ";
    for (int iIncrement = 0; iIncrement+1 < IncrementLevel; ++iIncrement) Message += "  ";
    Message += StringPrint(Format, Args...);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank_));
    fprintf(stdout, "%s\n", Message.c_str());
    fflush(stdout);
  }

}

template <typename... Ts> void logger::LogWarning(bool WriteCondition, const std::string &Format,
  const Ts &... Args) const {

  if (LoggingWarnings() && WriteCondition) {
    std::string Message = "ovk :: WARNING: " + StringPrint(Format, Args...);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank_));
    fprintf(stderr, "%s\n", Message.c_str());
    fflush(stderr);
  }

}

template <typename... Ts> void logger::LogError(bool WriteCondition, const std::string &Format,
  const Ts &... Args) const {

  if (LoggingErrors() && WriteCondition) {
    std::string Message = "ovk :: ERROR: " + StringPrint(Format, Args...);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank_));
    fprintf(stderr, "%s\n", Message.c_str());
    fflush(stderr);
  }

}

}}
