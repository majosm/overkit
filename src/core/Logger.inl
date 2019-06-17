// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename... Ts> void logger::LogError(bool WriteCondition, const std::string &Format,
  const Ts &... Args) const {

  if (LoggingErrors_ && WriteCondition) {
    std::string Prefix = "ovk :: [!] ERROR: ";
    std::string Message = StringPrint(Format, Args...);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank_));
    std::fprintf(stderr, "%s%s\n", Prefix.c_str(), Message.c_str());
    std::fflush(stderr);
  }

}

template <typename... Ts> void logger::LogWarning(bool WriteCondition, const std::string &Format,
  const Ts &... Args) const {

  if (LoggingWarnings_ && WriteCondition) {
    std::string Prefix = "ovk :: [!] WARNING: ";
    std::string Message = StringPrint(Format, Args...);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank_));
    std::fprintf(stderr, "%s%s\n", Prefix.c_str(), Message.c_str());
    std::fflush(stderr);
  }

}

template <typename... Ts> void logger::LogStatus(bool WriteCondition, int IncrementLevel, const
  std::string &Format, const Ts &... Args) const {

  if (LoggingStatus_ && WriteCondition) {
    std::string Prefix = "ovk :: [-] ";
    for (int iIncrement = 0; iIncrement+1 < IncrementLevel; ++iIncrement) Prefix += "  ";
    if (IncrementLevel > 0) Prefix += "* ";
    std::string Message = StringPrint(Format, Args...);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank_));
    std::fprintf(stdout, "%s%s\n", Prefix.c_str(), Message.c_str());
    std::fflush(stdout);
  }

}

template <typename... Ts> void logger::LogDebug(bool WriteCondition, int IncrementLevel, const
  std::string &Format, const Ts &... Args) const {

  if (LoggingDebug_ && WriteCondition) {
    std::string Prefix = "ovk :: [ ] ";
    for (int iIncrement = 0; iIncrement+1 < IncrementLevel; ++iIncrement) Prefix += "  ";
    if (IncrementLevel > 0) Prefix += "* ";
    std::string Message = StringPrint(Format, Args...);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank_));
    std::fprintf(stdout, "%s%s\n", Prefix.c_str(), Message.c_str());
    std::fflush(stdout);
  }

}

}}
