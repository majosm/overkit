// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

inline logger::logger(int Rank):
  Rank_(Rank)
{
  OVK_DEBUG_ASSERT(Rank_ >= 0, "Invalid log rank.");
}

inline logger &logger::EnableErrors() {

  LoggingErrors_ = true;

  return *this;

}

inline logger &logger::DisableErrors() {

  LoggingErrors_ = false;

  return *this;

}

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

inline logger &logger::EnableWarnings() {

  LoggingWarnings_ = true;

  return *this;

}

inline logger &logger::DisableWarnings() {

  LoggingWarnings_ = false;

  return *this;

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

inline logger &logger::SetStatusThreshold(int StatusThreshold) {

  StatusThreshold_ = StatusThreshold;

  return *this;

}

template <typename... Ts> void logger::LogStatus(bool WriteCondition, const std::string &Format,
  const Ts &... Args) const {

  if (StatusLevel_ <= StatusThreshold_ && WriteCondition) {
    std::string Prefix = "ovk :: [ ] ";
    for (int iIndent = 0; iIndent+1 < StatusIndent_; ++iIndent) Prefix += "  ";
    if (StatusIndent_ > 0) Prefix += "* ";
    std::string Message = StringPrint(Format, Args...);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank_));
    std::fprintf(stdout, "%s%s\n", Prefix.c_str(), Message.c_str());
    std::fflush(stdout);
  }

}

inline logger::status_level_handle logger::IncreaseStatusLevel(int Amount) {

  StatusLevel_ += Amount;

  return status_level_handle(StatusLevel_, Amount);

}

inline logger::status_indent_handle logger::IndentStatus(int Amount) {

  StatusIndent_ += Amount;

  return status_indent_handle(StatusIndent_, Amount);

}

inline logger::status_level_and_indent_handle logger::IncreaseStatusLevelAndIndent(int Amount) {

  StatusLevel_ += Amount;
  StatusIndent_ += Amount;

  return status_level_and_indent_handle(StatusLevel_, StatusIndent_, Amount);

}

inline logger &logger::EnableDebug() {

  LoggingDebug_ = true;

  return *this;

}

inline logger &logger::DisableDebug() {

  LoggingDebug_ = false;

  return *this;

}

template <typename... Ts> void logger::LogDebug(bool WriteCondition, const std::string &Format,
  const Ts &... Args) const {

  if (LoggingDebug_ && WriteCondition) {
    std::string Prefix = "ovk :: [&] ";
    std::string Message = StringPrint(Format, Args...);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank_));
    std::fprintf(stdout, "%s%s\n", Prefix.c_str(), Message.c_str());
    std::fflush(stdout);
  }

}

}}
