// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

namespace logger_internal {
void ReplaceRank(std::string &Message, int Rank);
}

template <typename... Ts> void LogStatus(logger &Logger, bool WriteCondition, int IncrementLevel,
  const std::string &Format, const Ts &... Args) {

  if (LoggingStatus(Logger) && WriteCondition) {
    std::string Message = "ovk :: ";
    for (int iIncrement = 0; iIncrement+1 < IncrementLevel; ++iIncrement) Message += "  ";
    Message += StringPrint(Format, Args...);
    logger_internal::ReplaceRank(Message, Logger.WriteRank_);
    fprintf(stdout, "%s\n", Message.c_str());
    fflush(stdout);
  }

}

template <typename... Ts> void LogWarning(logger &Logger, bool WriteCondition, const std::string
  &Format, const Ts &... Args) {

  if (LoggingWarnings(Logger) && WriteCondition) {
    std::string Message = "ovk :: WARNING: " + StringPrint(Format, Args...);
    logger_internal::ReplaceRank(Message, Logger.WriteRank_);
    fprintf(stderr, "%s\n", Message.c_str());
    fflush(stderr);
  }

}

template <typename... Ts> void LogError(logger &Logger, bool WriteCondition, const std::string
  &Format, const Ts &... Args) {

  if (LoggingErrors(Logger) && WriteCondition) {
    std::string Message = "ovk :: ERROR: " + StringPrint(Format, Args...);
    logger_internal::ReplaceRank(Message, Logger.WriteRank_);
    fprintf(stderr, "%s\n", Message.c_str());
    fflush(stderr);
  }

}

}}
