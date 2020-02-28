// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "support/CommandArgs.hpp"

#include <ovk/core/Array.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Set.hpp>

#include <cstdio>
#include <exception>
#include <string>
#include <utility>

namespace support {

command_args::command_args(bool WriteCondition):
  WriteCondition_(WriteCondition)
{}

command_args::command_args(ovk::map<std::string,std::string> Options, ovk::array<std::string>
  Arguments, bool WriteCondition):
  WriteCondition_(WriteCondition),
  Options_(std::move(Options)),
  Arguments_(std::move(Arguments))
{}

command_args_parser::command_args_parser(bool WriteCondition):
  WriteCondition_(WriteCondition)
{}

command_args command_args_parser::Parse(const ovk::array<std::string> &RawArguments) const {

  ovk::map<std::string,std::string> ParsedOptions;
  ovk::array<std::string> ParsedArguments;

  bool Help = false;

  ovk::map<std::string, option_data> AugmentedOptions = Options_;
  AugmentedOptions.Insert("help", 'h', value_type_flag::BOOL, "Display help information");

  ovk::map<char, std::string> AugmentedOptionShortNameToName = OptionShortNameToName_;
  AugmentedOptionShortNameToName.Insert('h', "help");

  auto ValidateOptionValue = [](const std::string &ValueString, value_type_flag
    ValueTypeFlag) -> std::string {
    std::string Result;
    bool Valid;
    switch (ValueTypeFlag) {
    case value_type_flag::BOOL:
      Valid = command_args_internal::Convert<bool>(ValueString).Present();
      break;
    case value_type_flag::INT:
      Valid = command_args_internal::Convert<int>(ValueString).Present();
      break;
    case value_type_flag::LONG_LONG:
      Valid = command_args_internal::Convert<long long>(ValueString).Present();
      break;
    case value_type_flag::DOUBLE:
      Valid = command_args_internal::Convert<double>(ValueString).Present();
      break;
    case value_type_flag::STRING:
      Valid = command_args_internal::Convert<std::string>(ValueString).Present();
      break;
    case value_type_flag::CHAR:
      Valid = command_args_internal::Convert<char>(ValueString).Present();
      break;
    default:
      Valid = false;
      break;
    }
    if (Valid) Result = ValueString;
    return Result;
  };

  int iArgument = 1;
  while (iArgument < RawArguments.Count()) {

    std::string RawArgument = RawArguments(iArgument);

    std::string NextRawArgument;
    if (iArgument+1 < RawArguments.Count()) {
      NextRawArgument = RawArguments(iArgument+1);
    }

    bool SkipNext = false;

    if (RawArgument[0] == '-') {

      std::string OptionName;
      std::string OptionValue;

      if (RawArgument[1] == '-') {

        std::size_t iEquals = RawArgument.find('=');
        if (iEquals == std::string::npos) iEquals = RawArgument.length();

        std::size_t NamePos = 2;
        std::size_t NameLength = iEquals-NamePos;

        OptionName = RawArgument.substr(NamePos, NameLength);

        auto OptionIter = AugmentedOptions.Find(OptionName);
        if (OptionIter == AugmentedOptions.End()) {
          if (WriteCondition_) {
            std::fprintf(stderr, "ERROR: Invalid option '%s'.\n", OptionName.c_str());
            std::fflush(stderr);
          }
          throw command_args_invalid_argument_error();
        }

        const option_data &OptionData = OptionIter->Value();

        std::size_t ValuePos = ovk::Min(iEquals+1, RawArgument.length());
        std::size_t ValueLength = RawArgument.length() - ValuePos;

        OptionValue = ValidateOptionValue(RawArgument.substr(ValuePos, ValueLength),
          OptionData.ValueTypeFlag);

        if (OptionValue == "" && OptionData.ValueTypeFlag == value_type_flag::BOOL) {
          // Don't allow equals sign with no value following
          if (iEquals == RawArgument.length()) {
            OptionValue = "true";
          }
        }

      } else {

        std::size_t iPos = 1;
        while (iPos < RawArgument.length()) {

          char OptionShortName = RawArgument[iPos];

          auto NameIter = AugmentedOptionShortNameToName.Find(OptionShortName);
          if (NameIter == AugmentedOptionShortNameToName.End()) {
            if (WriteCondition_) {
              std::fprintf(stderr, "ERROR: Invalid option '%c'.\n", OptionShortName);
              std::fflush(stderr);
            }
            throw command_args_invalid_argument_error();
          }

          OptionName = NameIter->Value();

          const option_data &OptionData = AugmentedOptions(OptionName);

          std::size_t ValuePos = iPos+1;
          std::size_t ValueLength = RawArgument.length()-ValuePos;

          if (ValueLength > 0) {
            OptionValue = ValidateOptionValue(RawArgument.substr(ValuePos, ValueLength),
              OptionData.ValueTypeFlag);
            if (OptionValue != "") iPos = RawArgument.length()-1;
          } else {
            OptionValue = ValidateOptionValue(NextRawArgument, OptionData.ValueTypeFlag);
            if (OptionValue != "") SkipNext = true;
          }

          if (OptionValue == "" && OptionData.ValueTypeFlag == value_type_flag::BOOL) {
            OptionValue = "true";
          }

          ++iPos;

        }

      }

      if (OptionName == "help") {
        Help = *command_args_internal::Convert<bool>(OptionValue);
      }

      ParsedOptions.Insert(OptionName, OptionValue);

    } else {

      ParsedArguments.Append(RawArgument);

    }

    iArgument += SkipNext ? 2 : 1;

  }

  if (Help && WriteCondition_) {
    if (HelpUsage_ != "") {
      std::printf("Usage: %s\n\n", HelpUsage_.c_str());
    }
    if (HelpDescription_ != "") {
      std::printf("Description: %s\n\n", HelpDescription_.c_str());
    }
    if (HelpLongDescription_ != "") {
      std::printf("%s\n\n", HelpLongDescription_.c_str());
    }
    std::printf("Options:\n");
    for (auto &Entry : AugmentedOptions) {
      const std::string &OptionName = Entry.Key();
      const option_data &OptionData = Entry.Value();
      std::printf("  --%s (-%c) -- %s\n", OptionName.c_str(), OptionData.ShortName,
        OptionData.Description.c_str());
    }
    std::fflush(stdout);
  }

  return {std::move(ParsedOptions), std::move(ParsedArguments), WriteCondition_};

}

ovk::optional<command_args> command_args_parser::Parse(const ovk::array<std::string> &RawArguments,
  captured_command_args_error &Error) const {

  Error.Reset();

  ovk::optional<command_args> MaybeCommandArgs;

  try {
    MaybeCommandArgs = Parse(RawArguments);
  } catch (const command_args_error &ParseError) {
    Error = ParseError.Capture();
  }

  return MaybeCommandArgs;

}

namespace command_args_internal {

namespace convert_internal {

namespace {

const ovk::set<std::string> TrueValues = {
  "true",
  "True",
  "TRUE"
};

const ovk::set<std::string> FalseValues = {
  "false",
  "False",
  "FALSE"
};

}

ovk::optional<bool> ConvertHelper(const std::string &ValueString, type_tag<bool>) {

  ovk::optional<bool> MaybeValue;

  if (TrueValues.Contains(ValueString)) {
    MaybeValue = true;
  } else if (FalseValues.Contains(ValueString)) {
    MaybeValue = false;
  } else {
    try {
      MaybeValue = bool(std::stoll(ValueString));
    } catch (const std::exception &) {}
  }

  return MaybeValue;

}

ovk::optional<int> ConvertHelper(const std::string &ValueString, type_tag<int>) {

  ovk::optional<int> MaybeValue;

  try {
    MaybeValue = std::stoi(ValueString);
  } catch (const std::exception &) {}

  return MaybeValue;

}

ovk::optional<long long> ConvertHelper(const std::string &ValueString, type_tag<long long>) {

  ovk::optional<long long> MaybeValue;

  try {
    MaybeValue = std::stoll(ValueString);
  } catch (const std::exception &) {}

  return MaybeValue;

}

ovk::optional<double> ConvertHelper(const std::string &ValueString, type_tag<double>) {

  ovk::optional<double> MaybeValue;

  try {
    MaybeValue = std::stoll(ValueString);
  } catch (const std::exception &) {}

  return MaybeValue;

}

ovk::optional<std::string> ConvertHelper(const std::string &ValueString, type_tag<std::string>) {

  return ValueString;

}

ovk::optional<char> ConvertHelper(const std::string &ValueString, type_tag<char>) {

  ovk::optional<char> MaybeValue;

  if (ValueString.length() == 1) {
    MaybeValue = ValueString[0];
  }

  return MaybeValue;

}

}

}

}
