// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_COMMAND_ARGS_HPP_LOADED
#define OVK_SUPPORT_COMMAND_ARGS_HPP_LOADED

#include <support/CommandArgs.h>

#include <ovk/core/Array.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Exception.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Optional.hpp>

#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>

namespace support {

enum class command_args_error_code : typename std::underlying_type<support_command_args_error>::type
  {
  NONE = SUPPORT_COMMAND_ARGS_ERROR_NONE,
  INVALID_ARGUMENT = SUPPORT_COMMAND_ARGS_ERROR_INVALID_ARGUMENT
};

class command_args_error;

using captured_command_args_error = ovk::core::captured_exception<command_args_error>;

}

namespace ovk {
namespace core {
template <> struct exception_traits<support::command_args_error> {
  using code_type = support::command_args_error_code;
  using underlying_type = typename std::underlying_type<support::command_args_error_code>::type;
  static constexpr support::command_args_error_code SuccessCode =
    support::command_args_error_code::NONE;
  static support::captured_command_args_error CaptureFromCode(support::command_args_error_code Code);
  static support::command_args_error_code GetCode(const support::command_args_error &Error);
  static void SyncAuxData(support::command_args_error &Error, int Root, comm_view Comm);
};
}}

namespace support {

class command_args_error : public std::runtime_error {
public:
  command_args_error(command_args_error_code Code, const char *ErrorString):
    runtime_error(ErrorString),
    Code_(Code)
  {}
  virtual ~command_args_error() noexcept {}
  command_args_error_code Code() const { return Code_; }
  virtual captured_command_args_error Capture() const = 0;
  virtual void SyncAuxData(int Root, ovk::comm_view Comm) {}
private:
  command_args_error_code Code_;
};

class command_args_invalid_argument_error : public command_args_error {
public:
  command_args_invalid_argument_error():
    command_args_error(command_args_error_code::INVALID_ARGUMENT,
      "support::command_args_invalid_argument_error")
  {}
  virtual captured_command_args_error Capture() const override { return *this; }
};

}

namespace ovk {
namespace core {
inline support::captured_command_args_error exception_traits<support::command_args_error>::
  CaptureFromCode(support::command_args_error_code Code) {
  switch (Code) {
  case support::command_args_error_code::INVALID_ARGUMENT:
    return support::command_args_invalid_argument_error();
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    return {};
  }
}
inline support::command_args_error_code exception_traits<support::command_args_error>::GetCode(const
  support::command_args_error &Error) {
  return Error.Code();
}
inline void exception_traits<support::command_args_error>::SyncAuxData(support::command_args_error
  &Error, int Root, comm_view Comm) {
  Error.SyncAuxData(Root, Comm);
}
}}

namespace support {

namespace command_args_internal {

enum class value_type_flag {
  BOOL,
  INT,
  LONG_LONG,
  DOUBLE,
  STRING,
  CHAR
};

namespace value_type_flag_internal {
template <typename T> struct helper {};
template <> struct helper<bool> {
  static constexpr value_type_flag GetFlag() { return value_type_flag::BOOL; }
};
template <> struct helper<int> {
  static constexpr value_type_flag GetFlag() { return value_type_flag::INT; }
};
template <> struct helper<long long> {
  static constexpr value_type_flag GetFlag() { return value_type_flag::LONG_LONG; }
};
template <> struct helper<double> {
  static constexpr value_type_flag GetFlag() { return value_type_flag::DOUBLE; }
};
template <> struct helper<std::string> {
  static constexpr value_type_flag GetFlag() { return value_type_flag::STRING; }
};
template <> struct helper<char> {
  static constexpr value_type_flag GetFlag() { return value_type_flag::CHAR; }
};
}
template <typename T> constexpr value_type_flag GetValueTypeFlag_() {
  return value_type_flag_internal::helper<T>::GetFlag();
}

namespace convert_internal {
template <typename T> struct type_tag {};
ovk::optional<bool> ConvertHelper(const std::string &ValueString, type_tag<bool>);
ovk::optional<int> ConvertHelper(const std::string &ValueString, type_tag<int>);
ovk::optional<long long> ConvertHelper(const std::string &ValueString, type_tag<long long>);
ovk::optional<double> ConvertHelper(const std::string &ValueString, type_tag<double>);
ovk::optional<std::string> ConvertHelper(const std::string &ValueString, type_tag<std::string>);
ovk::optional<char> ConvertHelper(const std::string &ValueString, type_tag<char>);
}
template <typename T> ovk::optional<T> Convert(const std::string &ValueString) {
  return convert_internal::ConvertHelper(ValueString, convert_internal::type_tag<T>());
}

}

class command_args {

public:

  command_args() = default;
  explicit command_args(bool WriteCondition);
  command_args(ovk::map<std::string,std::string> Options, ovk::array<std::string> Arguments, bool
    WriteCondition=true);

  bool WriteCondition() const { return WriteCondition_; }

  bool OptionIsPresent(const std::string &Name) const {
    return Options_.Contains(Name);
  }

  template <typename T> T GetOptionValue(const std::string &Name) const {
    OVK_DEBUG_ASSERT(Options_.Contains(Name), "Option not present.");
    T Value;
    const std::string &ValueString = Options_(Name);
    auto MaybeValue = command_args_internal::Convert<T>(ValueString);
    if (MaybeValue) {
      Value = MaybeValue.Release();
    } else {
      if (WriteCondition_) {
        std::fprintf(stderr, "ERROR: Invalid value '%s' for option '%s'.\n", ValueString.c_str(),
          Name.c_str());
        std::fflush(stderr);
      }
      throw command_args_invalid_argument_error();
    }
    return Value;
  }

  template <typename T> ovk::optional<T> GetOptionValue(const std::string &Name,
    captured_command_args_error &Error) const {
    Error.Reset();
    ovk::optional<T> MaybeValue;
    try {
      MaybeValue = GetOptionValue<T>(Name);
    } catch (const command_args_error &OptionError) {
      Error = OptionError.Capture();
    }
    return MaybeValue;
  }

  template <typename T> T GetOptionValue(const std::string &Name, T DefaultValue) const {
    T Value;
    auto Iter = Options_.Find(Name);
    if (Iter != Options_.End()) {
      const std::string &ValueString = Iter->Value();
      auto MaybeValue = command_args_internal::Convert<T>(ValueString);
      if (MaybeValue) {
        Value = MaybeValue.Release();
      } else {
        if (WriteCondition_) {
          std::fprintf(stderr, "ERROR: Invalid value '%s' for option '%s'.\n", ValueString.c_str(),
            Name.c_str());
          std::fflush(stderr);
        }
        throw command_args_invalid_argument_error();
      }
    } else {
      Value = std::move(DefaultValue);
    }
    return Value;
  }

  template <typename T> ovk::optional<T> GetOptionValue(const std::string &Name, T DefaultValue,
    captured_command_args_error &Error) const {
    Error.Reset();
    ovk::optional<T> MaybeValue;
    try {
      MaybeValue = GetOptionValue<T>(Name, std::move(DefaultValue));
    } catch (const command_args_error &OptionError) {
      Error = OptionError.Capture();
    }
    return MaybeValue;
  }

  int ArgumentCount() const { return Arguments_.Count(); }

  template <typename T> T GetArgument(int iArgument) const {
    OVK_DEBUG_ASSERT(iArgument < Arguments_.Count(), "Invalid argument index.");
    T Argument;
    auto MaybeArgument = command_args_internal::Convert<T>(Arguments_(iArgument));
    if (MaybeArgument) {
      Argument = MaybeArgument.Release();
    } else {
      if (WriteCondition_) {
        std::fprintf(stderr, "ERROR: Invalid argument '%s'.\n", Arguments_(iArgument).c_str());
        std::fflush(stderr);
      }
      throw command_args_invalid_argument_error();
    }
    return Argument;
  }

  template <typename T> ovk::optional<T> GetArgument(int iArgument, captured_command_args_error
    &Error) const {
    Error.Reset();
    ovk::optional<T> MaybeArgument;
    try {
      MaybeArgument = GetArgument<T>(iArgument);
    } catch (const command_args_error &ArgumentError) {
      Error = ArgumentError.Capture();
    }
    return MaybeArgument;
  }

private:

  bool WriteCondition_ = true;
  ovk::map<std::string,std::string> Options_;
  ovk::array<std::string> Arguments_;

};

class command_args_parser {

public:

  command_args_parser() = default;
  explicit command_args_parser(bool WriteCondition);

  bool WriteCondition() const { return WriteCondition_; }

  command_args_parser &SetHelpUsage(std::string HelpUsage) {
    HelpUsage_ = std::move(HelpUsage);
    return *this;
  }

  command_args_parser &SetHelpDescription(std::string HelpDescription) {
    HelpDescription_ = std::move(HelpDescription);
    return *this;
  }

  command_args_parser &SetHelpLongDescription(std::string HelpLongDescription) {
    HelpLongDescription_ = std::move(HelpLongDescription);
    return *this;
  }

  template <typename T> command_args_parser &AddOption(const std::string &Name, char ShortName,
    const std::string &Description) {
    Options_.Insert(Name, ShortName, command_args_internal::GetValueTypeFlag_<T>(), Description);
    OptionShortNameToName_.Insert(ShortName, Name);
    return *this;
  }

  command_args_parser &RemoveOption(const std::string &Name) {
    auto Iter = Options_.Find(Name);
    if (Iter != Options_.End()) {
      const option_data &OptionData = Iter->Value();
      OptionShortNameToName_.Erase(OptionData.ShortName);
      Options_.Erase(Iter);
    }
    return *this;
  }

  command_args_parser &ClearOptions() {
    Options_.Clear();
    OptionShortNameToName_.Clear();
    return *this;
  }

  command_args Parse(const ovk::array<std::string> &RawArguments) const;
  ovk::optional<command_args> Parse(const ovk::array<std::string> &RawArguments,
    captured_command_args_error &Error) const;

private:

  using value_type_flag = command_args_internal::value_type_flag;

  struct option_data {
    char ShortName;
    value_type_flag ValueTypeFlag;
    std::string Description;
    option_data(char ShortName_, value_type_flag ValueTypeFlag_, std::string Description_):
      ShortName(ShortName_),
      ValueTypeFlag(ValueTypeFlag_),
      Description(std::move(Description_))
    {}
  };

  bool WriteCondition_ = true;
  std::string HelpUsage_;
  std::string HelpDescription_;
  std::string HelpLongDescription_;
  ovk::map<std::string,option_data> Options_;
  ovk::map<char,std::string> OptionShortNameToName_;

};

}

#endif
