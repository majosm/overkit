// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "support/CommandArgsC.h"

#include "support/CommandArgs.hpp"

#include <ovk/core/Array.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Optional.hpp>

#include <cstring>
#include <string>
#include <utility>

#ifdef __cplusplus
extern "C" {
#endif

void support_CreateCommandArgs(support_command_args **CommandArgs, bool WriteCondition) {

  OVK_DEBUG_ASSERT(CommandArgs, "Invalid command args pointer.");

  auto CommandArgsCPPPtr = new support::command_args(WriteCondition);
  *CommandArgs = reinterpret_cast<support_command_args *>(CommandArgsCPPPtr);

}

void support_DestroyCommandArgs(support_command_args **CommandArgs) {

  OVK_DEBUG_ASSERT(CommandArgs, "Invalid command args pointer.");
  OVK_DEBUG_ASSERT(*CommandArgs, "Invalid command args pointer.");

  auto CommandArgsCPPPtr = reinterpret_cast<support::command_args *>(*CommandArgs);

  delete CommandArgsCPPPtr;

  *CommandArgs = nullptr;

}

void support_FillCommandArgs(support_command_args *CommandArgs, int NumOptions, char **OptionNames,
  char **OptionValues, int NumArguments, char **Arguments) {

  OVK_DEBUG_ASSERT(CommandArgs, "Invalid command args pointer.");
  OVK_DEBUG_ASSERT(NumOptions >= 0, "Invalid option count.");
  OVK_DEBUG_ASSERT(OptionNames || NumOptions == 0, "Invalid option names pointer.");
  OVK_DEBUG_ASSERT(OptionValues || NumOptions == 0, "Invalid option values pointer.");
  OVK_DEBUG_ASSERT(NumArguments >= 0, "Invalid argument count.");
  OVK_DEBUG_ASSERT(Arguments || NumArguments == 0, "Invalid arguments pointer.");

  ovk::map<std::string,std::string> OptionsCPP;
  OptionsCPP.Reserve(NumOptions);

  for (int iOption = 0; iOption < NumOptions; ++iOption) {
    OptionsCPP.Insert(OptionNames[iOption], OptionValues[iOption]);
  }

  ovk::array<std::string> ArgumentsCPP({NumArguments}, Arguments);

  auto &CommandArgsCPP = *reinterpret_cast<support::command_args *>(CommandArgs);
  CommandArgsCPP = support::command_args(std::move(OptionsCPP), std::move(ArgumentsCPP),
    CommandArgsCPP.WriteCondition());

}


bool support_CommandOptionIsPresent(const support_command_args *CommandArgs, const char *Name) {

  OVK_DEBUG_ASSERT(CommandArgs, "Invalid command args pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &CommandArgsCPP = *reinterpret_cast<const support::command_args *>(CommandArgs);
  return CommandArgsCPP.OptionIsPresent(Name);

}

void support_GetCommandOption(const support_command_args *CommandArgs, const char *Name,
  support_command_args_value_type ValueType, void *Value, support_command_args_error *Error) {

  OVK_DEBUG_ASSERT(CommandArgs, "Invalid command args pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");
  OVK_DEBUG_ASSERT(support_ValidCommandArgsValueType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Value, "Invalid value pointer.");
  OVK_DEBUG_ASSERT(Error, "Invalid error pointer.");

  auto &CommandArgsCPP = *reinterpret_cast<const support::command_args *>(CommandArgs);

  std::string NameCPP = Name;

  support::captured_command_args_error ErrorCPP;

  switch (ValueType) {
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_BOOL: {
    auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<bool>(NameCPP, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<bool *>(Value) = *MaybeValueCPP;
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_INT: {
    auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<int>(NameCPP, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<int *>(Value) = *MaybeValueCPP;
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_LONG_LONG: {
    auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<long long>(NameCPP, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<long long *>(Value) = *MaybeValueCPP;
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_DOUBLE: {
    auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<double>(NameCPP, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<double *>(Value) = *MaybeValueCPP;
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_STRING: {
    auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<std::string>(NameCPP, ErrorCPP);
    if (!ErrorCPP) {
      std::strcpy(static_cast<char *>(Value), (*MaybeValueCPP).c_str());
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_CHAR: {
    auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<char>(NameCPP, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<char *>(Value) = *MaybeValueCPP;
    }
    break;
  }
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  *Error = support_command_args_error(ErrorCPP.Code());

}

void support_GetCommandOptionIfPresent(const support_command_args *CommandArgs, const char *Name,
  support_command_args_value_type ValueType, void *Value, support_command_args_error *Error) {

  OVK_DEBUG_ASSERT(CommandArgs, "Invalid command args pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");
  OVK_DEBUG_ASSERT(support_ValidCommandArgsValueType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Value, "Invalid value pointer.");
  OVK_DEBUG_ASSERT(Error, "Invalid error pointer.");

  auto &CommandArgsCPP = *reinterpret_cast<const support::command_args *>(CommandArgs);

  std::string NameCPP = Name;

  support::captured_command_args_error ErrorCPP;

  if (CommandArgsCPP.OptionIsPresent(NameCPP)) {
    switch (ValueType) {
    case SUPPORT_COMMAND_ARGS_VALUE_TYPE_BOOL: {
      auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<bool>(NameCPP, ErrorCPP);
      if (!ErrorCPP) {
        *static_cast<bool *>(Value) = *MaybeValueCPP;
      }
      break;
    }
    case SUPPORT_COMMAND_ARGS_VALUE_TYPE_INT: {
      auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<int>(NameCPP, ErrorCPP);
      if (!ErrorCPP) {
        *static_cast<int *>(Value) = *MaybeValueCPP;
      }
      break;
    }
    case SUPPORT_COMMAND_ARGS_VALUE_TYPE_LONG_LONG: {
      auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<long long>(NameCPP, ErrorCPP);
      if (!ErrorCPP) {
        *static_cast<long long *>(Value) = *MaybeValueCPP;
      }
      break;
    }
    case SUPPORT_COMMAND_ARGS_VALUE_TYPE_DOUBLE: {
      auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<double>(NameCPP, ErrorCPP);
      if (!ErrorCPP) {
        *static_cast<double *>(Value) = *MaybeValueCPP;
      }
      break;
    }
    case SUPPORT_COMMAND_ARGS_VALUE_TYPE_STRING: {
      auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<std::string>(NameCPP, ErrorCPP);
      if (!ErrorCPP) {
        std::strcpy(static_cast<char *>(Value), (*MaybeValueCPP).c_str());
      }
      break;
    }
    case SUPPORT_COMMAND_ARGS_VALUE_TYPE_CHAR: {
      auto MaybeValueCPP = CommandArgsCPP.GetOptionValue<char>(NameCPP, ErrorCPP);
      if (!ErrorCPP) {
        *static_cast<char *>(Value) = *MaybeValueCPP;
      }
      break;
    }
    default:
      OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
      break;
    }
  }

  *Error = support_command_args_error(ErrorCPP.Code());

}

int support_GetCommandArgumentCount(const support_command_args *CommandArgs) {

  OVK_DEBUG_ASSERT(CommandArgs, "Invalid command args pointer.");

  auto &CommandArgsCPP = *reinterpret_cast<const support::command_args *>(CommandArgs);
  return CommandArgsCPP.ArgumentCount();

}

void support_GetCommandArgument(const support_command_args *CommandArgs, int iArgument,
  support_command_args_value_type ValueType, void *Argument, support_command_args_error *Error) {

  OVK_DEBUG_ASSERT(CommandArgs, "Invalid command args pointer.");
  OVK_DEBUG_ASSERT(iArgument >= 0, "Invalid argument index.");
  OVK_DEBUG_ASSERT(support_ValidCommandArgsValueType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Argument, "Invalid argument pointer.");
  OVK_DEBUG_ASSERT(Error, "Invalid error pointer.");

  auto &CommandArgsCPP = *reinterpret_cast<const support::command_args *>(CommandArgs);

  OVK_DEBUG_ASSERT(iArgument < CommandArgsCPP.ArgumentCount(), "Invalid argument index.");

  support::captured_command_args_error ErrorCPP;

  switch (ValueType) {
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_BOOL: {
    auto MaybeArgumentCPP = CommandArgsCPP.GetArgument<bool>(iArgument, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<bool *>(Argument) = *MaybeArgumentCPP;
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_INT: {
    auto MaybeArgumentCPP = CommandArgsCPP.GetArgument<int>(iArgument, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<int *>(Argument) = *MaybeArgumentCPP;
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_LONG_LONG: {
    auto MaybeArgumentCPP = CommandArgsCPP.GetArgument<long long>(iArgument, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<long long *>(Argument) = *MaybeArgumentCPP;
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_DOUBLE: {
    auto MaybeArgumentCPP = CommandArgsCPP.GetArgument<double>(iArgument, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<double *>(Argument) = *MaybeArgumentCPP;
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_STRING: {
    auto MaybeArgumentCPP = CommandArgsCPP.GetArgument<std::string>(iArgument, ErrorCPP);
    if (!ErrorCPP) {
      std::strcpy(static_cast<char *>(Argument), (*MaybeArgumentCPP).c_str());
    }
    break;
  }
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_CHAR: {
    auto MaybeArgumentCPP = CommandArgsCPP.GetArgument<char>(iArgument, ErrorCPP);
    if (!ErrorCPP) {
      *static_cast<char *>(Argument) = *MaybeArgumentCPP;
    }
    break;
  }
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  *Error = support_command_args_error(ErrorCPP.Code());

}

void support_CreateCommandArgsParser(support_command_args_parser **CommandArgsParser, bool
  WriteCondition) {

  OVK_DEBUG_ASSERT(CommandArgsParser, "Invalid command args parser pointer.");

  auto CommandArgsParserCPPPtr = new support::command_args_parser(WriteCondition);
  *CommandArgsParser = reinterpret_cast<support_command_args_parser *>(CommandArgsParserCPPPtr);

}

void support_DestroyCommandArgsParser(support_command_args_parser **CommandArgsParser) {

  OVK_DEBUG_ASSERT(CommandArgsParser, "Invalid command args parser pointer.");
  OVK_DEBUG_ASSERT(*CommandArgsParser, "Invalid command args parser pointer.");

  auto CommandArgsParserCPPPtr = reinterpret_cast<support::command_args_parser *>(
    *CommandArgsParser);

  delete CommandArgsParserCPPPtr;

  *CommandArgsParser = nullptr;

}

void support_SetCommandArgsParserHelpUsage(support_command_args_parser *CommandArgsParser, const
  char *HelpUsage) {

  OVK_DEBUG_ASSERT(CommandArgsParser, "Invalid command args parser pointer.");
  OVK_DEBUG_ASSERT(HelpUsage, "Invalid help usage pointer.");

  auto &CommandArgsParserCPP = *reinterpret_cast<support::command_args_parser *>(CommandArgsParser);
  CommandArgsParserCPP.SetHelpUsage(HelpUsage);

}

void support_SetCommandArgsParserHelpDescription(support_command_args_parser *CommandArgsParser,
  const char *HelpDescription) {

  OVK_DEBUG_ASSERT(CommandArgsParser, "Invalid command args parser pointer.");
  OVK_DEBUG_ASSERT(HelpDescription, "Invalid help description pointer.");

  auto &CommandArgsParserCPP = *reinterpret_cast<support::command_args_parser *>(CommandArgsParser);
  CommandArgsParserCPP.SetHelpDescription(HelpDescription);

}

void support_SetCommandArgsParserHelpLongDescription(support_command_args_parser *CommandArgsParser,
  const char *HelpLongDescription) {

  OVK_DEBUG_ASSERT(CommandArgsParser, "Invalid command args parser pointer.");
  OVK_DEBUG_ASSERT(HelpLongDescription, "Invalid help long description pointer.");

  auto &CommandArgsParserCPP = *reinterpret_cast<support::command_args_parser *>(CommandArgsParser);
  CommandArgsParserCPP.SetHelpLongDescription(HelpLongDescription);

}

void support_AddCommandArgsParserOption(support_command_args_parser *CommandArgsParser, const
  char *Name, char ShortName, support_command_args_value_type ValueType, const char *Description) {

  OVK_DEBUG_ASSERT(CommandArgsParser, "Invalid command args parser pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");
  OVK_DEBUG_ASSERT(support_ValidCommandArgsValueType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Description, "Invalid description pointer.");

  auto &CommandArgsParserCPP = *reinterpret_cast<support::command_args_parser *>(CommandArgsParser);

  switch (ValueType) {
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_BOOL:
    CommandArgsParserCPP.AddOption<bool>(Name, ShortName, Description);
    break;
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_INT:
    CommandArgsParserCPP.AddOption<int>(Name, ShortName, Description);
    break;
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_LONG_LONG:
    CommandArgsParserCPP.AddOption<long long>(Name, ShortName, Description);
    break;
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_DOUBLE:
    CommandArgsParserCPP.AddOption<double>(Name, ShortName, Description);
    break;
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_STRING:
    CommandArgsParserCPP.AddOption<std::string>(Name, ShortName, Description);
    break;
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_CHAR:
    CommandArgsParserCPP.AddOption<char>(Name, ShortName, Description);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

}

void support_RemoveCommandArgsParserOption(support_command_args_parser *CommandArgsParser, const
  char *Name) {

  OVK_DEBUG_ASSERT(CommandArgsParser, "Invalid command args parser pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &CommandArgsParserCPP = *reinterpret_cast<support::command_args_parser *>(CommandArgsParser);
  CommandArgsParserCPP.RemoveOption(Name);

}

void support_ClearCommandArgsParserOptions(support_command_args_parser *CommandArgsParser) {

  OVK_DEBUG_ASSERT(CommandArgsParser, "Invalid command args parser pointer.");

  auto &CommandArgsParserCPP = *reinterpret_cast<support::command_args_parser *>(CommandArgsParser);
  CommandArgsParserCPP.ClearOptions();

}

void support_ParseCommandArgs(const support_command_args_parser *CommandArgsParser, int
  NumRawArguments, char **RawArguments, support_command_args **CommandArgs,
  support_command_args_error *Error) {

  OVK_DEBUG_ASSERT(CommandArgsParser, "Invalid command args parser pointer.");
  OVK_DEBUG_ASSERT(NumRawArguments >= 0, "Invalid raw argument count.");
  OVK_DEBUG_ASSERT(RawArguments || NumRawArguments == 0, "Invalid raw arguments pointer.");
  OVK_DEBUG_ASSERT(CommandArgs, "Invalid command_args pointer.");
  OVK_DEBUG_ASSERT(Error, "Invalid error pointer.");

  auto &CommandArgsParserCPP = *reinterpret_cast<const support::command_args_parser *>(
    CommandArgsParser);

  ovk::array<std::string> RawArgumentsCPP({NumRawArguments}, RawArguments);

  support::captured_command_args_error ErrorCPP;
  auto MaybeCommandArgsCPP = CommandArgsParserCPP.Parse(RawArgumentsCPP, ErrorCPP);

  if (!ErrorCPP) {
    auto CommandArgsCPPPtr = new support::command_args(MaybeCommandArgsCPP.Release());
    *CommandArgs = reinterpret_cast<support_command_args *>(CommandArgsCPPPtr);
  } else {
    *CommandArgs = nullptr;
  }

  *Error = support_command_args_error(ErrorCPP.Code());

}

#ifdef __cplusplus
}
#endif
