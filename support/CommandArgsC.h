// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_COMMAND_ARGS_C_H_LOADED
#define OVK_SUPPORT_COMMAND_ARGS_C_H_LOADED

#include <support/CommandArgs.h>

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  SUPPORT_COMMAND_ARGS_VALUE_TYPE_BOOL = 0,
  SUPPORT_COMMAND_ARGS_VALUE_TYPE_INT,
  SUPPORT_COMMAND_ARGS_VALUE_TYPE_LONG_LONG,
  SUPPORT_COMMAND_ARGS_VALUE_TYPE_DOUBLE,
  SUPPORT_COMMAND_ARGS_VALUE_TYPE_STRING,
  SUPPORT_COMMAND_ARGS_VALUE_TYPE_CHAR
} support_command_args_value_type;

static inline bool support_ValidCommandArgsValueType(support_command_args_value_type ValueType) {

  switch (ValueType) {
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_BOOL:
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_INT:
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_LONG_LONG:
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_DOUBLE:
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_STRING:
  case SUPPORT_COMMAND_ARGS_VALUE_TYPE_CHAR:
    return true;
  default:
    return false;
  }

}

struct support_command_args;
typedef struct support_command_args support_command_args;

struct support_command_args_parser;
typedef struct support_command_args_parser support_command_args_parser;

void support_CreateCommandArgs(support_command_args **CommandArgs, bool WriteCondition);
void support_DestroyCommandArgs(support_command_args **CommandArgs);

void support_FillCommandArgs(support_command_args *CommandArgs, int NumOptions, char **OptionNames,
  char **OptionValues, int NumArguments, char **Arguments);

bool support_CommandOptionIsPresent(const support_command_args *CommandArgs, const char *Name);

void support_GetCommandOption(const support_command_args *CommandArgs, const char *Name,
  support_command_args_value_type ValueType, void *Value, support_command_args_error *Error);
void support_GetCommandOptionIfPresent(const support_command_args *CommandArgs, const char *Name,
  support_command_args_value_type ValueType, void *Value, support_command_args_error *Error);

int support_GetCommandArgumentCount(const support_command_args *CommandArgs);
void support_GetCommandArgument(const support_command_args *CommandArgs, int iArgument,
  support_command_args_value_type ValueType, void *Argument, support_command_args_error *Error);

void support_CreateCommandArgsParser(support_command_args_parser **CommandArgsParser, bool
  WriteCondition);
void support_DestroyCommandArgsParser(support_command_args_parser **CommandArgsParser);

void support_SetCommandArgsParserHelpUsage(support_command_args_parser *CommandArgsParser, const
  char *HelpUsage);
void support_SetCommandArgsParserHelpDescription(support_command_args_parser *CommandArgsParser,
  const char *HelpDescription);
void support_SetCommandArgsParserHelpLongDescription(support_command_args_parser *CommandArgsParser,
  const char *HelpLongDescription);
void support_AddCommandArgsParserOption(support_command_args_parser *CommandArgsParser, const
  char *Name, char ShortName, support_command_args_value_type ValueType, const char *Description);
void support_RemoveCommandArgsParserOption(support_command_args_parser *CommandArgsParser, const
  char *Name);
void support_ClearCommandArgsParserOptions(support_command_args_parser *CommandArgsParser);

void support_ParseCommandArgs(const support_command_args_parser *CommandArgsParser, int
  NumRawArguments, char **RawArguments, support_command_args **CommandArgs,
  support_command_args_error *Error);

#ifdef __cplusplus
}
#endif

#endif
