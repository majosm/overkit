// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXAMPLES_COMMON_H_LOADED
#define OVK_EXAMPLES_COMMON_H_LOADED

#include <support/CommandArgsC.h>
#include <support/Constants.h>
#include <support/Decomp.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EXAMPLES_PI SUPPORT_PI

#define examples_command_args_error support_command_args_error
#define EXAMPLES_COMMAND_ARGS_ERROR_NONE SUPPORT_COMMAND_ARGS_ERROR_NONE
#define EXAMPLES_COMMAND_ARGS_ERROR_INVALID_ARGUMENT SUPPORT_COMMAND_ARGS_ERROR_INVALID_ARGUMENT

#define examples_command_args_value_type support_command_args_value_type
#define EXAMPLES_COMMAND_ARGS_VALUE_TYPE_BOOL SUPPORT_COMMAND_ARGS_VALUE_TYPE_BOOL
#define EXAMPLES_COMMAND_ARGS_VALUE_TYPE_INT SUPPORT_COMMAND_ARGS_VALUE_TYPE_INT
#define EXAMPLES_COMMAND_ARGS_VALUE_TYPE_LONG_LONG SUPPORT_COMMAND_ARGS_VALUE_TYPE_LONG_LONG
#define EXAMPLES_COMMAND_ARGS_VALUE_TYPE_DOUBLE SUPPORT_COMMAND_ARGS_VALUE_TYPE_DOUBLE
#define EXAMPLES_COMMAND_ARGS_VALUE_TYPE_STRING SUPPORT_COMMAND_ARGS_VALUE_TYPE_STRING
#define EXAMPLES_COMMAND_ARGS_VALUE_TYPE_CHAR SUPPORT_COMMAND_ARGS_VALUE_TYPE_CHAR

#define examples_command_args support_command_args
#define examples_command_args_parser support_command_args_parser
#define examples_CreateCommandArgs support_CreateCommandArgs
#define examples_DestroyCommandArgs support_DestroyCommandArgs
#define examples_CommandOptionIsPresent support_CommandOptionIsPresent
#define examples_GetCommandOption support_GetCommandOption
#define examples_GetCommandOptionIfPresent support_GetCommandOptionIfPresent
#define examples_GetCommandArgumentCount support_GetCommandArgumentCount
#define examples_GetCommandArgument support_GetCommandArgument

#define examples_CreateCommandArgsParser support_CreateCommandArgsParser
#define examples_DestroyCommandArgsParser support_DestroyCommandArgsParser
#define examples_SetCommandArgsParserHelpUsage support_SetCommandArgsParserHelpUsage
#define examples_SetCommandArgsParserHelpDescription support_SetCommandArgsParserHelpDescription
#define examples_SetCommandArgsParserHelpLongDescription support_SetCommandArgsParserHelpLongDescription
#define examples_AddCommandArgsParserOption support_AddCommandArgsParserOption
#define examples_RemoveCommandArgsParserOption support_RemoveCommandArgsParserOption
#define examples_ClearCommandArgsParserOptions support_ClearCommandArgsParserOptions
#define examples_ParseCommandArgs support_ParseCommandArgs

void examples_DecomposeDomain(int NumGrids, const long long *NumPointsPerGrid, int NumProcs, int
  *GridProcRanges) {

  support_DecomposeDomain(NumGrids, NumPointsPerGrid, NumProcs, GridProcRanges);

}

void examples_CreateCartesianDecompDims(int Size, int NumDims, int *Dims) {

  support_CreateCartesianDecompDims(Size, NumDims, Dims);

}

void examples_CartesianDecomp(int NumDims, const int *Size, MPI_Comm CartComm, int *LocalRange) {

  int Zero[] = {0,0,0};

  support_CartesianDecomp(NumDims, Zero, Size, CartComm, LocalRange, LocalRange+3);

}

#ifdef __cplusplus
}
#endif

#endif
