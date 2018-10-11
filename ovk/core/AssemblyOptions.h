// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ASSEMBLY_OPTIONS_INCLUDED
#define OVK_CORE_ASSEMBLY_OPTIONS_INCLUDED

#include "ovk/core/ovkAssemblyOptions.h"

#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Logger.h"
#include "ovk/core/OrderedMap.h"

#ifdef __cplusplus
extern "C" {
#endif

void PRIVATE(CreateAssemblyOptions)(ovk_assembly_options **Options, int NumDims, int NumGrids,
  const t_ordered_map *GridIDs, t_logger *Logger, t_error_handler *ErrorHandler);
#define CreateAssemblyOptions(...) PRIVATE(CreateAssemblyOptions)(__VA_ARGS__)
void PRIVATE(DestroyAssemblyOptions)(ovk_assembly_options **Options);
#define DestroyAssemblyOptions(...) PRIVATE(DestroyAssemblyOptions)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
