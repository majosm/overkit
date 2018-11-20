// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DATA_TYPE_BASE_H_INCLUDED
#define OVK_CORE_DATA_TYPE_BASE_H_INCLUDED

#include <ovk/core/ConstantsBase.h>

#include <mpi.h>

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline int ovkDataTypeSize(ovk_data_type DataType);
static inline bool ovkDataTypeIsIntegral(ovk_data_type DataType);
static inline bool ovkDataTypeIsFloatingPoint(ovk_data_type DataType);
static inline bool ovkDataTypeIsSigned(ovk_data_type DataType);
static inline bool ovkDataTypeIsUnsigned(ovk_data_type DataType);
static inline MPI_Datatype ovkDataTypeToMPI(ovk_data_type DataType);

#ifdef __cplusplus
}
#endif

#include <ovk/core/DataTypeBase.inl>

#endif
