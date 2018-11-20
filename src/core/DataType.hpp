// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DATA_TYPE_HPP_INCLUDED
#define OVK_CORE_DATA_TYPE_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/DataTypeBase.h>
#include <ovk/core/Global.hpp>

#include <mpi.h>

namespace ovk {

inline int DataTypeSize(data_type DataType) { return ovkDataTypeSize(ovk_data_type(DataType)); }
inline bool DataTypeIsIntegral(data_type DataType) { return ovkDataTypeIsIntegral(ovk_data_type(DataType)); }
inline bool DataTypeIsFloatingPoint(data_type DataType) { return ovkDataTypeIsFloatingPoint(ovk_data_type(DataType)); }
inline bool DataTypeIsSigned(data_type DataType) { return ovkDataTypeIsSigned(ovk_data_type(DataType)); }
inline bool DataTypeIsUnsigned(data_type DataType) { return ovkDataTypeIsUnsigned(ovk_data_type(DataType)); }
inline MPI_Datatype DataTypeToMPI(data_type DataType) { return ovkDataTypeToMPI(ovk_data_type(DataType)); }

namespace core {

template <typename T> constexpr data_type GetDataType();
template <> constexpr data_type GetDataType<bool>() { return data_type::BOOL; }
template <> constexpr data_type GetDataType<unsigned char>() { return data_type::BYTE; }
template <> constexpr data_type GetDataType<int>() { return data_type::INT; }
template <> constexpr data_type GetDataType<long>() { return data_type::LONG; }
template <> constexpr data_type GetDataType<long long>() { return data_type::LONG_LONG; }
template <> constexpr data_type GetDataType<unsigned int>() { return data_type::UNSIGNED_INT; }
template <> constexpr data_type GetDataType<unsigned long>() { return data_type::UNSIGNED_LONG; }
template <> constexpr data_type GetDataType<unsigned long long>() { return data_type::UNSIGNED_LONG_LONG; }
template <> constexpr data_type GetDataType<float>() { return data_type::FLOAT; }
template <> constexpr data_type GetDataType<double>() { return data_type::DOUBLE; }

// Tried to make this constexpr, but apparently MPI types aren't always compile-time constants
template <typename T> inline MPI_Datatype GetMPIDataType();
template <> inline MPI_Datatype GetMPIDataType<char>() { return MPI_SIGNED_CHAR; }
template <> inline MPI_Datatype GetMPIDataType<short>() { return MPI_SHORT; }
template <> inline MPI_Datatype GetMPIDataType<int>() { return MPI_INT; }
template <> inline MPI_Datatype GetMPIDataType<long>() { return MPI_LONG; }
template <> inline MPI_Datatype GetMPIDataType<long long>() { return MPI_LONG_LONG; }
template <> inline MPI_Datatype GetMPIDataType<unsigned char>() { return MPI_UNSIGNED_CHAR; }
template <> inline MPI_Datatype GetMPIDataType<unsigned short>() { return MPI_UNSIGNED_SHORT; }
template <> inline MPI_Datatype GetMPIDataType<unsigned int>() { return MPI_UNSIGNED; }
template <> inline MPI_Datatype GetMPIDataType<unsigned long>() { return MPI_UNSIGNED_LONG; }
template <> inline MPI_Datatype GetMPIDataType<unsigned long long>() { return MPI_UNSIGNED_LONG_LONG; }
template <> inline MPI_Datatype GetMPIDataType<float>() { return MPI_FLOAT; }
template <> inline MPI_Datatype GetMPIDataType<double>() { return MPI_DOUBLE; }

}

}

#endif
