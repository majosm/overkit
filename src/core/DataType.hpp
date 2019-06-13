// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DATA_TYPE_HPP_INCLUDED
#define OVK_CORE_DATA_TYPE_HPP_INCLUDED

#include <ovk/core/DataType.h>
#include <ovk/core/Global.hpp>

#include <mpi.h>

#include <cstdint>

namespace ovk {

enum class data_type {
  BOOL = OVK_BOOL,
  BYTE = OVK_BYTE,
  INT = OVK_INT,
  LONG = OVK_LONG,
  LONG_LONG = OVK_LONG_LONG,
#ifdef INT32_MAX
  INT32 = OVK_INT32,
#endif
#ifdef INT64_MAX
  INT64 = OVK_INT64,
#endif
  UNSIGNED_INT = OVK_UNSIGNED_INT,
  UNSIGNED_LONG = OVK_UNSIGNED_LONG,
  UNSIGNED_LONG_LONG = OVK_UNSIGNED_LONG_LONG,
#ifdef UINT32_MAX
  UNSIGNED_INT32 = OVK_UNSIGNED_INT32,
#endif
#ifdef UINT64_MAX
  UNSIGNED_INT64 = OVK_UNSIGNED_INT64,
#endif
  FLOAT = OVK_FLOAT,
  DOUBLE = OVK_DOUBLE
};

inline bool ValidDataType(data_type DataType) {
  return ovkValidDataType(ovk_data_type(DataType));
}
inline int DataTypeSize(data_type DataType) {
  return ovkDataTypeSize(ovk_data_type(DataType));
}
inline bool DataTypeIsIntegral(data_type DataType) {
  return ovkDataTypeIsIntegral(ovk_data_type(DataType));
}
inline bool DataTypeIsFloatingPoint(data_type DataType) {
  return ovkDataTypeIsFloatingPoint(ovk_data_type(DataType));
}
inline bool DataTypeIsSigned(data_type DataType) {
  return ovkDataTypeIsSigned(ovk_data_type(DataType));
}
inline bool DataTypeIsUnsigned(data_type DataType) {
  return ovkDataTypeIsUnsigned(ovk_data_type(DataType));
}
inline MPI_Datatype DataTypeToMPI(data_type DataType) {
  return ovkDataTypeToMPI(ovk_data_type(DataType));
}

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

// Use unsigned char in place of bool for MPI sends/recvs
namespace data_type_internal {
template <typename T> struct mpi_compatible_type_helper { using type = T; };
template <> struct mpi_compatible_type_helper<bool> { using type = unsigned char; };
}
template <typename T> using mpi_compatible_type = typename data_type_internal::
  mpi_compatible_type_helper<T>::type;

}

}

#endif
