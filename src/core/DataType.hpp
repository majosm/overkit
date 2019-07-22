// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DATA_TYPE_HPP_INCLUDED
#define OVK_CORE_DATA_TYPE_HPP_INCLUDED

#include <ovk/core/DataType.h>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>

#include <mpi.h>

#include <cstdint>

namespace ovk {

enum class data_type : int {
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

template <typename T, typename=void> struct data_type_traits {
//   static constexpr data_type Type = ???;
//   using mpi_convert_type = ???;
// Tried to make this constexpr, but apparently MPI types aren't always compile-time constants
//   static const MPI_Datatype MPIType = ???;
};

namespace is_supported_data_type_internal {
template <typename T> using maybe_int = int;
template <typename T> constexpr std::true_type Test(maybe_int<decltype(data_type_traits<T>::Type)>)
  {
  return {};
}
template <typename T> constexpr std::false_type Test(...) { return {}; }
}
template <typename T> constexpr bool IsSupportedDataType() {
  return decltype(is_supported_data_type_internal::Test<T>(0))::value;
}

template <typename T, OVK_FUNCTION_REQUIRES(IsSupportedDataType<T>())> constexpr data_type
  GetDataType() {
  return data_type_traits<T>::Type;
}

template <typename T, OVK_FUNCTION_REQUIRES(IsSupportedDataType<T>())> MPI_Datatype GetMPIDataType()
  {
  return data_type_traits<T>::MPIType;
}

template <typename T> using mpi_compatible_type = typename data_type_traits<T>::mpi_convert_type;

template <> struct data_type_traits<bool> {
  static constexpr data_type Type = data_type::BOOL;
  using mpi_convert_type = unsigned char;
  static const MPI_Datatype MPIType = MPI_UNSIGNED_CHAR;
};

template <> struct data_type_traits<byte> {
  static constexpr data_type Type = data_type::BYTE;
  using mpi_convert_type = byte;
  static const MPI_Datatype MPIType = MPI_UNSIGNED_CHAR;
};

template <> struct data_type_traits<int> {
  static constexpr data_type Type = data_type::INT;
  using mpi_convert_type = int;
  static const MPI_Datatype MPIType = MPI_INT;
};

template <> struct data_type_traits<long> {
  static constexpr data_type Type = data_type::LONG;
  using mpi_convert_type = long;
  static const MPI_Datatype MPIType = MPI_LONG;
};

template <> struct data_type_traits<long long> {
  static constexpr data_type Type = data_type::LONG_LONG;
  using mpi_convert_type = long long;
  static const MPI_Datatype MPIType = MPI_LONG_LONG;
};

template <> struct data_type_traits<unsigned int> {
  static constexpr data_type Type = data_type::UNSIGNED_INT;
  using mpi_convert_type = unsigned int;
  static const MPI_Datatype MPIType = MPI_UNSIGNED;
};

template <> struct data_type_traits<unsigned long> {
  static constexpr data_type Type = data_type::UNSIGNED_LONG;
  using mpi_convert_type = unsigned long;
  static const MPI_Datatype MPIType = MPI_UNSIGNED_LONG;
};

template <> struct data_type_traits<unsigned long long> {
  static constexpr data_type Type = data_type::UNSIGNED_LONG_LONG;
  using mpi_convert_type = unsigned long long;
  static const MPI_Datatype MPIType = MPI_UNSIGNED_LONG_LONG;
};

template <> struct data_type_traits<float> {
  static constexpr data_type Type = data_type::FLOAT;
  using mpi_convert_type = float;
  static const MPI_Datatype MPIType = MPI_FLOAT;
};

template <> struct data_type_traits<double> {
  static constexpr data_type Type = data_type::DOUBLE;
  using mpi_convert_type = double;
  static const MPI_Datatype MPIType = MPI_DOUBLE;
};

}

}

#endif
