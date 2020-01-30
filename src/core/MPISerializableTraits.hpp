// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MPI_SERIALIZABLE_TRAITS_HPP_INCLUDED
#define OVK_CORE_MPI_SERIALIZABLE_TRAITS_HPP_INCLUDED

#include <ovk/core/DataTypeOps.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Handle.hpp>
#include <ovk/core/Interval.hpp>

#include <mpi.h>

#include <type_traits>

namespace ovk {
namespace core {

template <typename T, typename=void> struct mpi_serializable_traits {
/*
  using packed_type = ???;
  static handle<MPI_Datatype> CreateMPIType() { return ???; }
  static packed_type Pack(const T &Value) { return ???; }
  static T Unpack(const packed_type &PackedValue) { return ???; }
*/
};

namespace is_mpi_serializable_internal {
template <typename T> using maybe_int = int;
template <typename T> constexpr std::true_type Test(maybe_int<typename
  mpi_serializable_traits<T>::packed_type>) {
  return {};
}
template <typename T> constexpr std::false_type Test(...) { return {}; }
}
template <typename T> constexpr bool IsMPISerializable() {
  return decltype(is_mpi_serializable_internal::Test<T>(0))::value;
}

template <> struct mpi_serializable_traits<bool> {
  using packed_type = bool;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_C_BOOL, [](MPI_Datatype) {}};
  }
  static packed_type Pack(bool Value) { return Value; }
  static bool Unpack(packed_type PackedValue) { return PackedValue; }
};

template <> struct mpi_serializable_traits<byte> {
  using packed_type = byte;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_BYTE, [](MPI_Datatype) {}};
  }
  static packed_type Pack(byte Value) { return Value; }
  static byte Unpack(packed_type PackedValue) { return PackedValue; }
};

template <> struct mpi_serializable_traits<int> {
  using packed_type = int;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_INT, [](MPI_Datatype) {}};
  }
  static packed_type Pack(int Value) { return Value; }
  static int Unpack(packed_type PackedValue) { return PackedValue; }
};

template <> struct mpi_serializable_traits<long> {
  using packed_type = long;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_LONG, [](MPI_Datatype) {}};
  }
  static packed_type Pack(long Value) { return Value; }
  static long Unpack(packed_type PackedValue) { return PackedValue; }
};

template <> struct mpi_serializable_traits<long long> {
  using packed_type = long long;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_LONG_LONG, [](MPI_Datatype) {}};
  }
  static packed_type Pack(long long Value) { return Value; }
  static long long Unpack(packed_type PackedValue) { return PackedValue; }
};

template <> struct mpi_serializable_traits<unsigned int> {
  using packed_type = unsigned int;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_UNSIGNED, [](MPI_Datatype) {}};
  }
  static packed_type Pack(unsigned int Value) { return Value; }
  static unsigned int Unpack(packed_type PackedValue) { return PackedValue; }
};

template <> struct mpi_serializable_traits<unsigned long> {
  using packed_type = unsigned long;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_UNSIGNED_LONG, [](MPI_Datatype) {}};
  }
  static packed_type Pack(unsigned long Value) { return Value; }
  static unsigned long Unpack(packed_type PackedValue) { return PackedValue; }
};

template <> struct mpi_serializable_traits<unsigned long long> {
  using packed_type = unsigned long long;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_UNSIGNED_LONG_LONG, [](MPI_Datatype) {}};
  }
  static packed_type Pack(unsigned long long Value) { return Value; }
  static unsigned long long Unpack(packed_type PackedValue) { return PackedValue; }
};

template <> struct mpi_serializable_traits<float> {
  using packed_type = float;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_FLOAT, [](MPI_Datatype) {}};
  }
  static packed_type Pack(float Value) { return Value; }
  static float Unpack(packed_type PackedValue) { return PackedValue; }
};

template <> struct mpi_serializable_traits<double> {
  using packed_type = double;
  static handle<MPI_Datatype> CreateMPIType() {
    return {MPI_DOUBLE, [](MPI_Datatype) {}};
  }
  static packed_type Pack(double Value) { return Value; }
  static double Unpack(packed_type PackedValue) { return PackedValue; }
};

template <typename T, int N> struct mpi_serializable_traits<elem<T,N>, OVK_SPECIALIZATION_REQUIRES(
  IsMPISerializable<T>())> {
  using packed_type = elem<T,N>;
  static handle<MPI_Datatype> CreateMPIType() {
    auto ValueMPIType = mpi_serializable_traits<T>::CreateMPIType();
    return CreateMPIContiguousType(N, ValueMPIType);
  }
  static packed_type Pack(const elem<T,N> &Value) { return Value; }
  static elem<T,N> Unpack(const packed_type &PackedValue) { return PackedValue; }
};

template <typename T, int N> struct mpi_serializable_traits<interval<T,N>,
  OVK_SPECIALIZATION_REQUIRES(IsMPISerializable<T>())> {
  using packed_type = interval<T,N>;
  static handle<MPI_Datatype> CreateMPIType() {
    using tuple_type = typename interval<T,N>::tuple_type;
    auto TupleMPIType = mpi_serializable_traits<tuple_type>::CreateMPIType();
    return CreateMPIContiguousType(2, TupleMPIType);
  }
  static packed_type Pack(const interval<T,N> &Value) { return Value; }
  static interval<T,N> Unpack(const packed_type &PackedValue) { return PackedValue; }
};

}}

#endif
