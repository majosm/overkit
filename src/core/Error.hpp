// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ERROR_HPP_INCLUDED
#define OVK_CORE_ERROR_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/CommunicationOps.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Error.h>
#include <ovk/core/Exception.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>

#include <mpi.h>

#include <stdexcept>
#include <type_traits>

namespace ovk {

enum class error_code : typename std::underlying_type<ovk_error>::type {
  NONE = OVK_ERROR_NONE,
  MPI_NOT_INITIALIZED = OVK_ERROR_MPI_NOT_INITIALIZED,
  FILE_OPEN = OVK_ERROR_FILE_OPEN,
  FILE_READ = OVK_ERROR_FILE_READ,
  FILE_WRITE = OVK_ERROR_FILE_WRITE
};

class error;

namespace core {
template <> struct exception_traits<error> {
  using code_type = error_code;
  using code_underlying_type = typename std::underlying_type<error_code>::type;
  static constexpr error_code SuccessCode = error_code::NONE;
  static captured_exception<error> CaptureFromCode(error_code Code);
  static error_code GetCode(const error &Error);
  static void SyncAuxData(error &Error, int Root, comm_view Comm);
};
}

using captured_error = core::captured_exception<error>;

class error : public std::runtime_error {
public:
  error(error_code Code, const char *ErrorString):
    runtime_error(ErrorString),
    Code_(Code)
  {}
  virtual ~error() noexcept {}
  error_code Code() const { return Code_; }
  virtual captured_error Capture() const = 0;
  virtual void SyncAuxData(int Root, comm_view Comm) {}
protected:
  error_code Code_;
};

class mpi_error : public error {
public:
  mpi_error(error_code ErrorCode, const char *ErrorString):
    error(ErrorCode, ErrorString)
  {}
};

class mpi_not_initialized_error : public mpi_error {
public:
  mpi_not_initialized_error():
    mpi_error(error_code::MPI_NOT_INITIALIZED, "ovk::mpi_not_initialized_error")
  {}
  virtual captured_error Capture() const override { return *this; }
};

class io_error : public error {
public:
  io_error(error_code ErrorCode, const char *ErrorString):
    error(ErrorCode, ErrorString)
  {}
};

class file_open_error : public io_error {
public:
  file_open_error():
    io_error(error_code::FILE_OPEN, "ovk::file_open_error")
  {}
  explicit file_open_error(std::string FilePath):
    file_open_error()
  {
    FilePath_ = std::move(FilePath);
  }
  virtual captured_error Capture() const override { return *this; }
  virtual void SyncAuxData(int Root, comm_view Comm) override {
    core::BroadcastString(FilePath_, Root, Comm);
  }
protected:
  std::string FilePath_;
};

class file_read_error : public io_error {
public:
  file_read_error():
    io_error(error_code::FILE_READ, "ovk::file_read_error")
  {}
  explicit file_read_error(std::string FilePath):
    file_read_error()
  {
    FilePath_ = std::move(FilePath);
  }
  virtual captured_error Capture() const override { return *this; }
  virtual void SyncAuxData(int Root, comm_view Comm) override {
    core::BroadcastString(FilePath_, Root, Comm);
  }
protected:
  std::string FilePath_;
};

class file_write_error : public io_error {
public:
  file_write_error():
    io_error(error_code::FILE_WRITE, "ovk::file_write_error")
  {}
  explicit file_write_error(std::string FilePath):
    file_write_error()
  {
    FilePath_ = std::move(FilePath);
  }
  virtual captured_error Capture() const override { return *this; }
  virtual void SyncAuxData(int Root, comm_view Comm) override {
    core::BroadcastString(FilePath_, Root, Comm);
  }
protected:
  std::string FilePath_;
};

namespace core {
inline captured_exception<error> exception_traits<error>::CaptureFromCode(error_code Code) {
  switch (Code) {
  case error_code::MPI_NOT_INITIALIZED:
    return mpi_not_initialized_error();
  case error_code::FILE_OPEN:
    return file_open_error();
  case error_code::FILE_READ:
    return file_read_error();
  case error_code::FILE_WRITE:
    return file_write_error();
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    return {};
  }
}
inline error_code exception_traits<error>::GetCode(const error &Error) {
  return Error.Code();
}
inline void exception_traits<error>::SyncAuxData(error &Error, int Root, comm_view Comm) {
  Error.SyncAuxData(Root, Comm);
}
}

}

#endif
