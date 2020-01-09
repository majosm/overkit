// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_XDMF_HPP_LOADED
#define OVK_SUPPORT_XDMF_HPP_LOADED

#include <support/XDMF.h>

#ifdef OVK_HAVE_XDMF

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/CommunicationOps.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Exception.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Handle.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Range.hpp>

#include <hdf5.h>
#include <mpi.h>

#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>

namespace support {

enum class xdmf_error_code : typename std::underlying_type<support_xdmf_error>::type {
  NONE = SUPPORT_XDMF_ERROR_NONE,
  FILE_CREATE = SUPPORT_XDMF_ERROR_FILE_CREATE,
  FILE_OPEN = SUPPORT_XDMF_ERROR_FILE_OPEN,
  FILE_READ = SUPPORT_XDMF_ERROR_FILE_READ,
  FILE_WRITE = SUPPORT_XDMF_ERROR_FILE_WRITE
};

class xdmf_error;

using captured_xdmf_error = ovk::core::captured_exception<xdmf_error>;

}

namespace ovk {
namespace core {
template <> struct exception_traits<support::xdmf_error> {
  using code_type = support::xdmf_error_code;
  using underlying_type = typename std::underlying_type<support::xdmf_error_code>::type;
  static constexpr support::xdmf_error_code SuccessCode = support::xdmf_error_code::NONE;
  static support::captured_xdmf_error CaptureFromCode(support::xdmf_error_code Code);
  static support::xdmf_error_code GetCode(const support::xdmf_error &Error);
  static void SyncAuxData(support::xdmf_error &Error, int Root, ovk::comm_view Comm);
};
}}

namespace support {

class xdmf_error : public std::runtime_error {
public:
  xdmf_error(xdmf_error_code Code, const char *ErrorString):
    runtime_error(ErrorString),
    Code_(Code)
  {}
  virtual ~xdmf_error() noexcept {}
  xdmf_error_code Code() const { return Code_; }
  virtual captured_xdmf_error Capture() const = 0;
  virtual void SyncAuxData(int Root, ovk::comm_view Comm) {}
private:
  xdmf_error_code Code_;
};

class xdmf_file_create_error : public xdmf_error {
public:
  xdmf_file_create_error():
    xdmf_error(xdmf_error_code::FILE_CREATE, "support::xdmf_file_create_error")
  {}
  explicit xdmf_file_create_error(std::string FilePath):
    xdmf_file_create_error()
  {
    FilePath_ = std::move(FilePath);
  }
  virtual captured_xdmf_error Capture() const override { return *this; }
  virtual void SyncAuxData(int Root, ovk::comm_view Comm) override {
    ovk::core::BroadcastString(FilePath_, Root, Comm);
  }
private:
  std::string FilePath_;
};

class xdmf_file_open_error : public xdmf_error {
public:
  xdmf_file_open_error():
    xdmf_error(xdmf_error_code::FILE_OPEN, "support::xdmf_file_open_error")
  {}
  explicit xdmf_file_open_error(std::string FilePath):
    xdmf_file_open_error()
  {
    FilePath_ = std::move(FilePath);
  }
  virtual captured_xdmf_error Capture() const override { return *this; }
  virtual void SyncAuxData(int Root, ovk::comm_view Comm) override {
    ovk::core::BroadcastString(FilePath_, Root, Comm);
  }
private:
  std::string FilePath_;
};

class xdmf_file_read_error : public xdmf_error {
public:
  xdmf_file_read_error():
    xdmf_error(xdmf_error_code::FILE_READ, "support::xdmf_file_read_error")
  {}
  explicit xdmf_file_read_error(std::string FilePath):
    xdmf_file_read_error()
  {
    FilePath_ = std::move(FilePath);
  }
  virtual captured_xdmf_error Capture() const override { return *this; }
  virtual void SyncAuxData(int Root, ovk::comm_view Comm) override {
    ovk::core::BroadcastString(FilePath_, Root, Comm);
  }
private:
  std::string FilePath_;
};

class xdmf_file_write_error : public xdmf_error {
public:
  xdmf_file_write_error():
    xdmf_error(xdmf_error_code::FILE_WRITE, "support::xdmf_file_write_error")
  {}
  explicit xdmf_file_write_error(std::string FilePath):
    xdmf_file_write_error()
  {
    FilePath_ = std::move(FilePath);
  }
  virtual captured_xdmf_error Capture() const override { return *this; }
  virtual void SyncAuxData(int Root, ovk::comm_view Comm) override {
    ovk::core::BroadcastString(FilePath_, Root, Comm);
  }
private:
  std::string FilePath_;
};

}

namespace ovk {
namespace core {
inline support::captured_xdmf_error exception_traits<support::xdmf_error>::CaptureFromCode(
  support::xdmf_error_code Code) {
  switch (Code) {
  case support::xdmf_error_code::FILE_CREATE:
    return support::xdmf_file_create_error();
  case support::xdmf_error_code::FILE_OPEN:
    return support::xdmf_file_open_error();
  case support::xdmf_error_code::FILE_READ:
    return support::xdmf_file_read_error();
  case support::xdmf_error_code::FILE_WRITE:
    return support::xdmf_file_write_error();
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    return {};
  }
}
inline support::xdmf_error_code exception_traits<support::xdmf_error>::GetCode(const
  support::xdmf_error &Error) {
  return Error.Code();
}
inline void exception_traits<support::xdmf_error>::SyncAuxData(support::xdmf_error &Error, int Root,
  ovk::comm_view Comm) {
  Error.SyncAuxData(Root, Comm);
}
}}

namespace support {

class xdmf_grid_meta {
public:
  xdmf_grid_meta(std::string Name, const ovk::tuple<int> &Size):
    Name_(std::move(Name)),
    Size_(Size)
  {}
  const std::string &Name() const { return Name_; }
  const ovk::tuple<int> &Size() const { return Size_; }
private:
  std::string Name_;
  ovk::tuple<int> Size_;
};

enum class xdmf_attribute_type : typename std::underlying_type<support_xdmf_attribute_type>::type {
  INT = SUPPORT_XDMF_ATTRIBUTE_TYPE_INT,
  LONG_LONG = SUPPORT_XDMF_ATTRIBUTE_TYPE_LONG_LONG,
  DOUBLE = SUPPORT_XDMF_ATTRIBUTE_TYPE_DOUBLE,
};

}

namespace ovk {
namespace core {
template <> struct data_type_traits<support::xdmf_attribute_type> : data_type_traits<typename
  std::underlying_type<support::xdmf_attribute_type>::type> {};
}}

namespace support {

class xdmf_attribute_meta {
public:
  xdmf_attribute_meta(std::string Name, xdmf_attribute_type Type):
    Name_(std::move(Name)),
    Type_(Type)
  {}
  const std::string &Name() const { return Name_; }
  xdmf_attribute_type Type() const { return Type_; }
private:
  std::string Name_;
  xdmf_attribute_type Type_;
};

class xdmf {

public:

  xdmf() = default;

  explicit operator bool() const { return Path_ != ""; }

  xdmf &WriteGeometry(const std::string &GridName, int Dimension, ovk::field_view<const double>
    Data);
  xdmf &WriteGeometry(const std::string &GridName, int Dimension, ovk::field_view<const double>
    Data, const ovk::range &WriteRange);

  xdmf &WriteAttribute(const std::string &GridName, const std::string &AttributeName,
    ovk::field_view<const int> Data);
  xdmf &WriteAttribute(const std::string &GridName, const std::string &AttributeName,
    ovk::field_view<const int> Data, const ovk::range &WriteRange);
  xdmf &WriteAttribute(const std::string &GridName, const std::string &AttributeName,
    ovk::field_view<const long long> Data);
  xdmf &WriteAttribute(const std::string &GridName, const std::string &AttributeName,
    ovk::field_view<const long long> Data, const ovk::range &WriteRange);
  xdmf &WriteAttribute(const std::string &GridName, const std::string &AttributeName,
    ovk::field_view<const double> Data);
  xdmf &WriteAttribute(const std::string &GridName, const std::string &AttributeName,
    ovk::field_view<const double> Data, const ovk::range &WriteRange);

  xdmf &Close() {
    *this = xdmf();
    return *this;
  }

  const std::string &Path() const { return Path_; }
  int Dimension() const { return NumDims_; }
  const ovk::comm &Comm() const { return Comm_; }
  const ovk::array<xdmf_grid_meta> &Grids() const { return Grids_; }
  const ovk::array<xdmf_attribute_meta> &Attributes() const { return Attributes_; }

  static xdmf internal_Create(std::string &&Path, int NumDims, ovk::comm_view Comm,
    ovk::array<xdmf_grid_meta> &&Grids, ovk::array<xdmf_attribute_meta> &&Attributes);

  static xdmf internal_Open(std::string &&Path, ovk::comm_view Comm);

private:

  std::string Path_;
  int NumDims_;
  ovk::comm Comm_;
  ovk::array<xdmf_grid_meta> Grids_;
  ovk::array<xdmf_attribute_meta> Attributes_;
  ovk::core::handle<hid_t> HDF5File_;

  xdmf(std::string &&Path, int NumDims, ovk::comm &&Comm, ovk::array<xdmf_grid_meta> &&Grids,
    ovk::array<xdmf_attribute_meta> &&Attributes, ovk::core::handle<hid_t> &&HDF5File);

};

xdmf CreateXDMF(std::string Path, int NumDims, ovk::comm_view Comm, ovk::array<xdmf_grid_meta>
  Grids, ovk::array<xdmf_attribute_meta> Attributes);
ovk::optional<xdmf> CreateXDMF(std::string Path, int NumDims, ovk::comm_view Comm,
  ovk::array<xdmf_grid_meta> Grids, ovk::array<xdmf_attribute_meta> Attributes, captured_xdmf_error
  &Error);

xdmf OpenXDMF(std::string Path, ovk::comm_view Comm);
ovk::optional<xdmf> OpenXDMF(std::string Path, ovk::comm_view Comm, captured_xdmf_error &Error);

}

#endif

#endif
