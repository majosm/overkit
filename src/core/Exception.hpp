// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_EXCEPTION_HPP_INCLUDED
#define OVK_CORE_EXCEPTION_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <limits>
#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

template <typename T> struct exception_traits {
//   using code_type = ???;
//   using code_underlying_type = ???;
//   static constexpr code_type SuccessCode = ???;
//   static captured_exception<T> CaptureFromCode(code_type Code) { ??? }
//   static code_type GetCode(const T &Exception) { ??? }
//   static void SyncAuxData(T &Exception, int Root, comm_view Comm) { ??? }
};

namespace is_exception_internal {
template <typename T> using maybe_int = int;
template <typename T> constexpr std::true_type Test(maybe_int<typename
  exception_traits<T>::code_type>) {
  return {};
}
template <typename T> constexpr std::false_type Test(...) { return {}; }
}
template <typename T> constexpr bool IsException() {
  return decltype(is_exception_internal::Test<T>(0))::value;
}

namespace exception_code_type_internal {
template <typename T, typename=void> struct helper;
template <typename T> struct helper<T, OVK_SPECIALIZATION_REQUIRES(IsException<T>())> {
  using type = typename exception_traits<T>::code_type;
};
template <typename T> struct helper<T, OVK_SPECIALIZATION_REQUIRES(!IsException<T>())> {
  using type = std::false_type;
};
}
template <typename T> using exception_code_type = typename exception_code_type_internal::helper<T
  >::type;

namespace exception_code_underlying_type_internal {
template <typename T, typename=void> struct helper;
template <typename T> struct helper<T, OVK_SPECIALIZATION_REQUIRES(IsException<T>())> {
  using type = typename exception_traits<T>::code_underlying_type;
};
template <typename T> struct helper<T, OVK_SPECIALIZATION_REQUIRES(!IsException<T>())> {
  using type = std::false_type;
};
}
template <typename T> using exception_code_underlying_type = typename
  exception_code_underlying_type_internal::helper<T>::type;

template <typename ExceptionType, OVK_FUNCTION_REQUIRES(IsException<ExceptionType>())> constexpr
  exception_code_type<ExceptionType> ExceptionSuccessCode() {
  return exception_traits<ExceptionType>::SuccessCode;
}

template <typename ExceptionType> class captured_exception;

template <typename ExceptionType, OVK_FUNCTION_REQUIRES(IsException<ExceptionType>())>
  captured_exception<ExceptionType> CaptureExceptionFromCode(typename exception_traits<
  ExceptionType>::code_type Code) {
  return exception_traits<ExceptionType>::CaptureFromCode(Code);
}

template <typename ExceptionType, OVK_FUNCTION_REQUIRES(IsException<ExceptionType>())>
  exception_code_type<ExceptionType> GetExceptionCode(const ExceptionType &Exception) {
  return exception_traits<ExceptionType>::GetCode(Exception);
}

template <typename ExceptionType, OVK_FUNCTION_REQUIRES(IsException<ExceptionType>())>
  void SyncExceptionAuxData(ExceptionType &Exception, int Root, comm_view Comm) {
  exception_traits<ExceptionType>::SyncAuxData(Exception, Root, Comm);
}

template <typename ExceptionType> class captured_exception {

public:

  static_assert(IsException<ExceptionType>(), "Invalid exception type.");

  using exception_type = ExceptionType;
  using code_type = exception_code_type<exception_type>;

  captured_exception() = default;
  template <typename T, OVK_FUNCTION_REQUIRES(std::is_base_of<exception_type, core::remove_cvref<T>
    >::value)> captured_exception(T &&Exception):
    Exception_(new model<core::remove_cvref<T>>(std::forward<T>(Exception)))
  {}

  template <typename T, OVK_FUNCTION_REQUIRES(std::is_base_of<exception_type, core::remove_cvref<T>
    >::value)> captured_exception &operator=(T &&Exception) {
    Exception_.reset(new model<core::remove_cvref<T>>(std::forward<T>(Exception)));
    return *this;
  }

  explicit operator bool() const { return static_cast<bool>(Exception_); }

  bool Present() const { return static_cast<bool>(Exception_); }

  code_type Code() const {
    return Exception_ ? Exception_->Code() : ExceptionSuccessCode<exception_type>();
  }

  template <typename T> const T &Get() const {
    return static_cast<const model<T> *>(Exception_.get())->Exception_;
  }

  template <typename T> T &Get() {
    return static_cast<model<T> *>(Exception_.get())->Exception_;
  }

  captured_exception &Sync(comm_view Comm) {

    using code_type = exception_code_type<exception_type>;
    using code_underlying_type = exception_code_underlying_type<exception_type>;

    code_underlying_type MaxCodeValue = std::numeric_limits<code_underlying_type>::max();

    code_underlying_type CodeValue = Exception_ ? code_underlying_type(Exception_->Code()) :
      MaxCodeValue;
    code_underlying_type LowestCodeValue = CodeValue;
    MPI_Allreduce(MPI_IN_PLACE, &LowestCodeValue, 1, GetMPIDataType<code_underlying_type>(),
      MPI_MIN, Comm);

    if (LowestCodeValue != MaxCodeValue) {
      if (!Exception_) {
        *this = CaptureExceptionFromCode<exception_type>(code_type(LowestCodeValue));
      }
      int Root = Comm.Size();
      if (CodeValue == LowestCodeValue) {
        Root = Comm.Rank();
      }
      MPI_Allreduce(MPI_IN_PLACE, &Root, 1, MPI_INT, MPI_MIN, Comm);
      Exception_->SyncAuxData(Root, Comm);
    }

    return *this;

  }

  void Throw() const { Exception_->Throw(); }

  void Check() const { if (Exception_) Exception_->Throw(); }
  void Check(comm_view Comm) { Sync(Comm).Check(); }

  void Reset() { Exception_.reset(); }

  template <typename T> T Release() {
    T ReleasedException = std::move(static_cast<model<T> *>(Exception_.get())->Exception_);
    Exception_.reset();
    return ReleasedException;
  }

private:

  struct concept {
    code_type Code_;
    explicit concept(code_type Code): Code_(Code) {}
    virtual ~concept() noexcept {}
    code_type Code() const { return Code_; }
    virtual void Throw() = 0;
    virtual void SyncAuxData(int Root, comm_view Comm) = 0;
  };

  template <typename T> struct model final : concept {
    T Exception_;
    explicit model(T Exception):
      concept(GetExceptionCode<exception_type>(Exception)),
      Exception_(std::move(Exception))
    {}
    virtual void Throw() override { throw Exception_; }
    virtual void SyncAuxData(int Root, comm_view Comm) override {
      SyncExceptionAuxData<exception_type>(Exception_, Root, Comm);
    }
  };

  std::unique_ptr<concept> Exception_;

};

}}

#endif
