// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_STRING_WRAPPER_HPP_INCLUDED
#define OVK_CORE_STRING_WRAPPER_HPP_INCLUDED

#include <ovk/core/Global.hpp>

#include <memory>
#include <string>
#include <utility>

namespace ovk {
namespace core {

// std::string isn't noexcept movable until C++17; workaround
class string_wrapper {

public:

  string_wrapper():
    String_(new std::string())
  {}

  string_wrapper(std::string String):
    String_(new std::string(std::move(String)))
  {}

  string_wrapper(const char *String):
    String_(new std::string(String))
  {}

  string_wrapper(const string_wrapper &Other):
    String_(new std::string(*Other.String_))
  {}

  string_wrapper(string_wrapper &&Other) noexcept:
    String_(std::move(Other.String_))
  {}

  string_wrapper &operator=(const string_wrapper &Other) {
    if (String_) {
      *String_ = *Other.String_;
    } else {
      String_.reset(new std::string(*Other.String_));
    }
    return *this;
  }

  string_wrapper &operator=(string_wrapper &&Other) noexcept {
    String_ = std::move(Other.String_);
    return *this;
  }

  const std::string &operator*() const { return *String_; }
  std::string &operator*() { return *String_; }

  const std::string *operator->() const { return String_.get(); }
  std::string *operator->() { return String_.get(); }

  const std::string &Get() const { return *String_; }
  std::string &Get() { return *String_; }

  operator const std::string &() const { return *String_; }
  operator std::string &() { return *String_; }

private:

  std::unique_ptr<std::string> String_;

};

}}

#endif
