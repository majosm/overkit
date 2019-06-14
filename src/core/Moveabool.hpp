// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MOVEABOOL_HPP_INCLUDED
#define OVK_CORE_MOVEABOOL_HPP_INCLUDED

#include <ovk/core/Global.hpp>

namespace ovk {
namespace core {

// Determine whether an object is in a moved-from state or not in cases where there's no
// obvious resource being wrapped
class moveabool {

public:

  moveabool():
    Value_(false)
  {}

  moveabool(bool Value):
    Value_(Value)
  {}

  moveabool(const moveabool &Other):
    Value_(Other.Value_)
  {}
  moveabool(moveabool &&Other) noexcept:
    Value_(Other.Value_)
  {
    Other.Value_ = false;
  }

  moveabool &operator=(const moveabool &Other) {
    Value_ = Other.Value_;
    return *this;
  }
  moveabool &operator=(moveabool &&Other) noexcept {
    Value_ = Other.Value_;
    Other.Value_ = false;
    return *this;
  }

  operator bool() const { return Value_; }

  bool Value() const { return Value_; }

private:

  bool Value_;

};

}}

#endif
