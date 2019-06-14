// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_TESTS_MOCKS_NONDEFAULTCONSTRUCTIBLE_HPP_LOADED
#define OVK_TESTS_MOCKS_NONDEFAULTCONSTRUCTIBLE_HPP_LOADED

namespace tests {

// Non-default-constructible type
template <typename T> class nondefaultconstructible {
public:
  nondefaultconstructible(const T &Value):
    Value_(Value)
  {}
  const T &Value() const { return Value_; }
private:
  T Value_;
};

}

#endif
