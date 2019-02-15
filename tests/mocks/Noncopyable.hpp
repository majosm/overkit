// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_TESTS_MOCKS_NONCOPYABLE_HPP_LOADED
#define OVK_TESTS_MOCKS_NONCOPYABLE_HPP_LOADED

namespace tests {

// Non-copyable type
template <typename T> class noncopyable {
public:
  noncopyable() = default;
  noncopyable(const T &Value):
    Value_(Value)
  {}
  noncopyable(const noncopyable &Other) = delete;
  noncopyable(noncopyable &&Other) noexcept:
    Value_(Other.Value_)
  {
    Other.Value_ = T();
  }
  noncopyable &operator=(const noncopyable &Other) = delete;
  noncopyable &operator=(noncopyable &&Other) noexcept {
    Value_ = Other.Value_;
    Other.Value_ = T();
    return *this;
  }
  const T &Value() const { return Value_; }
private:
  T Value_;
};

}

#endif
