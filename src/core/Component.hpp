// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COMPONENT_HPP_INCLUDED
#define OVK_CORE_COMPONENT_HPP_INCLUDED

#include <ovk/core/DomainBase.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

class component {

public:

  component() = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, component>::value)>
    component(T &&Component):
    Component_(new model<T>(std::forward<T>(Component)))
  {}

  component(const component &Other) = delete;
  component(component &&Other) noexcept = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, component>::value)>
    component &operator=(T &&Component) {
    Component_.reset(new model<T>(std::forward<T>(Component)));
    return *this;
  }

  component &operator=(const component &Other) = delete;
  component &operator=(component &&Other) noexcept = default;

  explicit operator bool() { return static_cast<bool>(Component_); }

  bool Present() const { return static_cast<bool>(Component_); }

  template <typename T> const T &Get() const;
  template <typename T> T &Get();

  void Reset() { Component_.reset(); }

  template <typename T> T Release();

  void StartEdit() { return Component_->StartEdit(); }
  void EndEdit() { return Component_->EndEdit(); }

private:

  struct concept {
    virtual ~concept() noexcept {}
    virtual void StartEdit() = 0;
    virtual void EndEdit() = 0;
  };

  template <typename T> struct model final : concept {
    T Component_;
    explicit model(T Component):
      Component_(std::move(Component))
    {}
    virtual void StartEdit() override {
      return Component_.StartEdit();
    }
    virtual void EndEdit() override {
      return Component_.EndEdit();
    }
  };

  std::unique_ptr<concept> Component_;

};

}}

#include <ovk/core/Component.inl>

#endif
