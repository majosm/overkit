// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_EDITOR_HPP_INCLUDED
#define OVK_CORE_EDITOR_HPP_INCLUDED

#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <functional>
#include <memory>
#include <utility>

namespace ovk {

template <typename T> class edit_handle;

class editor {

public:

  editor() = default;

  editor(const editor &Other) = delete;
  editor(editor &&Other) noexcept;

  editor &operator=(const editor &Other) = delete;
  editor &operator=(editor &&Other) noexcept;

  bool Active() const { return static_cast<bool>(DeactivateFunc_); }

  template <typename F, OVK_FUNCDECL_REQUIRES(core::IsCallableWith<F>())> void Activate(F
    DeactivateFunc);

  template <typename T> edit_handle<T> Edit(T &Target);

  void Restore();

private:

  floating_ref_generator FloatingRefGenerator_;

  int RefCount_ = 0;
  // std::function is not noexcept movable until C++20
  std::unique_ptr<std::function<void()>> DeactivateFunc_;

  void Increment_();
  void Decrement_();

  template <typename T> friend class edit_handle;

};

template <typename T> class edit_handle {

public:

  edit_handle() = default;

  edit_handle(const edit_handle &Other);
  edit_handle(edit_handle &&Other) noexcept;

  edit_handle &operator=(edit_handle Other) noexcept;

  ~edit_handle() noexcept;

  explicit operator bool() const { return static_cast<bool>(Editor_); }

  T &operator*() const;

  T *operator->() const;

  T &Get() const;

  void Restore();

  T *Release();

private:

  floating_ref<editor> Editor_;
  T *Target_ = nullptr;

  edit_handle(editor &Editor, T &Target);

  friend class editor;

};

}

#include <ovk/core/Editor.inl>

#endif
