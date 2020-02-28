// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline editor::editor(editor &&Other) noexcept:
  FloatingRefGenerator_(std::move(Other.FloatingRefGenerator_)),
  RefCount_(Other.RefCount_),
  DeactivateFunc_(std::move(Other.DeactivateFunc_))
{
  Other.RefCount_ = 0;
}

inline editor &editor::operator=(editor &&Other) noexcept {

  FloatingRefGenerator_ = std::move(Other.FloatingRefGenerator_);
  RefCount_ = Other.RefCount_;
  DeactivateFunc_ = std::move(Other.DeactivateFunc_);

  Other.RefCount_ = 0;

  return *this;
}

template <typename F, OVK_FUNCDEF_REQUIRES(core::IsCallableWith<F>())> void editor::Activate(F
  DeactivateFunc) {

  DeactivateFunc_.reset(new std::function<void()>(std::move(DeactivateFunc)));

}

template <typename T> edit_handle<T> editor::Edit(T &Target) {

  return edit_handle<T>(*this, Target);

}

inline void editor::Restore() {

  Decrement_();

}

inline void editor::Increment_() {

  ++RefCount_;

}

inline void editor::Decrement_() {

  --RefCount_;
  if (RefCount_ == 0) {
    (*DeactivateFunc_)();
    DeactivateFunc_.reset();
  }

}

template <typename T> edit_handle<T>::edit_handle(editor &Editor, T &Target):
  Editor_(Editor.FloatingRefGenerator_.Generate(Editor)),
  Target_(&Target)
{
  Editor.Increment_();
}

template <typename T> edit_handle<T>::edit_handle(const edit_handle &Other):
  Editor_(Other.Editor_),
  Target_(Other.Target_)
{
  if (Editor_) {
    Editor_->Increment_();
  }
}

template <typename T> edit_handle<T>::edit_handle(edit_handle &&Other) noexcept:
  Editor_(std::move(Other.Editor_)),
  Target_(Other.Target_)
{
  Other.Target_ = nullptr;
}

template <typename T> edit_handle<T> &edit_handle<T>::operator=(edit_handle Other) noexcept {

  using std::swap;

  swap(Editor_, Other.Editor_);
  swap(Target_, Other.Target_);

  return *this;

}

template <typename T> edit_handle<T>::~edit_handle() noexcept {

  if (Editor_) {
    Editor_->Decrement_();
  }

}

template <typename T> T &edit_handle<T>::operator*() const {

  OVK_DEBUG_ASSERT(Editor_, "Invalid edit handle.");

  return *Target_;

}

template <typename T> T *edit_handle<T>::operator->() const {

  OVK_DEBUG_ASSERT(Editor_, "Invalid edit handle.");

  return Target_;

}

template <typename T> T &edit_handle<T>::Get() const {

  OVK_DEBUG_ASSERT(Editor_, "Invalid edit handle.");

  return *Target_;

}

template <typename T> void edit_handle<T>::Restore() {

  if (Editor_) {
    Editor_->Decrement_();
  }

}

template <typename T> T *edit_handle<T>::Release() {

  OVK_DEBUG_ASSERT(Editor_, "Invalid edit handle.");

  T *Target = Target_;

  Editor_.Reset();
  Target_ = nullptr;

  return Target;

}

}
