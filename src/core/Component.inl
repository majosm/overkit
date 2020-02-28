// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename T> const T &component::Get() const {

  return static_cast<const model<T> *>(Component_.get())->Component_;

}

template <typename T> T &component::Get() {

  return static_cast<model<T> *>(Component_.get())->Component_;

}

template <typename T> T component::Release() {

  T ReleasedComponent = std::move(static_cast<model<T> *>(Component_.get())->Component_);

  Component_.reset();

  return ReleasedComponent;

}

}}
