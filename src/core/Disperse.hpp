// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISPERSE_HPP_INCLUDED
#define OVK_CORE_DISPERSE_HPP_INCLUDED

#include <ovk/core/Context.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/DisperseMap.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

class disperse {

public:

  disperse() = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, disperse>::value)>
    disperse(T &&Disperse):
    Disperse_(new model<T>(std::forward<T>(Disperse)))
  {}

  disperse(const disperse &Other) = delete;
  disperse(disperse &&Other) noexcept = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, disperse>::value)>
    disperse &operator=(T &&Disperse) {
    Disperse_.reset(new model<T>(std::forward<T>(Disperse)));
    return *this;
  }

  disperse &operator=(const disperse &Other) = delete;
  disperse &operator=(disperse &&Other) noexcept = default;

  void Disperse(const void *PackedValues, void *FieldValues) {
    Disperse_->Disperse(PackedValues, FieldValues);
  }

private:

  class concept {
  public:
    virtual ~concept() noexcept {}
    virtual void Disperse(const void *PackedValues, void *FieldValues) = 0;
  };

  template <typename T> class model final : public concept {
  public:
    explicit model(T Disperse):
      Disperse_(std::move(Disperse))
    {}
    virtual void Disperse(const void *PackedValues, void *FieldValues) override {
      Disperse_.Disperse(PackedValues, FieldValues);
    }
  private:
    T Disperse_;
  };

  std::unique_ptr<concept> Disperse_;

};

disperse CreateDisperseOverwrite(std::shared_ptr<context> Context, const disperse_map &DisperseMap,
  data_type ValueType, int Count, const range &FieldValuesRange, array_layout FieldValuesLayout);
disperse CreateDisperseAppend(std::shared_ptr<context> Context, const disperse_map &DisperseMap,
  data_type ValueType, int Count, const range &FieldValuesRange, array_layout FieldValuesLayout);

}}

#endif
