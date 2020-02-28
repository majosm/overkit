// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_FIELD_OPS_HPP_INCLUDED
#define OVK_CORE_FIELD_OPS_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Requires.hpp>

#include <cstdio>

namespace ovk {
namespace core {

namespace print_field_internal {
inline std::string FormatValue(int Value, int Width) {
  char Format[16];
  std::sprintf(Format, "%%%ii", Width);
  char Chars[16];
  std::sprintf(Chars, Format, Value);
  return Chars;
}
inline std::string FormatValue(long long Value, int Width) {
  char Format[16];
  std::sprintf(Format, "%%%illd", Width);
  char Chars[16];
  std::sprintf(Chars, Format, Value);
  return Chars;
}
inline std::string FormatValue(double Value, int Width) {
  char Format[16];
  std::sprintf(Format, "%%%i.%ie", Width, Width-6);
  char Chars[16];
  std::sprintf(Chars, Format, Value);
  return Chars;
}
inline std::string FormatValue(bool Value, int) {
  return Value ? "X" : "-";
}
}

template <typename T> void PrintField(field_view<const T> Field, const range &Range, int Width=1) {

  for (int k = Range.Begin(2); k < Range.End(2); ++k) {
    if (Range.Size(2) > 1) {
      std::printf("k = %i:\n", k); std::fflush(stdout);
    }
    for (int j = Range.End(1)-1; j >= Range.Begin(1); --j) {
      for (int i = Range.Begin(0); i < Range.End(0); ++i) {
        std::printf(" %s ", print_field_internal::FormatValue(Field(i,j,k), Width).c_str());
      }
      std::printf("\n"); std::fflush(stdout);
    }
    std::printf("\n"); std::fflush(stdout);
  }

}

template <typename T> void PrintField(field_view<const T> Field, int Width=1) {
  PrintField(Field, Field.Extents(), Width);
}

template <typename FieldType, OVK_FUNCTION_REQUIRES(IsField<FieldType>())> void PrintField(const
  FieldType &Field, const range &Range, int Width=1) {
  field_view<const array_value_type<FieldType>> View(Field);
  PrintField(View, Range, Width);
}

template <typename FieldType, OVK_FUNCTION_REQUIRES(IsField<FieldType>())> void PrintField(const
  FieldType &Field, int Width=1) {
  field_view<const array_value_type<FieldType>> View(Field);
  PrintField(View, Width);
}

}}

#endif
