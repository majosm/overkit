// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/TextProcessing.hpp"

#include "ovk/core/Global.hpp"

#include <string>

namespace ovk {
namespace core {

std::string StringReplace(const std::string &String, const std::string &Substring, const std::string
  &Replacement) {

  std::string NewString;

  size_t PrevPos = 0;
  size_t Pos;
  while ((Pos = String.find(Substring, PrevPos)) != std::string::npos) {
    NewString.append(String, PrevPos, Pos-PrevPos);
    NewString += Replacement;
    PrevPos = Pos + Substring.length();
  }
  NewString.append(String, PrevPos, Pos);

  return NewString;

}

}}
