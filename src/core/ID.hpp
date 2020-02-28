// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ID_HPP_INCLUDED
#define OVK_CORE_ID_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/Set.hpp>

namespace ovk {

inline int NextAvailableID(const set<int> &IDs) {

  int PrevID = -1;

  auto Iter = IDs.Begin();
  while (Iter != IDs.End()) {
    if (*Iter - PrevID > 1) {
      break;
    }
    PrevID = *Iter;
    ++Iter;
  }

  return PrevID + 1;

}

}

#endif
