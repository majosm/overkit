// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "MiscUtils.h"

#include "Global.h"

#include <algorithm>
#include <numeric>

namespace {

template <typename T> struct sort_permutation_compare {
  explicit sort_permutation_compare(const T *Array): array(Array) {}
  bool operator()(const T &Left, const T &Right) const { return array[Left] < array[Right]; }
  const T *array;
};

}

extern "C" {

void PRIVATE(SortPermutation_size_t)(size_t N, const size_t *Array, size_t *Permutation) {

  for (size_t i = 0; i < N; ++i) {
    Permutation[i] = i;
  }

  std::sort(Permutation, Permutation+N, sort_permutation_compare<size_t>(Array));

}

}
