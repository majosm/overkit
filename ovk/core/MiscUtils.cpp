// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/MiscUtils.h"

#include "ovk/core/Global.h"

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

void PRIVATE(SortPermutation_long_long)(long long N, const long long *Array, long long *Permutation) {

  for (long long i = 0; i < N; ++i) {
    Permutation[i] = i;
  }

  std::sort(Permutation, Permutation+N, sort_permutation_compare<long long>(Array));

}

}
