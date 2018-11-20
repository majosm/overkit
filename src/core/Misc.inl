// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename T> void SortPermutation(long long N, const T *Array, long long *Permutation) {

  for (long long i = 0; i < N; ++i) {
    Permutation[i] = i;
  }

  auto Compare = [&Array](long long Left, long long Right) -> bool {
    return Array[Left] < Array[Right];
  };

  std::sort(Permutation, Permutation+N, Compare);

}

}}
