// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename ArrayType, OVK_FUNCDEF_REQUIRES(IsArray<ArrayType>() && ArrayRank<ArrayType>()
  == 1)> void SortPermutation(const ArrayType &Array, array_view<long long> Permutation) {

  for (long long i = 0; i < Permutation.Count(); ++i) {
    Permutation[i] = i;
  }

  auto Compare = [&Array](long long Left, long long Right) -> bool {
    return Array[Left] < Array[Right];
  };

  std::sort(Permutation.Begin(), Permutation.End(), Compare);

}

}}
