// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MISC_UTILS_INCLUDED
#define OVK_CORE_MISC_UTILS_INCLUDED

#include "ovk/core/Global.h"

#ifdef __cplusplus
extern "C" {
#endif

void PRIVATE(SortPermutation_long_long)(long long N, const long long *Array, long long *Permutation);
#define SortPermutation_long_long(...) PRIVATE(SortPermutation_long_long)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
