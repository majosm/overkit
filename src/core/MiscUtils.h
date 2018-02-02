// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MISC_UTILS_INCLUDED
#define OVK_CORE_MISC_UTILS_INCLUDED

#include "Global.h"

#ifdef __cplusplus
extern "C" {
#endif

void PRIVATE(SortPermutation_size_t)(size_t N, const size_t *Array, size_t *Permutation);
#define SortPermutation_size_t(...) PRIVATE(SortPermutation_size_t)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
