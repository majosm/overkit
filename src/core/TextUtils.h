// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_TEXT_UTILS_INCLUDED
#define OVK_CORE_TEXT_UTILS_INCLUDED

#include "Global.h"

#include <stdio.h>
#include <string.h>

static inline void SizeToString(size_t N, char *NString) {

  int i, j;

  char UnformattedNString[32];

  sprintf(UnformattedNString, "%zu", N);

  int NumDigits = strlen(UnformattedNString);
  int NumBeforeComma = ((NumDigits-1) % 3) + 1;

  strncpy(NString, UnformattedNString, NumBeforeComma);

  j = NumBeforeComma;
  for (i = NumBeforeComma; i < NumDigits; ++i) {
    if ((i-NumBeforeComma) % 3 == 0) {
      NString[j] = ',';
      ++j;
    }
    NString[j] = UnformattedNString[i];
    ++j;
  }

  NString[j] = '\0';

}

static inline void IntToString(int N, char *NString) {

  if (N >= 0) {
    SizeToString(N, NString);
  } else {
    NString[0] = '-';
    SizeToString(abs(N), NString+1);
  }

}

static inline void PluralizeLabel(size_t Count, const char *PluralLabel, const char *SingularLabel,
  char *String) {

  char CountString[32];
  SizeToString(Count, CountString);

  if (Count > 1 || Count == 0) {
    sprintf(String, "%s %s", CountString, PluralLabel);
  } else {
    sprintf(String, "%s %s", CountString, SingularLabel);
  }

}

#endif
