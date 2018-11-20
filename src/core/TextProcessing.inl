// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

inline std::string FormatNumber(size_t N) {

  int iInput;

  std::string NString;

  char UnformattedNStringChars[32];
  std::sprintf(UnformattedNStringChars, "%zu", N);

  int NumDigits = strlen(UnformattedNStringChars);
  int NumBeforeComma = ((NumDigits-1) % 3) + 1;

  for (iInput = 0; iInput < NumBeforeComma; ++iInput) {
    NString += UnformattedNStringChars[iInput];
  }

  for (iInput = NumBeforeComma; iInput < NumDigits; ++iInput) {
    if ((iInput-NumBeforeComma) % 3 == 0) NString += ',';
    NString += UnformattedNStringChars[iInput];
  }

  return NString;

}

inline std::string FormatNumber(long long N) {

  std::string NString;

  if (N >= 0) {
    NString = FormatNumber((size_t)N);
  } else {
    NString = '-';
    NString += FormatNumber((size_t)llabs(N));
  }

  return NString;

}

inline std::string FormatNumber(int N) {

  std::string NString;

  if (N >= 0) {
    NString = FormatNumber((size_t)N);
  } else {
    NString = '-';
    NString += FormatNumber((size_t)abs(N));
  }

  return NString;

}

template <typename IntegerType> std::string FormatNumber(IntegerType N, const std::string
  &PluralLabel, const std::string &SingularLabel) {

  std::string NString;

  if (N > 1 || N == 0) {
    NString = FormatNumber(N) + " " + PluralLabel;
  } else {
    NString = FormatNumber(N) + " " + SingularLabel;
  }

  return NString;

}

namespace text_processing_internal {

template <typename T> void StringPrintAppend(std::string &String, const std::string &Format,
  const T &Arg) {
  std::vector<char> Buffer(1);
  int WriteSize = snprintf(Buffer.data(), 1, Format.c_str(), Arg);
  Buffer.resize(WriteSize+1);
  std::sprintf(Buffer.data(), Format.c_str(), Arg);
  String.append(Buffer.begin(), Buffer.begin()+WriteSize);
}

// Specialization for std::string to avoid explicit use of .c_str() everywhere
inline void StringPrintAppend(std::string &String, const std::string &Format, const std::string
  &Arg) {
  std::vector<char> Buffer(1);
  int WriteSize = snprintf(Buffer.data(), 1, Format.c_str(), Arg.c_str());
  Buffer.resize(WriteSize+1);
  std::sprintf(Buffer.data(), Format.c_str(), Arg.c_str());
  String.append(Buffer.begin(), Buffer.begin()+WriteSize);
}

template <typename T1, typename T2, typename... Ts> void StringPrintAppend(std::string &String,
  const std::string &Format, const T1 &Arg1, const T2 &Arg2, const Ts &... RemainingArgs) {

  size_t iSplitPos = 0;
  // C++11 requires String[String.length()] to be well-defined, so accessing it is OK
  for (size_t iPos = 1; iPos < Format.length(); ++iPos) {
    if (Format[iPos] == '%' && Format[iPos+1] != '%') {
      iSplitPos = iPos;
      break;
    }
  }

  StringPrintAppend(String, Format.substr(0, iSplitPos), Arg1);
  StringPrintAppend(String, Format.substr(iSplitPos, std::string::npos), Arg2, RemainingArgs...);

}

}

template <typename T, typename... Ts> std::string StringPrint(const std::string &Format, const T
  &Arg, const Ts &... RemainingArgs) {

  std::string String;

  size_t NumSplits = 0;
  // C++11 requires String[String.length()] to be well-defined, so accessing it is OK
  for (size_t iPos = 0; iPos < Format.length(); ++iPos) {
    if (Format[iPos] == '%' && Format[iPos+1] != '%') {
      ++NumSplits;
    }
  }

  // Can't use C++ debug assert because it depends on StringPrint
  OVK_DEBUG_ASSERT_C(NumSplits == 1+sizeof...(RemainingArgs), "Incorrect number of arguments.");

  size_t iSplitPos = Format.length();
  // C++11 requires String[String.length()] to be well-defined, so accessing it is OK
  for (size_t iPos = 0; iPos < Format.length(); ++iPos) {
    if (Format[iPos] == '%' && Format[iPos+1] != '%') {
      iSplitPos = iPos;
      break;
    }
  }

  String.append(Format, 0, iSplitPos);

  text_processing_internal::StringPrintAppend(String, Format.substr(iSplitPos, std::string::npos),
    Arg, RemainingArgs...);

  return String;

}

inline std::string StringPrint(const std::string &Format) {

  return {Format};

}

}}
