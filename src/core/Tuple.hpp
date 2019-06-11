// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_TUPLE_HPP_INCLUDED
#define OVK_CORE_TUPLE_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>

namespace ovk {

template <typename T> using tuple = elem<T,MAX_DIMS>;

template <typename T> inline tuple<T> MakeUniformTuple(int NumDims, T Value, T PadValue=T()) {

  tuple<T> Tuple;

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    Tuple(iDim) = Value;
  }
  for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
    Tuple(iDim) = PadValue;
  }

  return Tuple;
}

}

#endif
