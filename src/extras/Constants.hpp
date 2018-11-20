// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_CONSTANTS_HPP_INCLUDED
#define OVK_EXTRAS_CONSTANTS_HPP_INCLUDED

#include <ovk/extras/ConstantsBase.h>
#include <ovk/extras/Global.hpp>

namespace ovk {

enum class endian {
  LITTLE = OVK_LITTLE_ENDIAN,
  BIG = OVK_BIG_ENDIAN
};

#ifdef OVK_XPACC
enum class xintout_format {
  STANDARD = OVK_XINTOUT_STANDARD,
  EXTENDED = OVK_XINTOUT_EXTENDED
};
#endif

inline bool ValidEndian(endian Endian) { return ovkValidEndian(ovk_endian(Endian)); }
#ifdef OVK_XPACC
inline bool ValidXINTOUTFormat(xintout_format Format) { return ovkValidXINTOUTFormat(ovk_xintout_format(Format)); }
#endif

}

#endif
