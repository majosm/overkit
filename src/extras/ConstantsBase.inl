// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifdef __cplusplus
extern "C" {
#endif

static inline bool ovkValidEndian(ovk_endian Endian) {

  switch (Endian) {
  case OVK_LITTLE_ENDIAN:
  case OVK_BIG_ENDIAN:
    return true;
  default:
    return false;
  }

}

#ifdef OVK_XPACC
static inline bool ovkValidXINTOUTFormat(ovk_xintout_format Format) {

  switch (Format) {
  case OVK_XINTOUT_STANDARD:
  case OVK_XINTOUT_EXTENDED:
    return true;
  default:
    return false;
  }

}
#endif

#ifdef __cplusplus
}
#endif
