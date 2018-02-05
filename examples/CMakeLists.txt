# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#==================
# Overkit examples
#==================

if(EXAMPLES)

# Link to overkit and overkit-support
include_directories("${CMAKE_SOURCE_DIR}/ovk/core")
# include_directories("${CMAKE_SOURCE_DIR}/support")
# link_libraries(overkit-support overkit)
link_libraries(overkit overkit-extras)

add_executable(Test Test.c)
# add_executable(Debug Debug.c)
if(XPACC)
add_executable(ReadInterp ReadInterp.c)
endif()

endif()
