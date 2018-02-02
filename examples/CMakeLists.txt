# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#==================
# Overkit examples
#==================

if(EXAMPLES)

# Link to overkit and overkit-support
include_directories("${CMAKE_SOURCE_DIR}/src/core")
# include_directories("${CMAKE_SOURCE_DIR}/src/support")
# link_libraries(overkit-support overkit)
link_libraries(overkit)

add_executable(Test Test.c)
# add_executable(Debug Debug.c)

endif()
