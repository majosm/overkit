# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

execute_process(
  COMMAND perl "${SOURCE_DIR}/misc/coverage-report.pl" .
  WORKING_DIRECTORY "${BINARY_DIR}/coverage"
)
