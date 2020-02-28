// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "tests/MPIPrinter.hpp"

#include "gtest/gtest.h"

#include <ovk/core/TextProcessing.hpp>

#include <mpi.h>

#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

using ovk::core::StringPrint;

int main(int argc, char *argv[]) {

  // Fork a child process to check for deadlock (thanks Alex)
  pid_t ChildProcessID = fork();

  if (ChildProcessID != 0) {

    testing::InitGoogleTest(&argc, argv);

    MPI_Init(&argc, &argv);

    auto& Listeners = testing::UnitTest::GetInstance()->listeners();

    auto *DefaultPrinter = Listeners.default_result_printer();
    Listeners.Release(DefaultPrinter);
    delete DefaultPrinter;

    auto *MPIPrinter = new tests::mpi_printer();
    Listeners.Append(MPIPrinter);

    int Result = RUN_ALL_TESTS();

    // No deadlock detected; kill child process before it kills us
    kill(ChildProcessID, SIGKILL);

    Listeners.Release(MPIPrinter);
    delete MPIPrinter;

    MPI_Finalize();

    return Result;

  } else {

    char *TimeLimitString = getenv("OVK_TEST_TIME_LIMIT");

    int TimeLimit;
    if (TimeLimitString) {
      TimeLimit = std::stoi(TimeLimitString);
    } else {
      // Default to 5 minutes
      TimeLimit = 300;
    }

    sleep(TimeLimit);

    // Deadlock (maybe) detected; kill parent process
    std::cout << "[  FAILED  ] Test time limit exceeded" << std::endl << std::flush;
    kill(getppid(), SIGKILL);

    return 1;

  }

}
