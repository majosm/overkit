// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

// Parts of code below are adapted from Google Test's source code.
// Google Test license:

// Copyright 2008, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef TESTS_MPI_PRINTER_HPP_INCLUDED
#define TESTS_MPI_PRINTER_HPP_INCLUDED

#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>
#include <ovk/core/TextProcessing.hpp>

#include <mpi.h>

#include <iostream>
#include <numeric>
#include <vector>

namespace tests {

using ovk::core::FormatNumber;
using ovk::core::StringPrint;
using ovk::core::comm;

class mpi_printer : public testing::EmptyTestEventListener {

public:

  mpi_printer():
    testing::EmptyTestEventListener(),
    Comm_(MPI_COMM_WORLD)
  {}

  // Fired before each iteration of tests starts.  There may be more than
  // one iteration if GTEST_FLAG(repeat) is set. iteration is the iteration
  // index, starting from 0.
  virtual void OnTestIterationStart(const testing::UnitTest& unit_test, int iteration) override {

    if (Comm_.Rank() == 0) {
      std::string NumTestsString = FormatNumber(unit_test.test_to_run_count(), "tests", "test");
      std::string NumTestSuitesString = FormatNumber(unit_test.test_case_to_run_count(),
        "test suites", "test suite");
      std::string NumProcsString = FormatNumber(Comm_.Size(), "processes", "process");
      std::cout << StringPrint("[==========] Running %s from %s on %s", NumTestsString,
        NumTestSuitesString, NumProcsString) << std::endl << std::flush;
    }

    MPI_Barrier(Comm_);

  }

  // Fired before environment set-up for each iteration of tests starts.
  virtual void OnEnvironmentsSetUpStart(const testing::UnitTest& unit_test) override {

    if (Comm_.Rank() == 0) {
      std::cout << "[----------] Global test environment set-up" << std::endl << std::flush;
    }

    MPI_Barrier(Comm_);

  }

  // Fired before the test case starts.
  virtual void OnTestCaseStart(const testing::TestCase& test_case) override {

    if (Comm_.Rank() == 0) {
      std::string NumTestsString = FormatNumber(test_case.test_to_run_count(), "tests", "test");
      std::cout << StringPrint("[----------] %s from %s", NumTestsString, test_case.name());
      if (!test_case.type_param()) {
        std::cout << std::endl;
      } else {
        std::cout << StringPrint(", where TypeParam = %s\n", test_case.type_param());
      }
      std::cout << std::flush;
    }

    MPI_Barrier(Comm_);

  }

  // Fired before the test starts.
  virtual void OnTestStart(const testing::TestInfo& test_info) override {

    if (Comm_.Rank() == 0) {
      std::cout << StringPrint("[ RUN      ] %s.%s", test_info.test_case_name(), test_info.name());
      std::cout << std::endl;
      std::cout << std::flush;
    }

    MPI_Barrier(Comm_);

  }

  // Fired after a failed assertion or a SUCCEED() invocation.
  // If you want to throw an exception from this function to skip to the next
  // TEST, it must be AssertionException defined above, or inherited from it.
  virtual void OnTestPartResult(const testing::TestPartResult& test_part_result) override {

    // If the test part succeeded, we don't need to do anything.
    if (test_part_result.type() == testing::TestPartResult::kSuccess)
      return;

    TestResults_.push_back(test_part_result);

  }

  // Fired after the test ends.
  virtual void OnTestEnd(const testing::TestInfo& test_info) override {

    char TestPassedLocal = test_info.result()->Passed();
    std::vector<char> TestPassedAll(Comm_.Size());
    MPI_Gather(&TestPassedLocal, 1, MPI_CHAR, TestPassedAll.data(), 1, MPI_CHAR, 0, Comm_);

    bool TestPassed = true;
    for (int iRank = 0; iRank < Comm_.Size(); ++iRank) {
      TestPassed = TestPassed && TestPassedAll[iRank];
    }

    if (!TestPassed) {

      std::string LocalResultString;
      if (!TestPassedLocal) {
        for (int iResult = 0; iResult < int(TestResults_.size()); ++iResult) {
          const testing::TestPartResult &test_part_result = TestResults_[iResult];
          LocalResultString += PrintTestPartResultToString(test_part_result) + '\n';
        }
      }

      if (Comm_.Rank() > 0) {
        if (!TestPassedLocal) {
          int NumChars = int(LocalResultString.length());
          MPI_Send(&NumChars, 1, MPI_INT, 0, 0, Comm_);
          MPI_Send(LocalResultString.c_str(), NumChars, MPI_CHAR, 0, 0, Comm_);
        }
      } else {
        std::string ResultString;
        if (!TestPassedLocal) {
          ResultString += "[  FAILED  ] Rank 0 produced the following error(s):\n";
          ResultString += LocalResultString;
        }
        for (int iRank = 1; iRank < Comm_.Size(); ++iRank) {
          if (!TestPassedAll[iRank]) {
            std::vector<char> RemoteResultChars;
            int NumChars;
            MPI_Recv(&NumChars, 1, MPI_INT, iRank, 0, Comm_, MPI_STATUS_IGNORE);
            RemoteResultChars.resize(NumChars);
            MPI_Recv(RemoteResultChars.data(), NumChars, MPI_CHAR, iRank, 0, Comm_,
              MPI_STATUS_IGNORE);
            ResultString += StringPrint("[  FAILED  ] Rank %i produced the following error(s):\n",
              iRank);
            ResultString.append(RemoteResultChars.begin(), RemoteResultChars.end());
          }
        }
        PrintTestPartResult(ResultString);
      }

    }

    MPI_Barrier(Comm_);

    if (Comm_.Rank() == 0) {
      if (TestPassed) {
        std::cout << "[       OK ] ";
      } else {
        std::cout << "[  FAILED  ] ";
      }
      std::cout << StringPrint("%s.%s", test_info.test_case_name(), test_info.name());
      if (!TestPassed) {
        PrintFullTestCommentIfPresent(test_info);
      }
    }

    // Don't know how to query this externally
    // if (GTEST_FLAG(print_time)) {
    if (true) {
      int ElapsedTimeMs = int(test_info.result()->elapsed_time());
      int MaxElapsedTimeMs;
      MPI_Reduce(&ElapsedTimeMs, &MaxElapsedTimeMs, 1, MPI_INT, MPI_MAX, 0, Comm_);
      if (Comm_.Rank() == 0) {
        std::string MaxElapsedTimeString = StringPrint("%i ms", MaxElapsedTimeMs);
        std::cout << StringPrint(" (%s)", MaxElapsedTimeString);
        std::cout << std::endl;
      }
    } else {
      if (Comm_.Rank() == 0) {
        std::cout << std::endl;
      }
    }

    if (Comm_.Rank() == 0) {
      std::cout << std::flush;
    }

    TestResults_.clear();

    MPI_Barrier(Comm_);

  }

  // Fired after the test case ends.
  virtual void OnTestCaseEnd(const testing::TestCase& test_case) override {

    int ElapsedTimeMs = int(test_case.elapsed_time());

    int MaxElapsedTimeMs;
    MPI_Reduce(&ElapsedTimeMs, &MaxElapsedTimeMs, 1, MPI_INT, MPI_MAX, 0, Comm_);

    if (Comm_.Rank() == 0) {
      std::string NumTestsString = FormatNumber(test_case.test_to_run_count(), "tests", "test");
      std::string MaxElapsedTimeString = StringPrint("%i ms", MaxElapsedTimeMs);
      std::cout << StringPrint("[----------] %s from %s (%s total)", NumTestsString,
        test_case.name(), MaxElapsedTimeString) << std::endl << std::flush;
    }

    MPI_Barrier(Comm_);

  }

  // Fired before environment tear-down for each iteration of tests starts.
  virtual void OnEnvironmentsTearDownStart(const testing::UnitTest& unit_test) override {

    if (Comm_.Rank() == 0) {
      std::cout << "[----------] Global test environment tear-down" << std::endl << std::flush;
    }

    MPI_Barrier(Comm_);

  }

  // Fired after each iteration of tests finishes.
  virtual void OnTestIterationEnd(const testing::UnitTest& unit_test,
                                  int iteration) override {

    if (Comm_.Rank() == 0) {
      std::string NumTestsString = FormatNumber(unit_test.test_to_run_count(), "tests", "test");
      std::string NumTestSuitesString = FormatNumber(unit_test.test_case_to_run_count(), "test suites",
        "test suite");
      std::cout << StringPrint("[==========] %s from %s ran", NumTestsString, NumTestSuitesString);
    }

    // Don't know how to query this externally
    // if (GTEST_FLAG(print_time)) {
    if (true) {
      int ElapsedTimeMs = int(unit_test.elapsed_time());
      int MaxElapsedTimeMs;
      MPI_Reduce(&ElapsedTimeMs, &MaxElapsedTimeMs, 1, MPI_INT, MPI_MAX, 0, Comm_);
      if (Comm_.Rank() == 0) {
        std::string MaxElapsedTimeString = StringPrint("%i ms", MaxElapsedTimeMs);
        std::cout << StringPrint(" (%s total)", MaxElapsedTimeString);
        std::cout << std::endl;
      }
    } else {
      if (Comm_.Rank() == 0) {
        std::cout << std::endl;
      }
    }

    int NumTestSuites = unit_test.total_test_case_count();

    std::vector<char> LocalTestResults;
    for (int i = 0; i < NumTestSuites; ++i) {
      const testing::TestCase& test_case = *unit_test.GetTestCase(i);
      if (!test_case.should_run()) continue;
      int NumTests = test_case.total_test_count();
      for (int j = 0; j < NumTests; ++j) {
        const testing::TestInfo& test_info = *test_case.GetTestInfo(j);
        if (!test_info.should_run()) continue;
        LocalTestResults.push_back(test_info.result()->Passed());
      }
    }

    int TotalTests = LocalTestResults.size();

    std::vector<char> TestResults;

    if (Comm_.Rank() > 0) {
      MPI_Reduce(LocalTestResults.data(), NULL, TotalTests, MPI_CHAR, MPI_LAND, 0, Comm_);
    } else {
      TestResults.resize(TotalTests);
      MPI_Reduce(LocalTestResults.data(), TestResults.data(), TotalTests, MPI_CHAR, MPI_LAND, 0,
        Comm_);
    }

    if (Comm_.Rank() == 0) {

      int NumPassedTests = std::accumulate(TestResults.begin(), TestResults.end(), 0);
      int NumFailedTests = TotalTests - NumPassedTests;

      std::string NumPassedTestsString = FormatNumber(NumPassedTests, "tests", "test");
      std::cout << "[  PASSED  ] " << NumPassedTestsString << std::endl;

      if (NumFailedTests > 0) {
        std::string NumFailedTestsString1 = FormatNumber(NumFailedTests, "tests", "test");
        std::cout << "[  FAILED  ] " << NumFailedTestsString1 << ", listed below:";
        std::cout << std::endl;
        int iNextTest = 0;
        for (int i = 0; i < NumTestSuites; ++i) {
          const testing::TestCase& test_case = *unit_test.GetTestCase(i);
          if (!test_case.should_run()) continue;
          int NumTests = test_case.total_test_count();
          for (int j = 0; j < NumTests; ++j) {
            const testing::TestInfo& test_info = *test_case.GetTestInfo(j);
            if (!test_info.should_run()) continue;
            if (!TestResults[iNextTest]) {
              std::cout << StringPrint("[  FAILED  ] %s.%s", test_case.name(), test_info.name());
              PrintFullTestCommentIfPresent(test_info);
              std::cout << std::endl;
            }
            ++iNextTest;
          }
        }
        std::cout << std::endl;
        std::string NumFailedTestsString2 = FormatNumber(NumFailedTests, "FAILED TESTS",
          "FAILED TEST");
        std::cout << NumFailedTestsString2 << std::endl;
      }
      int NumDisabledTests = unit_test.reportable_disabled_test_count();
      // Don't know how to query this externally
      // if (NumDisabledTests > 0 && !GTEST_FLAG(also_run_disabled_tests)) {
      if (NumDisabledTests > 0) {
        if (NumFailedTests == 0) {
          std::cout << std::endl;  // Add a spacer if no FAILURE banner is displayed.
        }
        std::string NumDisabledTestsString = FormatNumber(NumDisabledTests, "DISABLED TESTS",
          "DISABLED TEST");
        std::cout << "  YOU HAVE " << NumDisabledTestsString << std::endl << std::endl;
      }
      std::cout << std::flush;

    }

    MPI_Barrier(Comm_);

  }

private:

  comm Comm_;
  std::vector<testing::TestPartResult> TestResults_;

  void PrintFullTestCommentIfPresent(const testing::TestInfo& test_info) {
    const char* const type_param = test_info.type_param();
    const char* const value_param = test_info.value_param();

    if (type_param != NULL || value_param != NULL) {
      std::cout << ", where ";
      if (type_param != NULL) {
        std::cout << StringPrint("TypeParam = %s", type_param);
        if (value_param != NULL)
          std::cout << " and ";
      }
      if (value_param != NULL) {
        std::cout << StringPrint("GetParam() = %s", value_param);
      }
    }
  }

  // Formats a source file path and a line number as they would appear
  // in an error message from the compiler used to compile this code.
  std::string FormatFileLocation(const char* file, int line) {
    static const char kUnknownFile[] = "unknown file";
    const std::string file_name(file == NULL ? kUnknownFile : file);
    if (line < 0) {
      return file_name + ":";
    }
  #ifdef _MSC_VER
    return file_name + "(" + StringPrint("%i", line) + "):";
  #else
    return file_name + ":" + StringPrint("%i", line) + ":";
  #endif  // _MSC_VER
  }

  // Converts a TestPartResult::Type enum to human-friendly string
  // representation.  Both kNonFatalFailure and kFatalFailure are translated
  // to "Failure", as the user usually doesn't care about the difference
  // between the two when viewing the test result.
  const char * TestPartResultTypeToString(testing::TestPartResult::Type type) {
    switch (type) {
      case testing::TestPartResult::kSuccess:
        return "Success";
      case testing::TestPartResult::kNonFatalFailure:
      case testing::TestPartResult::kFatalFailure:
#ifdef _MSC_VER
        return "error: ";
#else
        return "Failure\n";
#endif
      default:
        return "Unknown result type";
    }
  }

  // Prints a TestPartResult to an std::string.
  std::string PrintTestPartResultToString(const testing::TestPartResult& test_part_result) {
    return (testing::Message()
            << FormatFileLocation(test_part_result.file_name(), test_part_result.line_number())
            << " " << TestPartResultTypeToString(test_part_result.type())
            << test_part_result.message()).GetString();
  }

  // Prints a TestPartResult.
  void PrintTestPartResult(const std::string &ResultString) {
    std::cout << ResultString << std::flush;
    // If the test program runs in Visual Studio or a debugger, the
    // following statements add the test part result message to the Output
    // window such that the user can double-click on it to jump to the
    // corresponding source code location; otherwise they do nothing.
#if GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE
    // We don't call OutputDebugString*() on Windows Mobile, as printing
    // to stdout is done by OutputDebugString() there already - we don't
    // want the same message printed twice.
    OutputDebugStringA(ResultString.c_str());
#endif
  }

};

class mpi_test : public testing::Test {

public:

  mpi_test() = default;

  virtual void SetUp() override {
    Comm_ = comm(MPI_COMM_WORLD);
  }

  virtual void TearDown() override {
    Comm_.Reset();
  }

  const comm &TestComm() { return Comm_; }

private:

  comm Comm_;

};

}

#endif
