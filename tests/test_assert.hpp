#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

inline void test_fail(const std::string& message) {
  std::cerr << "TEST FAILURE: " << message << "\n";
  std::exit(1);
}

#define TEST_CHECK(cond)                 \
  do {                                   \
    if (!(cond)) {                       \
      test_fail("CHECK failed: " #cond); \
    }                                    \
  } while (false)

#define TEST_NEAR(a, b, tol)                                               \
  do {                                                                      \
    if (std::fabs(static_cast<double>(a) - static_cast<double>(b)) >        \
        static_cast<double>(tol)) {                                         \
      test_fail(std::string("NEAR failed: ") + #a + " vs " + #b);           \
    }                                                                       \
  } while (false)
