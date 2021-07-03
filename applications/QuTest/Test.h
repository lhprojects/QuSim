#pragma once

inline int& qu_num_failed_tests()
{
    static int num_failed_tests;
    return num_failed_tests;
}

inline void qu_test_assert(char const* x,
    char const *file,
    int line, bool v)
{
    if (!v) {
        printf("failed (%s:%d): %s\n", file, line, x);
        qu_num_failed_tests()++;
    }
}
inline void qu_test_less(char const* x,
    char const* file,
    int line, double v, double tol)
{
    if (!(v < tol)) {
        printf("failed (%s:%d): %s(%8.2E) < %8.2E\n", file, line, x, v, tol);
        qu_num_failed_tests()++;
    }
}

#define QU_TEST_ASSERT(x) qu_test_assert(#x, __FILE__, __LINE__, (x))
#define QU_TEST_LESS(x, tol) qu_test_less(#x, __FILE__, __LINE__, (x), tol)
