#define _CRT_SECURE_NO_WARNINGS

#include <QuSim.h>
#include "../src/Cal.h"
#include "../QuTest/Test.h"
void testFFT();

int main()
{
    testCal();
    testFFT();
    printf("%d tests failed\n", qu_num_failed_tests());
    return 0;
}

