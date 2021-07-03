#ifndef QUSIM_BENCHMARK
#define QUSIM_BENCHMARK

void begin_section(std::string sec)
{
    printf("%*s\n", (int)std::max(0, 20 + (int)sec.length() / 2), sec.c_str());
    printf("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
}

void end_section()
{
    printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
    printf("                                                 \n");
}

#endif
