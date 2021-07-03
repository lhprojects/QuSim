
#include <QuSim.h>
#include <chrono>

struct Test {
    SolverMethod met;
    char const *name;
    Options opts;
};

int main()
{
    Options smal_round_err;
    smal_round_err.SmallRoundError(true);

    Test tests[] = {
        { SolverMethod::ImplicitMidpointMethod , "midpoint", Options() },
    { SolverMethod::ImplicitMidpointMethod , "midpoint-", smal_round_err },
    { SolverMethod::ExplicitRungeKuttaO4Classical , "rk4", Options() },
    { SolverMethod::ExplicitRungeKuttaO4Classical , "rk4-", smal_round_err },
    { SolverMethod::GaussLegendreO4 , "glo4", Options() },
    { SolverMethod::GaussLegendreO4 , "glo4-", smal_round_err },
    { SolverMethod::ExplicitRungeKuttaO6Luther1967 , "rko6", Options() },
    { SolverMethod::ExplicitRungeKuttaO6Luther1967 , "rko6-", smal_round_err },
    };

    int dims[] = {
        2000,
        10000,
        20000,
        100000,
        200000,
        1000000,
    };

    printf("%-10s ", "");
    printf("%15s ", "Dim");
    for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
        printf("%6d ", dims[j]);
    }
    printf("\n");

    for (int i = 0; i < sizeof(tests) / sizeof(Test); ++i) {
        printf("%-10s ", tests[i].name);
        printf("%15s ", "Dim/Time [M/s]");
        int n = sizeof(dims) / sizeof(int);
        for (int j = 0; j < n; ++j) {
            int dim = dims[j];

            Solver1D syst;
            syst.init(FunctorWrapper("exp(-x*x)"), -10, 10, dim, 0.5, 1, I, tests[i].met, 1, 1, tests[i].opts);


            auto t0 = std::chrono::system_clock::now();
            for(int i = 0; i < 50; ++i)
                syst.Compute();
            auto t1 = std::chrono::system_clock::now();

            auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0) / 10;
            printf("%6.1f ", 50.0*dim / d.count());
        }
        printf("\n");

    }

    return 0;
}

