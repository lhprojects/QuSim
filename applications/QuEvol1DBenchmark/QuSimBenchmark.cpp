
#include "QuSim.h"
#include <chrono>
#include <algorithm>

struct Test {
    SolverMethod met;
    char const *name;
    Options opts;
};

int main()
{

    Options space_O2;
    space_O2.SpaceOrder(2);

    Options space_O4;
    space_O2.SpaceOrder(4);

    Options cuda;
    cuda.Cuda().Batch(100);


    Test tests[] = {
#ifdef QUSIM_USE_CUDA
    { SolverMethod::SplittingMethodO2 , "splitO2+cuda+batch=10", cuda },
    { SolverMethod::SplittingMethodO4 , "splitO4+cuda+batch=10", cuda },
#endif
    { SolverMethod::SplittingMethodO2 , "splitO2+kiss", Options() },
    { SolverMethod::SplittingMethodO4 , "splitO4+kiss", Options() },
    { SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", Options() },
    { SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
    { SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", Options() },
    { SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
    { SolverMethod::GaussLegendreO6 , "gaussO6+spaceO4", Options() },
    { SolverMethod::Eigen , "eigen", Options() },
    };

    int dims[] = {
        100,
        200,
        400,
        1000,
        10000,
        100000,
        //200000,
    };

    printf("%-30s ", "");
    printf("%20s ", "Dim");
    for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
        printf("%6d ", dims[j]);
    }
    printf("\n");

    for (int i = 0; i < sizeof(tests)/sizeof(Test); ++i) {
        printf("%-30s ", tests[i].name);
        printf("%20s ", "Iters*Dim/Time [M/s]");
        int n = sizeof(dims) / sizeof(int);
        if (tests[i].met == SolverMethod::ImplicitMidpointMethod) n = std::min(5, n);
        if (tests[i].met == SolverMethod::GaussLegendreO4) n = std::min(5, n);
        if (tests[i].met == SolverMethod::GaussLegendreO6) n = std::min(5, n);
        if (tests[i].met == SolverMethod::Eigen) n = std::min(4, n);

        using microseconds = std::chrono::duration<double, std::micro>;
        using seconds = std::chrono::duration<double, std::ratio<1,1>>;
        for (int j = 0; j < n; ++j) {
            Evolver1D syst;
            auto vf = [](Real x) -> Real { return exp(-(x - 2)*(x - 2)); };
            syst.init([](Real x) -> Complex { return exp(I*x)*exp(-x * x); },
                true, 0.01, false,
                vf,
                -10, 10, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
                tests[i].opts);

            int dim = dims[j];

            auto t0 = std::chrono::system_clock::now();
            auto t1 = std::chrono::system_clock::now();
            for (;;) {
                int iter = std::min(100*20000 / dim, 1);
                for (int it = 0; it < iter; ++it) {
                    syst.step();
                }
                t1 = std::chrono::system_clock::now();
                if (std::chrono::duration_cast<seconds>(t1 - t0).count() > 1) {
                    break;
                }
            }

            double walltime = std::chrono::duration_cast<seconds>(t1 - t0).count();
            double iters = syst.Time() / 0.01;

            if (strstr(tests[i].name, "cuda")) {
                iters *= 100;
            }
            printf("%6.2f ", (dim * iters) / walltime/1E6);
        }
        printf("\n");

    }

    return 0;
}

