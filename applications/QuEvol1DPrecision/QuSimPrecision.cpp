
#include <QuSim.h>
#include <chrono>
#include "../QuBenchmark/Benchmark.h"

struct Test {
    SolverMethod met;
    char const *name;
    Options opts;
};

void test_tunneling()
{

    printf("Test the tunneling probability\n");
    printf("Inital wave function: C*gauss(x, -20, 5)*exp(I*x)\n");
    printf("Potential: exp(-x*x)\n");

    Options space_O2;
    space_O2.SpaceOrder(2);

    Options space_O4;
    space_O2.SpaceOrder(4);

    Options space_O6;
    space_O6.SpaceOrder(6);

    Test tests[] = {
    { SolverMethod::SplittingMethodO2 , "splitO2", Options() },
    { SolverMethod::SplittingMethodO4 , "splitO4", Options() },
    { SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", space_O2 },
    { SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
    { SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO6", space_O6 },
    { SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
    { SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", space_O4 },
    { SolverMethod::GaussLegendreO4 , "gaussO4+spaceO6", space_O6 },
    { SolverMethod::GaussLegendreO6 , "gaussO6+spaceO2", space_O2 },
    { SolverMethod::GaussLegendreO6 , "gaussO6+spaceO4", space_O4 },
    { SolverMethod::GaussLegendreO6 , "gaussO6+spaceO6", space_O6 },
    { SolverMethod::Eigen , "eigen", Options() },
    };

    int dims[] = {
        1000,
        2000,
        3000,
        4000,
        //5000,
        //6000,
    };

    double times[] = {
        1E-2,
        2E-2,
        3E-2,
        4E-2,
    };


    if (0) {
        Evolver1D syst;
        syst.init(FunctorWrapper("gauss(x, -20, 5)*exp(I*x)"), true, 5E-3, false, FunctorWrapper("exp(-x*x)"),
            -100, 100, 8000, BoundaryCondition::Period, SolverMethod::SplittingMethodO4, 1, 1,
            Options());
        for (; syst.Time() < 64;) {
            syst.step();
        }
        printf("ref 1 %.16f\n", syst.Norm2());
    }

    if (0) {
        Evolver1D syst;
        syst.init(FunctorWrapper("gauss(x, -20, 5)*exp(I*x)"), true, 5E-3, false, FunctorWrapper("exp(-x*x)"),
            -100, 100, 8000, BoundaryCondition::Period, SolverMethod::GaussLegendreO6, 1, 1,
            Options());
        for (; syst.Time() < 64;) {
            syst.step();
        }
        printf("ref 1 %.16f\n", syst.Norm2());

    }
    double ref = (0.1593774409 + 0.1593774408) / 2;

    for (int i = 0; i < sizeof(tests) / sizeof(Test); ++i) {
        begin_section(tests[i].name);



        printf("%-10s ", "Dim\\Dt");
        for (int k = 0; k < sizeof(times) / sizeof(double); ++k) {
            printf("%10.3f ", times[k]);
        }
        printf("\n");


        int n = sizeof(dims) / sizeof(int);
        if (tests[i].met == SolverMethod::Eigen) n = 1;
        for (int j = 0; j < n; ++j) {
            int dim = dims[j];

            printf("%-10d ", dim);
            for (int k = 0; k < sizeof(times) / sizeof(double); ++k) {
                double time = times[k];

                Evolver1D syst;
                syst.init(FunctorWrapper("gauss(x, -20, 5)*exp(I*x)"), true, time, false, FunctorWrapper("exp(-x*x)"),
                    -100, 100, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
                    tests[i].opts);

                auto t0 = std::chrono::system_clock::now();
                for (; syst.Time() < 64;) {
                    syst.step();
                }
                auto t1 = std::chrono::system_clock::now();

                auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

                auto x0 = syst.GetX0();
                auto x1 = syst.GetX1();
                printf("%+10.1E ", syst.Norm2(0.5*(x0+x1), x1) - ref);
            }
            printf("\n");
        }
        printf("\n");

        end_section();
    }

}

void test_speed()
{
    printf("Test the Speed\n");
    printf("Inital wave function: C*gauss(x, 0, 5)*exp(I*x)\n");
    printf("Potential: 0\n");

    Options space_O2;
    space_O2.SpaceOrder(2);

    Options space_O4;
    space_O4.SpaceOrder(4);

    Options space_O6;
    space_O6.SpaceOrder(6);

    Options cuda = Options().Cuda();

    Test tests[] = {
    
#ifdef QUSIM_USE_CUDA
    { SolverMethod::SplittingMethodO2 , "splitO2+cuda", cuda },
    { SolverMethod::SplittingMethodO4 , "splitO4+cuda", cuda },
#endif
    { SolverMethod::SplittingMethodO2 , "splitO2+kiss", Options() },
    { SolverMethod::SplittingMethodO4 , "splitO4+kiss", Options() },
    { SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", space_O2 },
    { SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
    { SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO6", space_O6 },
    { SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
    { SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", space_O4 },
    { SolverMethod::GaussLegendreO4 , "gaussO4+spaceO6", space_O6 },
    { SolverMethod::GaussLegendreO6 , "gaussO6+spaceO2", space_O2 },
    { SolverMethod::GaussLegendreO6 , "gaussO6+spaceO4", space_O4 },
    { SolverMethod::GaussLegendreO6 , "gaussO6+spaceO6", space_O6 },
    { SolverMethod::Eigen , "eigen", Options() },
    };

    int dims[] = {
        1000,
        2000,
        3000,
        4000,
        //5000,
        //6000,
    };

    double times[] = {
        1E-2,
        2E-2,
        5E-2,
        1E-1,
    };

    for (int i = 0; i < sizeof(tests) / sizeof(Test); ++i) {
        begin_section(tests[i].name);



        printf("%-10s ", "Dim\\Dt");
        for (int k = 0; k < sizeof(times) / sizeof(double); ++k) {
            printf("%10.3f ", times[k]);
        }
        printf("\n");


        int n = sizeof(dims) / sizeof(int);
        if (tests[i].met == SolverMethod::Eigen) n = 1;
        for (int j = 0; j < n; ++j) {
            int dim = dims[j];

            printf("%-10d ", dim);
            for (int k = 0; k < sizeof(times) / sizeof(double); ++k) {
                double time = times[k];

                Evolver1D syst;
                syst.init(FunctorWrapper("gauss(x, 0, 5)*exp(I*x)"), true, time, false, FunctorWrapper("0"),
                    -100, 100, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
                    tests[i].opts);

                double x0 = syst.Xavg();
                for (; ;) {

                    syst.step();
                    if (syst.Time() > 0.1) break;
                }
                double x1 = syst.Xavg();
                double vel = (x1 - x0) / syst.Time();
                printf("%+10.1E ", vel - 1.);
            }
            printf("\n");
        }
        printf("\n");
        end_section();
    }

}

int main()
{

    test_speed();
    test_tunneling();
    return 0;
}

