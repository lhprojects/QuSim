#define _CRT_SECURE_NO_WARNINGS
#include <QuSim.h>
#include "../QuBenchmark/Benchmark.h"

void TestInverseMatrix(double p)
{
    begin_section("Test Reflaction Ratio fo Inverse Matrix Method");
    printf("V(x)=%lf exp(-x*x)\n", p);

    auto f = [&](Real x) { return p * exp(-x * x); };
    Solver1D solver;

    solver.init(f, -10, 10, 20000, 0.5, 1, I,
        SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1,
        Options().SmallRoundError(true));
    solver.Compute();

    if (p == 1) {
        double R = 1 - 0.1355284179587569045304922;
        printf("solver error %.2E\n", solver.GetR() - R);
    }

    printf("%6s %5s | %10s %10s %10s\n",
        "method", "bins",
        "SpaceO2", "SpaceO4", "SpaceO6");

    for (int i = 0; i < 10; ++i) {

        double x0 = -150;
        double x1 = 150;
        size_t n = 500 + 500 * size_t(i);


        Options space_o2 = Options().SpaceOrder(2);
        Options space_o4 = Options().SpaceOrder(4);
        Options space_o6 = Options().SpaceOrder(6);

        QuScatteringInverseMatrix1D inv1;
        inv1.init(f, x0, x1, n, 0.5, 1,
            SolverMethod::MatrixInverse, 1, 1, space_o2);
        inv1.Compute();

        QuScatteringInverseMatrix1D inv2;
        inv2.init(f, x0, x1, n, 0.5, 1,
            SolverMethod::MatrixInverse, 1, 1, space_o4);
        inv2.Compute();

        QuScatteringInverseMatrix1D inv3;
        inv3.init(f, x0, x1, n, 0.5, 1,
            SolverMethod::MatrixInverse, 1, 1, space_o6);
        inv3.Compute();

        auto geti = [&](double x) -> size_t {  return (size_t)((x - x0) / (x1 - x0) * n); };
        printf("%6s %5d | %10.2E %10.2E %10.2E\n",
            "|psi|",
            (int)n,
            abs2(inv1.GetPsi()[geti(-10)]) - solver.GetR(),
            abs2(inv2.GetPsi()[geti(-10)]) - solver.GetR(),
            abs2(inv3.GetPsi()[geti(-10)]) - solver.GetR());
        printf("%6s %5d | %10.2E %10.2E %10.2E\n",
            "GetR()",
            (int)n,
            inv1.GetR() - solver.GetR(),
            inv2.GetR() - solver.GetR(),
            inv3.GetR() - solver.GetR());
    }
    end_section();
}

void TestMethodOfInverseMatrix()
{

    begin_section("Test Method of Inversing Matrix");
    double R = 1 - 0.1355284179587569045304922;

    printf("%6s | %10s %10s %10s %10s\n",
        "Grids", "LU", "DiagonalPre", "IncompleteLUTPre", "IdentityPre");

    for (int i = 0; i < 10; ++i) {

        double x0 = -150;
        double x1 = 150;
        size_t n = 500 + 500 * size_t(i);

        Options opt1;
        opt1.SpaceOrder(6).MatrixSolverLU();

        Options opt2;
        opt2.SpaceOrder(6).MatrixSolverBiCGSTAB().DiagonalPreconditioner();

        Options opt3;
        opt3.SpaceOrder(6).MatrixSolverBiCGSTAB().IncompleteLUTPreconditioner();

        Options opt4;
        opt3.SpaceOrder(6).MatrixSolverBiCGSTAB().IdentityPreconditioner();

        QuScatteringInverseMatrix1D inv1;
        inv1.init(FunctorWrapper("1*exp(-x*x)"), x0, x1, n, 0.5, 1,
            SolverMethod::MatrixInverse, 1, 1, opt1);
        inv1.Compute();

        QuScatteringInverseMatrix1D inv2;
        inv2.init(FunctorWrapper("1*exp(-x*x)"), x0, x1, n, 0.5, 1,
            SolverMethod::MatrixInverse, 1, 1, opt2);
        inv2.Compute();

        QuScatteringInverseMatrix1D inv3;
        inv3.init(FunctorWrapper("1*exp(-x*x)"), x0, x1, n, 0.5, 1,
            SolverMethod::MatrixInverse, 1, 1, opt3);
        inv3.Compute();

        QuScatteringInverseMatrix1D inv4;
        inv4.init(FunctorWrapper("1*exp(-x*x)"), x0, x1, n, 0.5, 1,
            SolverMethod::MatrixInverse, 1, 1, opt4);
        inv4.Compute();

        printf("%6d | % 10.2E % 10.2E % 10.2E % 10.2E\n",
            (int)n,
            abs(inv1.GetR() - R),
            abs(inv2.GetR() - R),
            abs(inv3.GetR() - R),
            abs(inv4.GetR() - R));
    }

    end_section();
}

void TestNaiveBornSerise(double p)
{

    begin_section("Test NaiveBornSerise");
    printf("V(x) = %.2f * exp(-x*x)\n", p);

    auto vfunc = [=](double x) {
        return p * exp(-x * x);
    };

    Solver1D solver;
    solver.init(vfunc, -10, 10, 20000, 0.5, 1, I,
        SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1, Options());
    solver.Compute();


    QuPerturbation1D per;
    Options opts;
    opts.Order(1);
    double epsilon = 0.01;
    per.init(vfunc, -2000, 2000, 20000, 0.5, epsilon,
        1, SolverMethod::BornSerise, 1, 1, opts);

    for (int i = 0; i < 20; ++i) {


        if (i == 0) {

            printf("  Energy Gap                %lf\n", per.GetEnergyGap());
            printf("  Epsilon                   %lf\n", per.GetEpsilon());
            printf("  Energy Gap / Epsilon      %lf\n", per.GetEnergyGap() / per.GetEpsilon());
            printf("  Epsilon Decay Length      %lf\n", per.GetEpsilonDecayLength());
            printf("  Epsilon Boundary Error    %lf\n", per.GetEpsilonBoundaryError());
            printf("%5s\n", "Grids");
        }

        per.Compute();
        printf("%5d %10.3E (Exact %10.3E)\n", 1 + i, per.GetR(), solver.GetR());

    }

    end_section();
}


void testPerburbation()
{
    begin_section("Test Perturbation for various potential");
    printf("%10s | %10s %10s %10s %10s %10s %10s | %10s\n",
            "V0", "NaivO1", "NaivO2", "NaivO10", "NaivO2Sp2", "PreO5", "PreO20", "Exact");
    for (int i = 0; i < 10; ++i) {

        Real v0 = 0.05 + 0.1*i;
        auto vfunc = [&](Real x) { return v0 * exp(-x * x); };
        Solver1D solver;
        solver.init(vfunc, -10, 10, 20000, 0.5, 1, I,
            SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1, Options());
        solver.Compute();

        QuPerturbation1D per1;
        per1.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
            1, SolverMethod::BornSerise, 1, 1, Options());
        per1.Compute();


        QuPerturbation1D per2;
        Options o2 = Options().Order(2);
        per2.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
            1, SolverMethod::BornSerise, 1, 1, o2);
        per2.Compute();

        QuPerturbation1D per3;
        Options o3 = Options().Order(10);
        per3.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
            1, SolverMethod::BornSerise, 1, 1, o3);
        per3.Compute();

        QuPerturbation1D per4;
        Options o4 = Options().Order(2).SplitN(3);
        per4.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
            1, SolverMethod::BornSerise, 1, 1, o4);
        per4.Compute();

        QuPerturbation1D per5;
        Options o5 = Options().Order(5).Preconditional(true);
        per5.init(vfunc, -100, 100, 2000, 0.5, 0.51,
            1, SolverMethod::BornSerise, 1, 1, o5);
        per5.Compute();

        QuPerturbation1D per6;
        Options o6 = Options().Order(20).Preconditional(true);;
        per6.init(vfunc, -100, 100, 2000, 0.5, 0.51,
            1, SolverMethod::BornSerise, 1, 1, o6);
        per6.Compute();

        printf("%10lf | %10lf %10lf %10f %10lf %10lf %10lf | %10lf\n",
            v0, per1.GetR(), per2.GetR(), per3.GetR(), per4.GetR(), per5.GetR(), per6.GetR(),
            solver.GetR());
    }
    end_section();
}

void testPerburbativeConditioner(char const *name, Real p, int n = 100)
{
    printf("Test %f*exp(-x^2)\n", p);
    printf("Preconditioner %s\n", name);
    printf("%10s | ", "Ord\\Slow");


    Real v0 = p;
    auto vfunc = [&](Real x) { return v0 * exp(-x * x); };

    Solver1D solver;
    {

        Options opts = Options().SmallRoundError(true);
        solver.init(vfunc, -10, 10, 2000, 0.5, 1, I,
            SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1, opts);
        solver.Compute();
    }

    QuPerturbation1D per3[9];
    {
        for (int i = 0; i < 9; ++i) {


            Options opts = Options().Preconditional(true).Order(10).Slow(0.2 + 0.1*i);
            opts.SetString("preconditioner", name);
            per3[i].init(vfunc, -150, 150, 10000, 0.5, 0,
                1, SolverMethod::BornSerise, 1, 1, opts);
            printf(" %8.1f", 0.2 + 0.1*i);
        }
        printf("\n");
    }

    for (int i = 0; i < n; ++i) {

        printf("%10d | ", i * 10);

        for (int i = 0; i < 9; ++i) {
            per3[i].Compute();
            printf(" %8.1E", per3[i].GetR() - solver.GetR());
        }
        printf("\n");

    }

}


int main()
{


    TestInverseMatrix(0.1);
    TestInverseMatrix(1);

    TestMethodOfInverseMatrix();

    TestNaiveBornSerise(0.1);
    TestNaiveBornSerise(0.4);

    testPerburbation();

    testPerburbativeConditioner("Vellekoop", 0.5);
    testPerburbativeConditioner("Hao1", 0.5);
    testPerburbativeConditioner("Hao2", 0.5);

    testPerburbativeConditioner("Vellekoop", 1);
    testPerburbativeConditioner("Hao1", 1);
    testPerburbativeConditioner("Hao2", 1);


    testPerburbativeConditioner("Vellekoop", 5, 500);
    testPerburbativeConditioner("Hao1", 5, 500);
    testPerburbativeConditioner("Hao2", 5, 500);


}

