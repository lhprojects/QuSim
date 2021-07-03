
#include <QuSim.h>
#include "../QuBenchmark/Benchmark.h"

double born(double mass, double hbar, double p, double v0, double alpha, double theta)
{
    // V = v0 exp(-alpha(x*x + y*y))
    // m^2 / (2 pi p hbar^3) | exp(I delta k r) V(r) dr|^2
    double kx = (1 - cos(theta))*p / hbar;
    double ky = sin(theta) * p / hbar;
    double mat = v0 * Pi * exp(-(kx*kx + ky * ky) / (4 * alpha)) / alpha;
    return mass * mass / (2 * Pi * p * hbar*hbar*hbar) * mat*mat;

}

void testBorn()
{
    begin_section("Test Born @ first order");
    printf("V(x) = 0.001 exp(-x^2 - y^2)\n");
    Options opts;
    QuScatteringInverseMatrix2D solver;
    solver.init([](Real x, Real y) { return 0.001*exp(-x * x - y * y); }, -100, 100, 300, -100, 100, 300,
        0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
    solver.Compute();

    printf("%12s | %12s %12s | %12s\n", "theta", "inv.mat.", "ana.for.", "ratio");
    for (int i = 0; i < 10; ++i) {
        Real theta = i * 2 * Pi / 10;
        Real cosx = cos(theta);
        Real cosy = sin(theta);
        Real a = solver.ComputeXSection(cosx, cosy);
        Real b = born(1, 1, 1, 0.001, 1, theta);
        printf("%12lf | %12.8E %12.8E | %12.8E\n", theta, a, b, a / b);
    }
    end_section();
}

void testInverseAndPerburbation()
{
    Real v0 = 0.5;
    begin_section("Inverse And Perburbation");
    printf("V0 = %lf\n", v0);
    auto vf = [&](Real x, Real y) { return v0*exp(-x * x - y * y); };

    QuScatteringInverseMatrix2D solver;
    {
        Options opts;
        solver.init(vf, -100, 100, 1000, -100, 100, 1000,
            0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
        solver.Compute();
    }

    QuPerturbation2D per1;
    {
        Options opts;
        opts.Order(1);
        per1.init(vf, -100, 100, 1000, -100, 100, 1000,
            0.5, 0.1, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
        per1.Compute();
    }

    QuPerturbation2D per2;
    {
        Options opts;
        opts.Order(2);
        per2.init(vf, -100, 100, 1000, -100, 100, 1000,
            0.5, 0.1, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
        per2.Compute();
    }

    QuPerturbation2D per10;
    {
        Options opts;
        opts.Order(20);
        per10.init(vf, -100, 100, 1000, -100, 100, 1000,
            0.5, 0.1, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
        per10.Compute();
    }

    QuPerturbation2D perp10;
    {
        Options opts;
        opts.Order(20).Preconditional(true).VellekoopPreconditioner();
        perp10.init(vf, -100, 100, 1000, -100, 100, 1000,
            0.5, 0.0, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
        perp10.Compute();
    }


    printf("%12s | %12s %12s %12s %12s %12s %12s\n", "theta", "born",
        "per.1", "per.2", "per.20", "perp.20", "inv.mat.");
    for (int i = 0; i < 10; ++i) {
        Real theta = i * 2 * Pi / 10;
        Real cosx = cos(theta);
        Real cosy = sin(theta);
        printf("%12lf | %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
            theta,
            born(1, 1, 1, v0, 1, 0),
            per1.ComputeXSection(cosx, cosy),
            per2.ComputeXSection(cosx, cosy),
            per10.ComputeXSection(cosx, cosy),
            perp10.ComputeXSection(cosx, cosy),
            solver.ComputeXSection(cosx, cosy)
            );
    }
    end_section();
}


void testPerburbationConverge(Real v0, int n = 100)
{
    begin_section("Test Perburbation Convege");
    printf("V0 = %lf\n", v0);
    auto vf = [&](Real x, Real y) { return v0*exp(-x * x - y * y); };

    Real ref;
    QuScatteringInverseMatrix2D solver;
    {
        Options opts;
        opts.SpaceOrder(2);
        solver.init(vf, -100, 100, 400, -100, 100, 400,
            0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
        solver.Compute();
        ref = solver.ComputeXSection(1, 0);
        //ref = born(1, 1, 1, v0, 1, 0);
    }

    QuPerturbation2D perp0;
    {
        Options opts;
        opts.Preconditional(false).Order(1);
        opts.Slow(0.5);
        perp0.init(vf, -150, 150, 1000, -150, 150, 1000,
            0.5, 0.005, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
    }

    QuPerturbation2D perp1;
    {
        Options opts;
        opts.Preconditional(true).VellekoopPreconditioner().Order(1);
        opts.Slow(0.5);
        perp1.init(vf, -150, 150, 1000, -150, 150, 1000,
            0.5, 0.0, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
    }

    QuPerturbation2D perp2;
    {
        Options opts;
        opts.Preconditional(true).Hao1Preconditioner().Order(1);
        opts.Slow(0.5);
        perp2.init(vf, -150, 150, 1000, -150, 150, 1000,
            0.5, 0.0, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
    }

    QuPerturbation2D perp3;
    {
        Options opts;
        opts.Preconditional(true).BornIdentityPreconditioner().Order(1);
        opts.Slow(0.5);
        if (v0 > 1.5) opts.Slow(0.1);
        perp3.init(vf, -150, 150, 1000, -150, 150, 1000,
            0.5, 0.0, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
    }


    printf("%15s | %16s %16s %16s %16s %16s\n",
        "O\\Xsect\\Met", "inv.mat.", "Naive", "Vellekoop", "Hao1", "Hao2");
    for (int i = 0; i < n; ++i) {
        Real theta = 0;
        Real cosx = cos(theta);
        Real cosy = sin(theta);
        printf("%15d | %16.10E %16.10E %16.10E %16.10E %16.10E | %6.2E %6.2E %6.2E %6.2E\n",
            i,
            ref,
            perp0.ComputeXSection(cosx, cosy),
            perp1.ComputeXSection(cosx, cosy),
            perp2.ComputeXSection(cosx, cosy),
            perp3.ComputeXSection(cosx, cosy),
            perp0.GetDeltaPsiNorm(),
            perp1.GetDeltaPsiNorm(),
            perp2.GetDeltaPsiNorm(),
            perp3.GetDeltaPsiNorm()

        );
        if(v0 <= 1) perp0.Compute();
        perp1.Compute();
        perp2.Compute();
        perp3.Compute();
    }
    end_section();
}

void testInvMatVsBornAnalForDifferentPotential()
{

    begin_section("Compare inverse matrix x section with born.anal. for different potential");
    printf("%10s | %10s %10s\n", "v0", "inv.mat", "born");
    for (int i = 0; i < 10; ++i) {
        Real v0 = 0.1 + 0.1 * i;
        Options opts;
        QuScatteringInverseMatrix2D solver;
        solver.init([&](Real x, Real y) { return v0 *exp(-x * x - y * y); }, -100, 100, 300, -100, 100, 300,
            0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
        solver.Compute();

        Real a = solver.ComputeXSection(1, 0);
        Real b = born(1, 1, 1, v0, 1, 0);
        printf("%10lf| %10.5E %10.5E\n", v0, a , b);
    }
    end_section();
}

void TestTotalXSectionConvergeWithNumberOfSamplingPoints()
{
    begin_section("Test total x section converge with number of smapling points");
    Options opts;
    QuScatteringInverseMatrix2D solver;
    solver.init([](Real x, Real y) { return 1 * exp(-x * x - y * y); }, -100, 100, 300, -100, 100, 300,
        0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
    solver.Compute();

    printf("%8s %8s", "Grids", "Err.\n");
    Real ref = solver.ComputeTotalXSection(1000);
    for (int i = 0; i < 20; ++i) {
        int n = 1 + i;
        Real xs = solver.ComputeTotalXSection(n);
        printf("%8d %8.1E\n", n, xs - ref);
    }
    end_section();
}

int main()
{
    testPerburbationConverge(1);

    testInverseAndPerburbation();
    testBorn();
    testInvMatVsBornAnalForDifferentPotential();

    TestTotalXSectionConvergeWithNumberOfSamplingPoints();

    testPerburbationConverge(0.2);
    testPerburbationConverge(0.6);
    testPerburbationConverge(2);
}
