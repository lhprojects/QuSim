#include <QuSim.h>
#include "../QuBenchmark/Benchmark.h"

#ifdef QUSIM_USE_CUDA
static const bool QuUseCuda = true;
#else
static const bool QuUseCuda = false;
#endif

// x-section scatter by born
double born(double mass, double hbar, double p,
	double v0, double alpha,
	double cosx, double cosy, double cosz)
{
	double kx = cosx *p / hbar;
	double ky = cosy * p / hbar;
	double kz = (cosz-1) * p / hbar;
	double mat = v0 * pow(Pi/alpha, 1.5) * exp(-(kx*kx + ky*ky + kz*kz) / (4 * alpha));
	return mass * mass / (2*Pi*2*Pi*hbar*hbar*hbar*hbar) * mat*mat;

}

void testXsectionCal()
{
    begin_section("test XsectionCal");
    printf("Diff X-section\n");
	Options opts;
	QuScatteringInverseMatrix3D solver;
	Real p = 1;
	Real e = p * p / 2;
	solver.init([](Real x, Real y, Real z) { return 0.01 * exp(-(x * x + y * y + z * z)); },
		-10, 10, 100, -10, 10, 100, -10, 10, 100,
		e, 0, 0, 1, SolverMethod::MatrixInverse, 1, 1, opts);

	for (int i = 0; i < 11; ++i) {
		Real theta = (0.0 + i * 0.1) * Pi;
        Real cosx = sin(theta) * cos(0.1234);
        Real cosy = sin(theta) * sin(0.1234);
		Real cosz = cos(theta);
		Real a = solver.ComputeXSection(cosx, cosy, cosz);
		Real b = born(1, 1, p, 0.01, 1, cosx, cosy, cosz);
		if (i == 0) {
			printf("%6s | %12s %12s %12s\n", "theta", "NaivePer", "AnaBorn", "RelDiff");
		}
		printf("%6.3E | %12.5E %12.5E %12.5E\n", theta, a, b, (a - b) / b);
	}

	end_section();

}

#include <thread>
std::vector<double> fil(QuScatteringProblemSolver3D &solver)
{
	std::vector<double> as;
	for (int i = 0; i < 11; ++i) {
        Real theta = (0.0 + i * 0.1) * Pi;
        Real cosx = sin(theta);
        Real cosy = 0;
        Real cosz = cos(theta);
        as.push_back(solver.ComputeXSection(cosx, cosy, cosz));
	}
	return as;
}

void testInverseMatrix(bool cuda = true)
{
	if (!QuUseCuda) cuda = false;

	begin_section("test Inverse Matrix");
	double v0 = 1.0;
	Real p = 1;
	Real e = p * p / 2;

	std::vector<double> mat_x;
	std::vector<double> per1_x;
	std::vector<double> per10_x;
	std::vector<double> per20_x;

	std::vector<double> perp10_x;
	std::vector<double> perp20_x;

	{
		QuScatteringInverseMatrix3D mat;
		Options opts;
		opts.MatrixSolverBiCGSTAB().IdentityPreconditioner();
		opts.MatrixSolverBiCGSTABMaxIters(1E-4);
		opts.CpuParVec();
        mat.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
            -20, 20, 200, -20, 20, 200, -20, 20, 200,
            e, 0, 0, 1, SolverMethod::MatrixInverse, 1, 1, opts);

		mat.Compute();
		mat_x = fil(mat);
	}

	{
		QuScatteringInverseMatrix3D mat;
		Options opts;
		opts.SpaceOrder(4);
		opts.MatrixSolverBiCGSTAB().IdentityPreconditioner();
		opts.MatrixSolverBiCGSTABMaxIters(1E-4);
		opts.CpuParVec();
		mat.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-20, 20, 200, -20, 20, 200, -20, 20, 200,
			e, 0, 0, 1, SolverMethod::MatrixInverse, 1, 1, opts);

		mat.Compute();
		per1_x = fil(mat);
	}

	{
		QuScatteringInverseMatrix3D mat;
		Options opts;
		opts.MatrixSolverBiCGSTAB().IdentityPreconditioner();
		opts.MatrixSolverBiCGSTABMaxIters(1E-4);
		opts.CpuParVec();
		mat.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-20, 20, 200, -20, 20, 200, -20, 20, 200,
			e, 0, 0, 1, SolverMethod::MatrixInverse, 1, 1, opts);

		mat.Compute();
		per10_x = fil(mat);
	}
	{
		QuScatteringInverseMatrix3D mat;
		Options opts;
		opts.SpaceOrder(4);
		opts.MatrixSolverBiCGSTAB().IdentityPreconditioner();
		opts.MatrixSolverBiCGSTABMaxIters(1E-4);
		opts.CpuParVec();
		mat.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-20, 20, 200, -20, 20, 200, -20, 20, 200,
			e, 0, 0, 1, SolverMethod::MatrixInverse, 1, 1, opts);

		mat.Compute();
		per20_x = fil(mat);
	}

	{
		QuPerturbation3D perp20;
		Options per_opts;
		per_opts.Preconditional(true).VellekoopPreconditioner();
		per_opts.Order(256);
		per_opts.Slow(0.5);
		if (cuda) per_opts.Cuda();
		perp20.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-20, 20, 400, -20, 20, 400, -20, 20, 400,
			e, 0.1,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, per_opts);
		perp20.Compute();
		perp10_x = fil(perp20);
	}
	{
		QuPerturbation3D perp20;
		Options per_opts;
		per_opts.Preconditional(true).VellekoopPreconditioner();
		per_opts.Order(512);
		per_opts.Slow(0.5);
		if (cuda) per_opts.Cuda();
		perp20.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-20, 20, 400, -20, 20, 400, -20, 20, 400,
			e, 0.1,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, per_opts);
		perp20.Compute();
		perp20_x = fil(perp20);
	}

	for (int i = 0; i < 11; ++i) {
		Real theta = (0.0 + i * 0.1) * Pi;
		Real cosx = sin(theta);
		Real cosy = 0;
		Real cosz = cos(theta);
		Real b = born(1, 1, p, v0, 1, cosx, cosy, cosz);

		if (i == 0) {
			printf("%10s | %12s %12s %12s %12s %12s %12s\n",
				"theta", "inv.mat", "inv.mat.mi", "inv.mat.md", "inv.mat.mid",
				"perpO20", "perpO40");
		}
        printf("%10.3E | %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
            theta, mat_x[i], per1_x[i], per10_x[i], per20_x[i],
            perp10_x[i], perp20_x[i]);
	}

	end_section();
}

void testNaiveBorn(double width, int grids, bool cuda = true)
{
	double v0 = 1.0;
	begin_section("test Naive Born");
    printf("width = %7.1f, n = %4d, v0 = %lf, cuda = %d\n", width, grids, v0, (int)cuda);

	Real p = 1;
	Real e = p * p / 2;

	std::vector<double> per1_x;
	std::vector<double> per10_x;
	std::vector<double> per20_x;

	std::vector<double> perp20_x;
	std::vector<double> perp40_x;


    {
        QuPerturbation3D per1;
        Options per_opts;
        per_opts.Order(1);
        if (cuda) per_opts.Cuda();
        per1.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
            -width, width, grids, -width, width, grids, -width, width, grids,
            e, 0.05,
            0, 0, 1, SolverMethod::BornSerise, 1, 1, per_opts);
        per1.Compute();
        per1_x = fil(per1);
    }

    {
		QuPerturbation3D per10;
		Options per_opts;
        per_opts.Order(10);
		if (cuda) per_opts.Cuda();
        per10.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-width, width, grids, -width, width, grids, -width, width, grids,
			e, 0.05,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, per_opts);
        per10.Compute();
		per10_x = fil(per10);
	}

	{
		QuPerturbation3D per20;
		Options per_opts;
		per_opts.Order(20);
		if (cuda) per_opts.Cuda();
		per20.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-width, width, grids, -width, width, grids, -width, width, grids,
			e, 0.05,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, per_opts);
		per20.Compute();
		per20_x = fil(per20);
	}

    {
        QuPerturbation3D perp20;
        Options per_opts;
        per_opts.Preconditional(true).VellekoopPreconditioner();
        per_opts.Order(256);
		per_opts.Slow(0.5);
		if (cuda) per_opts.Cuda();
        perp20.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-width, width, grids, -width, width, grids, -width, width, grids,
			e, 0.05,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, per_opts);
        perp20.Compute();
        perp20_x = fil(perp20);
    }
    {
        QuPerturbation3D perp20;
        Options per_opts;
        per_opts.Preconditional(true).VellekoopPreconditioner();
        per_opts.Order(512);
		per_opts.Slow(0.5);
        if (cuda) per_opts.Cuda();
        perp20.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-width, width, grids, -width, width, grids, -width, width, grids,
			e, 0.1,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, per_opts);
        perp20.Compute();
        perp40_x = fil(perp20);
    }

	for (int i = 0; i < 11; ++i) {
        Real theta = (0.0 + i * 0.1) * Pi;

		if (i == 0) {
            printf("%10s | %12s %12s %12s %12s %12s\n",
                "theta", "per1", "per10", "per20", "perp20", "perp40");
        }
        printf("%10.3E | %12.5E %12.5E %12.5E %12.5E %12.5E\n",
            theta, per1_x[i], per10_x[i], per20_x[i], 
            perp20_x[i], perp40_x[i]);
	}

	end_section();
}

std::vector<double> fill(QuPerturbation3D &solver, size_t n)
{
	std::vector<double> as;
	for (int i = 0; i < n; ++i) {
        as.push_back(solver.ComputeXSection(0, 0, 1));
        solver.Compute();
	}
	return as;
}

void testPerburbationConverge(Real v0, int n = 100, bool cuda = true)
{
#ifndef QUSIM_USE_CUDA
	return;
#endif

	begin_section("Test Perburbation Convege");
	printf("V0 = %lf\n", v0);
	auto const e = 0.5;

	std::vector<double> per_nav;
	std::vector<double> per_vel;
	std::vector<double> per_hao1;
	std::vector<double> per_hao2;
	int order = 5;


	{
		QuPerturbation3D solver;
		Options opts;
		opts.Preconditional(true).Hao1Preconditioner().Order(order);
		opts.Slow(0.5);
		if (cuda) opts.Cuda();
		solver.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-20, 20, 400, -20, 20, 400, -20, 20, 400,
			e, 0.05,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, opts);
		per_hao1 = fill(solver, n);
	}

	{
		QuPerturbation3D solver;
		Options opts;
		opts.Preconditional(true).Hao2Preconditioner().Order(order);
		opts.Slow(0.5);
		if (cuda) opts.Cuda();
		solver.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-20, 20, 400, -20, 20, 400, -20, 20, 400,
			e, 0.05,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, opts);
		per_hao2 = fill(solver, n);
	}


	{
		QuPerturbation3D solver;
		Options opts;
		opts.Preconditional(false).Order(order);
		opts.Slow(0.5);
		if (cuda) opts.Cuda();
		solver.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-30, 30, 400, -30, 30, 400, -30, 30, 400,
			e, 0.05,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, opts);
		per_nav = fill(solver, n);
	}

	{
		QuPerturbation3D solver;
		Options opts;
		opts.Preconditional(true).VellekoopPreconditioner().Order(order);
		opts.Slow(0.5);
		if (cuda) opts.Cuda();
		solver.init([=](Real x, Real y, Real z) { return v0 * exp(-(x * x + y * y + z * z)); },
			-20, 20, 400, -20, 20, 400, -20, 20, 400,
			e, 0.05,
			0, 0, 1, SolverMethod::BornSerise, 1, 1, opts);
		per_vel = fill(solver, n);
	}


	printf("%15s | %16s %16s %16s %16s\n",
		"O\\Xsect\\Met", "Naive", "Vellekoop", "Hao1", "Hao2");

	for (int i = 0; i < n; ++i) {
		Real theta = 0;
		Real cosx = cos(theta);
		Real cosy = sin(theta);
		printf("%15d | %16.10E %16.10E %16.10E %16.10E\n",
			i,
			per_nav[i],
			per_vel[i],
			per_hao1[i],
			per_hao2[i]
		);
	}
	end_section();
}

int main()
{
	try {
		testXsectionCal();


		testInverseMatrix();
		testPerburbationConverge(1);

#if defined(QUSIM_USE_CUDA)
        testNaiveBorn(20, 200, true);
        testNaiveBorn(20, 300, true);
        testNaiveBorn(30, 300, true);
#endif


	}
	catch (std::runtime_error& re) {
        printf("exception: %s\n", re.what());
        exit(1);
	}
}

