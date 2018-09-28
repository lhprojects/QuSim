#include "../qsim/QuSim.h"

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

void testPerturbationConverge()
{
	std::map<std::string, std::string> opts;
	QuPerturbation3D solver;

	Real p = 1;
	Real e = p * p / 2;

	opts["order"] = "1";
	solver.init([](Real x, Real y, Real z) { return 0.01*exp(-(x * x + y * y + z * z)); },
		-100, 100, 400, -100, 100, 400, -100, 100, 400,
		e, 0,
		0, 0, 1, SolverMethod::BornSerise, 1, 1, opts);


	for (int i = 0; i < 11; ++i) {
		Real theta = 0;
		Real cosx = sin(theta);
		Real cosy = 0;
		Real cosz = cos(theta);
		Real a = solver.ComputeXSection(cosx, cosy, cosz);
		Real b = born(1, 1, p, 0.01, 1, cosx, cosy, cosz);
		printf("%d | %10.5E %10.5E %10.5E\n", i, a, b, (a - b) / b);

		solver.Compute();
	}

}

void test0()
{
	std::map<std::string, std::string> opts;
	QuScatteringInverseMatrix3D solver;
	Real p = 1;
	Real e = p * p / 2;
	solver.init([](Real x, Real y, Real z) { return 0.01*exp(-(x*x + y * y + z * z)); },
		-10, 10, 100, -10, 10, 100, -10, 10, 100,
		e, 0, 0, 1, SolverMethod::MatrixInverse, 1, 1, opts);

	//solver.Compute();

	for (int i = 0; i < 11; ++i) {
		Real theta = (0.0 + i * 0.1) * Pi;
		Real cosx = sin(theta);
		Real cosy = 0;
		Real cosz = cos(theta);
		Real a = solver.ComputeXSection(cosx, cosy, cosz);
		Real b = born(1, 1, p, 0.01, 1, cosx, cosy, cosz);
		printf("%5.3E | %10.5E %10.5E %10.5E\n", theta, a, b, (a-b)/b);
	}

}


void test1()
{
	std::map<std::string, std::string> opts;
	opts["matrix_solver"] = "BiCGSTAB";
	//opts["preconditioner"] = "DiagonalPreconditioner";
	opts["preconditioner"] = "IdentityPreconditioner";
	QuScatteringInverseMatrix3D solver;
	Real p = 1;
	Real e = p * p / 2;
	solver.init([](Real x, Real y, Real z) { return 0.01*exp(-(x*x + y * y + z * z)); },
		-25, 25, 100, -25, 25, 100, -25, 25, 100,
		e, 0, 0, 1, SolverMethod::MatrixInverse, 1, 1, opts);

	solver.Compute();

	for (int i = 0; i < 11; ++i) {
		Real theta = (0.0 + i * 0.1) * Pi;
		Real cosx = sin(theta);
		Real cosy = 0;
		Real cosz = cos(theta);
		Real a = solver.ComputeXSection(cosx, cosy, cosz);
		Real b = born(1, 1, p, 0.01, 1, cosx, cosy, cosz);
		printf("%5.3E | %10.5E %10.5E %10.5E\n", theta, a, b, (a - b) / b);
	}

}

int main()
{
	testPerturbationConverge();
	test0();
	test1();
}
