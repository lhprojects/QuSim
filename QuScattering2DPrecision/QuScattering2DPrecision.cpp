
#include "../qsim/QuSim.h"

double born(double mass, double hbar, double p, double v0, double alpha, double theta)
{
	// V = v0 exp(-alpha(x*x + y*y))
	// m^2 / (2 pi p hbar^3) | exp(I delta k r) V(r) dr|^2
	double kx = (1 - cos(theta))*p / hbar;
	double ky = sin(theta) * p / hbar;
	double mat = v0 * Pi * exp(-(kx*kx + ky * ky) / (4 * alpha)) / alpha;
	return mass * mass / (2 * Pi * p * hbar*hbar*hbar) * mat*mat;

}

void test0()
{

	std::map<std::string, std::string> opts;
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

}

void test1()
{


	for (int i = 0; i < 10; ++i) {
		Real v0 = 0.1 + 0.1 * i;
		std::map<std::string, std::string> opts;
		QuScatteringInverseMatrix2D solver;
		solver.init([&](Real x, Real y) { return v0 *exp(-x * x - y * y); }, -100, 100, 300, -100, 100, 300,
			0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
		solver.Compute();

		Real a = solver.ComputeXSection(1, 0);
		Real b = born(1, 1, 1, v0, 1, 0);
		printf("%lf %lf\n", a , b);
	}

}

void test2()
{

	std::map<std::string, std::string> opts;
	QuScatteringInverseMatrix2D solver;
	solver.init([](Real x, Real y) { return 10*exp(-x * x - y * y); }, -100, 100, 300, -100, 100, 300,
		0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
	solver.Compute();

	printf("%12s | %12s\n", "theta", "ana.for.");
	for (int i = 0; i < 10; ++i) {
		Real theta = i * 2 * Pi / 10;
		Real cosx = cos(theta);
		Real cosy = sin(theta);
		Real a = solver.ComputeXSection(cosx, cosy);
		printf("%12lf | %12.8E\n", theta, a);
	}


}


void test3()
{

	std::map<std::string, std::string> opts;
	QuScatteringInverseMatrix2D solver;
	solver.init([](Real x, Real y) { return 1 * exp(-x * x - y * y); }, -100, 100, 300, -100, 100, 300,
		0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
	solver.Compute();

	Real ref = solver.ComputeTotalXSection(1000);
	for (int i = 0; i < 20; ++i) {
		int n = 1 + i;
		Real xs = solver.ComputeTotalXSection(n);
		printf("%d %8.1E\n", n, xs - ref);
	}


}

int main()
{
	test3();
	test2();
	test1();
	test0();

}
