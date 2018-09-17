
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


}

int main()
{


	std::map<std::string, std::string> opts;
	QuScatteringInverseMatrix2D solver;
	solver.init([](Real x, Real y) { return 0.01*exp(-x * x - y * y); }, -100, 100, 200, -100, 100, 200,
		0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
	solver.Compute();

	printf("%12s | %12s %12s | %12s\n", "theta", "inv.mat.", "ana.for.", "ratio");
	for (int i = 0; i < 10; ++i) {
		Real theta = i * 2 * Pi / 10;
		Real cosx = cos(theta);
		Real cosy = sin(theta);
		Real a = solver.ComputeXSection(cosx, cosy);
		Real b = born(1, 1, 1, 0.01, 1, theta);
		printf("%12lf | %12.8E %12.8E | %12.8E\n", theta, a, b, a/b);
	}

}
