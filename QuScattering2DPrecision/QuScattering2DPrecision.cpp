
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

void testBorn()
{
	printf("test Born\n");
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

void testInverseAndPerburbation()
{
	Real v0 = 0.2;
	printf("test Inverse And Perburbation v0 = %lf\n", v0);
	auto vf = [&](Real x, Real y) { return v0*exp(-x * x - y * y); };

	QuScatteringInverseMatrix2D solver;
	{
		std::map<std::string, std::string> opts;
		solver.init(vf, -100, 100, 300, -100, 100, 300,
			0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
		solver.Compute();
	}

	QuPerturbation2D per1;
	{
		std::map<std::string, std::string> opts;
		opts["order"] = "1";
		per1.init(vf, -100, 100, 300, -100, 100, 300,
			0.5, 0.1, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
		per1.Compute();
	}

	QuPerturbation2D per2;
	{
		std::map<std::string, std::string> opts;
		opts["order"] = "2";
		per2.init(vf, -100, 100, 300, -100, 100, 300,
			0.5, 0.1, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
		per2.Compute();
	}

	QuPerturbation2D per10;
	{
		std::map<std::string, std::string> opts;
		opts["order"] = "10";
		per10.init(vf, -100, 100, 300, -100, 100, 300,
			0.5, 0.1, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
		per10.Compute();
	}

	QuPerturbation2D perp10;
	{
		std::map<std::string, std::string> opts;
		opts["preconditional"] = "1";
		opts["preconditioner"] = "Vellekoop";
		opts["order"] = "100";
		perp10.init(vf, -100, 100, 400, -100, 100, 400,
			0.5, 0.0, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
		perp10.Compute();
	}


	printf("%12s | %12s %12s %12s %12s %12s %12s\n", "theta", "born", "inv.mat.",
		"per.1", "per.2", "per.10", "perp.10");
	for (int i = 0; i < 10; ++i) {
		Real theta = i * 2 * Pi / 10;
		Real cosx = cos(theta);
		Real cosy = sin(theta);
		printf("%12lf | %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
			theta,
			born(1, 1, 1, v0, 1, 0),
			solver.ComputeXSection(cosx, cosy),
			per1.ComputeXSection(cosx, cosy),
			per2.ComputeXSection(cosx, cosy),
			per10.ComputeXSection(cosx, cosy),
			perp10.ComputeXSection(cosx, cosy)
			);
	}

}


void testPerburbationConverge(Real v0, int n = 100)
{
	printf("test Perburbation Convege v0 = %lf\n", v0);
	auto vf = [&](Real x, Real y) { return v0*exp(-x * x - y * y); };

	QuScatteringInverseMatrix2D solver;
	{
		std::map<std::string, std::string> opts;
		opts["space_order"] = "4";
		solver.init(vf, -100, 100, 500, -100, 100, 500,
			0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
		//solver.Compute();
	}

	QuPerturbation2D perp0;
	{
		std::map<std::string, std::string> opts;
		opts["preconditional"] = "0";
		opts["order"] = "1";
		//opts["slow"] = "0.5";
		//opts["fft_lib"] = "FFTW";
		perp0.init(vf, -150, 150, 1000, -150, 150, 1000,
			0.5, 0.005, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
	}

	QuPerturbation2D perp1;
	{
		std::map<std::string, std::string> opts;
		opts["preconditional"] = "1";
		opts["preconditioner"] = "Vellekoop";
		opts["order"] = "1";
		//opts["slow"] = "0.5";
		//opts["fft_lib"] = "FFTW";
		perp1.init(vf, -150, 150, 1000, -150, 150, 1000,
			0.5, 0.0, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
	}

	QuPerturbation2D perp2;
	{
		std::map<std::string, std::string> opts;
		opts["preconditional"] = "1";
		opts["preconditioner"] = "Hao1";
		opts["order"] = "1";
		//opts["slow"] = "0.5";
		//opts["fft_lib"] = "FFTW";
		perp2.init(vf, -150, 150, 1000, -150, 150, 1000,
			0.5, 0.0, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
	}

	QuPerturbation2D perp3;
	{
		std::map<std::string, std::string> opts;
		opts["preconditional"] = "1";
		opts["preconditioner"] = "Hao2";
		opts["order"] = "1";
		//opts["slow"] = "0.5";
		//opts["fft_lib"] = "FFTW";
		perp3.init(vf, -150, 150, 1000, -150, 150, 1000,
			0.5, 0.0, 1, 0, SolverMethod::BornSerise, 1, 1, opts);
	}


	printf("%5s | %16s %16s %16s %16s\n", "order", "inv.mat.", "born", "Vellekoop", "Hao1", "Hao2");
	for (int i = 0; i < n; ++i) {
		Real theta = 0;
		Real cosx = cos(theta);
		Real cosy = sin(theta);
		printf("%5d | %16.10E %16.10E %16.10E %16.10E %16.10E | %6.2E %6.2E %6.2E %6.2E\n",
			i,
			5.92704E-2,
			perp0.ComputeXSection(cosx, cosy),
			perp1.ComputeXSection(cosx, cosy),
			perp2.ComputeXSection(cosx, cosy),
			perp3.ComputeXSection(cosx, cosy),
			perp0.GetDeltaPsiNorm(),
			perp1.GetDeltaPsiNorm(),
			perp2.GetDeltaPsiNorm(),
			perp3.GetDeltaPsiNorm()

		);
		perp0.Compute();
		perp1.Compute();
		perp2.Compute();
		perp3.Compute();
	}

}

void testInvMatVsBorn()
{

	printf("Compare inverse matrix x section with born\n");
	printf("%10s | %10s %10s\n", "v0", "inv.mat", "born");
	for (int i = 0; i < 10; ++i) {
		Real v0 = 0.1 + 0.1 * i;
		std::map<std::string, std::string> opts;
		QuScatteringInverseMatrix2D solver;
		solver.init([&](Real x, Real y) { return v0 *exp(-x * x - y * y); }, -100, 100, 300, -100, 100, 300,
			0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
		solver.Compute();

		Real a = solver.ComputeXSection(1, 0);
		Real b = born(1, 1, 1, v0, 1, 0);
		printf("%10lf| %10.5E %10.5E\n", v0, a , b);
	}

}

void testTotalXSection()
{
	printf("test total x section converge");
	std::map<std::string, std::string> opts;
	QuScatteringInverseMatrix2D solver;
	solver.init([](Real x, Real y) { return 1 * exp(-x * x - y * y); }, -100, 100, 300, -100, 100, 300,
		0.5, 1, 0, SolverMethod::MatrixInverse, 1, 1, opts);
	solver.Compute();

	printf("Test Total X section\n");
	printf("%2s %8s", "N", "Err.\n");
	Real ref = solver.ComputeTotalXSection(1000);
	for (int i = 0; i < 20; ++i) {
		int n = 1 + i;
		Real xs = solver.ComputeTotalXSection(n);
		printf("%2d %8.1E\n", n, xs - ref);
	}


}

int main()
{
	testPerburbationConverge(0.2, 1000);
	testPerburbationConverge(0.6, 2000);
	testInverseAndPerburbation();
	testInvMatVsBorn();
	testBorn();


	testTotalXSection();

}
