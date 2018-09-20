#define _CRT_SECURE_NO_WARNINGS
#include "../qsim/QuSim.h"


void test0(double p)
{
	printf("Test inverse matrix method V(x)=%lf exp(-x*x)\n", p);
	auto f = [&](Real x) { return p * exp(-x * x); };
	Solver1D solver;
	std::map<std::string, std::string> small;
	small["small_round_error"] = "1";
	solver.init(f, -10, 10, 20000, 0.5, 1, I,
		SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1, std::map<std::string, std::string>());
	solver.Compute();

	if (p == 1) {
		double R = 1 - 0.1355284179587569045304922;
		printf("solver error %.2E\n", solver.GetR() - R);
	}

	for (int i = 0; i < 10; ++i) {

		double x0 = -150;
		double x1 = 150;
		size_t n = 500 + 500 * i;

		std::map<std::string, std::string> space_o2;
		space_o2["space_order"] = "2";
		std::map<std::string, std::string> space_o4;
		space_o4["space_order"] = "4";
		std::map<std::string, std::string> space_o6;
		space_o6["space_order"] = "6";

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
		printf("%5d | %10.2E %10.2E %10.2E\n", (int)n,
			abs2(inv1.GetPsi()[geti(-10)]) - solver.GetR(),
			abs2(inv2.GetPsi()[geti(-10)]) - solver.GetR(),
			abs2(inv3.GetPsi()[geti(-10)]) - solver.GetR());
		printf("%5d | %10.2E %10.2E %10.2E\n", (int)n,
			inv1.GetR() - solver.GetR(),
			inv2.GetR() - solver.GetR(),
			inv3.GetR() - solver.GetR());
	}
}

void test10()
{

	printf("Test inverse matrix method\n");
	double R = 1 - 0.1355284179587569045304922;

	for (int i = 0; i < 10; ++i) {

		double x0 = -150;
		double x1 = 150;
		size_t n = 500 + 500 * i;

		std::map<std::string, std::string> opt1;
		opt1["space_order"] = "6";
		opt1["matrix_solver"] = "LU";

		std::map<std::string, std::string> opt2;
		opt2["space_order"] = "6";
		opt2["matrix_solver"] = "BiCGSTAB";
		opt2["preconditioner"] = "DiagonalPreconditioner";

		std::map<std::string, std::string> opt3;
		opt3["space_order"] = "6";
		opt3["matrix_solver"] = "BiCGSTAB";
		opt3["preconditioner"] = "IdentityPreconditioner";

		std::map<std::string, std::string> opt4;
		opt4["space_order"] = "6";
		opt4["matrix_solver"] = "BiCGSTAB";
		opt4["preconditioner"] = "IncompleteLUT";

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


		auto geti = [&](double x) -> size_t {  return (size_t)((x - x0) / (x1 - x0) * n); };
		printf("%10.2E %10.2E %10.2E %10.2E\n",
			abs2(inv1.GetPsi()[geti(-10)]) - R,
			abs2(inv2.GetPsi()[geti(-10)]) - R,
			abs2(inv3.GetPsi()[geti(-10)]) - R,
			abs2(inv4.GetPsi()[geti(-10)]) - R);
	}
}


void test1()
{

	printf("Test R T calculation\n");
	Solver1D solver;
	solver.init(FunctorWrapper("0.000001*exp(-x*x)"), -10, 10, 20000, 0.5, 1, I,
		SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1, std::map<std::string, std::string>());
	solver.Compute();

	QuPerturbation1D per;
	per.init(FunctorWrapper("exp(-x*x)"), -50, 50, 20000, 0.5, 0.001,
		1, SolverMethod::BornSerise, 1, 1, std::map<std::string, std::string>());

	per.Compute();
	printf("    R %g (Exact %g)\n", per.GetR(), solver.GetR()/(0.000001*0.000001));

}

void test1d5()
{

	printf("Test perburbation 1 order correction\n");
	Solver1D solver;
	solver.init(FunctorWrapper("0.1*exp(-x*x)"), -10, 10, 20000, 0.5, 1, I,
		SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1, std::map<std::string, std::string>());
	solver.Compute();

	for (int i = 0; i < 10; ++i) {
		QuPerturbation1D per;
		char b[10];
		sprintf(b, "%d", i);

		std::map<std::string, std::string> opts;
		opts["order"] = b;
		per.init(FunctorWrapper("0.1*exp(-x*x)"), -1000, 1000, 20000, 0.5, 0.01,
			1, SolverMethod::BornSerise, 1, 1, opts);

		if (i == 0) {
			printf("  Max Momentum %lf\n", per.GetMaxMomentum());
			printf("  Max Energy %lf\n", per.GetMaxEnergy());

			printf("  Momentum Gap %lf\n", per.GetMomentumGap());
			printf("  Energy Gap %lf\n", per.GetEnergyGap());

			printf("  Epsilon %lf\n", per.GetEpsilon());
			printf("  Epsilon Momentum Width %lf\n", per.GetEpsilonMomentumWidth());

			printf("  Momentum Gap / Epsilon Momentum Width %lf\n", per.GetMomentumGap() / per.GetEpsilonMomentumWidth());
			printf("  Energy Gap / Epsilon %lf\n", per.GetEnergyGap() / per.GetEpsilon());
			printf("  Epsilon Boundary Error %lf\n", per.GetEpsilonBoundaryError());

		}

		per.Compute();
		printf("  R O%d %g (Exact %g)\n", 1 + i, per.GetR(),solver.GetR());

	}

}

void testPerburbation()
{
	printf("%10s | %10s %10s %10s %10s %10s %10s | %10s\n",
		    "V0", "RO1", "RO2", "RO10", "RO3N3", "RP5", "RP20", "Exact");
	for (int i = 0; i < 10; ++i) {

		Real v0 = 0.05 + 0.1*i;
		auto vfunc = [&](Real x) { return v0 * exp(-x * x); };
		Solver1D solver;
		solver.init(vfunc, -10, 10, 20000, 0.5, 1, I,
			SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1, std::map<std::string, std::string>());
		solver.Compute();

		QuPerturbation1D per1;
		per1.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
			1, SolverMethod::BornSerise, 1, 1, std::map<std::string, std::string>());
		per1.Compute();

		QuPerturbation1D per2;
		std::map<std::string, std::string> o2;
		o2["order"] = "2";
		per2.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
			1, SolverMethod::BornSerise, 1, 1, o2);
		per2.Compute();

		QuPerturbation1D per3;
		std::map<std::string, std::string> o3;
		o3["order"] = "10";
		per3.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
			1, SolverMethod::BornSerise, 1, 1, o3);
		per3.Compute();

		QuPerturbation1D per4;
		std::map<std::string, std::string> o4;
		o4["order"] = "3";
		o4["split_n"] = "3";
		per4.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
			1, SolverMethod::BornSerise, 1, 1, o4);
		per4.Compute();

		QuPerturbation1D per5;
		std::map<std::string, std::string> o5;
		o5["order"] = "5";
		o5["preconditional"] = "1";
		per5.init(vfunc, -100, 100, 2000, 0.5, 0.51,
			1, SolverMethod::BornSerise, 1, 1, o5);
		per5.Compute();

		QuPerturbation1D per6;
		std::map<std::string, std::string> o6;
		o6["order"] = "20";
		o6["preconditional"] = "1";
		per6.init(vfunc, -100, 100, 2000, 0.5, 0.51,
			1, SolverMethod::BornSerise, 1, 1, o6);
		per6.Compute();

		printf("%10lf | %10lf %10lf %10f %10lf %10lf %10lf| %10lf\n",
			v0, per1.GetR(), per2.GetR(), per3.GetR(), per4.GetR(), per5.GetR(), per6.GetR(),
			solver.GetR());
	}

}

void testPerburbativeConditioner()
{
	printf("%6s | %8s %8s %8s\n",
		"Order", "Vellekoop", "Hao1", "Hao2");

	char order[10];
	sprintf(order, "%d", 1);

	Real v0 = 1;
	auto vfunc = [&](Real x) { return v0 * exp(-x * x); };

	Solver1D solver;
	{

		std::map<std::string, std::string> opts;
		opts["small_round_error"] = "1";
		solver.init(vfunc, -10, 10, 2000, 0.5, 1, I,
			SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1, opts);
		solver.Compute();
	}

	QuPerturbation1D per1;
	{
		std::map<std::string, std::string> opts;
		opts["preconditional"] = "1";
		opts["order"] = order;
		opts["preconditioner"] = "Vellekoop";
		per1.init(vfunc, -150, 150, 10000, 0.5, 0,
			1, SolverMethod::BornSerise, 1, 1, opts);
	}

	QuPerturbation1D per2;
	{
		std::map<std::string, std::string> opts;
		opts["preconditional"] = "1";
		opts["order"] = order;
		opts["preconditioner"] = "Hao1";
		per2.init(vfunc, -150, 150, 10000, 0.5, 0,
			1, SolverMethod::BornSerise, 1, 1, opts);
	}

	QuPerturbation1D per3[10];
	{
		for (int i = 0; i < 10; ++i) {
			char b[10];
			sprintf(b, "%lf", 0.1 + 0.1*i);
			std::map<std::string, std::string> opts;
			opts["preconditional"] = "1";
			opts["order"] = order;
			opts["preconditioner"] = "Hao2";
			opts["slow"] = b;
			per3[i].init(vfunc, -150, 150, 10000, 0.5, 0,
				1, SolverMethod::BornSerise, 1, 1, opts);
		}
	}

	for (int i = 0; i < 2000; ++i) {

		per1.Compute();
		per2.Compute();

		printf("%6d | %8.1E %8.1E",
			i, per1.GetR() - solver.GetR(),
			per2.GetR() - solver.GetR());

		for (int i = 0; i < 10; ++i) {
			per3[i].Compute();
			printf(" %8.1E", per3[i].GetR() - solver.GetR());
		}
		printf("\n");

	}

}


int main()
{

	testPerburbativeConditioner();
	testPerburbation();
	test0(0.1);
	test0(1);
	test10();
	test1();
	test1d5();
}

