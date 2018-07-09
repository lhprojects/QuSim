
#include "../qsim/System.h"
#include <chrono>

struct Test {
	SolverMethod met;
	char const *name;
	std::map<std::string, std::string> opts;
};

int main()
{

	std::map<std::string, std::string> split;
	split["split_time_2"] = "1";

	Test tests[] = {
	{ SolverMethod::ImplicitMidpointMethod , "midpoint", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+split", split },
	{ SolverMethod::GaussLegendreO4 , "gaussO4", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+split", split },
	{ SolverMethod::GaussLegendreO6 , "gaussO6", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO6 , "gaussO6+split", split },
	{ SolverMethod::SplittingMethodO2 , "splitO2", std::map<std::string, std::string>() },
	{ SolverMethod::SplittingMethodO4 , "splitO4", std::map<std::string, std::string>() },
	{ SolverMethod::Eigen , "eigen", std::map<std::string, std::string>() },
	};

	int dims[] = {
		100,
		200,
		1000,
		2000,
		10000,
		20000,
		100000,
		200000,
	};

	printf("%-20s ", "");
	printf("%20s ", "Dim");
	for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
		printf("%8d ", dims[j]);
	}
	printf("\n");

	for (int i = 0; i < sizeof(tests)/sizeof(Test); ++i) {
		printf("%-20s ", tests[i].name);
		printf("%10s ", "Iters*Dim/Time [M/s]");
		int n = sizeof(dims) / sizeof(int);
		if (tests[i].met == SolverMethod::ImplicitMidpointMethod) n = 5;
		if (tests[i].met == SolverMethod::GaussLegendreO4) n = 5;
		if (tests[i].met == SolverMethod::GaussLegendreO6) n = 5;
		if (tests[i].met == SolverMethod::Eigen) n = 4;
		for (int j = 0; j < n; ++j) {
			System1D syst;
			syst.init("exp(I*x)*exp(-x*x)", true, 0.01, false, "exp(-(x-2)*(x-2))",
				-10, 10, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
				tests[i].opts);

			int dim = dims[j];
			int iter = 500000 * 10 / dim;

			auto t0 = std::chrono::system_clock::now();
			for (int it = 0; it < iter; ++it) {
				syst.step();
			}
			auto t1 = std::chrono::system_clock::now();

			auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
			printf("%8.3f ", (dim*iter)/(1E-6*d.count())*1E-6);
		}
		printf("\n");

	}

    return 0;
}

