
#include "../qsim/QuSim.h"
#include <chrono>

struct Test {
	SolverMethod met;
	char const *name;
	std::map<std::string, std::string> opts;
};

int main()
{

	Test tests[] = {
		{ SolverMethod::ImplicitMidpointMethod , "midpoint", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO4 , "gaussO4", std::map<std::string, std::string>() },
	{ SolverMethod::ExplicitRungeKuttaO4Classical , "RK4", std::map<std::string, std::string>() },
	};

	int dims[] = {
		2000,
		10000,
		20000,
		100000,
		200000,
		1000000,
	};

	printf("%-10s ", "");
	printf("%15s ", "Dim");
	for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
		printf("%6d ", dims[j]);
	}
	printf("\n");

	for (int i = 0; i < sizeof(tests) / sizeof(Test); ++i) {
		printf("%-10s ", tests[i].name);
		printf("%15s ", "Dim/Time [M/s]");
		int n = sizeof(dims) / sizeof(int);
		for (int j = 0; j < n; ++j) {
			int dim = dims[j];

			Solver1D syst;
			syst.init(FunctorWrapper("exp(-x*x)"), -10, 10, dim, 0.5, 1, I, tests[i].met, 1, 1);


			auto t0 = std::chrono::system_clock::now();
			for(int i = 0; i < 10; ++i)
				syst.Compute();
			auto t1 = std::chrono::system_clock::now();

			auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0) / 10;
			printf("%6.2f ", dim / (1E-6*d.count())*1E-6);
		}
		printf("\n");

	}

	return 0;
}

