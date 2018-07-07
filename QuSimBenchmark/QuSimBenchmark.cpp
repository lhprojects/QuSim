
#include "../qsim/System.h"
#include <chrono>

int main()
{

	SolverMethod solver[] = {
		SolverMethod::SplittingMethodO2,
		SolverMethod::SplittingMethodO4,
		SolverMethod::ImplicitMidpointMethod,
		SolverMethod::Eigen,
	};
	char const *name[] = {
		"SplitO2",
		"SplitO4",
		"Midpoint",
		"Eigen",
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

	printf("%10s ", "");
	printf("%20s ", "Dim");
	for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
		printf("%8d ", dims[j]);
	}
	printf("\n");

	for (int i = 0; i < 4; ++i) {
		printf("%10s ", name[i]);
		printf("%10s ", "Iters*Dim/Time [M/s]");
		int n = sizeof(dims) / sizeof(int);
		if (solver[i] == SolverMethod::ImplicitMidpointMethod) n = 6;
		if (solver[i] == SolverMethod::Eigen) n = 4;
		for (int j = 0; j < n; ++j) {
			System1D syst;
			syst.init("exp(I*x)*exp(-x*x)", true, 0.01, false, "exp(-(x-2)*(x-2))",
				-10, 10, dims[j], BoundaryCondition::Period, solver[i], 1, 1);
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

