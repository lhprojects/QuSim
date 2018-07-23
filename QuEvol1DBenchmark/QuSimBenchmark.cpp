
#include "../qsim/QuSim.h"
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

	std::map<std::string, std::string> space_O2;
	space_O2["space_O2"] = "1";

	std::map<std::string, std::string> space_O4;
	space_O4["space_O2"] = "0";

	std::map<std::string, std::string> fftw;
	fftw["fft_lib"] = "FFTW";

	std::map<std::string, std::string> cuda;
	cuda["fft_lib"] = "cuda";

	Test tests[] = {
#ifdef USE_CUDA
	{ SolverMethod::SplittingMethodO2 , "splitO2+cuda", cuda },
#endif
	{ SolverMethod::SplittingMethodO2 , "splitO2+kiss", std::map<std::string, std::string>() },
	{ SolverMethod::SplittingMethodO2 , "splitO2+fftw", fftw },
	{ SolverMethod::SplittingMethodO4 , "splitO4", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+split+spaceO2", split },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+split+spaceO4", split },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
	{ SolverMethod::GaussLegendreO6 , "gaussO6+spaceO4", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO6 , "gaussO6+split+spaceO4", split },
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

	printf("%-30s ", "");
	printf("%20s ", "Dim");
	for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
		printf("%6d ", dims[j]);
	}
	printf("\n");

	for (int i = 0; i < sizeof(tests)/sizeof(Test); ++i) {
		printf("%-30s ", tests[i].name);
		printf("%20s ", "Iters*Dim/Time [M/s]");
		int n = sizeof(dims) / sizeof(int);
		if (tests[i].met == SolverMethod::ImplicitMidpointMethod) n = 5;
		if (tests[i].met == SolverMethod::GaussLegendreO4) n = 5;
		if (tests[i].met == SolverMethod::GaussLegendreO6) n = 5;
		if (tests[i].met == SolverMethod::Eigen) n = 4;
		for (int j = 0; j < n; ++j) {
			Evolver1D syst;
			syst.init(FunctorWrapper("exp(I*x)*exp(-x*x)"), true, 0.01, false, FunctorWrapper("exp(-(x-2)*(x-2))"),
				-10, 10, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
				tests[i].opts);

			int dim = dims[j];
			int iter = 400000 * 10 / dim;

			auto t0 = std::chrono::system_clock::now();
			for (int it = 0; it < iter; ++it) {
				syst.step();
			}
			auto t1 = std::chrono::system_clock::now();

			auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
			printf("%6.2f ", (dim*iter)/(1E-6*d.count())*1E-6);
		}
		printf("\n");

	}

    return 0;
}

