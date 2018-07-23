
#define _CRT_SECURE_NO_WARNINGS
#include "../qsim/QuSim.h"
#include <chrono>

struct Test {
	SolverMethod met;
	char const *name;
	std::map<std::string, std::string> opts;
};

int main()
{

	std::map<std::string, std::string> space_O2;
	space_O2["space_O2"] = "1";

	std::map<std::string, std::string> space_O4;
	space_O4["space_O2"] = "0";

	std::map<std::string, std::string> fftw;
	fftw["fft_lib"] = "FFTW";

	std::map<std::string, std::string> cuda;
	cuda["fft_lib"] = "cuda";

	Test tests[] = {
		{ SolverMethod::SplittingMethodO2 , "splitO2+cuda", cuda },
	{ SolverMethod::SplittingMethodO2 , "splitO2+fftw", fftw },
	{ SolverMethod::SplittingMethodO2 , "splitO2", std::map<std::string, std::string>() },
	{ SolverMethod::SplittingMethodO4 , "splitO4", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
	{ SolverMethod::GaussLegendreO6 , "gaussO6+spaceO4", std::map<std::string, std::string>() },
	};

	int dims[] = {
		100,
		200,
		400,
		800,
	};

	printf("%-30s ", "");
	printf("%25s ", "Dim");
	for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
		char b[20];
		sprintf(b, "%d*%d", dims[j], dims[j]);
		printf("%13s ", b);
	}
	printf("\n");

	for (int i = 0; i < sizeof(tests) / sizeof(Test); ++i) {
		printf("%-30s ", tests[i].name);
		printf("%25s ", "Iters*Dim/Time [M/s]");
		int n = sizeof(dims) / sizeof(int);
		for (int j = 0; j < n; ++j) {
			Evolver2D syst;
			syst.init(Functor2DWrapper("exp(I*x)*exp(-x*x)*exp(-y*y)"), true, 0.01, false, Functor2DWrapper("exp(-x*x)"),
				-10, 10, dims[j],-10, 10, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
				tests[i].opts);

			int dim = dims[j];
			int iter = 5000000 * 10 / (dim*dim);

			auto t0 = std::chrono::system_clock::now();
			for (int it = 0; it < iter; ++it) {
				syst.step();
			}
			auto t1 = std::chrono::system_clock::now();

			auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
			printf("%13.2f ", (dim*dim*iter) / (1E-6*d.count())*1E-6);
		}
		printf("\n");

	}

	return 0;
}

