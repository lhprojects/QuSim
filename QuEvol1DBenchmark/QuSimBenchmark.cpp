
#include "../qsim/QuSim.h"
#include <chrono>

struct Test {
	SolverMethod met;
	char const *name;
	Options opts;
};

int main()
{

	Options space_O2;
	space_O2.SetBool("space_O2", true);

	Options space_O4;
	space_O4.SetBool("space_O2", false);

	Options fftw;
	fftw.SetString("fft_lib", "FFTW");

	Options cuda;
	cuda.SetString("fft_lib", "cuda");

	Options cuda_10;
	cuda_10.SetString("fft_lib", "cuda");
	cuda_10.SetInt("batch", 10);

	Options cuda_10_single;
	cuda_10_single.SetString("fft_lib", "cuda");
	cuda_10_single.SetInt("batch", 100);
	cuda_10_single.SetString("cuda_precision", "single");

	Test tests[] = {
#ifdef USE_CUDA
	//{ SolverMethod::SplittingMethodO2 , "splitO2+cuda", cuda },
	{ SolverMethod::SplittingMethodO2 , "splitO2+cuda_10", cuda_10_single },
	//{ SolverMethod::SplittingMethodO2 , "splitO2+cuda_10_single", cuda_10_single },
#endif
	{ SolverMethod::SplittingMethodO2 , "splitO2+kiss", Options() },
	{ SolverMethod::SplittingMethodO2 , "splitO2+fftw", fftw },
	{ SolverMethod::SplittingMethodO4 , "splitO4", Options() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", Options() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", Options() },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
	{ SolverMethod::GaussLegendreO6 , "gaussO6+spaceO4", Options() },
	{ SolverMethod::Eigen , "eigen", Options() },
	};

	int dims[] = {
		1000,
		5000,
		10000,
		//20000,
		//100000,
		//200000,
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
		if (tests[i].met == SolverMethod::ImplicitMidpointMethod && n >= 5) n = 5;
		if (tests[i].met == SolverMethod::GaussLegendreO4 && n >= 5) n = 5;
		if (tests[i].met == SolverMethod::GaussLegendreO6 && n >= 5) n = 5;
		if (tests[i].met == SolverMethod::Eigen && n >= 4) n = 4;
		for (int j = 0; j < n; ++j) {
			Evolver1D syst;
			auto vf = [](Real x) -> Real { return exp(-(x - 2)*(x - 2)); };
			syst.init([](Real x) -> Complex { return exp(I*x)*exp(-x * x); },
				true, 0.01, false,
				vf,
				-10, 10, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
				tests[i].opts);

			int dim = dims[j];
			int iter = 100000 * 10 / dim;

			auto t0 = std::chrono::system_clock::now();
			for (int it = 0; it < iter; ++it) {
				syst.step();
			}
			auto t1 = std::chrono::system_clock::now();

			auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
			printf("%6.2f ", (dim*syst.Time()/0.01)/(1E-6*d.count())*1E-6);
		}
		printf("\n");

	}

    return 0;
}

