
#define _CRT_SECURE_NO_WARNINGS
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
	space_O2.SetBool("space_O2", false);

	Options fftw;
	fftw.FFTW();

	Options cuda;
	cuda.Cuda();

	Options cuda_10;
	cuda_10.Cuda().Batch(10);

	Options cuda_10_single;
	cuda_10_single.Cuda().Batch(10).CudaPrecisionSingle();

	Test tests[] = {
#ifdef USE_CUDA
	{ SolverMethod::SplittingMethodO2 , "splitO2+cuda_batch10_single", cuda_10_single },
	{ SolverMethod::SplittingMethodO2 , "splitO2+cuda_batch10", cuda_10 },
	{ SolverMethod::SplittingMethodO2 , "splitO2+cuda", cuda },
#endif
#if 1
	{ SolverMethod::SplittingMethodO2 , "splitO2+fftw", fftw },
	{ SolverMethod::SplittingMethodO2 , "splitO2", Options() },
	{ SolverMethod::SplittingMethodO4 , "splitO4", Options() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", Options() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", Options() },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
	{ SolverMethod::GaussLegendreO6 , "gaussO6+spaceO4", Options() },
#endif
	};

	int dims[] = {
		200,
		400,
	};

	printf("%-30s ", "");
	printf("%25s ", "Dim");
	for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
		char b[20];
		sprintf(b, "%d*%d", dims[j], dims[j]);
		printf("%13s ", b);
	}
	printf("\n");

	int ntest = sizeof(tests) / sizeof(Test);
	for (int i = 0; i < ntest; ++i) {
		printf("%-30s ", tests[i].name);
		printf("%25s ", "Iters*Dim/Time [M/s]");
		int n = sizeof(dims) / sizeof(int);
		for (int j = 0; j < n; ++j) {
			Evolver2D syst;
			syst.init([](Real x, Real y) { return exp(I*x)*exp(-x * x)*exp(-y * y); }, true, 0.01, false,
				[](Real x, Real y) { return exp(-x * x); },
				-10, 10, dims[j], -10, 10, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
				tests[i].opts);

			int dim = dims[j];
			int iter = 1;

			//int batch = i < 2 ? 10 : 1;
			//iter /= batch;

			//if (iter < 2) iter = 2;

			auto t0 = std::chrono::system_clock::now();
			auto t1 = std::chrono::system_clock::now();
			for (;;) {
				for (int it = 0; it < iter; ++it) {
					syst.step();
				}
				t1 = std::chrono::system_clock::now();
				if (std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count() > 1) {
					break;
				}
			}

			auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
			printf("%13.2f ", (dim*dim*syst.Time()/0.01) / (1E-6*d.count())*1E-6);
		}
		printf("\n");

	}

	return 0;
}

