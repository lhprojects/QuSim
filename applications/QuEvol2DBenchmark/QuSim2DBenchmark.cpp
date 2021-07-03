
#define _CRT_SECURE_NO_WARNINGS
#include <QuSim.h>
#include <chrono>

struct Test {
	SolverMethod met;
	char const *name;
	Options opts;
};

int main()
{

	Options space_O2;
	space_O2.SpaceOrder(2);

	Options space_O4;
	space_O4.SpaceOrder(4);

	Options space_O6;
	space_O6.SpaceOrder(6);

	Options cuda;
	cuda.Cuda().Batch(100);

	Test tests[] = {
#ifdef QUSIM_USE_CUDA
	{ SolverMethod::SplittingMethodO2 , "splitO2+cuda", cuda },
	{ SolverMethod::SplittingMethodO4 , "splitO2+cuda", cuda },
#endif

#if 1
	{ SolverMethod::SplittingMethodO2 , "splitO2+kiss", Options() },
	{ SolverMethod::SplittingMethodO4 , "splitO4+kiss", Options() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", space_O2 },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO6", space_O6 },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", space_O4 },
	{ SolverMethod::GaussLegendreO6 , "gaussO6+spaceO6", space_O6 },
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
		printf("%25s ", "Iters*Dim^2/Time [M/s]");
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

			using seconds = std::chrono::duration<double, std::ratio<1> >;
			auto t0 = std::chrono::system_clock::now();
			auto t1 = std::chrono::system_clock::now();
			for (;;) {
				for (int it = 0; it < iter; ++it) {
					syst.step();
				}
				t1 = std::chrono::system_clock::now();
				if (std::chrono::duration_cast<seconds>(t1 - t0).count() > 1) {
					break;
				}
			}

			double time = std::chrono::duration_cast<seconds>(t1 - t0).count();
			double iters = syst.Time() / 0.01;
			if (strstr(tests[i].name, "cuda") == 0) {
				iters *= 100;
			}
            printf("%13.2f ", dim * dim * iters / time * 1E-6);
		}
		printf("\n");

	}

	return 0;
}

