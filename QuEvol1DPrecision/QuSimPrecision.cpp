
#include "../qsim/QuSim.h"

#include <chrono>

struct Test {
	SolverMethod met;
	char const *name;
	std::map<std::string, std::string> opts;
};

void test_tunneling()
{

	printf("Test the tunneling probability\n");
	printf("Inital wave function: C*gauss(x, -20, 5)*exp(I*x)\n");
	printf("Potential: exp(-x*x)\n");

	std::map<std::string, std::string> space_O2;
	space_O2["space_O2"] = "1";

	std::map<std::string, std::string> space_O4;
	space_O4["space_O2"] = "0";

	Test tests[] = {
		{ SolverMethod::SplittingMethodO2 , "splitO2", std::map<std::string, std::string>() },
	{ SolverMethod::SplittingMethodO4 , "splitO4", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", std::map<std::string, std::string>() },
	{ SolverMethod::Eigen , "eigen", std::map<std::string, std::string>() },
	};

	int dims[] = {
		1000,
		2000,
		3000,
		4000,
		5000,
		6000,
	};

	double times[] = {
		1E-2,
		2E-2,
		3E-2,
		4E-2,
	};


	if (0) {
		Evolver1D syst;
		syst.init(FunctorWrapper("gauss(x, -20, 5)*exp(I*x)"), true, 5E-3, false, FunctorWrapper("exp(-x*x)"),
			-100, 100, 8000, BoundaryCondition::Period, SolverMethod::SplittingMethodO4, 1, 1,
			std::map<std::string, std::string>());
		for (; syst.Time() < 64;) {
			syst.step();
		}
		printf("ref 1 %.16f\n", syst.NormRight());
	}

	if (0) {
		Evolver1D syst;
		syst.init(FunctorWrapper("gauss(x, -20, 5)*exp(I*x)"), true, 5E-3, false, FunctorWrapper("exp(-x*x)"),
			-100, 100, 8000, BoundaryCondition::Period, SolverMethod::GaussLegendreO6, 1, 1,
			std::map<std::string, std::string>());
		for (; syst.Time() < 64;) {
			syst.step();
		}
		printf("ref 1 %.16f\n", syst.NormRight());

	}
	double ref = (0.1593774409 + 0.1593774408) / 2;

	for (int i = 0; i < sizeof(tests) / sizeof(Test); ++i) {
		printf("         Error for %-30s\n", tests[i].name);



		printf("%-10s ", "Dim\\Dt");
		for (int k = 0; k < sizeof(times) / sizeof(double); ++k) {
			printf("%10.3f ", times[k]);
		}
		printf("\n");


		for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
			int dim = dims[j];

			printf("%-10d ", dim);
			for (int k = 0; k < sizeof(times) / sizeof(double); ++k) {
				double time = times[k];

				Evolver1D syst;
				syst.init(FunctorWrapper("gauss(x, -20, 5)*exp(I*x)"), true, time, false, FunctorWrapper("exp(-x*x)"),
					-100, 100, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
					tests[i].opts);

				auto t0 = std::chrono::system_clock::now();
				for (; syst.Time() < 64;) {
					syst.step();
				}
				auto t1 = std::chrono::system_clock::now();

				auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

				printf("%+10.1E ", syst.NormRight() - ref);
			}
			printf("\n");
		}
		printf("\n");

	}

}

void test_speed()
{

	printf("Test the Speed\n");
	printf("Inital wave function: C*gauss(x, 0, 5)*exp(I*x)\n");
	printf("Potential: 0\n");

	std::map<std::string, std::string> space_O2;
	space_O2["space_O2"] = "1";

	std::map<std::string, std::string> space_O4;
	space_O4["space_O2"] = "0";

	std::map<std::string, std::string> fftw;
	fftw["fft_lib"] = "FFTW";

	std::map<std::string, std::string> cuda;
	cuda["fft_lib"] = "cuda";

	Test tests[] = {
	
	{ SolverMethod::SplittingMethodO2 , "splitO2+fftw", fftw },
	{ SolverMethod::SplittingMethodO2 , "splitO2+cuda", cuda },
	{ SolverMethod::SplittingMethodO2 , "splitO2+kiss", std::map<std::string, std::string>() },
	{ SolverMethod::SplittingMethodO4 , "splitO4", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO2", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint+spaceO4", space_O4 },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO2", space_O2 },
	{ SolverMethod::GaussLegendreO4 , "gaussO4+spaceO4", std::map<std::string, std::string>() },
	{ SolverMethod::Eigen , "eigen", std::map<std::string, std::string>() },
	};

	int dims[] = {
		1000,
		2000,
		3000,
		4000,
		5000,
		6000,
	};

	double times[] = {
		1E-2,
		2E-2,
		5E-2,
		1E-1,
	};

	double ref = 40;

	for (int i = 0; i < sizeof(tests) / sizeof(Test); ++i) {
		printf("         Error for %-30s\n", tests[i].name);



		printf("%-10s ", "Dim\\Dt");
		for (int k = 0; k < sizeof(times) / sizeof(double); ++k) {
			printf("%10.3f ", times[k]);
		}
		printf("\n");


		for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
			int dim = dims[j];

			printf("%-10d ", dim);
			for (int k = 0; k < sizeof(times) / sizeof(double); ++k) {
				double time = times[k];

				Evolver1D syst;
				syst.init(FunctorWrapper("gauss(x, 0, 5)*exp(I*x)"), true, time, false, FunctorWrapper("0"),
					-100, 100, dims[j], BoundaryCondition::Period, tests[i].met, 1, 1,
					tests[i].opts);

				auto t0 = std::chrono::system_clock::now();
				for (; ;) {
					syst.step();
					if (syst.Time() > 40 - 1E-6) break;
				}
				auto t1 = std::chrono::system_clock::now();

				auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
				//printf("%lf\n", syst.Xavg());
				//printf("%lf\n", syst.Time());
				printf("%+10.1E ", syst.Xavg()/ref - 1);
			}
			printf("\n");
		}
		printf("\n");

	}

}

int main()
{

	test_speed();
	test_tunneling();
	return 0;
}

