
#include "../qsim/QuSim.h"

struct Test {
	SolverMethod met;
	char const *name;
	std::map<std::string, std::string> opts;
};

void test_tunneling2()
{

	printf("Test the tunneling probability static schrodinger eq\n");
	printf("Potential: exp(-x*x)\n");


	Test tests[] = {
		{ SolverMethod::ImplicitMidpointMethod , "midpoint", std::map<std::string, std::string>() },
	{ SolverMethod::ExplicitRungeKuttaO4Classical , "rk4", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO4 , "glo4", std::map<std::string, std::string>() },
	};

	int dims[] = {
		1000,
		2000,
		4000,
		8000,
	};

	double ref = 0.1355284179587569045304922;
	if (1) {
		Solver1D syst;
		syst.init(FunctorWrapper("exp(-x*x)"), -10, 10, 32000*10, 0.5, 1, -I, SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1);
		syst.Compute();
		ref = syst.GetT();
		printf("%.20E\n", syst.GetT());
	}


	printf("%-20s ", "Dim");
	for (int k = 0; k < sizeof(dims) / sizeof(int); ++k) {
		printf("%10.3d ", dims[k]);
	}
	printf("\n");

	for (int i = 0; i < sizeof(tests) / sizeof(Test); ++i) {
		printf("%-20s ", tests[i].name);


		for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
			int dim = dims[j];

			Solver1D syst;
			syst.init(FunctorWrapper("exp(-x*x)"), -10, 10, dim, 0.5, 1, -I, tests[i].met, 1, 1);
			syst.Compute();

			printf("%+10.1E ", syst.GetT()/ref - 1);

		}
		printf("\n");

	}

}

void test_tunneling3()
{

	printf("Test the tunneling probability static schrodinger eq\n");
	printf("Potential: 10*exp(-x*x)\n");


	Test tests[] = {
		{ SolverMethod::ImplicitMidpointMethod , "midpoint", std::map<std::string, std::string>() },
	{ SolverMethod::ExplicitRungeKuttaO4Classical , "rk4", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO4 , "glo4", std::map<std::string, std::string>() },
	};

	int dims[] = {
		1000,
		2000,
		4000,
		8000,
	};

	double ref = 8.666893019856854462914111E-9;

	printf("%-20s ", "Dim");
	for (int k = 0; k < sizeof(dims) / sizeof(int); ++k) {
		printf("%10.3d ", dims[k]);
	}
	printf("\n");

	for (int i = 0; i < sizeof(tests) / sizeof(Test); ++i) {
		printf("%-20s ", tests[i].name);


		for (int j = 0; j < sizeof(dims) / sizeof(int); ++j) {
			int dim = dims[j];

			Solver1D syst;
			syst.init(FunctorWrapper("10*exp(-x*x)"), -10, 10, dim, 0.5, 1, -I, tests[i].met, 1, 1);
			syst.Compute();

			printf("%+10.1E ", syst.GetT()/ref - 1);

		}
		printf("\n");

	}

}

int main()
{

	test_tunneling2();
	test_tunneling3();
	return 0;
}

