
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

	std::map<std::string, std::string> smal_round_err;
	smal_round_err["small_round_error"] = "0";

	Test tests[] = {
		{ SolverMethod::ImplicitMidpointMethod , "midpoint", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint-", smal_round_err },
	{ SolverMethod::ExplicitRungeKuttaO4Classical , "rk4", std::map<std::string, std::string>() },
	{ SolverMethod::ExplicitRungeKuttaO4Classical , "rk4-", smal_round_err },
	{ SolverMethod::GaussLegendreO4 , "glo4", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO4 , "glo4-", smal_round_err },
	{ SolverMethod::ExplicitRungeKuttaO6Luther1967 , "rko6", std::map<std::string, std::string>() },
	{ SolverMethod::ExplicitRungeKuttaO6Luther1967 , "rko6-", smal_round_err },
	};

	int dims[] = {
		1000,
		2000,
		4000,
		8000,
		16000,
		32000,
		64000,
	};


	double ref = 0.1355284179587569045304922;


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
			syst.init(FunctorWrapper("exp(-x*x)"), -10, 10, dim, 0.5, 1, -I, tests[i].met, 1, 1, tests[i].opts);
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

	std::map<std::string, std::string> smal_round_err;
	smal_round_err["small_round_error"] = "0";

	Test tests[] = {
		{ SolverMethod::ImplicitMidpointMethod , "midpoint", std::map<std::string, std::string>() },
	{ SolverMethod::ImplicitMidpointMethod , "midpoint-", smal_round_err },
	{ SolverMethod::ExplicitRungeKuttaO4Classical , "rk4", std::map<std::string, std::string>() },
	{ SolverMethod::ExplicitRungeKuttaO4Classical , "rk4-", smal_round_err },
	{ SolverMethod::GaussLegendreO4 , "glo4", std::map<std::string, std::string>() },
	{ SolverMethod::GaussLegendreO4 , "glo4-", smal_round_err },
	{ SolverMethod::ExplicitRungeKuttaO6Luther1967 , "rko6", std::map<std::string, std::string>() },
	{ SolverMethod::ExplicitRungeKuttaO6Luther1967 , "rko6-", smal_round_err },
	};

	int dims[] = {
		1000,
		2000,
		4000,
		8000,
		16000,
		32000,
		64000,
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
			syst.init(FunctorWrapper("10*exp(-x*x)"), -10, 10, dim, 0.5, 1, -I, tests[i].met, 1, 1, tests[i].opts);
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

