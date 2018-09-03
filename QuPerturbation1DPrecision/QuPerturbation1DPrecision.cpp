#include "../qsim/QuSim.h"

int main()
{

	QuPerturbation1D per;
	per.init(FunctorWrapper("exp(-x*x)"), -40000, 40000, 1000000, 0.5, 0.001,
		1, SolverMethod::Eigen, 1, 1, std::map<std::string, std::string>());
	
	per.Compute();
	printf("R %lf\n", per.GetR());
	printf("T %lf\n", per.GetT());
	printf("Momentum Gap %lf\n", per.GetMomentumGap());
	printf("Energy Gap %lf\n", per.GetEnergyGap());

	printf("Epsilon %lf\n", per.GetEpsilon());
	printf("Epsilon Momentum Width %lf\n", per.GetEpsilonMomentumWidth());
	printf("Momentum Gap / Epsilon Momentum Width %lf\n", per.GetMomentumGap() / per.GetEpsilonMomentumWidth());
	printf("Epsilon boundary error %lf", per.GetEpsilonBoundaryError());
}

