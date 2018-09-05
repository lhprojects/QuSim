#include "../qsim/QuSim.h"

int main()
{

	QuPerturbation1D per;
	per.init(FunctorWrapper("0.0001*exp(-x*x)"), -5000, 5000, 200000, 0.5, 0.001,
		1, SolverMethod::BornSerise, 1, 1, std::map<std::string, std::string>());

	printf("Max MOmentum %lf\n", per.GetMaxMomentum());
	printf("Max Energy %lf\n", per.GetMaxEnergy());

	printf("Momentum Gap %lf\n", per.GetMomentumGap());
	printf("Energy Gap %lf\n", per.GetEnergyGap());

	printf("Epsilon %lf\n", per.GetEpsilon());
	printf("Epsilon Momentum Width %lf\n", per.GetEpsilonMomentumWidth());

	printf("Momentum Gap / Epsilon Momentum Width %lf\n", per.GetMomentumGap() / per.GetEpsilonMomentumWidth());
	printf("Energy Gap / Epsilon %lf\n", per.GetEnergyGap() / per.GetEpsilon());
	printf("Epsilon Boundary Error %lf\n", per.GetEpsilonBoundaryError());

	per.Compute();
	printf("R %lf\n", per.GetR());


	printf("%10s | %10s %10s %10s %10s | %10s\n", "V0", "RO1", "RO2", "RO10", "RO3N4", "Exact");
	for (int i = 0; i < 50; ++i) {
		
		Real v0 = 0.01 + 0.05*i;
		auto vfunc = [&](Real x) { return v0*exp(-x * x); };
		Solver1D solver;
		solver.init(vfunc, -10, 10, 20000, 0.5, 1, I,
			SolverMethod::ExplicitRungeKuttaO4Classical, 1, 1, std::map<std::string, std::string>());
		solver.Compute();

		QuPerturbation1D per1;
		per1.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
			1, SolverMethod::BornSerise, 1, 1, std::map<std::string, std::string>());
		per1.Compute();

		QuPerturbation1D per2;
		std::map<std::string, std::string> o2;
		o2["order"] = "2";
		per2.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
			1, SolverMethod::BornSerise, 1, 1, o2);
		per2.Compute();

		QuPerturbation1D per3;
		std::map<std::string, std::string> o3;
		o3["order"] = "10";
		per3.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
			1, SolverMethod::BornSerise, 1, 1, o3);
		per3.Compute();

		QuPerturbation1D per4;
		std::map<std::string, std::string> o4;
		o4["order"] = "3";
		o4["split_n"] = "3";
		per4.init(vfunc, -5000, 5000, 200000, 0.5, 0.002,
			1, SolverMethod::BornSerise, 1, 1, o4);
		per4.Compute();

		printf("%10lf | %10lf %10lf %10f %10lf | %10lf\n", 
			v0, per1.GetR(), per2.GetR(), per3.GetR(), per4.GetR(), solver.GetR());
	}

}

