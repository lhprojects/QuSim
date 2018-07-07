#define _CRT_SECURE_NO_WARNINGS

#include "../qsim/System.h"
#include "../qsim/Cal.h"
#include "../qsim/SplittingMethod.h"

int main()
{
	testCal();
	SplittingMethod sm;
	int n = 10;
	double t = 1E-7;
	sm.initSystem1D("1", true, t, false, "0", -10, 10, n,
		BoundaryCondition::Period, SolverMethod::SplittingMethodO2, 1, 1);

	double dx = 20. / n;
	PsiVector psi;
	PsiVector tpsi;
	psi.resize(n);
	tpsi.resize(n);
	psi[n / 2] = 1;
	sm.ExpT(tpsi, psi, 1);

	tpsi[n / 2] -= 1;
	int i = 0;
	for (auto x : tpsi) {
		printf("%d %+10f, %+10f\n", i, x.real()/t*dx*dx, x.imag()/t*dx*dx);
		++i;
	}

	return 0;
}

