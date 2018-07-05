#include "SystemImpl.h"
#include "Cal.h"


void SystemImpl::init(char const * psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	char const *vs, BoundaryCondition b,
	SolverMethod solver,Real mass, Real hbar)
{
	fStep = 0;

	fFN = force_normalization;
	fPsiStr = psi;

	fVStr = vs;
	fBoundaryCondition = b;

	fDt = dt;
	fFNES = force_normalization_each_step;
	fMass = mass;
	fHbar = hbar;

	fSolverMethod = solver;

}


Real SystemImpl1D::CalPotEn()
{
	Real norm2 = 0;
	for (size_t i = 0; i < fN; ++i) {
		norm2 += abs2(fPsi[i])*fV[i] * fDx;
	}
	return norm2 / Norm2();
}

Real SystemImpl1D::CalKinEn()
{
	return 0;
}

Real SystemImpl1D::PotEn()
{
	return CalPotEn();
}

Real SystemImpl1D::KinEn()
{
	return CalKinEn();
}

Real SystemImpl1D::EnPartialT()
{
	Complex en = 0;
	for (size_t i = 0; i < fN; ++i) {
		en += I * fHbar*conj(fPsi[i])*(fPsi[i] - fLastLastPsi[i]) / (2.0*fDt)*fDx;
	}
	return en.real();
}

void SystemImpl1D::init(char const *psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	char const *vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar)
{

	SystemImpl::init(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, b, solver,
		mass, hbar);

	fX0 = x0;
	fDx = (x1 - x0) / n;
	fN = n;

	initPsi();
	fLastLastPsi = fPsi;
	fLastPsi = fPsi;

	initPotential();

}


void SystemImpl1D::initPotential()
{
	Cal cal(fVStr.c_str());
	fV.resize(fN);

	Complex f = -1.0 / fHbar * fDt * 0.5;

	for (size_t i = 0; i < fN; ++i) {
		Real x = getX(i);
		cal.SetVarVal("x", Complex(x));
		Complex com = cal.Val();
		fV[i] = com.real();

	}
	//dump(fV, "V.txt");
}

void SystemImpl1D::initPsi()
{
	fPsi.resize(fN);
	Cal cal(fPsiStr.c_str());
	for (size_t i = 0; i < fN; ++i) {
		Real x = getX(i);
		cal.SetVarVal("x", Complex(x));
		Complex com = cal.Val();
		fPsi[i] = com;
	}

	if (fFN) {
		Scale(fPsi, 1 / sqrt(Norm2()));
	}
	//dump(fPsi, "psi.txt");
}
