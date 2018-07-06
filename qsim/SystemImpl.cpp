#include "SystemImpl.h"
#include "Cal.h"


void SystemImpl::initSystem(char const * psi, bool force_normalization,
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
	return en.real() / Norm2();
}

void SystemImpl1D::step()
{
	fLastLastPsi = fLastPsi;
	fLastPsi = fPsi;

	update_psi();

	if (fFNES) {
		Scale(fPsi, 1.0 / sqrt(Norm2()));
	}
	++fStep;

}

Real SystemImpl1D::Norm2()
{
	return Norm2(fPsi);
}

void SystemImpl1D::initSystem1D(char const *psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	char const *vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar)
{

	initSystem(psi, force_normalization,
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















void SystemImpl2D::initSystem2D(char const * psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	char const * vs, Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny,
	BoundaryCondition b,
	SolverMethod solver, Real mass, Real hbar)
{
	initSystem(psi, force_normalization, dt, force_normalization_each_step,
		vs, b, solver, mass, hbar);

	fX0 = x0;
	fDx = (x1 - x0) / nx;
	fNx = nx;
	fY0 = y0;
	fDy = (y1 - y0) / ny;
	fNy = ny;

	initPsi();
	fLastPsi = fPsi;
	fLastLastPsi = fLastPsi;

	initPotential();
	//Real x = Norm2();

}

void SystemImpl2D::initPsi()
{
	fPsi.resize(fNy, fNx);
	Cal cal(fPsiStr.c_str());
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			Real x = getX(i);
			Real y = getY(j);
			cal.SetVarVal("x", Complex(x));
			cal.SetVarVal("y", Complex(y));
			Complex com = cal.Val();
			fPsi(j, i) = com;
		}
	}

	if (fFN) {
		fPsi *= 1.0 / sqrt(Norm2());
	}
	double x = Norm2();
	//dump(fPsi, "psi.txt");

}

void SystemImpl2D::initPotential()
{
	fV.resize(fNy, fNx);
	Cal cal(fVStr.c_str());
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			Real x = getX(i);
			Real y = getY(j);
			cal.SetVarVal("x", Complex(x));
			cal.SetVarVal("y", Complex(y));
			Complex com = cal.Val();
			fV(j, i) = com.real();
		}
	}
	//dump_matrix_real(fV, "pot.txt");
	Real x = Norm2();

}


void SystemImpl2D::step()
{
	Real y = Norm2();
	fLastLastPsi = fLastPsi;
	fLastPsi = fPsi;

	//Real x = Norm2();
	update_psi();

	if (fFNES) {
		fPsi *= 1.0 / sqrt(Norm2());
	}
	++fStep;

}

Real SystemImpl2D::EnPartialT()
{

	Real pot = 0;
	Real norm2 = 0;
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			Complex dpsidt = (fPsi(j, i) - fLastLastPsi(j, i)) / (2.0*fDt);
			pot += (I * fHbar * conj(fPsi(j, i)) * dpsidt).real();
			norm2 += abs2(fPsi(j, i));
		}
	}
	return pot / norm2;

}


Real SystemImpl2D::Norm2()
{
	return fPsi.squaredNorm()*fDx*fDy;
}

Real SystemImpl2D::CalPotEn() const
{
	Real pot = 0;
	Real norm2 = 0;
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			pot += abs2(fPsi(j, i)) * fV(j, i);
			norm2 += abs2(fPsi(j, i));
		}
	}
	return pot / norm2;
}

Real SystemImpl2D::CalKinEn() const
{
	return 0;
}

