#define _CRT_SECURE_NO_WARNINGS

#include "System.h"
#include "SystemImpl.h"
#include "SystemHalfVTHalfV.h"
#include "SystemEigen.h"

using std::abs;

void dump(PsiVector const &v, char const *fn)
{
	FILE *f = fopen(fn," w");
	for (auto &x : v) {
		fprintf(f, "(% .20lf, % .20lf)\n", x.real(), x.imag());
	}
	fclose(f);
};

void dump(std::vector<Real> const &v, char const *fn)
{
	FILE *f = fopen(fn, " w");
	for (auto &x : v) {
		fprintf(f, "% .20lf\n", x);
	}
	fclose(f);
};

void System::init(char const *psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	char const *vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar)
{
	if (solver == SolverMethod::HalfVTHalfV) {
		fImpl = new SystemHalfVTHalfV();
	} else {
		fImpl = new SystemEigen();
	}
	fImpl->init(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar);
}

void System::step()
{
	fImpl->step();
}

PsiVector const & System::GetPsi()
{
	return fImpl->fPsi;
}

std::vector<Real> const & System::GetV()
{
	return fImpl->fV;
}

Real System::Xavg()
{
	return fImpl->Xavg();
}

Real System::Norm2()
{
	return fImpl->Norm2();
}

Real System::Time()
{
	return fImpl->Time();
}

UInt System::GetN()
{
	return fImpl->fN;
}

Real System::NormLeft()
{
	return fImpl->NormLeft();
}

Real System::NormRight()
{
	return fImpl->NormRight();
}

Real System::PotEn()
{
	return fImpl->PotEn();
}

Real System::KinEn()
{
	return fImpl->KinEn();
}

Real System::EnPartialT()
{
	return fImpl->EnPartialT();
}
