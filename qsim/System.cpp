#define _CRT_SECURE_NO_WARNINGS

#include "System.h"
#include "SystemImpl.h"
#include "SystemHalfVTHalfV.h"
#include "SystemEigen.h"

using std::abs;


void System::step()
{
	fImpl->step();
}

Real System::Norm2()
{
	return fImpl->Norm2();
}

Real System::Time()
{
	return fImpl->Time();
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














void System1D::init(char const *psi, bool force_normalization,
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
	System::fImpl = fImpl;
	fImpl->init(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar);
}


PsiVector const & System1D::GetPsi()
{
	return fImpl->fPsi;
}

std::vector<Real> const & System1D::GetV()
{
	return fImpl->fV;
}

Real System1D::Xavg()
{
	return fImpl->Xavg();
}

UInt System1D::GetN()
{
	return fImpl->fN;
}

Real System1D::NormLeft()
{
	return fImpl->NormLeft();
}

Real System1D::NormRight()
{
	return fImpl->NormRight();
}
