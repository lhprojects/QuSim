
#include "System.h"
#include "SystemImpl.h"
#include "SplittingMethod.h"
#include "EigenMethod.h"
#include "GaussLegendreMethod.h"
#include "SplittingMethod2D.h"

using std::abs;

System::System() {
	fImpl = nullptr;
}

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
	Real mass, Real hbar,
	std::map<std::string, std::string> const &opts)
{
	if (solver == SolverMethod::SplittingMethodO2) {
		fImpl = new SplittingMethod();
	} else if (solver == SolverMethod::SplittingMethodO4) {
		fImpl = new SplittingMethod();
	} else if(solver == SolverMethod::Eigen) {
		fImpl = new EigenMethod();
	} else if (solver == SolverMethod::ImplicitMidpointMethod) {
		fImpl = new GaussLegendreMethod();
	} else if (solver == SolverMethod::GaussLegendreO4) {
		fImpl = new GaussLegendreMethod();
	} else if (solver == SolverMethod::GaussLegendreO6) {
		fImpl = new GaussLegendreMethod();
	} else {
		throw std::runtime_error("unspported solver");
	}
	System::fImpl = fImpl;
	fImpl->initSystem1D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar, opts);
}

System1D::System1D() {
	fImpl = nullptr;
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















System2D::System2D()
{
	fImpl = nullptr;
}

void System2D::init(char const * psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	char const * vs, Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar,
	std::map<std::string, std::string> const &opts)
{
	if (solver == SolverMethod::SplittingMethodO2) {
		fImpl = new SplittingMethod2D();
	} else if (solver == SolverMethod::SplittingMethodO4) {
		fImpl = new SplittingMethod2D();
	} else {
		throw std::runtime_error("unsupported solver");
	}
	System::fImpl = fImpl;
	fImpl->initSystem2D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, nx,
		y0, y1, ny,
		b, solver,
		mass, hbar, std::map<std::string, std::string>());

}

Eigen::MatrixXcd const & System2D::GetPsi()
{
	//double x = Norm2();
	return fImpl->fPsi;
}

Eigen::MatrixXd const & System2D::GetV()
{
	return fImpl->fV;
}

UInt System2D::GetNx()
{
	return fImpl->fNx;
}

UInt System2D::GetNy()
{
	return fImpl->fNy;
}
