
#include "System.h"
#include "SystemImpl.h"
#include "SplittingMethod.h"
#include "EigenMethod.h"
#include "GaussLegendreMethod.h"
#include "SplittingMethod2D.h"
#include "Cal.h"

using std::abs;

System::System() {
	fImpl = nullptr;
}

System::~System()
{
	
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














void System1D::init(std::function<Complex(Real)> const &psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar,
	std::map<std::string, std::string> const &opts)
{
	if (solver == SolverMethod::SplittingMethodO2) {
		fImpl.reset(new SplittingMethod());
	} else if (solver == SolverMethod::SplittingMethodO4) {
		fImpl.reset(new SplittingMethod());
	} else if(solver == SolverMethod::Eigen) {
		fImpl.reset(new EigenMethod());
	} else if (solver == SolverMethod::ImplicitMidpointMethod) {
		fImpl.reset(new GaussLegendreMethod());
	} else if (solver == SolverMethod::GaussLegendreO4) {
		fImpl.reset(new GaussLegendreMethod());
	} else if (solver == SolverMethod::GaussLegendreO6) {
		fImpl.reset(new GaussLegendreMethod());
	} else {
		throw std::runtime_error("unspported solver");
	}

	((SystemImpl1D*)fImpl.get())->initSystem1D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar, opts);
}

System1D::System1D() {
	fImpl = nullptr;
}


PsiVector const & System1D::GetPsi()
{
	return ((SystemImpl1D*)fImpl.get())->fPsi;
}

std::vector<Real> const & System1D::GetV()
{
	return ((SystemImpl1D*)fImpl.get())->fV;
}

Real System1D::Xavg()
{
	return ((SystemImpl1D*)fImpl.get())->Xavg();
}

UInt System1D::GetN()
{
	return ((SystemImpl1D*)fImpl.get())->fN;
}

Real System1D::NormLeft()
{
	return ((SystemImpl1D*)fImpl.get())->NormLeft();
}

Real System1D::NormRight()
{
	return ((SystemImpl1D*)fImpl.get())->NormRight();
}















System2D::System2D()
{
	fImpl = nullptr;
}

void System2D::init(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real, Real)> const &vs, Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar,
	std::map<std::string, std::string> const &opts)
{
	if (solver == SolverMethod::SplittingMethodO2) {
		fImpl.reset(new SplittingMethod2D());
	} else if (solver == SolverMethod::SplittingMethodO4) {
		fImpl.reset(new SplittingMethod2D());
	} else {
		throw std::runtime_error("unsupported solver");
	}

	((SystemImpl2D*)fImpl.get())->initSystem2D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, nx,
		y0, y1, ny,
		b, solver,
		mass, hbar, std::map<std::string, std::string>());

}

Eigen::MatrixXcd const & System2D::GetPsi()
{
	//double x = Norm2();
	return ((SystemImpl2D*)fImpl.get())->fPsi;
}

Eigen::MatrixXd const & System2D::GetV()
{
	return ((SystemImpl2D*)fImpl.get())->fV;
}

UInt System2D::GetNx()
{
	return ((SystemImpl2D*)fImpl.get())->fNx;
}

UInt System2D::GetNy()
{
	return ((SystemImpl2D*)fImpl.get())->fNy;
}










FunctorWrapper::FunctorWrapper(char const *str) : fCal(std::make_shared<Cal>(str))
{
}

Complex FunctorWrapper::operator()(Real x)
{
	fCal->SetVarVal("x", x);
	return fCal->Val();
}

FunctorWrapper::~FunctorWrapper()
{
}

Functor2DWrapper::Functor2DWrapper(char const *str) : fCal(std::make_shared<Cal>(str))
{
}

Complex Functor2DWrapper::operator()(Real x, Real y)
{
	fCal->SetVarVal("x", x);
	fCal->SetVarVal("y", y);
	return fCal->Val();
}

Functor2DWrapper::~Functor2DWrapper()
{
}
