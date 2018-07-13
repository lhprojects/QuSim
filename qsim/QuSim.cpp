
#include "QuSim.h"
#include "EvolverImpl.h"
#include "SplittingMethod.h"
#include "EigenMethod.h"
#include "GaussLegendreMethod.h"
#include "SplittingMethod2D.h"
#include "Cal.h"

using std::abs;

Evolver::Evolver() {
	fImpl = nullptr;
}

Evolver::~Evolver()
{
	
}

void Evolver::step()
{
	fImpl->step();
}

Real Evolver::Norm2()
{
	return fImpl->Norm2();
}

Real Evolver::Time()
{
	return fImpl->Time();
}


Real Evolver::PotEn()
{
	return fImpl->PotEn();
}

Real Evolver::KinEn()
{
	return fImpl->KinEn();
}

Real Evolver::EnPartialT()
{
	return fImpl->EnPartialT();
}














void Evolver1D::init(std::function<Complex(Real)> const &psi, bool force_normalization,
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

	((EvolverImpl1D*)fImpl.get())->initSystem1D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar, opts);
}

Evolver1D::Evolver1D() {
	fImpl = nullptr;
}


PsiVector const & Evolver1D::GetPsi()
{
	return ((EvolverImpl1D*)fImpl.get())->fPsi;
}

std::vector<Real> const & Evolver1D::GetV()
{
	return ((EvolverImpl1D*)fImpl.get())->fV;
}

Real Evolver1D::Xavg()
{
	return ((EvolverImpl1D*)fImpl.get())->Xavg();
}

UInt Evolver1D::GetN()
{
	return ((EvolverImpl1D*)fImpl.get())->fN;
}

Real Evolver1D::NormLeft()
{
	return ((EvolverImpl1D*)fImpl.get())->NormLeft();
}

Real Evolver1D::NormRight()
{
	return ((EvolverImpl1D*)fImpl.get())->NormRight();
}















Evolver2D::Evolver2D()
{
	fImpl = nullptr;
}

void Evolver2D::init(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
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

	((EvolverImpl2D*)fImpl.get())->initSystem2D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, nx,
		y0, y1, ny,
		b, solver,
		mass, hbar, std::map<std::string, std::string>());

}

Eigen::MatrixXcd const & Evolver2D::GetPsi()
{
	//double x = Norm2();
	return ((EvolverImpl2D*)fImpl.get())->fPsi;
}

Eigen::MatrixXd const & Evolver2D::GetV()
{
	return ((EvolverImpl2D*)fImpl.get())->fV;
}

UInt Evolver2D::GetNx()
{
	return ((EvolverImpl2D*)fImpl.get())->fNx;
}

UInt Evolver2D::GetNy()
{
	return ((EvolverImpl2D*)fImpl.get())->fNy;
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
