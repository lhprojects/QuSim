
#include "QuSim.h"
#include "EvolverImpl.h"
#include "SplittingMethod.h"
#include "SplittingMethod1DCUDA.h"
#include "EigenMethod.h"
#include "GaussLegendreMethod.h"
#include "SplittingMethod2D.h"
#include "SplittingMethod2DCUDA.h"
#include "GaussLegendreMethod2D.h"

#include "IVPSolverImpl.h"
#include "SolverImpl.h"

#include "Perturbation1DImpl.h"
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
	if (solver == SolverMethod::SplittingMethodO2 || solver == SolverMethod::SplittingMethodO4) {
		auto it = opts.find("fft_lib");
		if (it != opts.end() && it->second == "cuda") {
			fImpl.reset(CreateSplittingMethod1DCUDA(opts));
		} else {
			fImpl.reset(new SplittingMethod());
		}
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

size_t Evolver1D::GetN()
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
		auto it = opts.find("fft_lib");
		if (it != opts.end() && it->second == "cuda") {
			fImpl.reset(CreateSplittingMethod2DCUDA(opts));
		} else {
			fImpl.reset(new SplittingMethod2D());
		}
	} else if (solver == SolverMethod::SplittingMethodO4) {
		auto it = opts.find("fft_lib");
		if (it != opts.end() && it->second == "cuda") {
			fImpl.reset(CreateSplittingMethod2DCUDA(opts));
		} else {
			fImpl.reset(new SplittingMethod2D());
		}
	} else if (solver == SolverMethod::ImplicitMidpointMethod
		|| solver == SolverMethod::GaussLegendreO4
		|| solver == SolverMethod::GaussLegendreO6) {
		fImpl.reset(new GaussLegendreMethod2D());
	} else {
		throw std::runtime_error("unsupported solver");
	}

	((EvolverImpl2D*)fImpl.get())->initSystem2D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, nx,
		y0, y1, ny,
		b, solver,
		mass, hbar, opts);

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

size_t Evolver2D::GetNx()
{
	return ((EvolverImpl2D*)fImpl.get())->fNx;
}

size_t Evolver2D::GetNy()
{
	return ((EvolverImpl2D*)fImpl.get())->fNy;
}










FunctorWrapper::FunctorWrapper(char const *str) : fCal(std::make_shared<Cal>(str))
{
	fCal->SetVarVal("x", 0);
	fX = &fCal->GetVarVal("x");
	fCal->GenPseudoCode();
}

Complex FunctorWrapper::operator()(Real x)
{
	*fX = x;
	return fCal->RunPseudoCode();
}

FunctorWrapper::~FunctorWrapper()
{
}

Functor2DWrapper::Functor2DWrapper(char const *str) : fCal(std::make_shared<Cal>(str))
{
	fCal->SetVarVal("x", 0);
	fCal->SetVarVal("y", 0);
	fX = &fCal->GetVarVal("x");
	fY = &fCal->GetVarVal("y");
	fCal->GenPseudoCode();
}

Complex Functor2DWrapper::operator()(Real x, Real y)
{
	*fX = x;
	*fY = y;
	return fCal->RunPseudoCode();
}

Functor2DWrapper::~Functor2DWrapper()
{
}


void Solver::Compute()
{
	return fImpl->Compute();
}

Solver::Solver()
{
	fImpl = nullptr;
}

Solver::~Solver()
{
}


void Solver1D::init(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n,
	Real en,
	Complex initPsi,
	Complex initPsiPrime, SolverMethod met,
	Real mass, Real hbar,
	std::map<std::string, std::string> const &opts)
{
	fImpl.reset(new SolverImpl1D());
	((SolverImpl1D*)fImpl.get())->initSystem1D(v, x0, x1, n, en, initPsi, initPsiPrime, met, mass, hbar, opts);
}

PsiVector const & Solver1D::GetPsi()
{
	return ((SolverImpl1D*)fImpl.get())->fPsi;
}

std::vector<Real> const & Solver1D::GetV()
{
	return ((SolverImpl1D*)fImpl.get())->fV;
}

size_t Solver1D::GetNPoints()
{
	return ((SolverImpl1D*)fImpl.get())->fNPoints;
}

Eigen::Matrix2cd Solver1D::GetTMat()
{
	return ((SolverImpl1D*)fImpl.get())->fTMat;
}

Real Solver1D::GetT()
{
	return ((SolverImpl1D*)fImpl.get())->fT;
}

Real Solver1D::GetR()
{
	return ((SolverImpl1D*)fImpl.get())->fR;
}

Real Solver1D::GetEnergy()
{
	return ((SolverImpl1D*)fImpl.get())->fE;
}

Complex Solver1D::InitPsi()
{
	return ((SolverImpl1D*)fImpl.get())->fPsi[0];
}

Complex Solver1D::InitPsiPrime()
{
	return ((SolverImpl1D*)fImpl.get())->fPsiPrime[0];
}

Complex Solver1D::FinalPsi()
{
	return ((SolverImpl1D*)fImpl.get())->fPsi[((SolverImpl1D*)fImpl.get())->fNPoints - 1];
}

Complex Solver1D::FinalPsiPrime()
{
	return ((SolverImpl1D*)fImpl.get())->fPsiPrime[((SolverImpl1D*)fImpl.get())->fNPoints - 1];
}

QuPerturbation::QuPerturbation()
{
}

QuPerturbation::~QuPerturbation()
{
}

QuPerturbation1D::QuPerturbation1D()
{
}

void QuPerturbation1D::init(std::function<Complex(Real)> const & v, Real x0, Real x1, size_t n, Real en, Real epsilon,
	Real direction, SolverMethod met, Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{
	fImpl.reset(new QuPerturbation1DImpl());
	((QuPerturbation1DImpl*)fImpl.get())->InitQuPerturbation1D(v, x0, x1, n, en, epsilon, direction, met, mass, hbar, opts);
}

PsiVector const & QuPerturbation1D::GetPsi()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->fPsi;
}

std::vector<Real> const & QuPerturbation1D::GetV()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->fV;
}

size_t QuPerturbation1D::GetNPoints()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->fNx;
}

Real QuPerturbation1D::GetT()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->fT;
}

void QuPerturbation1D::Compute()
{
	static_cast<QuPerturbation1DImpl*>(fImpl.get())->Compute();
}

Real QuPerturbation1D::GetR()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->fR;
}


Real QuPerturbation1D::GetEnergy()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->fE;
}

Real QuPerturbation1D::GetMomentum()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->GetMomentum();
}

Real QuPerturbation1D::GetMaxEnergy()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->GetMaxEnergy();
}

Real QuPerturbation1D::GetMaxMomentum()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->GetMaxMomentum();
}

Real QuPerturbation1D::GetMomentumGap()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->GetMomentumGap();
}

Real QuPerturbation1D::GetEpsilonMomentumWidth()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->GetEpsilonMomentumWidth();
}

Real QuPerturbation1D::GetEnergyGap()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->GetEnergyGap();
}


Real QuPerturbation1D::GetEpsilon()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->fEpsilon;
}

Real QuPerturbation1D::GetEpsilonBoundaryError()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl.get())->GetEpsilonBoundaryError();
}
