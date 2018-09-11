
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
#include "ComplexPotentialIVPSolver1DImpl.h"

#include "Perturbation1DImpl.h"
#include "ScatteringProblemSolverInverseMatrix.h"
#include "Cal.h"
#include "View.h"

using std::abs;

#define DELETE(x) if(x) { delete x; x = nullptr; }

Evolver::Evolver() {
	fImpl = nullptr;
}

Evolver::~Evolver()
{
	DELETE(fImpl);
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
			fImpl = CreateSplittingMethod1DCUDA(opts);
		} else {
			fImpl = new SplittingMethod();
		}
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

	static_cast<EvolverImpl1D*>(fImpl)->initSystem1D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar, opts);
}

Evolver1D::Evolver1D() {
	fImpl = nullptr;
}


VectorView<Complex> Evolver1D::GetPsi()
{
	return View(static_cast<EvolverImpl1D*>(fImpl)->fPsi);
}

VectorView<Real> Evolver1D::GetV()
{
	return View(static_cast<EvolverImpl1D*>(fImpl)->fV);
}

Real Evolver1D::Xavg()
{
	return static_cast<EvolverImpl1D*>(fImpl)->Xavg();
}

size_t Evolver1D::GetN()
{
	return static_cast<EvolverImpl1D*>(fImpl)->fN;
}

Real Evolver1D::NormLeft()
{
	return static_cast<EvolverImpl1D*>(fImpl)->NormLeft();
}

Real Evolver1D::NormRight()
{
	return static_cast<EvolverImpl1D*>(fImpl)->NormRight();
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
			fImpl = CreateSplittingMethod2DCUDA(opts);
		} else {
			fImpl = new SplittingMethod2D();
		}
	} else if (solver == SolverMethod::SplittingMethodO4) {
		auto it = opts.find("fft_lib");
		if (it != opts.end() && it->second == "cuda") {
			fImpl = CreateSplittingMethod2DCUDA(opts);
		} else {
			fImpl = new SplittingMethod2D();
		}
	} else if (solver == SolverMethod::ImplicitMidpointMethod
		|| solver == SolverMethod::GaussLegendreO4
		|| solver == SolverMethod::GaussLegendreO6) {
		fImpl = new GaussLegendreMethod2D();
	} else {
		throw std::runtime_error("unsupported solver");
	}

	static_cast<EvolverImpl2D*>(fImpl)->initSystem2D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, nx,
		y0, y1, ny,
		b, solver,
		mass, hbar, opts);

}

MatrixView<Complex> Evolver2D::GetPsi()
{
	//double x = Norm2();
	return View(static_cast<EvolverImpl2D*>(fImpl)->fPsi);
}

MatrixView<Real> Evolver2D::GetV()
{
	return View(static_cast<EvolverImpl2D*>(fImpl)->fV);
}

size_t Evolver2D::GetNx()
{
	return static_cast<EvolverImpl2D*>(fImpl)->fNx;
}

size_t Evolver2D::GetNy()
{
	return static_cast<EvolverImpl2D*>(fImpl)->fNy;
}










FunctorWrapper::FunctorWrapper(char const *str) : fCal(new Cal(str)), fRef(new int(1))
{
	fCal->SetVarVal("x", 0);
	fX = &fCal->GetVarVal("x");
	fCal->GenPseudoCode();
}

FunctorWrapper::FunctorWrapper(FunctorWrapper const &r) : fCal(r.fCal), fRef(r.fRef)
{
	(*fRef)++;
	fX = r.fX;
}

FunctorWrapper &FunctorWrapper::operator=(FunctorWrapper const &r)
{
	(*r.fRef)++;
	(*fRef)--;
	if (*fRef == 0) {
		DELETE(fCal);
		DELETE(fRef);
	}
	fCal = r.fCal;
	fRef = r.fRef;
	fX = r.fX;
	return *this;
}

Complex FunctorWrapper::operator()(Real x)
{
	*fX = x;
	return fCal->RunPseudoCode();
}

FunctorWrapper::~FunctorWrapper()
{
	(*fRef)--;
	if (*fRef == 0) {
		DELETE(fCal);
		DELETE(fRef);
	}
	fCal = nullptr;
	fRef = nullptr;
}

Functor2DWrapper::Functor2DWrapper(char const *str) : fCal(new Cal(str)), fRef(new int(1))
{
	fCal->SetVarVal("x", 0);
	fCal->SetVarVal("y", 0);
	fX = &fCal->GetVarVal("x");
	fY = &fCal->GetVarVal("y");
	fCal->GenPseudoCode();
}

Functor2DWrapper::Functor2DWrapper(Functor2DWrapper const &r) : fCal(r.fCal), fRef(r.fRef)
{
	(*fRef)++;
	fX = r.fX;
	fY = r.fY;
}

Functor2DWrapper &Functor2DWrapper::operator=(Functor2DWrapper const &r)
{
	(*r.fRef)++;
	(*fRef)--;
	if (*fRef == 0) {
		DELETE(fCal);
		DELETE(fRef);
	}
	fCal = r.fCal;
	fRef = r.fRef;
	fX = r.fX;
	fY = r.fY;
	return *this;
}

Complex Functor2DWrapper::operator()(Real x, Real y)
{
	*fX = x;
	*fY = y;
	return fCal->RunPseudoCode();
}

Functor2DWrapper::~Functor2DWrapper()
{
	(*fRef)--;
	if (*fRef == 0) {
		DELETE(fCal);
		DELETE(fRef);
	}
	fCal = nullptr;
	fRef = nullptr;
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
	DELETE(fImpl);
}


void Solver1D::init(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n,
	Real en,
	Complex initPsi,
	Complex initPsiPrime, SolverMethod met,
	Real mass, Real hbar,
	std::map<std::string, std::string> const &opts)
{
	if (opts.find("complex_potential") != opts.end() && opts.find("complex_potential")->second != "0")
		fImpl = new ComplexPotentialIVPSolver1DImpl();
	else
		fImpl = new SolverImpl1D();
	static_cast<SolverImpl1D*>(fImpl)->initSystem1D(v, x0, x1, n, en, initPsi, initPsiPrime, met, mass, hbar, opts);
}

VectorView<Complex> Solver1D::GetPsi()
{
	return View(static_cast<SolverImpl1D*>(fImpl)->fPsi);
}

VectorView<Real> Solver1D::GetV()
{
	return View(static_cast<SolverImpl1D*>(fImpl)->fV);
}

size_t Solver1D::GetNPoints()
{
	return static_cast<SolverImpl1D*>(fImpl)->fNPoints;
}

MatrixView<Real> Solver1D::GetTMat()
{
	return View(static_cast<SolverImpl1D*>(fImpl)->fTMat);
}

Real Solver1D::GetT()
{
	return static_cast<SolverImpl1D*>(fImpl)->fT;
}

Real Solver1D::GetR()
{
	return static_cast<SolverImpl1D*>(fImpl)->fR;
}

Real Solver1D::GetEnergy()
{
	return static_cast<SolverImpl1D*>(fImpl)->fE;
}

Complex Solver1D::InitPsi()
{
	return static_cast<SolverImpl1D*>(fImpl)->fPsi[0];
}

Complex Solver1D::InitPsiPrime()
{
	return static_cast<SolverImpl1D*>(fImpl)->fPsiPrime[0];
}

Complex Solver1D::FinalPsi()
{
	return static_cast<SolverImpl1D*>(fImpl)->fPsi[static_cast<SolverImpl1D*>(fImpl)->fNPoints - 1];
}

Complex Solver1D::FinalPsiPrime()
{
	return static_cast<SolverImpl1D*>(fImpl)->fPsiPrime[static_cast<SolverImpl1D*>(fImpl)->fNPoints - 1];
}





QuScatteringProblemSolver::QuScatteringProblemSolver() : fImpl(nullptr)
{
}

QuScatteringProblemSolver::~QuScatteringProblemSolver()
{
	DELETE(fImpl);
}

void QuScatteringProblemSolver::Compute()
{
	static_cast<ScatteringSolverImpl*>(fImpl)->Compute();
}






VectorView<Complex> QuScatteringProblemSolver1D::GetPsi()
{
	return View(static_cast<ScatteringSolver1DImpl*>(fImpl)->fPsiX);
}

VectorView<Real> QuScatteringProblemSolver1D::GetV()
{
	return View(static_cast<ScatteringSolver1DImpl*>(fImpl)->fV);
}

size_t QuScatteringProblemSolver1D::GetNPoints()
{
	return static_cast<ScatteringSolver1DImpl*>(fImpl)->fNx;
}

Real QuScatteringProblemSolver1D::GetT()
{
	return static_cast<ScatteringSolver1DImpl*>(fImpl)->fT;
}


Real QuScatteringProblemSolver1D::GetR()
{
	return static_cast<ScatteringSolver1DImpl*>(fImpl)->fR;
}


Real QuScatteringProblemSolver1D::GetEnergy()
{
	return static_cast<ScatteringSolver1DImpl*>(fImpl)->fE;
}

Real QuScatteringProblemSolver1D::GetMomentum()
{
	return static_cast<ScatteringSolver1DImpl*>(fImpl)->GetMomentum();
}




// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/*                    QuScatteringInverseMatrix1D                          */
void QuScatteringInverseMatrix1D::init(std::function<Complex(Real)> const & v, Real x0, Real x1, size_t n, Real en, 
	Real direction, SolverMethod met, Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{

	fImpl = new ScatteringProblemSolverInverseMatrix1D();
	static_cast<ScatteringProblemSolverInverseMatrix1D*>(fImpl)->InitScatteringSolver1D(v, x0, x1, n, en, direction, met, mass, hbar, opts);

}
/*                    QuScatteringInverseMatrix1D                          */
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/*                    QuPerturbation1D                          */

QuPerturbation1D::QuPerturbation1D()
{
}

void QuPerturbation1D::init(std::function<Complex(Real)> const & v, Real x0, Real x1, size_t n, Real en, Real epsilon,
	Real direction, SolverMethod met, Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{
	fImpl = new QuPerturbation1DImpl();
	static_cast<QuPerturbation1DImpl*>(fImpl)->InitPerturbation1D(v, x0, x1, n, en, epsilon, direction, met, mass, hbar, opts);
}


Real QuPerturbation1D::GetMaxEnergy()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl)->GetMaxEnergy();
}

Real QuPerturbation1D::GetMaxMomentum()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl)->GetMaxMomentum();
}

Real QuPerturbation1D::GetMomentumGap()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl)->GetMomentumGap();
}

Real QuPerturbation1D::GetEpsilonMomentumWidth()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl)->GetEpsilonMomentumWidth();
}

Real QuPerturbation1D::GetEnergyGap()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl)->GetEnergyGap();
}


Real QuPerturbation1D::GetEpsilon()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl)->fEpsilon;
}

Real QuPerturbation1D::GetEpsilonBoundaryError()
{
	return static_cast<QuPerturbation1DImpl*>(fImpl)->GetEpsilonBoundaryError();
}
/*                    QuPerturbation1D                          */
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Calculator::Calculator(char const * expr) : fImpl(nullptr)
{
	fImpl = new Cal(expr);
}

Calculator::~Calculator()
{
	DELETE(fImpl);
}

Complex & Calculator::SetVaraible(char const *str, Complex v)
{
	fImpl->SetVarVal(str, v);
	return fImpl->GetVarVal(str);
}

Complex & Calculator::GetVaraible(char const *str)
{
	return fImpl->GetVarVal(str);
}

Complex Calculator::Evaluate()
{
	return fImpl->Val();
}
