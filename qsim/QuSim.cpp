
#include "QuSim.h"
#include "EvolverImpl.h"
#include "SplittingMethod1D.h"
#include "SplittingMethod2D.h"
#include "EigenMethod.h"
#include "SplittingMethod1DCUDA.h"
#include "SplittingMethod2DCUDA.h"
#include "GaussLegendreMethod1D.h"
#include "GaussLegendreMethod2D.h"

#include "IVPSolverImpl.h"
#include "SolverImpl.h"
#include "ComplexPotentialIVPSolver1DImpl.h"

#include "Perturbation1DImpl.h"
#include "Perturbation2DImpl.h"
#include "Perturbation3DImpl.h"
#include "ScatteringProblemSolverInverseMatrix1D.h"
#include "ScatteringProblemSolverInverseMatrix2D.h"
#include "ScatteringProblemSolverInverseMatrix3D.h"

#include "Cal.h"
#include "View.h"

#include "OptionsImpl.h"

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

SolverMethod Evolver::GetMethod()
{
	return fImpl->fSolverMethod;
}














void Evolver1D::init(std::function<Complex(Real)> const &psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar,
	Options const &opts)
{
	if (solver == SolverMethod::SplittingMethodO2 || solver == SolverMethod::SplittingMethodO4) {

		if (opts.fOpts->GetString("fft_lib", "") == "cuda") {
			fImpl = CreateSplittingMethod1DCUDA(*opts.fOpts);
		} else {
			fImpl = new SplittingMethod1D();
		}
	} else if (solver == SolverMethod::Eigen) {
		fImpl = new EigenMethod();
	} else if (solver == SolverMethod::ImplicitMidpointMethod
		|| solver == SolverMethod::GaussLegendreO4
		|| solver == SolverMethod::GaussLegendreO6) {
		fImpl = new GaussLegendreMethod1D();
	} else {
		throw std::runtime_error("unspported solver");
	}

	static_cast<EvolverImpl1D*>(fImpl)->initSystem1D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar, *opts.fOpts);
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
	Options const &opts)
{
	if (solver == SolverMethod::SplittingMethodO2) {

		if (opts.fOpts->GetString("fft_lib", "") == "cuda") {
			fImpl = CreateSplittingMethod2DCUDA(*opts.fOpts);
		} else {
			fImpl = new SplittingMethod2D();
		}
	} else if (solver == SolverMethod::SplittingMethodO4) {
		if (opts.fOpts->GetString("fft_lib", "") == "cuda") {
			fImpl = CreateSplittingMethod2DCUDA(*opts.fOpts);
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
		mass, hbar, *opts.fOpts);

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

SolverMethod Solver::GetMethod()
{
	return fImpl->fMethod;
}



void Solver1D::init(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n,
	Real en,
	Complex initPsi,
	Complex initPsiPrime, SolverMethod met,
	Real mass, Real hbar,
	Options const &opts)
{
	if (opts.fOpts->GetBool("complex_potential", false))
		fImpl = new ComplexPotentialIVPSolver1DImpl();
	else
		fImpl = new SolverImpl1D();
	static_cast<SolverImpl1D*>(fImpl)->initSystem1D(v, x0, x1, n, en, initPsi,
		initPsiPrime, met, mass, hbar, *opts.fOpts);
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

Real QuScatteringProblemSolver::GetEnergy()
{
	return static_cast<ScatteringSolverImpl*>(fImpl)->fE;
}

Real QuScatteringProblemSolver::GetMomentum()
{
	return static_cast<ScatteringSolverImpl*>(fImpl)->GetMomentum();
}

SolverMethod QuScatteringProblemSolver::GetMethod()
{
	return fImpl->fMet;
}






VectorView<Complex> QuScatteringProblemSolver1D::GetPsi()
{
	return View(static_cast<ScatteringSolver1DImpl*>(fImpl)->fPsiX);
}

VectorView<Complex> QuScatteringProblemSolver1D::GetPsi0()
{
	return View(static_cast<ScatteringSolver1DImpl*>(fImpl)->fPsi0X);
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






// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/*                    QuScatteringInverseMatrix1D                          */
void QuScatteringInverseMatrix1D::init(std::function<Complex(Real)> const & v, Real x0, Real x1, size_t n, Real en, 
	Real direction, SolverMethod met, Real mass, Real hbar, Options const & opts)
{

	fImpl = new ScatteringProblemSolverInverseMatrix1D();
	static_cast<ScatteringProblemSolverInverseMatrix1D*>(fImpl)->InitScatteringSolver1D(v, x0, x1,
		n, en, direction, met, mass, hbar, *opts.fOpts);

}
/*                    QuScatteringInverseMatrix1D                          */
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/*                    QuPerturbation1D                          */

QuPerturbation1D::QuPerturbation1D()
{
}

void QuPerturbation1D::init(std::function<Complex(Real)> const & v, Real x0, Real x1, size_t n, Real en, Real epsilon,
	Real direction, SolverMethod met, Real mass, Real hbar, Options const & opts)
{
	fImpl = new QuPerturbation1DImpl();
	static_cast<QuPerturbation1DImpl*>(fImpl)->InitPerturbation1D(v, x0, x1, n,
		en, epsilon, direction, met, mass, hbar, *opts.fOpts);
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

MatrixView<Complex> QuScatteringProblemSolver2D::GetPsi0()
{
	return View(static_cast<ScatteringSolver2DImpl*>(fImpl)->fPsi0X);
}

MatrixView<Complex> QuScatteringProblemSolver2D::GetPsi()
{
	return View(static_cast<ScatteringSolver2DImpl*>(fImpl)->fPsiX);
}

MatrixView<Real> QuScatteringProblemSolver2D::GetV()
{
	return View(static_cast<ScatteringSolver2DImpl*>(fImpl)->fV);
}

Real QuScatteringProblemSolver2D::ComputeXSection(Real cosx, Real cosy)
{
	return static_cast<ScatteringSolver2DImpl*>(fImpl)->ComputeXSection(cosx, cosy);
}

Real QuScatteringProblemSolver2D::ComputeTotalXSection(Int n)
{
	return static_cast<ScatteringSolver2DImpl*>(fImpl)->ComputeTotalXSection(n);
}

void QuScatteringInverseMatrix2D::init(std::function<Complex(Real, Real)> const & v, Real x0,
	Real x1, size_t nx, Real y0, Real y1, size_t ny,
	Real en, Real directionx, Real directiony, SolverMethod met, Real mass, 
	Real hbar, Options const & opts)
{
	fImpl = new ScatteringProblemSolverInverseMatrix2D();
	static_cast<ScatteringProblemSolverInverseMatrix2D*>(fImpl)->InitScatteringSolver2D(v, x0, x1,
		nx, y0, y1, ny,
		en, directionx, directiony, met, mass, hbar, *opts.fOpts);
}


void QuPerturbation2D::init(std::function<Complex(Real, Real)> const & v, Real x0, Real x1, size_t nx, 
	Real y0, Real y1, size_t ny, Real en, Real epsilon,
	Real directionx, Real directiony, SolverMethod met, Real mass,
	Real hbar, Options const & opts)
{
	fImpl = new QuPerturbation2DImpl();
	static_cast<QuPerturbation2DImpl*>(fImpl)->InitPerturbation2D(v, x0, x1,
		nx, y0, y1, ny,
		en, epsilon,
		directionx, directiony, met, mass, hbar, *opts.fOpts);
}

MatrixView<Real> QuPerturbation2D::GetVabsb()
{
	return View(static_cast<QuPerturbation2DImpl*>(fImpl)->fVasb.data(),
		static_cast<QuPerturbation2DImpl*>(fImpl)->fNx,
		static_cast<QuPerturbation2DImpl*>(fImpl)->fNy);
}

Real QuPerturbation2D::GetDeltaPsiNorm()
{
	return static_cast<QuPerturbation2DImpl*>(fImpl)->fNormDeltaPsi;
}


Tensor3View<Complex> QuScatteringProblemSolver3D::GetPsi()
{
	auto impl = static_cast<ScatteringProblemSolverInverseMatrix3D*>(fImpl);
	return View(impl->fPsiX.data(), impl->fNx, impl->fNy, impl->fNz);
}

Tensor3View<Real> QuScatteringProblemSolver3D::GetV()
{
	auto impl = static_cast<ScatteringProblemSolverInverseMatrix3D*>(fImpl);
	return View(impl->fV.data(), impl->fNx, impl->fNy, impl->fNz);
}

Real QuScatteringProblemSolver3D::ComputeXSection(Real cosx, Real cosy, Real cosz)
{
	return static_cast<ScatteringProblemSolverInverseMatrix3D*>(fImpl)->ComputeXSection(cosx, cosy, cosz);
}

Real QuScatteringProblemSolver3D::ComputeTotalXSection(Int npsi, Int ntheta)
{
	return static_cast<ScatteringProblemSolverInverseMatrix3D*>(fImpl)->ComputeTotalXSection(npsi, ntheta);
}

void QuScatteringInverseMatrix3D::init(std::function<Complex(Real, Real, Real)> const & v,
	Real x0, Real x1, size_t nx, Real y0, Real y1, size_t ny, Real z0, Real z1, size_t nz,
	Real en, Real directionx, Real directiony, Real directionz,
	SolverMethod met, Real mass, Real hbar, Options const & opts)
{
	fImpl = new ScatteringProblemSolverInverseMatrix3D();
	static_cast<ScatteringProblemSolverInverseMatrix3D*>(fImpl)->InitScatteringSolver3D(v,
		x0, x1, nx,
		y0, y1, ny,
		z0, z1, nz,
		en, directionx, directiony, directionz,
		SolverMethod::MatrixInverse, mass, hbar, *opts.fOpts);

}

void QuPerturbation3D::init(std::function<Complex(Real, Real, Real)> const & v,
	Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny,
	Real z0, Real z1, size_t nz,
	Real en, Real epsilon,
	Real directionx, Real directiony, Real directionz,
	SolverMethod met, Real mass, Real hbar, Options const & opts)
{
	fImpl = new QuPerturbation3DImpl();
	static_cast<QuPerturbation3DImpl*>(fImpl)->InitPerturbation3D(v,
		x0, x1, nx,
		y0, y1, ny,
		z0, z1, nz,
		en, epsilon,
		directionx, directiony, directionz,
		SolverMethod::BornSerise, mass, hbar, *opts.fOpts);
}





Options::Options()
{
	fOpts = new OptionsImpl();
}

Options::Options(Options const &r)
{
	fOpts = new OptionsImpl(*r.fOpts);
}

Options & Options::operator=(Options const &r)
{
	if(fOpts) delete fOpts;
	fOpts = new OptionsImpl(*r.fOpts);
	return *this;
}

Options::~Options()
{
	if(fOpts) delete fOpts;
}

Real Options::GetReal(char const *key)
{
	return fOpts->GetReal(key);
}

Complex Options::GetComplex(char const *key)
{
	return fOpts->GetComplex(key);
}

Int Options::GetInt(char const *key)
{
	return fOpts->GetInt(key);
}

char const * Options::GetString(char const *key)
{
	return fOpts->GetString(key).c_str();
}

void Options::SetBool(char const *key, bool v)
{
	fOpts->SetBool(key, v);
}

void Options::SetReal(char const *key, Real v)
{
	fOpts->SetReal(key, v);
}

void Options::SetComplex(char const *key, Complex v)
{
	fOpts->SetComplex(key, v);
}

void Options::SetInt(char const *key, Int v)
{
	fOpts->SetInt(key, v);
}

void Options::SetString(char const *key, char const * v)
{
	fOpts->SetString(key, v);
}
