#pragma once

#include "QuSim.h"
#include "Linear.h"
#include "eigen/Eigen/Dense"

struct ScatteringSolverImpl {

	ScatteringSolverImpl() : fHbar(0), fMass(0), fE(0), fMet(), fOpts() { }
	~ScatteringSolverImpl() { }

	void InitScatteringSolver(Real en, SolverMethod met, Real mass, Real hbar, std::map<std::string, std::string> const &opts)
	{
		
		const_cast<Real&>(fE) = en;
		const_cast<SolverMethod&>(fMet) = met;
		const_cast<Real&>(fMass) = mass;
		const_cast<Real&>(fHbar) = hbar;
		const_cast<std::map<std::string, std::string>&>(fOpts) = opts;
	}
	virtual void Compute() = 0;

	Real GetMomentum() { return sqrt(2 * fMass*fE); }

	Real const fE;
	SolverMethod const fMet;
	Real const fMass;
	Real const fHbar;
	std::map<std::string, std::string> const fOpts;
};

struct ScatteringSolver1DImpl : ScatteringSolverImpl {

	ScatteringSolver1DImpl() : fNx(0), fVFunc(), fX0(), fDx(), fK0(), fV(), fPsi0X() { }

	virtual void InitScatteringSolver1D(
		std::function<Complex(Real)> const & v,
		Real x0,
		Real x1,
		size_t n,
		Real en,
		Real direction,
		SolverMethod met,
		Real mass,
		Real hbar,
		std::map<std::string, std::string> const &opts);

	Real GetX(ptrdiff_t i) const { return fX0 + fDx * i; }


	size_t const fNx;
	Real const fX0;
	Real const fDx;
	std::function<Complex(Real)> const fVFunc;
	Real const fK0;
	std::vector<Real> const fV;
	PsiVector const fPsi0X;

	PsiVector fPsiX;
	Real fT;
	Real fR;
private:
	void InitPotential();
};

struct ScatteringSolver2DImpl : ScatteringSolverImpl {

	ScatteringSolver2DImpl() : fNx(0), fNy(0), fVFunc(), fX0(), fY0(), fDx(), fDy(),
		fK0Y(), fK0X(), fV(), fPsi0X(), fK0() {}

	virtual void InitScatteringSolver2D(
		std::function<Complex(Real, Real)> const & v,
		Real x0,
		Real x1,
		size_t nx,
		Real y0,
		Real y1,
		size_t ny,
		Real en,
		Real directionx,
		Real directiony,
		SolverMethod met,
		Real mass,
		Real hbar,
		std::map<std::string, std::string> const &opts);

	Real GetX(ptrdiff_t i) const { return fX0 + fDx * i; }
	Real GetY(ptrdiff_t i) const { return fY0 + fDy * i; }

	virtual Real ComputeXSection(Real cosx, Real cosy) = 0;
	virtual Real ComputeTotalXSection() = 0;

	Real GetMomentum()
	{
		return sqrt(2 * fMass*fE);
	}

	size_t const fNx;
	size_t const fNy;
	Real const fX0;
	Real const fY0;
	Real const fDx;
	Real const fDy;
	std::function<Complex(Real, Real)> const fVFunc;
	Real const fK0;
	Real const fK0X;
	Real const fK0Y;
	Eigen::MatrixXd const fV;
	Eigen::MatrixXcd const fPsi0X;
	Eigen::MatrixXcd fPsiX;
private:
	void InitPotential();
};
