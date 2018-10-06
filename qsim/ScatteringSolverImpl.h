#pragma once

#include "QuSim.h"
#include "Linear.h"
#include "OptionsImpl.h"
#include "eigen/Eigen/Dense"

struct ScatteringSolverImpl {

	ScatteringSolverImpl() : fHbar(0), fMass(0), fE(0), fMet(), fOpts(), fK0() { }
	virtual ~ScatteringSolverImpl() { }

	void InitScatteringSolver(Real en, SolverMethod met, Real mass, Real hbar, OptionsImpl const &opts)
	{
		
		const_cast<Real&>(fE) = en;
		const_cast<SolverMethod&>(fMet) = met;
		const_cast<Real&>(fMass) = mass;
		const_cast<Real&>(fHbar) = hbar;
		const_cast<OptionsImpl&>(fOpts) = opts;
		const_cast<Real&>(fK0) = sqrt(2 * fMass * fE) / fHbar;
	}
	virtual void Compute() = 0;

	Real GetMomentum() { return sqrt(2 * fMass*fE); }

	Real const fE;
	Real const fK0;
	SolverMethod const fMet;
	Real const fMass;
	Real const fHbar;
	OptionsImpl const fOpts;
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
		OptionsImpl const &opts);

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
		fK0Y(), fK0X(), fV(), fPsi0X() {}

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
		OptionsImpl const &opts);

	Real GetX(ptrdiff_t i) const { return fX0 + fDx * i; }
	Real GetY(ptrdiff_t i) const { return fY0 + fDy * i; }

	ptrdiff_t Idx(ptrdiff_t iy, ptrdiff_t ix) const { return ix * fNy + iy; }

	ptrdiff_t IdxFold(ptrdiff_t iy, ptrdiff_t ix) const
	{
		if (iy < 0) iy += fNy;
		else if (iy >= (ptrdiff_t)fNy) iy -= fNy;

		if (ix < 0) ix += fNx;
		else if (ix >= (ptrdiff_t)fNx) ix -= fNx;

		return Idx(iy, ix);
	}

	virtual Real ComputeXSection(Real cosx, Real cosy);
	virtual Real ComputeTotalXSection(Int n);

	size_t const fNx;
	size_t const fNy;
	Real const fX0;
	Real const fY0;
	Real const fDx;
	Real const fDy;
	std::function<Complex(Real, Real)> const fVFunc;
	Real const fK0X;
	Real const fK0Y;
	Eigen::MatrixXd const fV;
	Eigen::MatrixXcd const fPsi0X;
	Eigen::MatrixXcd fPsiX;
	Eigen::MatrixXcd flastPsiX;
	Real fNormDeltaPsi;
private:
	void InitPotential();
};

struct ScatteringSolver3DImpl : ScatteringSolverImpl {

	ScatteringSolver3DImpl() : fNx(0), fNy(0), fNz(0), fVFunc(), fX0(), fY0(), fZ0(), fDx(), fDy(), fDz(),
		fK0X(), fK0Y(), fK0Z(), fV(), fPsi0X() { }

	virtual void InitScatteringSolver3D(
		std::function<Complex(Real, Real, Real)> const & v,
		Real x0,
		Real x1,
		size_t nx,
		Real y0,
		Real y1,
		size_t ny,
		Real z0,
		Real z1,
		size_t nz,
		Real en,
		Real directionx,
		Real directiony,
		Real directionz,
		SolverMethod met,
		Real mass,
		Real hbar,
		OptionsImpl const &opts);

	Real GetX(ptrdiff_t i) const { return fX0 + fDx * i; }
	Real GetY(ptrdiff_t i) const { return fY0 + fDy * i; }
	Real GetZ(ptrdiff_t i) const { return fZ0 + fDz * i; }

	ptrdiff_t Idx(ptrdiff_t iz, ptrdiff_t iy, ptrdiff_t ix) const { return (ix * fNy + iy) * fNz + iz; }

	ptrdiff_t IdxFold(ptrdiff_t iz, ptrdiff_t iy, ptrdiff_t ix) const {
		if (iz < 0) iz += fNz;
		else if (iz >= (ptrdiff_t)fNz) iz -= fNz;

		if (iy < 0) iy += fNy;
		else if (iy >= (ptrdiff_t)fNy) iy -= fNy;

		if (ix < 0) ix += fNx;
		else if (ix >= (ptrdiff_t)fNx) ix -= fNx;

		return Idx(iz, iy, ix);
	}

	// cosx = cos(psi)
	// cosy = sin(psi)
	// cosz = cos(theta)
	virtual Real ComputeXSection(Real cosx, Real cosy, Real cosz);
	// npsi = number of sampling points for psi
	// ntheta = number of sampling points for theta
	virtual Real ComputeTotalXSection(Int npsi, Int ntheta);

	size_t const fNx;
	size_t const fNy;
	size_t const fNz;
	Real const fX0;
	Real const fY0;
	Real const fZ0;
	Real const fDx;
	Real const fDy;
	Real const fDz;
	std::function<Complex(Real, Real, Real)> const fVFunc;
	Real const fK0X;
	Real const fK0Y;
	Real const fK0Z;

	Complex const FV;

	Eigen::VectorXd const fV;
	bool fVComputeInTime;
	Eigen::VectorXcd const fPsi0X;
	bool fPsi0XComputeInTime;
	Eigen::VectorXcd fPsiX;
	static const bool fLastPsiXRecord = false;
	std::vector<Complex> fLastPsiX;
private:
	void InitPotential();
};
