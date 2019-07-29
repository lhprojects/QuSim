#include "EvolverImpl.h"


void EvolverImpl::initSystem(bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	BoundaryCondition b,
	SolverMethod solver,Real mass, Real hbar, OptionsImpl const &opts)
{
	fStep = 0;

	fFN = force_normalization;

	fBoundaryCondition = b;

	fDt = dt;
	fFNES = force_normalization_each_step;
	fMass = mass;
	fHbar = hbar;

	fSolverMethod = solver;
	const_cast< OptionsImpl &>(fOpts) = opts;
}












Real EvolverImpl1D::CalPotEn()
{
	Real norm2 = 0;
	for (size_t i = 0; i < fN; ++i) {
		norm2 += abs2(fPsi[i])*fV[i] * fDx;
	}
	return norm2 / Norm2();
}

Real EvolverImpl1D::CalKinEn()
{
	return 0;
}

Real EvolverImpl1D::PotEn()
{
	return CalPotEn();
}

Real EvolverImpl1D::KinEn()
{
	return CalKinEn();
}

Real EvolverImpl1D::EnPartialT()
{
	Complex en = 0;
	for (size_t i = 0; i < fN; ++i) {
		en += I * fHbar*conj(fPsi[i])*(fPsi[i] - fLastLastPsi[i]) / (2.0*fDt)*fDx;
	}
	return en.real() / Norm2();
}

void EvolverImpl1D::step()
{
	fLastLastPsi = fLastPsi;
	fLastPsi = fPsi;

	update_psi();

	if (fFNES) {
		Scale(fPsi, 1.0 / sqrt(Norm2()));
	}
	++fStep;

}

Real EvolverImpl1D::Norm2()
{
	return Norm2(fPsi);
}

void EvolverImpl1D::initSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real)> const &v, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar, OptionsImpl const &opts)
{

	initSystem(force_normalization,
		dt, force_normalization_each_step,
		b, solver,
		mass, hbar, opts);

	fPsi0Func = psi;
	fVFunc = v;

	fX0 = x0;
	fDx = (x1 - x0) / n;
	fN = n;

	initPsi();
	fLastLastPsi = fPsi;
	fLastPsi = fPsi;

	initPotential();

}

void EvolverImpl1D::initPotential()
{
	fV.resize(fN);

	for (size_t i = 0; i < fN; ++i) {
		Real x = getX(i);
		Complex com = fVFunc(x);
		fV[i] = com.real();

	}
	//dump(fV, "V.txt");
}

void EvolverImpl1D::initPsi()
{
	fPsi.resize(fN);
	for (size_t i = 0; i < fN; ++i) {
		Real x = getX(i);
		Complex com = fPsi0Func(x);
		fPsi[i] = com;
	}

	if (fFN) {
		Scale(fPsi, 1 / sqrt(Norm2()));
	}
	//dump(fPsi, "psi.txt");
}















void EvolverImpl2D::initSystem2D(std::function<Complex(Real, Real)> const &psi,
	bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real, Real)> const &v, Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny,
	BoundaryCondition b,
	SolverMethod solver, Real mass, Real hbar, OptionsImpl const &opts)
{
	initSystem(force_normalization, dt, force_normalization_each_step,
		b, solver, mass, hbar, opts);

	fVFunc = v;
	fPsi0Func = psi;

	fX0 = x0;
	fDx = (x1 - x0) / nx;
	fNx = nx;
	fY0 = y0;
	fDy = (y1 - y0) / ny;
	fNy = ny;

	initPsi();
	fLastStep = fStep;
	fLastPsi = fPsi;
	//fLastLastPsi = fLastPsi;

	initPotential();
	//Real x = Norm2();

}

void EvolverImpl2D::initPsi()
{
	fPsi.resize(fNy, fNx);
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			Real x = getX(i);
			Real y = getY(j);
			Complex com = fPsi0Func(x, y);
			fPsi(j, i) = com;
		}
	}

	if (fFN) {
		fPsi *= 1.0 / sqrt(Norm2());
	}
	//double x = Norm2();

}

void EvolverImpl2D::initPotential()
{
	fV.resize(fNy, fNx);
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			Real x = getX(i);
			Real y = getY(j);
			Complex com = fVFunc(x, y);
			fV(j, i) = com.real();
		}
	}
	//dump_matrix_real(fV, "pot.txt");
	//Real x = Norm2();

}


void EvolverImpl2D::step()
{
	//Real y = Norm2();
	//fLastLastPsi = fLastPsi;
	fLastStep = fStep;
	fLastPsi = fPsi;

	//Real x = Norm2();
	update_psi();

	if (fFNES) {
		fPsi *= 1.0 / sqrt(Norm2());
	}
	++fStep;

}

Real EvolverImpl2D::EnPartialT()
{

	Real pot = 0;
	Real norm2 = 0;
	UInt step = fStep - fLastStep;
	if (step == 0) return Real(0);
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			Complex dpsidt = (fPsi(j, i) - fLastPsi(j, i)) / (fDt*Real(step));
			pot += (I * fHbar * 0.5*conj(fPsi(j, i) + fLastPsi(j, i)) * dpsidt).real();
			norm2 += abs2(0.5*(fLastPsi(j, i) + fPsi(j, i)));
		}
	}
	return pot / norm2;

}


Real EvolverImpl2D::Norm2()
{
	return fPsi.squaredNorm()*fDx*fDy;
}

Real EvolverImpl2D::CalPotEn() const
{
	Real pot = 0;
	Real norm2 = 0;
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			pot += abs2(fPsi(j, i)) * fV(j, i);
			norm2 += abs2(fPsi(j, i));
		}
	}
	return pot / norm2;
}

Real EvolverImpl2D::CalKinEn() const
{
	return 0;
}

