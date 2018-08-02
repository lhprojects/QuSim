#define _SCL_SECURE_NO_WARNINGS

#include "GaussLegendreMethod.h"

void GaussLegendreMethod::initSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization, Complex dt,
	bool force_normalization_each_step, std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver, Real mass, Real hbar,
	std::map<std::string, std::string> const &opts)
{
	EvolverImpl1D::initSystem1D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar, opts);

	if (b != BoundaryCondition::Period) {
		throw std::runtime_error("not supported boundary condition");
	}

	//            1 - 1/2 i H Dt /hbar - 1/12 (H Dt /hbar)^2 + ... 
	//  psi  ->  --------------------------------------------------   psi
	//            1 + 1/2 i H Dt /hbar - 1/12 (H Dt /hbar)^2 + ...

	// in case of split time

	//            1 - 1/2 i H Dt/2 /hbar - 1/12 (H Dt/2 /hbar)^2 + ...     1 - 1/2 i H Dt/2 /hbar - 1/12 (H Dt/2 /hbar)^2 + ... 
	//  psi  ->  ------------------------------------------------------ . ------------------------------------------------------   psi
	//            1 + 1/2 i H Dt/2 /hbar - 1/12 (H Dt/2 /hbar)^2 + ...     1 + 1/2 i H Dt/2 /hbar - 1/12 (H Dt/2 /hbar)^2 + ...

	fh.resize(fN, fN); // = H Dt / hbar

	bool space_O2 = false;
	if (fOpts.find("space_O2") != fOpts.end()) {
		space_O2 = fOpts.find("space_O2")->second != "0";
	} else if (SolverMethod::ImplicitMidpointMethod == fSolverMethod) {
		space_O2 = true;
	}

	Complex hdt = fDt;

	if (space_O2) {
		for (int i = 0; i < fN; ++i) {
			fh.insert(i, i) = (fV[i] + hbar * hbar / (2 * fMass) * 2 / (fDx*fDx))*hdt / hbar;
			fh.insert(i, i + 1 >= fN ? 0 : i + 1) = (hbar * hbar / (2 * fMass) * (-1) / (fDx*fDx))*hdt / hbar;
			fh.insert(i, i - 1 < 0 ? fN - 1 : i - 1) = (hbar * hbar / (2 * fMass) * (-1) / (fDx*fDx))*hdt / hbar;
		}
	} else {
		Complex f = -(hbar * hbar / (2 * fMass) * 1 / (fDx*fDx))*hdt / hbar;
		for (int i = 0; i < fN; ++i) {

			fh.insert(i, i - 2 < 0 ? fN + i - 2 : i - 2) = -1. / 12 * f;
			fh.insert(i, i - 1 < 0 ? fN + i - 1 : i - 1) = 4. / 3 * f;
			fh.insert(i, i + 0) = fV[i] * hdt / hbar - 5. / 2 * f;
			fh.insert(i, i + 1 > fN - 1 ? i + 1 - fN : i + 1) = 4. / 3 * f;
			fh.insert(i, i + 2 > fN - 1 ? i + 2 - fN : i + 2) = -1. / 12 * f;

		}
	}


	Eigen::SparseMatrix<Complex> id(fN, fN);
	id.setIdentity();

	if (SolverMethod::GaussLegendreO6 == fSolverMethod) {
		fn = id - 0.5 * fh * I - 1 / 10.0 * (fh*fh) + 1 / 120.0*I*(fh*fh*fh);
		fd = id + 0.5 * fh * I - 1 / 10.0 * (fh*fh) - 1 / 120.0*I*(fh*fh*fh);
	} else if (SolverMethod::GaussLegendreO4 == fSolverMethod) {
		fn = id - 0.5 * fh * I - 1 / 12.0 * (fh*fh);
		fd = id + 0.5 * fh * I - 1 / 12.0 * (fh*fh);
	} else {
		fn = id - 0.5 * fh * I;
		fd = id + 0.5 * fh * I;
	}

	fLU.compute(fd);

	ffnPsi.resize(fN);
}

void GaussLegendreMethod::update_psi()
{

	for (int i = 0; i < fN; ++i) {
		ffnPsi[i] = fPsi[i];
	}

	ffnPsi = fn * ffnPsi;
	ffnPsi = fLU.solve(ffnPsi);

	for (int i = 0; i < fN; ++i) {
		fPsi[i] = ffnPsi[i];
	}


}
