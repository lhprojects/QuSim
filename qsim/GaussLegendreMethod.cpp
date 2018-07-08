
#include "GaussLegendreMethod.h"

void GaussLegendreMethod::initSystem1D(char const * psi, bool force_normalization, Complex dt,
	bool force_normalization_each_step, char const * vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver, Real mass, Real hbar)
{
	SystemImpl1D::initSystem1D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar);

	if (b != BoundaryCondition::Period) {
		throw std::runtime_error("not supported boundary condition");
	}

	//            1 - 1/2 i H Dt /hbar - 1/12 (H Dt /hbar)^2 
	//  psi  ->  --------------------------------------------   psi
	//            1 + 1/2 i H Dt /hbar - 1/12 (H Dt /hbar)^2

	fh.resize(fN, fN); // = H Dt / hbar
	for (int i = 0; i < fN; ++i) {
		fh.insert(i, i) = (fV[i] + hbar * hbar / (2 * fMass) * 2 / (fDx*fDx))*fDt / hbar;
		fh.insert(i, i + 1 >= fN ? 0 : i + 1) = (hbar * hbar / (2 * fMass) * (-1) / (fDx*fDx))*fDt / hbar;
		fh.insert(i, i - 1 < 0 ? fN - 1 : i - 1) = (hbar * hbar / (2 * fMass) * (-1) / (fDx*fDx))*fDt / hbar;
	}


	Eigen::SparseMatrix<Complex> id(fN, fN);
	id.setIdentity();

	if (SolverMethod::GaussLegendreO4 == fSolverMethod) {
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
