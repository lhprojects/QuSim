#include "EigenMethod.h"

EigenMethod::EigenMethod()
{
}

void EigenMethod::initSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar,
	std::map<std::string, std::string> const &opts)
{
	EvolverImpl1D::initSystem1D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar, opts);


	if (b == BoundaryCondition::InfiniteWall) {
		throw std::runtime_error("not supported boundary conditions!");
	} else if(b == BoundaryCondition::Period) {

		fH.resize(fN, fN);
		fH.setZero();
		Real f = -fHbar * fHbar / (fDx*fDx * 2 * fMass);

		if (0) {
			fH(0, fN - 1) = f;
			fH(0, 0) = -2 * f;
			fH(0, 1) = f;
			fH(0, 0) += fV[0];
			for (ptrdiff_t i = 1; i < (ptrdiff_t)fN - 1; ++i) {
				fH(i, i - 1) = f;
				fH(i, i + 1) = f;
				fH(i, i + 0) = -2 * f;
				fH(i, i + 0) += fV[i];
			}
			fH(fN - 1, fN - 1) = -2 * f;
			fH(fN - 1, fN - 2) = f;
			fH(fN - 1, 0) = f;
			fH(fN - 1, fN - 1) += fV[fN - 1];
		} else {
			for (ptrdiff_t i = 0; i < (ptrdiff_t)fN; ++i) {
				fH(i, i - 2 < 0 ? fN + i - 2: i - 2) = -1./12 * f;
				fH(i, i - 1 < 0 ? fN + i - 1 : i - 1) = 4./3 * f;
				fH(i, i + 0) = -5./2 * f;
				fH(i, i + 1 > (ptrdiff_t)fN - 1 ? i + 1 - fN : i + 1) = 4./3 * f;
				fH(i, i + 2 > (ptrdiff_t)fN - 1 ? i + 2 - fN : i + 2) = -1./12 * f;
				fH(i, i + 0) += fV[i];
			}
		}
		fSolver.compute(fH);

		Eigen::VectorXd const &eigenvalues = fSolver.eigenvalues();
		expDt.resize(fN);

		//dump_matrix_real(fH, "H.txt");
		//dump_real(eigenvalues, "eigen.txt");


		for (ptrdiff_t i = 0; i < (ptrdiff_t)eigenvalues.size(); ++i) {
			expDt[i] = exp(-I * eigenvalues[i] * fDt / fHbar);
		}
		//dump_comp(expDt, "expDt.txt");

		psi0.resize(fN);

		for (ptrdiff_t i = 0; i < (ptrdiff_t)eigenvalues.size(); ++i) {
			psi0(i) = fPsi[i];
		}

		psi0_eigenspace = (psi0.transpose() * fSolver.eigenvectors()).transpose();
		//dump_matrix_comp(psi0_eigenspace, "psi0_eigenspace.txt");

		this->psi = psi0;
		this->psi_eigenspace = psi0_eigenspace;

	} else  {
		throw std::runtime_error("not supported boundary conditions!");
	}


}

void EigenMethod::update_psi()
{
	for (ptrdiff_t i = 0; i < (ptrdiff_t)psi.size(); ++i) {
		psi_eigenspace[i] *= expDt(i);
	}
	psi = fSolver.eigenvectors() * psi_eigenspace;
	//dump_comp(psi_eigenspace, "psi_eigenspace.txt");
	//dump_comp(psi, "psi.txt");
	for (ptrdiff_t i = 0; i < (ptrdiff_t)psi.size(); ++i) {
		fPsi[i] = psi(i);
	}
	//dump_comp(fPsi, "fPsi.txt");
}
