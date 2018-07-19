#define _SCL_SECURE_NO_WARNINGS
#include "GaussLegendreMethod2D.h"

void GaussLegendreMethod2D::initSystem2D(std::function<Complex(Real, Real)> const & psi,
	bool force_normalization, Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real, Real)> const & vs, Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny, BoundaryCondition b,
	SolverMethod solver, Real mass, Real hbar,
	std::map<std::string, std::string> const & opts)
{
	EvolverImpl2D::initSystem2D(psi, force_normalization, dt, force_normalization_each_step,
		vs, x0, x1, nx,
		y0, y1, ny, b,
		solver, mass, hbar,
		opts);

	if (b != BoundaryCondition::Period) {
		throw std::runtime_error("not supported boundary condition");
	}

	fh.resize(fNx*fNy, fNx*fNy); // = H Dt / hbar

	bool space_O2 = false;
	if (fOpts.find("space_O2") != fOpts.end()) {
		space_O2 = fOpts.find("space_O2")->second != "0";
	} else if (SolverMethod::ImplicitMidpointMethod == fSolverMethod) {
		space_O2 = true;
	}

	Complex hdt = fDt;

	Complex fv = hdt / hbar;
	Complex fkx = -(hbar * hbar / (2 * fMass) * 1 / (fDx*fDx))*hdt / hbar;
	Complex fky = -(hbar * hbar / (2 * fMass) * 1 / (fDy*fDy))*hdt / hbar;

	std::vector<Eigen::Triplet<Complex> > elems;
	if (space_O2) {

		for (ptrdiff_t i = 0; i < (ptrdiff_t)fNx; ++i) {
			for (ptrdiff_t j = 0; j < (ptrdiff_t)fNy; ++j) {

				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i, j), -2.*fkx - 2.*fky + fV(j, i)*fv));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i, j+1), fky));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i, j-1), fky));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i+1, j), fkx));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i-1, j), fkx));
			}
		}

	} else {

		for (ptrdiff_t i = 0; i < (ptrdiff_t)fNx; ++i) {
			for (ptrdiff_t j = 0; j < (ptrdiff_t)fNy; ++j) {

				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i, j), -5. / 2 * fkx - 5. / 2 * fky + fV(j, i) * fv));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i - 1, j), 4. / 3 * fkx));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i + 1, j), 4. / 3 * fkx));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i - 2, j), -1. / 12 * fkx));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i + 2, j), -1. / 12 * fkx));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i, j - 1), 4. / 3 * fky));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i, j + 1), 4. / 3 * fky));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i, j - 2), -1. / 12 * fky));
				elems.push_back(Eigen::Triplet<Complex>((int)Index(i, j), (int)IndexP(i, j + 2), -1. / 12 * fky));
			}
		}

	}
	fh.setFromTriplets(elems.begin(), elems.end());

	Eigen::SparseMatrix<Complex> id(fNx*fNy, fNx*fNy);
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

	ffnPsi1.resize(fNx*fNy);
	ffnPsi2.resize(fNx*fNy);


}

void GaussLegendreMethod2D::update_psi()
{
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			ffnPsi1(Index(i, j)) = fPsi(j, i);
		}
	}

	ffnPsi2 = fn * ffnPsi1;
	ffnPsi1 = fLU.solve(ffnPsi2);

	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			fPsi(j, i) = ffnPsi1(Index(i, j));
		}
	}

}

Real GaussLegendreMethod2D::CalKinEn() const
{
	return Real();
}
