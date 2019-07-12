#define  _CRT_SECURE_NO_WARNINGS
#include "ScatteringProblemSolverInverseMatrix1D.h"
void ScatteringProblemSolverInverseMatrix1D::InitScatteringSolver1D(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n, Real en, Real direction, SolverMethod met,
	Real mass, Real hbar, OptionsImpl const & opts)
{
	ScatteringSolver1DImpl::InitScatteringSolver1D(v, x0, x1, n, en, direction, met, mass, hbar, opts);

	const_cast<int&>(fSpaceOrder) = (int)opts.GetInt("space_order", 2);
	const_cast<bool&>(fPreferPreciseSmallWaveFunction) = opts.GetBool("PreferPreciseSmallWaveFunction", false);
	const_cast<std::vector<Real>&>(fVabs).resize(fNx);

	fMatrixSolver.Init(opts);

	fEMinusH.resize(fNx, fNx);

	{
		auto fold = [&](ptrdiff_t i) -> ptrdiff_t {
			if (i < 0) i += fNx;
			else if (i >= (ptrdiff_t)fNx) i -= fNx;
			return i;
		};
		std::vector<Eigen::Triplet<Complex> > elems;
		for (ptrdiff_t i = 0; i < (ptrdiff_t)fNx; ++i) {

			Complex asb = 0;

			if (true) {
				Real lambda = 2 * Pi / sqrt(2 * fE*fMass) * fHbar;
				Real x = GetX(i);

				if (4 * lambda * 4 * 2 > fNx *fDx) {
					throw std::runtime_error("too small size of to fill absorbtion layer");
				}
				Real xx;
				if (i < (ptrdiff_t)fNx / 2) xx = (x - fX0) / (4 * lambda);
				else  xx = ((fX0 + fDx * fNx) - x) / (4 * lambda);
				asb = -I * 1.5 *fE* exp(-xx * xx);
			}
			const_cast<std::vector<Real>&>(fVabs)[i] = asb.imag();

			Real t = fHbar * fHbar / (2 * fMass) * 1 / (fDx*fDx);
			Complex dig = fE - (fV[i] + asb);
			if (fSpaceOrder <= 2) {
				Eigen::Triplet<Complex> tr1((int)i, (int)fold(i - 1), t);
				elems.push_back(tr1);

				Eigen::Triplet<Complex> tr0((int)i, (int)i, -2 * t + dig);
				elems.push_back(tr0);

				Eigen::Triplet<Complex> tr2((int)i, (int)fold(i + 1), t);
				elems.push_back(tr2);
			} else if (fSpaceOrder <= 4) {
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i - 2), -1. / 12 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i - 1), 16. / 12 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 0), -30. / 12 * t + dig));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 1), 16. / 12 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 2), -1. / 12 * t));
			} else if (fSpaceOrder <= 6) {
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i - 3), 2. / 180 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i - 2), -27. / 180 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i - 1), 270. / 180 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 0), -490. / 180 * t + dig));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 1), 270. / 180 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 2), -27. / 180 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 3), 2. / 180 * t));
			} else {
				throw std::runtime_error("too higher order");
			}

		}
		fEMinusH.setFromTriplets(elems.begin(), elems.end());
	}
	fEMinusH.makeCompressed();

}

void ScatteringProblemSolverInverseMatrix1D::Compute()
{

	Eigen::VectorXcd v;
	v.resize(fNx);
	for (size_t i = 0; i < fNx; ++i) v(i) = fPsi0X[i];

	if (fPreferPreciseSmallWaveFunction) {
		for (size_t i = 0; i < fNx; ++i) {
			v(i) *= -I*fVabs[i];
		}
	} else {
		for (size_t i = 0; i < fNx; ++i) {
			v(i) *= fV[i];
		}
	}

	Eigen::VectorXcd v1;
	fMatrixSolver.Solve(fEMinusH, v, v1);

	if (fPreferPreciseSmallWaveFunction) {
		// v1 is the psi
		for (size_t i = 0; i < fNx; ++i) v1(i) = v1(i) - fPsi0X[i];
	}

	for (size_t i = 0; i < fNx; ++i) fPsiX[i] = v1(i);

	ComputeRT();

}
