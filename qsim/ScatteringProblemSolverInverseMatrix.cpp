#define  _CRT_SECURE_NO_WARNINGS
#include "ScatteringProblemSolverInverseMatrix.h"
void ScatteringProblemSolverInverseMatrix1D::InitScatteringSolver1D(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n, Real en, Real direction, SolverMethod met,
	Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{
	ScatteringSolver1DImpl::InitScatteringSolver1D(v, x0, x1, n, en, direction, met, mass, hbar, opts);
	
	{
		int order = 2;
		auto it = opts.find("space_order");
		if (it != opts.end()) {
			if (sscanf(it->second.c_str(), "%d", &order) < 1) {
				throw std::runtime_error("can't parse space order");
			}
		}
		const_cast<int&>(fOrder) = order;
	}

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

			Real lambda = 2 * Pi / sqrt(2 * fE*fMass) * fHbar;
			Real x = GetX(i);

			Complex asb = 0;

			if (true) {

				if (4 * lambda * 4 * 2 > fNx *fDx) {
					throw std::runtime_error("too small size of to fill absorbtion layer");
				}
				Real xx;
				if (i < (ptrdiff_t)fNx / 2) xx = (x - fX0) / (4 * lambda);
				else  xx = ((fX0 + fDx * fNx) - x) / (4 * lambda);
				asb = -I * 1.5 *fE* exp(-xx * xx);
			} else {
				Real xx = 0;
				Real L = 100 / fK0;
				if (2 * L > fNx *fDx) {
					throw std::runtime_error("too small size of to fill absorbtion layer");
				}
				if (x - fX0 < L) xx = fX0 + L - x;
				else if (x > fX0 + fNx * fDx - L) xx = x - (fX0 + fNx * fDx - L);

				if (xx > 0) {
					Real alpha = fK0 / 3;
					Real alx = alpha * xx;
					asb = -0.5*fMass*fHbar*fHbar*alpha * alpha*(6 - alx + 2.*I*fK0*xx)*pow(alx, 6 - 1);
					asb /= 1 + alx * (1 + alx / 2 * (1 + alx / 3 * (1 + alx / 4 * (1 + alx / 5 * (1 + alx / 6)))));
					asb /= (6 * 5 * 4 * 3 * 2 * 1);
				}
			}

			Real t = fHbar * fHbar / (2 * fMass) * 1 / (fDx*fDx);
			Complex dig = fE - (fV[i] + asb);
			if (fOrder <= 2) {
				Eigen::Triplet<Complex> tr1((int)i, (int)fold(i - 1), t);
				elems.push_back(tr1);

				Eigen::Triplet<Complex> tr0((int)i, (int)i, -2 * t + dig);
				elems.push_back(tr0);

				Eigen::Triplet<Complex> tr2((int)i, (int)fold(i + 1), t);
				elems.push_back(tr2);
			} else if (fOrder <= 4) {
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i - 2), -1. / 12 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i - 1), 16. / 12 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 0), -30. / 12 * t + dig));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 1), 16. / 12 * t));
				elems.push_back(Eigen::Triplet<Complex>((int)i, (int)fold(i + 2), -1. / 12 * t));
			} else if (fOrder <= 6) {
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

	for (size_t i = 0; i < fNx; ++i) {
		v(i) *= fV[i];
	}

	Eigen::VectorXcd v1;
	fMatrixSolver.Solve(fEMinusH, v, v1);
	
	for (size_t i = 0; i < fNx; ++i) fPsiX[i] = v1(i);

	if (0) {
		fR = abs2(fPsiX[fNx / 2 - fNx / 8]);
		fT = abs2(fPsiX[fNx / 2 + fNx / 8]);

	} else {
		Complex r = 0;
		Complex t = 0;
		for (size_t i = 0; i < fNx; ++i) {
			r += (fPsi0X[i] + fPsiX[i])*fV[i] * exp(+I * fK0*GetX(i));
			t += (fPsi0X[i] + fPsiX[i])*fV[i] * exp(-I * fK0*GetX(i));
		}
		fR = abs(r / (2 * fK0 * I));
		fT = abs(t / (2 * fK0 * I));

	}

}
