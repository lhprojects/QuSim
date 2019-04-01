#include "ScatteringProblemSolverInverseMatrix2D.h"
#include "View.h"

void ScatteringProblemSolverInverseMatrix2D::InitScatteringSolver2D(std::function<Complex(Real, Real)> const & v,
	Real x0, Real x1, size_t nx, Real y0, Real y1, size_t ny, Real en,
	Real directionx, Real directiony, SolverMethod met, Real mass,
	Real hbar, OptionsImpl const & opts)
{

	ScatteringSolver2DImpl::InitScatteringSolver2D(v, x0, x1, nx, y0, y1, ny, en,
		directionx, directiony, met, mass, hbar, opts);

	const_cast<int&>(fOrder) = (int)opts.GetInt("space_order", 2);

	fMatrixSolver.Init(opts);

	fEMinusH.resize(fNx*fNy, fNx*fNy);

	{
		auto foldx = [&](ptrdiff_t i) -> ptrdiff_t {
			if (i < 0) i += fNx;
			else if (i >= (ptrdiff_t)fNx) i -= fNx;
			return i;
		};

		auto foldy = [&](ptrdiff_t i) -> ptrdiff_t {
			if (i < 0) i += fNy;
			else if (i >= (ptrdiff_t)fNy) i -= fNy;
			return i;
		};

		auto global = [&](ptrdiff_t y, ptrdiff_t x)->int {
			return (int)(foldy(y) + foldx(x) * fNy);
		};

		std::vector<Eigen::Triplet<Complex> > elems;
		for (ptrdiff_t i = 0; i < (ptrdiff_t)fNx; ++i) {
			for (ptrdiff_t j = 0; j < (ptrdiff_t)fNy; ++j) {

				Real const lambda = 2 * Pi / sqrt(2 * fE*fMass) * fHbar;
				Real const x = GetX(i);
				Real const y = GetY(j);

				Complex asb = 0;
				{

					if (3 * lambda * 3 * 2 > fNx *fDx || 3 * lambda * 3 * 2 > fNy *fDy) {
						char buf[1024];
						sprintf(buf, "%s: lambda: %f; absorbtion layer thickness*2: %f; space size: %f*%f\n",
							"no enough room to fill absorbtion layer", lambda, 18 * lambda, fNx*fDx, fNy*fDy);
						throw std::runtime_error(buf);
					}
					Real xx;
					Real yy;
					if (i < (ptrdiff_t)fNx / 2) xx = (x - fX0) / (3 * lambda);
					else  xx = ((fX0 + fDx * fNx) - x) / (3 * lambda);

					if (j < (ptrdiff_t)fNy / 2) yy = (y - fY0) / (3 * lambda);
					else  yy = ((fY0 + fDy * fNy) - y) / (3 * lambda);

					asb = -I * 1.5 *fE* (exp(-yy * yy) + exp(-xx * xx));
					//const_cast<Eigen::MatrixXd&>(fV)(j, i) += 1.5 *fE* (exp(-yy * yy) + exp(-xx * xx));
				}

				Complex const dig = fE - (fV(j, i) + asb);
				Real const tx = 1. / (2 * fMass)*fHbar*fHbar / (fDx*fDx);
				Real const ty = 1. / (2 * fMass)*fHbar*fHbar / (fDy*fDy);
				if (fOrder <= 2) {

					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i), -2 * tx - 2 * ty + dig));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i - 1), tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i + 1), tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j - 1, i), ty));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j + 1, i), ty));

				} else if (fOrder <= 4) {

					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i), -30. / 12 * tx - 30. / 12 * ty + dig));

					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i - 1), 16. / 12 * tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i + 1), 16. / 12 * tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i - 2), -1. / 12 * tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i + 2), -1. / 12 * tx));

					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j - 1, i), 16. / 12 * ty));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j + 1, i), 16. / 12 * ty));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j - 2, i), -1. / 12 * ty));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j + 2, i), -1. / 12 * ty));
				} else if (fOrder <= 6) {

					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i), -490. / 180 * tx - 490. / 180 * ty + dig));

					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i - 1), 270. / 180 * tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i + 1), 270. / 180 * tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i - 2), -27. / 180 * tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i + 2), -27. / 180 * tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i - 3), 2. / 180 * tx));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j, i + 3), 2. / 180 * tx));

					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j - 1, i), 270. / 180 * ty));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j + 1, i), 270. / 180 * ty));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j - 2, i), -27. / 180 * ty));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j + 2, i), -27. / 180 * ty));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j - 3, i), 2. / 180 * ty));
					elems.push_back(Eigen::Triplet<Complex>(global(j, i), global(j + 3, i), 2. / 180 * ty));

				} else {
					throw std::runtime_error("too higher order");
				}


				

			}
		}
		fEMinusH.setFromTriplets(elems.begin(), elems.end());
	}
	fEMinusH.makeCompressed();

}


void ScatteringProblemSolverInverseMatrix2D::Compute()
{

	Eigen::VectorXcd v;
	v.resize(fNx*fNy);
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			v(j + i * fNy) = fPsi0X(j, i) * fV(j, i);
		}
	}

	Eigen::VectorXcd v1;
	fMatrixSolver.Solve(fEMinusH, v, v1);

	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			fPsiX(j, i) = v1(j + i * fNy);
		}
	}

}

