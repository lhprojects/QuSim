#include "ScatteringProblemSolverInverseMatrix3D.h"

void ScatteringProblemSolverInverseMatrix3D::InitScatteringSolver3D(std::function<Complex(Real, Real, Real)> const & v,
	Real x0, Real x1, size_t nx, Real y0, Real y1, size_t ny, Real z0, Real z1, size_t nz,
	Real en, Real directionx, Real directiony, Real directionz,
	SolverMethod met, Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{

	ScatteringSolver3DImpl::InitScatteringSolver3D(v, x0, x1, nx, y0, y1, ny, z0, z1, nz, en,
		directionx, directiony, directionz, met, mass, hbar, opts);


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

	fEMinusH.resize(fNx*fNy*fNz, fNx*fNy*fNz);


}

void ScatteringProblemSolverInverseMatrix3D::Compute()
{
	InitEMinusH();
	fMatrixSolver.Solve(fEMinusH, fPsi0X, fPsiX);

}

Real ScatteringProblemSolverInverseMatrix3D::ComputeXSection(Real cosx, Real cosy, Real cosz)
{
	Real const nm = 1 / sqrt(cosx * cosx + cosy * cosy + cosz * cosz);
	cosx *= nm;
	cosy *= nm;
	cosz *= nm;

	Complex r = 0;

	// Prolbem: (nabla^2+k0^2) delta psi = 2m/hbar^2 v psi
	// Green's: (nabla^2+k0^2) G = delta(x)
	//     => : delta psi  = \int G (2m / hbar^2) v psi dx dy dz
	// G = 1/(4 Pi) exp(i k z) / z
	//     => : delta psi  = \int psi v exp(i k r) / r   1/(4 Pi) * [(2m / hbar^2) dx dy dz

	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			for (size_t k = 0; k < fNz; ++k) {
				r += (fPsi0X(Idx(k, j, i)) + fPsiX(Idx(k, j, i)))
					*fV(Idx(k, j, i))
					* exp(-I * (fK0 * cosx * GetX(i) + fK0 * cosy * GetY(j) + fK0 * cosz * GetZ(k)));
			}
		}
	}

	Real dX_dOmega = abs2(r*(2 * fMass) / (fHbar*fHbar)*fDx*fDy*fDz/(4 * Pi));

	return dX_dOmega;
}

Real ScatteringProblemSolverInverseMatrix3D::ComputeTotalXSection(Int npsi, Int ntheta)
{
	Real t = 0;
	for (Int j = 0; j < ntheta; ++j) {
		for (Int i = 0; i < npsi; ++i) {
			Real phi = i * 2 * Pi / npsi;
			Real theta = i * Pi / ntheta;
			Real cosx = sin(theta)*cos(phi);
			Real cosy = sin(theta)*sin(phi);
			Real cosz = cos(theta);
			Real xsec = ComputeXSection(cosx, cosy, cosz)*sin(theta);
			t += xsec;
		}
	}
	t *= 4 * Pi / (ntheta * npsi);
	return t;
}

void ScatteringProblemSolverInverseMatrix3D::InitEMinusH()
{
	{

		std::vector<Eigen::Triplet<Complex> > elems;
		for (ptrdiff_t i = 0; i < (ptrdiff_t)fNx; ++i) {
			for (ptrdiff_t j = 0; j < (ptrdiff_t)fNy; ++j) {
				for (ptrdiff_t k = 0; k < (ptrdiff_t)fNz; ++k) {

					Real const lambda = 2 * Pi / sqrt(2 * fE*fMass) * fHbar;
					Real const x = GetX(i);
					Real const y = GetY(j);
					Real const z = GetY(k);

					Complex asb = 0;
					{

						if (3 * lambda * 3 * 2 > fNx *fDx || 3 * lambda * 3 * 2 > fNy *fDy || 3 * lambda * 3 * 2 > fNz * fDz) {
							throw std::runtime_error("too small size of to fill absorbtion layer");
						}
						Real xx;
						Real yy;
						Real zz;
						if (i < (ptrdiff_t)fNx / 2) xx = (x - fX0) / (4 * lambda);
						else  xx = ((fX0 + fDx * fNx) - x) / (4 * lambda);

						if (j < (ptrdiff_t)fNy / 2) yy = (y - fY0) / (4 * lambda);
						else  yy = ((fY0 + fDy * fNy) - y) / (4 * lambda);

						if (k < (ptrdiff_t)fNz / 2) zz = (z - fZ0) / (4 * lambda);
						else  zz = ((fZ0 + fDz * fNz) - z) / (4 * lambda);

						asb = -I * 1.5 *fE* (exp(-xx * xx) + exp(-yy * yy) + exp(-zz * zz));
					}

					Complex const dig = fE - (fV(Idx(k, j, i)) + asb);
					Real const tx = 1. / (2 * fMass)*fHbar*fHbar / (fDx*fDx);
					Real const ty = 1. / (2 * fMass)*fHbar*fHbar / (fDy*fDy);
					Real const tz = 1. / (2 * fMass)*fHbar*fHbar / (fDz*fDz);
					if (fOrder <= 2) {

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i), -2 * (tx + ty + tz) + dig));

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i - 1), tx));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i + 1), tx));

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j - 1, i), ty));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j + 1, i), ty));

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k - 1, j, i), tz));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k + 1, j, i), tz));

					} else if (fOrder <= 4) {

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i), -30. / 12 * (tx + ty + tz) + dig));

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i - 1), 16. / 12 * tx));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i + 1), 16. / 12 * tx));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i - 2), -1. / 12 * tx));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i + 2), -1. / 12 * tx));

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j - 1, i), 16. / 12 * ty));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j + 1, i), 16. / 12 * ty));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j - 2, i), -1. / 12 * ty));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j + 2, i), -1. / 12 * ty));

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k - 1, j, i), 16. / 12 * tz));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k + 1, j, i), 16. / 12 * tz));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k - 2, j, i), -1. / 12 * tz));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k + 2, j, i), -1. / 12 * tz));
					} else if (fOrder <= 6) {

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i), -490. / 180 * (tx + ty + tz) + dig));

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i - 1), 270. / 180 * tx));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i + 1), 270. / 180 * tx));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i - 2), -27. / 180 * tx));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i + 2), -27. / 180 * tx));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i - 3), 2. / 180 * tx));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j, i + 3), 2. / 180 * tx));

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j - 1, i), 270. / 180 * ty));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j + 1, i), 270. / 180 * ty));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j - 2, i), -27. / 180 * ty));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j + 2, i), -27. / 180 * ty));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j - 3, i), 2. / 180 * ty));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k, j + 3, i), 2. / 180 * ty));

						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k - 1, j, i), 270. / 180 * tz));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k + 1, j, i), 270. / 180 * tz));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k - 2, j, i), -27. / 180 * tz));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k + 2, j, i), -27. / 180 * tz));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k - 3, j, i), 2. / 180 * tz));
						elems.push_back(Eigen::Triplet<Complex>(Idx(k, j, i), IdxFold(k + 3, j, i), 2. / 180 * tz));

					} else {
						throw std::runtime_error("too higher order");
					}


				}
			}
		}
		fEMinusH.setFromTriplets(elems.begin(), elems.end());
	}
	fEMinusH.makeCompressed();
}
