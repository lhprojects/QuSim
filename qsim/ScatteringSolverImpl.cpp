#include "ScatteringSolverImpl.h"

void ScatteringSolver1DImpl::InitScatteringSolver1D(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n, Real en, Real direction,
	SolverMethod met, Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{
	InitScatteringSolver(en, met, mass, hbar, opts);

	const_cast<size_t&>(fNx) = n;
	const_cast<Real&>(fX0) = x0;
	const_cast<Real&>(fDx) = (x1 - x0) / n;
	const_cast<std::function<Complex(Real)>&>(fVFunc) = v;
	const_cast<Real&>(fK0) = sqrt(2 * fMass * fE) / fHbar;

	InitPotential();

	const_cast<PsiVector&>(fPsi0X).resize(fNx);
	for (size_t i = 0; i < fNx; ++i) {
		const_cast<PsiVector&>(fPsi0X)[i] = exp(fK0 * GetX(i) * I);
	}
	fPsiX.resize(fNx);
}

void ScatteringSolver1DImpl::InitPotential()
{
	const_cast<std::vector<Real>&>(fV).resize(fNx);

	for (size_t i = 0; i < fNx; ++i) {
		Real x = GetX(i);
		Complex com = fVFunc(x);
		const_cast<std::vector<Real>&>(fV)[i] = com.real();

	}
}



void ScatteringSolver2DImpl::InitScatteringSolver2D(std::function<Complex(Real, Real)> const & v,
	Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny,
	Real en,
	Real directionx,
	Real directiony,
	SolverMethod met, Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{
	InitScatteringSolver(en, met, mass, hbar, opts);

	const_cast<size_t&>(fNx) = nx;
	const_cast<size_t&>(fNy) = ny;
	const_cast<Real&>(fX0) = x0;
	const_cast<Real&>(fY0) = y0;
	const_cast<Real&>(fDx) = (x1 - x0) / nx;
	const_cast<Real&>(fDy) = (y1 - y0) / ny;
	const_cast<std::function<Complex(Real, Real)>&>(fVFunc) = v;
	{
		Real const nm = 1 / sqrt(directionx*directionx + directiony * directiony);
		directionx *= nm;
		directiony *= nm;
	}
	const_cast<Real&>(fK0X) = fK0 * directionx;
	const_cast<Real&>(fK0Y) = fK0 * directiony;

	InitPotential();

	const_cast<Eigen::MatrixXcd&>(fPsi0X).resize(fNy, fNx);
	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			const_cast<Eigen::MatrixXcd&>(fPsi0X)(j, i) = exp((fK0X * GetX(i) + fK0Y * GetY(j)) * I);
		}
	}
	fPsiX.resize(fNy, fNx);
	flastPsiX.resize(fNy, fNx);
	fNormDeltaPsi = 0;
}

void ScatteringSolver2DImpl::InitPotential()
{
	const_cast<Eigen::MatrixXd&>(fV).resize(fNy, fNx);

	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			Real x = GetX(i);
			Real y = GetX(j);
			Complex com = fVFunc(x, y);
			const_cast<Eigen::MatrixXd&>(fV)(j, i) = com.real();
		}
	}
}

Real ScatteringSolver2DImpl::ComputeXSection(Real cosx, Real cosy)
{
	Real const nm = 1 / sqrt(cosx*cosx + cosy * cosy);
	cosx *= nm;
	cosy *= nm;

	Complex r = 0;

	// Prolbem: (nabla^2+k0^2) delta psi = 2m/hbar^2 v psi
	// Green's: (nabla^2+k0^2) G = delta(x)
	//     => : delta psi  = \int G (2m / hbar^2) v psi dx dy
	// G = 1/4 (1st kind Hankel, alpha = 0) (k r)
	// (1st kind Hankel, alpha = 0)(z) ~ sqrt(2/(Pi z)) exp(i k z)
	//     => : delta psi  = \int psi v exp(i k r)  [(2m / hbar^2)  1/4 sqrt(2/(Pi k r))] dx dy

	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			r += (fPsi0X(j, i) + fPsiX(j, i))
				*fV(j, i)
				* exp(-I * (fK0 * cosx * GetX(i) + fK0 * cosy * GetY(j)));
		}
	}

	Real dX_dTheta = abs2(r*(2 * fMass) / (fHbar*fHbar)*fDx*fDy*(1. / 4))*(2 / (Pi * fK0));

	return dX_dTheta;

}

Real ScatteringSolver2DImpl::ComputeTotalXSection(Int n)
{
	Real xsec_t = 0;
	for (Int i = 0; i < n; ++i) {
		Real theta = i * 2 * Pi / n;
		Real cosx = cos(theta);
		Real cosy = sin(theta);
		Real xsec = ComputeXSection(cosx, cosy);
		xsec_t += xsec;
	}
	xsec_t *= 2 * Pi / n;

	return xsec_t;
}




void ScatteringSolver3DImpl::InitScatteringSolver3D(std::function<Complex(Real, Real, Real)> const & v, Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny, Real z0, Real z1, size_t nz,
	Real en,
	Real directionx, Real directiony, Real directionz,
	SolverMethod met, Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{
	InitScatteringSolver(en, met, mass, hbar, opts);

	const_cast<size_t&>(fNx) = nx;
	const_cast<size_t&>(fNy) = ny;
	const_cast<size_t&>(fNz) = nz;
	const_cast<Real&>(fX0) = x0;
	const_cast<Real&>(fY0) = y0;
	const_cast<Real&>(fZ0) = z0;
	const_cast<Real&>(fDx) = (x1 - x0) / nx;
	const_cast<Real&>(fDy) = (y1 - y0) / ny;
	const_cast<Real&>(fDz) = (z1 - z0) / nz;
	const_cast<std::function<Complex(Real, Real, Real)>&>(fVFunc) = v;
	{
		Real const nm = 1 / sqrt(directionx * directionx + directiony * directiony + directionz * directionz);
		directionx *= nm;
		directiony *= nm;
		directionz *= nm;
	}
	const_cast<Real&>(fK0X) = fK0 * directionx;
	const_cast<Real&>(fK0Y) = fK0 * directiony;
	const_cast<Real&>(fK0Z) = fK0 * directionz;

	InitPotential();

	const_cast<Eigen::VectorXcd&>(fPsi0X).resize(fNz*fNy*fNx);
	const_cast<Eigen::VectorXcd&>(fPsiX).resize(fNz*fNy*fNx);

	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			for (size_t k = 0; k < fNz; ++k) {
				const_cast<Eigen::VectorXcd&>(fPsi0X)(Idx(k, j, i))
					= exp((fK0X * GetX(i) + fK0Y * GetY(j) + fK0Z * GetZ(k)) * I);
				const_cast<Eigen::VectorXcd&>(fPsiX)(Idx(k, j, i)) = 0;
			}
		}
	}
}

void ScatteringSolver3DImpl::InitPotential()
{
	const_cast<Eigen::VectorXd& >(fV).resize(fNx*fNy*fNz);

	for (size_t i = 0; i < fNx; ++i) {
		for (size_t j = 0; j < fNy; ++j) {
			for (size_t k = 0; k < fNz; ++k) {
				Real x = GetX(i);
				Real y = GetX(j);
				Real z = GetX(k);
				ptrdiff_t idx = Idx(k, j, i);
				Complex com = fVFunc(x, y, z);
				const_cast<Eigen::VectorXd& >(fV)(idx) = com.real();
			}
		}
	}
}
