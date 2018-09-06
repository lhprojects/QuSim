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

