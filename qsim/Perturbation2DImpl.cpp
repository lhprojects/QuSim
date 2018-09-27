#define _CRT_SECURE_NO_WARNINGS
#include "Perturbation2DImpl.h"

void QuPerturbation2DImpl::InitPerturbation2D(std::function<Complex(Real, Real)> const & v, Real x0,
	Real x1, size_t nx, Real y0, Real y1, size_t ny, Real en, Real epsilon,
	Real directionx, Real directiony, SolverMethod met,
	Real mass, Real hbar, std::map<std::string, std::string> const & opts)
{
	ScatteringSolver2DImpl::InitScatteringSolver2D(v, x0, x1, nx, y0, y1, ny,
		en, directionx, directiony,
		met, mass, hbar, opts);


	fFourierTransformOptions.Init(opts);
	fFFT.reset(FourierTransform2D::Create(fNy, fNx, false, fFourierTransformOptions.fLib));
	fInvFFT.reset(FourierTransform2D::Create(fNy, fNx, true, fFourierTransformOptions.fLib));
	const_cast<Real&>(fEpsilon) = epsilon;

	fPsiK.resize(fNy*fNx);

	fPsiX.setZero();

	if (fMet == SolverMethod::BornSerise) {

		{
			int the_order = 0;
			auto it = opts.find("order");
			if (it != opts.end()) {
				auto &order = it->second;
				if (sscanf(order.c_str(), "%d", &the_order) < 1) {
					throw std::runtime_error("not valid order");
				}
			}
			const_cast<int&>(fOrder) = the_order;
		}

		fPerturbationOptions.Init(opts);
		if (fPerturbationOptions.fPreconditional) {

			{
				fVasb.resize(fNx*fNy);
				Real lambda = 2 * Pi / (sqrt(2 * fMass*fE) / hbar);
				if (fNx*fDx < 24 * lambda) {
					throw std::runtime_error("too small size");
				}
				if (fNy*fDy < 24 * lambda) {
					throw std::runtime_error("too small size");
				}
				PerburbationUtility::GaussAsbLayer2D(fNx, fNy, fDx, fDy,
					fVasb.data(),
					fHbar, mass, fE, 4.0);

			}

			PreconditionalBornSerise pbs;
			Real minEpsilon = pbs.GetMinEpsilon(fNx*fNy, fV.data(), fVasb.data());
			const_cast<Real&>(fEpsilon) = (epsilon < minEpsilon ? minEpsilon : epsilon);

			ftmp1.resize(fNx*fNy);
		}

	}



}

void QuPerturbation2DImpl::Compute()
{
	auto X2K = [this](Complex const *psix, Complex *psik) {
		fFFT->Transform(psix, psik);
	};

	auto K2X = [this](Complex const *psik, Complex *psix) {
		fInvFFT->Transform(psik, psix);
		for (size_t i = 0; i < fNx; ++i) {
			for (size_t j = 0; j < fNy; ++j) {
				psix[Idx(j, i)] *= 1. / (fNx*fNy);
			}
		}
	};

	std::copy(fPsiX.data(), fPsiX.data() + fNx*fNy, flastPsiX.data());

	if (fPerturbationOptions.fPreconditional) { // Preconditional Born serise

		PreconditionalBornSerise pbs;
		for (int i = 0; i < fOrder; ++i) {
			pbs.Update2D(fNx, fNy, fPsi0X.data(), fPsiX.data(), fPsiK.data(),
				fV.data(), fVasb.data(), ftmp1.data(), fEpsilon, fE,
				fMass, fHbar, fDx, fDy,
				X2K, K2X,
				fPerturbationOptions.fSlow,
				fPerturbationOptions.fPreconditioner);
		}

	} else { // naive born serise
		BornSerise bs;

		for (int i = 0; i < fOrder; ++i) {
			bs.Update2D(fNx, fNy, fPsi0X.data(), fPsiX.data(),
				fPsiK.data(), fV.data(), fEpsilon, fE, fMass, fHbar, fDx, fDy,
				X2K, K2X);
		}
	}

	Real norm = 0;
	for (size_t i = 0; i < fNx * fNy; ++i) {
		norm += abs2(fPsiX.data()[i] - flastPsiX.data()[i])*fDx*fDy;
	}
	norm = sqrt(norm);
	fNormDeltaPsi = norm;
}
