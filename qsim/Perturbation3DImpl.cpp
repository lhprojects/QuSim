#define _CRT_SECURE_NO_WARNINGS
#include "Perturbation3DImpl.h"

void QuPerturbation3DImpl::InitPerturbation3D(std::function<Complex(Real, Real, Real)> const & v,
	Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny,
	Real z0, Real z1, size_t nz,
	Real en, Real epsilon,
	Real directionx, Real directiony, Real directionz,
	SolverMethod met,
	Real mass, Real hbar, OptionsImpl const & opts)
{
	ScatteringSolver3DImpl::InitScatteringSolver3D(v,
		x0, x1, nx,
		y0, y1, ny,
		z0, z1, nz,
		en,
		directionx, directiony, directionz,
		met, mass, hbar, opts);

	fNormDeltaPsi = 0;
	fFourierTransformOptions.Init(opts);
	fFFT.reset(FourierTransform3D::Create(fNz, fNy, fNx, false, fFourierTransformOptions.fLib));
	fInvFFT.reset(FourierTransform3D::Create(fNz, fNy, fNx, true, fFourierTransformOptions.fLib));
	const_cast<Real&>(fEpsilon) = epsilon;

	fPsiK.resize(fNz*fNy*fNx);

	fPsiX.setZero();

	if (fMet == SolverMethod::BornSerise) {

		const_cast<int&>(fOrder) = (int)opts.GetInt("order", 0);

		fPerturbationOptions.Init(opts);
		if (fPerturbationOptions.fPreconditional) {

			{
				fVasb.resize(fNx*fNy*fNz);
				Real lambda = 2 * Pi / (sqrt(2 * fMass*fE) / hbar);
				if (fNx*fDx < 24 * lambda) {
					throw std::runtime_error("too small size");
				}
				if (fNy*fDy < 24 * lambda) {
					throw std::runtime_error("too small size");
				}
				if (fNz*fDz < 24 * lambda) {
					throw std::runtime_error("too small size");
				}
				PerburbationUtility::GaussAsbLayer3D(fNx, fNy, fNz,
					fDx, fDy, fDz,
					fVasb.data(),
					fHbar, mass, fE, 4.0);

			}

			PreconditionalBornSerise pbs;
			pbs.GetEpsilon(epsilon, fPerturbationOptions.fPreconditioner,
				fNx*fNy, fV.data(), fVasb.data());
			const_cast<Real&>(fEpsilon) = epsilon;

			ftmp1.resize(fNx*fNy*fNz);
		}

	} else {
		throw std::runtime_error("not support method!");
	}



}

void QuPerturbation3DImpl::Compute()
{
	auto X2K = [this](Complex const *psix, Complex *psik) {
		fFFT->Transform(psix, psik);
	};

	auto K2X = [this](Complex const *psik, Complex *psix) {
		fInvFFT->Transform(psik, psix);
		for (size_t i = 0; i < fNx; ++i) {
			for (size_t j = 0; j < fNy; ++j) {
				for (size_t k = 0; k < fNz; ++k) {
					psix[Idx(k, j, i)] *= 1. / (fNx*fNy*fNz);
				}
			}
		}
	};

	if (fLastPsiXRecord) {
		for (size_t i = 0; i < fNx*fNy*fNz; ++i) {
			fLastPsiX.data()[i] = fPsiX.data()[i];
		}
	}

	if (fPerturbationOptions.fPreconditional) { // Preconditional Born serise

		PreconditionalBornSerise pbs;
		for (int i = 0; i < fOrder; ++i) {
			pbs.Update3D(fNx, fNy, fNz, fPsi0X.data(), fPsiX.data(), fPsiK.data(),
				fV.data(), fVasb.data(), ftmp1.data(), fEpsilon, fE,
				fMass, fHbar, fDx, fDy, fDz,
				X2K, K2X,
				fPerturbationOptions.fSlow,
				fPerturbationOptions.fPreconditioner);
		}

	} else { // naive born serise
		BornSerise bs;

		for (int i = 0; i < fOrder; ++i) {
			bs.Update3D(fNx, fNy, fNz, fPsi0X.data(), fPsiX.data(),
				fPsiK.data(), fV.data(), fEpsilon, fE, fMass, fHbar, fDx, fDy, fDz,
				X2K, K2X);
		}
	}

	Real norm = 0;
	if (fLastPsiXRecord) {
		for (size_t i = 0; i < fNx * fNy; ++i) {
			norm += abs2(fPsiX.data()[i] - fLastPsiX.data()[i])*fDx*fDy*fDz;
		}
		norm = sqrt(norm);
	}
	fNormDeltaPsi = norm;
}
