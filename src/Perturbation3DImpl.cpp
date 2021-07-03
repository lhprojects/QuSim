#define _CRT_SECURE_NO_WARNINGS
#include "Perturbation3DImpl.h"
#include "Utils.h"

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
	fFourierTransformOptions.Init(opts, fDeviceType);
    fFFT.reset(FourierTransform3D::Create(fNx, fNy, fNz, false, fFourierTransformOptions.fLib));
    fInvFFT.reset(FourierTransform3D::Create(fNx, fNy, fNz, true, fFourierTransformOptions.fLib));

    InitPerturbationCommon(opts, epsilon, fN, fDevice.get());

	if (fMet == SolverMethod::BornSerise) {

		if (fPreconditional) {

			Real lambda = 20/ (sqrt(2 * fMass * fE) / hbar);
            if (fNx * fDx < lambda) {
                throw std::runtime_error(std::string("too small size:")
                    + " lambda: " + std::to_string(lambda)
                    + " world_size_x: " + std::to_string(fNx * fDx)
                );
            }
            if (fNy * fDy < lambda) {
                throw std::runtime_error(std::string("too small size:")
                    + " lambda: " + std::to_string(lambda)
                    + " world_size_y: " + std::to_string(fNy * fDy)
                );
            }
            if (fNz * fDz < lambda) {
                throw std::runtime_error(std::string("too small size:")
                    + " lambda: " + std::to_string(lambda)
                    + " world_size_z: " + std::to_string(fNz * fDz)
                );
            }
			PerburbationUtility::GaussMAsbLayer3D(fNx, fNy, fNz,
				fDx, fDy, fDz,
				mutable_ptr_cast(fVHost),
				fHbar, mass, fE, 4.0);
			if (!fDevice->OnMainMem()) {
				fDevice->ToDevice(mutable_ptr_cast(fV), fVHost, fN);
			}

			PreconditionalBornSeries pbs;
			pbs.GetEpsilon(epsilon, fPreconditioner, fNx * fNy * fNz, fV, fDevice.get());
			mutable_cast(fEpsilon) = epsilon;
		}

	} else {
		throw std::runtime_error("not support method!");
	}
}

void QuPerturbation3DImpl::Compute()
{

	if (fLastPsiXRecord) {
		fDevice->Copy(fLastPsiX, fPsiX, fN);
	}

	if (fPreconditional) { // Preconditional Born serise

		PreconditionalBornSeries pbs;
		for (int i = 0; i < fOrder; ++i) {
#if 1
			pbs.Update3D(fNx, fNy, fNz, fPsi0X, fPsiX, fPsiK,
				fV, fTmpPsi, fEpsilon, fE,
				fMass, fHbar, fDx, fDy, fDz,
				fSlow,
				fPreconditioner,
				fFFT.get(), fInvFFT.get(), fDevice.get());
#else
			pbs.Update3D(fNx, fNy, fNz, fPsi0X, fPsiX, fPsiK,
				fV, fTmpPsi, fEpsilon, fE,
				fMass, fHbar, fDx, fDy, fDz,
				fSlow,
				fPreconditioner,
				fFFT.get(), fInvFFT.get(), fDevice.get());
#endif
		}

	} else { // naive born serise
		BornSeries bs;

		for (int i = 0; i < fOrder; ++i) {
#if 1
			bs.Update3D(fNx, fNy, fNz, fPsi0X, fPsiX,
				fPsiK, fV, fEpsilon, fE, fMass, fHbar, fDx, fDy, fDz,
				fFFT.get(), fInvFFT.get(), fDevice.get());
#else
			bs.Update3D(fNx, fNy, fNz, fPsi0X, fPsiX,
				fPsiK, fV, fEpsilon, fE, fMass, fHbar, fDx, fDy, fDz,
				fFFT.get(), fInvFFT.get(), fDevice.get());
#endif
		}
	}

    if (!fDevice->OnMainMem()) {
        fDevice->ToHost(fPsiXHost, fPsiX, fN);
    }

	Real norm = 0;
	if (fLastPsiXRecord) {
		Real norm = fDevice->MinusNorm2(fPsiX, fLastPsiX, fN);
		Real dfnorm = fDevice->Norm2(fPsiX, fN);
		fNormDeltaPsi = norm / dfnorm;
	}
	fNormDeltaPsi = norm;
}
