#define _CRT_SECURE_NO_WARNINGS
#include "Perturbation2DImpl.h"
#include "Utils.h"

void QuPerturbation2DImpl::InitPerturbation2D(std::function<Complex(Real, Real)> const & v, Real x0,
    Real x1, size_t nx, Real y0, Real y1, size_t ny, Real en, Real epsilon,
    Real directionx, Real directiony, SolverMethod met,
    Real mass, Real hbar, OptionsImpl const & opts)
{
    auto const k0 = QuCalMomentum(en, mass);
    ScatteringSolver2DImpl::InitScatteringSolver2D(v, x0, x1, nx, y0, y1, ny,
        en, k0, directionx, directiony,
        met, mass, hbar, opts);

    fFourierTransformOptions.Init(opts, fDeviceType);
    fFFT.reset(FourierTransform2D::Create(fNx, fNy, false, fFourierTransformOptions.fLib));
    fInvFFT.reset(FourierTransform2D::Create(fNx, fNy, true, fFourierTransformOptions.fLib));
    InitPerturbationCommon(opts, epsilon, fN, fDevice.get());

    if (fMet == SolverMethod::BornSerise) {

        if (fPreconditional) {

            {
                Real lambda = QuCalLambda(fMass, fE, fHbar);
                if (fNx * fDx < 24 * lambda) {
                    throw std::runtime_error("too small size");
                }
                if (fNy * fDy < 24 * lambda) {
                    throw std::runtime_error("too small size");
                }

                PerburbationUtility::GaussMAsbLayer2D(fNx, fNy, fDx, fDy,
                    mutable_ptr_cast(fVHost),
                    fHbar, mass, fE, 4.0);
                if (!fDevice->OnMainMem()) {
                    fDevice->ToDevice(mutable_ptr_cast(fV), fVHost, fN);
                }

            }

            PreconditionalBornSeries pbs;
            pbs.GetEpsilon(epsilon, fPreconditioner, fNx * fNy, fV, fDevice.get());
            mutable_cast(fEpsilon) = epsilon;
        }

    } else {
        throw std::runtime_error("not support method!");
    }
}

void QuPerturbation2DImpl::Compute()
{
    fDevice->Copy(fLastPsiX, fPsiX, fN);

    if (fPreconditional) { // Preconditional Born serise

        PreconditionalBornSeries pbs;
        for (int i = 0; i < fOrder; ++i) {
            pbs.Update2D(fNx, fNy, fPsi0X, fPsiX, fPsiK,
                fV, fTmpPsi, fEpsilon, fE,
                fMass, fHbar, fDx, fDy,
                fSlow,
                fPreconditioner,
                fFFT.get(), fInvFFT.get(), fDevice.get());
        }
        if (!fDevice->OnMainMem()) {
            fDevice->ToHost(fPsiXHost, fPsiX, fN);
        }

    } else { // naive born serise
        BornSeries bs;
        for (int i = 0; i < fOrder; ++i) {
            bs.Update2D(fNx, fNy, fPsi0X, fPsiX,
                fPsiK, fV, fEpsilon, fE, fMass, fHbar, fDx, fDy,
                fFFT.get(), fInvFFT.get(), fDevice.get());
        }
    }

    Real norm = fDevice->MinusNorm2(fPsiX, fLastPsiX, fN);
    Real dfnorm = fDevice->Norm2(fPsiX, fN);
    fNormDeltaPsi = sqrt(norm / dfnorm);
}
