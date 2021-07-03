#define _CRT_SECURE_NO_WARNINGS
#include "Perturbation1DImpl.h"
#include "Utils.h"
#include "QuSim.h"
#include "Perburbation.h"

void QuPerturbation1DImpl::InitPerturbation1D(std::function<Complex(Real)> const & v,
    Real x0, Real x1, size_t n, Real en, Real epsilon, Real direction,
    SolverMethod met, Real mass, Real hbar,
    OptionsImpl const & opts)
{
    InitScatteringSolver1D(v, x0, x1, n, en, direction, met, mass, hbar, opts);
    fFourierTransformOptions.Init(opts, fDeviceType);
    fFFT.reset(FourierTransform1D::Create(fNx, false, fFourierTransformOptions.fLib));
    fInvFFT.reset(FourierTransform1D::Create(fNx, true, fFourierTransformOptions.fLib));
    InitPerturbationCommon(opts, epsilon, fN, fDevice.get());

    fDevice->SetZero(fPsiX, fN);

    if (fMet == SolverMethod::BornSerise) {

        if (fPreconditional) {

            Real lambda = QuCalLambda(fMass, fE, fHbar);
            if (fNx * fDx < 24 * lambda) {
                throw std::runtime_error("too small world size to fit 24 wave length");
            }
            PerburbationUtility::GaussMAsbLayer1D(fNx, fDx, mutable_ptr_cast(fV),
                fHbar, mass, fE, 4.0);

            PreconditionalBornSeries pbs;
            pbs.GetEpsilon(epsilon, fPreconditioner, n, fV, fDevice.get());
            mutable_cast(fEpsilon) = epsilon;

        }

    } else {
        throw std::exception("Unknown method!\n");
    }
}


void QuPerturbation1DImpl::Compute()
{
    
    if (fPreconditional) { // Preconditional Born serise

        PreconditionalBornSeries pbs;
        for (int i = 0; i < fOrder; ++i) {
            pbs.Update1D(fNx, fPsi0X, fPsiX, fPsiK,
                fV, fTmpPsi, fEpsilon, fE,
                fMass, fHbar, fDx,
                fSlow,
                fPreconditioner, 
                fFFT.get(), fInvFFT.get(), fDevice.get());
        }

    } else { // naive born serise
        BornSeries bs;

        for (int i = 0; i < fOrder; ++i) {
            bs.Update1D(fNx, fPsi0X, fPsiX,
                fPsiK, fV, fEpsilon, fE, fMass, fHbar, fDx,
                fFFT.get(), fInvFFT.get(), fDevice.get());
        }
    }

    ComputeRT();
}

