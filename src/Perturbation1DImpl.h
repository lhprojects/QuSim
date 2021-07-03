#pragma once

#include "ScatteringSolverImpl.h"
#include "FourierTransform.h"
#include "FourierTransformOptions.h"
#include "Linear.h"
#include "Perburbation.h"
#include "PerturbationOptions.h"
#include "Utils.h"
#include <memory>

struct QuPerturbation1DImpl : ScatteringSolver1DImpl, PerturbationCommon
{

    virtual void InitPerturbation1D(
        std::function<Complex(Real)> const & v,
        Real x0,
        Real x1,
        size_t n,
        Real en,
        Real epsilon,
        Real direction,
        SolverMethod met,
        Real mass,
        Real hbar,
        OptionsImpl const &opts);

    void Compute() override;

    Real GetMomentum()
    {
        return sqrt(2 * fMass * fE);
    }

    Real GetMaxEnergy()
    {
        return 0.5 * QuSqr(GetMaxMomentum()) / fMass;
    }

    Real GetMaxMomentum()
    {
        return 2 * Pi / fDx * fHbar;
    }

    Real GetMomentumGap()
    {
        return 2 * Pi / (fDx * fNx) * fHbar;
    }

    Real GetEpsilonMomentumWidth()
    {
        return fMass * abs(fEpsilon) / sqrt(2 * fMass*fE);
    }

    Real GetEnergyGap()
    {
        return sqrt(2 * fMass * fE) / fMass * GetMomentumGap();
    }

    Real GetEpsilonDecayLength()
    {
        return 2 * Pi / (GetMomentum() / fHbar * fEpsilon / (2 * fE));
    }

    Real GetEpsilonBoundaryError()
    {
        return exp(-(fNx * fDx / 2) * (GetMomentum() / fHbar) * fEpsilon / (2 * fE));
    }


    FourierTransformOptions fFourierTransformOptions;
    std::shared_ptr<FourierTransform1D> fFFT;
    std::shared_ptr<FourierTransform1D> fInvFFT;


};
