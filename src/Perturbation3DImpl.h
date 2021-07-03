#pragma once

#include "ScatteringSolverImpl.h"
#include "FourierTransform.h"
#include "FourierTransformOptions.h"
#include "Linear.h"
#include "Perburbation.h"
#include "PerturbationOptions.h"

#include <memory>

struct QuPerturbation3DImpl : ScatteringSolver3DImpl, PerturbationCommon {

    virtual void InitPerturbation3D(
        std::function<Complex(Real, Real, Real)> const & v,
        Real x0,
        Real x1,
        size_t nx,
        Real y0,
        Real y1,
        size_t ny,
        Real z0,
        Real z1,
        size_t nz,
        Real en,
        Real epsilon,
        Real directionx,
        Real directiony,
        Real directionz,
        SolverMethod met,
        Real mass,
        Real hbar,
        OptionsImpl const &opts);


    void Compute() override;

    FourierTransformOptions fFourierTransformOptions;
    std::shared_ptr<FourierTransform3D> fFFT;
    std::shared_ptr<FourierTransform3D> fInvFFT;


    Real fNormDeltaPsi;

};
