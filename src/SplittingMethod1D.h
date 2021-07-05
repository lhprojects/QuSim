#pragma once
#include "EvolverImpl.h"
#include "FourierTransform.h"
#include "FourierTransformOptions.h"
#include "Device.h"
#include <memory>

struct SplittingMethod1D : QuEvolver1DImpl
{

    FourierTransformOptions fFourierTransformOptions;

    // period only
    mutable std::shared_ptr<FourierTransform1D> fft_N;
    mutable std::shared_ptr<FourierTransform1D> inv_fft_N;
    ComplexType *fFTPsi = nullptr;


    // infinite wall
    mutable std::shared_ptr<FourierTransform1D> inv_fft_2N;
    ComplexType * const fIWPsi = nullptr;
    ComplexType * const fIWKPsi = nullptr;


    ~SplittingMethod1D()
    {
        if (fDevice) {
            fDevice->SafeFree(fFTPsi);
            fDevice->SafeFree(mutable_cast(fIWPsi));
            fDevice->SafeFree(mutable_cast(fIWKPsi));
        }
    }

    void InitSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization,
        Complex dt, bool force_normalization_each_step,
        std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
        BoundaryCondition b, SolverMethod solver,
        Real mass, Real hbar, OptionsImpl const &opts) override;

    void UpdatePsi() override;
    Real CalKinEn() const override;

    void InitExpV();
    // vpsi = exp(-i/hbar V Dt) psi
    void ExpV(Complex *psi, Real t) const;
    void ExpT(Complex *psi, Real t) const;


};
