#pragma once

#include "EvolverImpl.h"
#include "eigen/Eigen/Dense"
#include <memory>
#include "FourierTransform.h"
#include "FourierTransformOptions.h"

struct SplittingMethod2D : QuEvolver2DImpl
{
    ComplexType const* const fExpVDt_0D5 = nullptr;     // exp(-I V dt)
    ComplexType const* const fExpVDt_C1 = nullptr;      // exp(-I V dt)
    ComplexType const* const fExpVDt_C2 = nullptr;      // exp(-I V dt)

    ComplexType const* const fExpTDt = nullptr;         // exp(-I T dt)
    ComplexType const* const fExpTDt_D1 = nullptr;      // exp(-I T dt)
    ComplexType const* const fExpTDt_D2 = nullptr;      // exp(-I T dt)

    bool fCacheExp = true;

    // period only
    mutable Complex *fFTPsi;
    mutable Eigen::VectorXcd fPsiYIn;
    mutable Eigen::VectorXcd fPsiYOut;
    FourierTransformOptions fFourierTransformOptions;
    std::shared_ptr<FourierTransform2D > fft;
    std::shared_ptr<FourierTransform2D > inv_fft;

    void InitSystem2D(std::function<Complex(Real, Real)> const &psi,
        bool force_normalization,
        Complex dt, bool force_normalization_each_step,
        std::function<Complex(Real, Real)> const &vs, Real x0, Real x1, size_t nx,
        Real y0, Real y1, size_t ny,
        BoundaryCondition b, SolverMethod solver,
        Real mass, Real hbar, OptionsImpl const &opts) override;

    void UpdatePsi() override;
    Real CalKinEn() const override;

private:
    void InitExpV();
    void InitExpT();
public:
    void ExpV(ComplexType *psi, RealType tt) const;
    void ExpT(ComplexType* psi, RealType tt) const;


};
