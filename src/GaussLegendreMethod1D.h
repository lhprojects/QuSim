#pragma once
#include "EvolverImpl.h"
#include "eigen/Eigen/Sparse"
#include "GaussLegendreMethod.h"


struct GaussLegendreMethod1D : QuEvolver1DImpl {

    void InitSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization,
        Complex dt, bool force_normalization_each_step,
        std::function<Complex(Real)> const &v, Real x0, Real x1, size_t n,
        BoundaryCondition b, SolverMethod solver,
        Real mass, Real hbar,
        OptionsImpl const &opts) override;

    // update fPsi
    void UpdatePsi() override;
    Real CalKinEn() const override;


    QuGaussLegendreMethod fGaussLegendreMethod;
    int const fSpaceOrder = 2;

};
