#pragma once
#include "EvolverImpl.h"
#include "eigen/Eigen/Sparse"
#include "GaussLegendreMethod.h"

struct GaussLegendreMethod2D : QuEvolver2DImpl {


    Int const fSpaceOrder = 2;
    QuGaussLegendreMethod fGaussLegendreMethod;

	void InitSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &vs, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, OptionsImpl const &opts) override;

	void UpdatePsi() override;
	Real CalKinEn() const override;



};

