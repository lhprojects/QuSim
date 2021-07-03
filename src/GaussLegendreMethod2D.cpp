#define _SCL_SECURE_NO_WARNINGS
#include "GaussLegendreMethod2D.h"
#include "NumericalDiff.h"

void GaussLegendreMethod2D::InitSystem2D(std::function<Complex(Real, Real)> const & psi,
	bool force_normalization, Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real, Real)> const & vs, Real x0, Real x1, size_t nx,
	Real y0, Real y1, size_t ny, BoundaryCondition b,
	SolverMethod solver, Real mass, Real hbar,
	OptionsImpl const & opts)
{
	QuEvolver2DImpl::InitSystem2D(psi, force_normalization, dt, force_normalization_each_step,
		vs, x0, x1, nx,
		y0, y1, ny, b,
		solver, mass, hbar,
		opts);

	if (b != BoundaryCondition::Period) {
		throw std::runtime_error("not supported boundary condition");
	}

	mutable_cast(fSpaceOrder) = GetSpaceOrder(fOpts, fSolverMethod);

	std::vector<Eigen::Triplet<Complex> > elems;

    FillElems2D(elems, fNx, fNy, fSpaceOrder,
        Zero(), fVHost, fDt / fHbar, fHbar, fMass, fDx, fDy);

	fGaussLegendreMethod.Init(opts, fN);
	fGaussLegendreMethod.Seth(std::move(elems));
	fGaussLegendreMethod.ComputeLU(fSolverMethod);

}

void GaussLegendreMethod2D::UpdatePsi()
{	
	fGaussLegendreMethod.UpdatePsi(fPsi);
	fStep += 1;
}

Real GaussLegendreMethod2D::CalKinEn() const
{
	return Real();
}

