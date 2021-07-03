#define _SCL_SECURE_NO_WARNINGS

#include "GaussLegendreMethod1D.h"
#include "GaussLegendreMethod.h"
#include "NumericalDiff.h"


void GaussLegendreMethod1D::InitSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization, Complex dt,
    bool force_normalization_each_step, std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
    BoundaryCondition b, SolverMethod solver, Real mass, Real hbar,
    OptionsImpl const &opts)
{
    QuEvolver1DImpl::InitSystem1D(psi, force_normalization,
        dt, force_normalization_each_step,
        vs, x0, x1, n, b, solver,
        mass, hbar, opts);

    if (b != BoundaryCondition::Period) {
        throw std::runtime_error("not supported boundary condition");
    }

    mutable_cast(fSpaceOrder) = GetSpaceOrder(opts, fSolverMethod);
    std::vector<Eigen::Triplet<Complex> > elems;
    FillElems1D(elems, fNx, fSpaceOrder, Zero(), fV, fDt / fHbar, fHbar, fMass, fDx);

    fGaussLegendreMethod.Init(opts, fN);
    fGaussLegendreMethod.Seth(std::move(elems));
    fGaussLegendreMethod.ComputeLU(fSolverMethod);
}

void GaussLegendreMethod1D::UpdatePsi()
{
    fGaussLegendreMethod.UpdatePsi(fPsi);
    fStep += 1;
}

Real GaussLegendreMethod1D::CalKinEn() const
{
    return CalKinEnergy(fPsi, fNx, fSpaceOrder, fHbar, fMass, fDx, fDevice.get());
}