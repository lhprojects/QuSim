#define  _CRT_SECURE_NO_WARNINGS

#include "ScatteringProblemSolverInverseMatrix1D.h"
#include "Utils.h"
#include "Perburbation.h"
#include "GaussLegendreMethod.h"
#include "NumericalDiff.h"

void ScatteringProblemSolverInverseMatrix1D::Initialize(std::function<Complex(Real)> const & v,
    Real x0, Real x1, size_t n, Real en, Real direction, SolverMethod met,
    Real mass, Real hbar, OptionsImpl const & opts)
{
    auto const so = (int)opts.GetInt("space_order", 2);
    auto const dx = (x1 - x0) / n;
    auto const k0 = QuCalMomentum(en, mass) / hbar;
    
    auto const p0 = ComputeK01D(en, k0, dx, so, hbar, mass) * hbar;

    ScatteringSolver1DImpl::InitScatteringSolver1D(v, x0, x1, n, en, p0,
        direction, met, mass, hbar, opts);

    InitCommon(fDevice.get(), fN, opts);	

    {		
        PerburbationUtility::GaussMAsbLayer1D(fNx, fDx, mutable_ptr_cast(fV), fHbar, fMass, fE, 4.);

        std::vector<Eigen::Triplet<Complex> > elems;
        FillElems1D(elems, fNx, fSpaceOrder, fE, fV, Minus1(), fHbar, fMass, fDx);

        fEMinusH.setFromTriplets(elems.begin(), elems.end());
        fEMinusH.makeCompressed();
    }

}

void ScatteringProblemSolverInverseMatrix1D::Compute()
{

    InverseMatrixCompute(this);
    ComputeRT();

}
