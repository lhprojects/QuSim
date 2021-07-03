#include "ScatteringProblemSolverInverseMatrix2D.h"
#include "View.h"
#include "Utils.h"
#include "Perburbation.h"
#include "NumericalDiff.h"

void ScatteringProblemSolverInverseMatrix2D::InitScatteringSolver2D(std::function<Complex(Real, Real)> const & v,
	Real x0, Real x1, size_t nx, Real y0, Real y1, size_t ny, Real en,
	Real directionx, Real directiony, SolverMethod met, Real mass,
	Real hbar, OptionsImpl const & opts)
{

	ScatteringSolver2DImpl::InitScatteringSolver2D(v, x0, x1, nx, y0, y1, ny, en,
		directionx, directiony, met, mass, hbar, opts);

	InitCommon(fDevice.get(), fN, opts);

	PerburbationUtility::GaussMAsbLayer2D(fNx, fNy, fDx, fDy, mutable_ptr_cast(fV), fHbar, fMass, fE, 4.);

    std::vector<Eigen::Triplet<Complex> > elems;
    FillElems2D(elems, fNx, fNy, fSpaceOrder, fE, fV, Minus1(), fHbar, fMass, fDx, fDy);
    fEMinusH.setFromTriplets(elems.begin(), elems.end());
    fEMinusH.makeCompressed();

}


void ScatteringProblemSolverInverseMatrix2D::Compute()
{
	InverseMatrixCompute(this);
}

