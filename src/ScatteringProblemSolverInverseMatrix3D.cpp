#include "ScatteringProblemSolverInverseMatrix3D.h"
#include "ScatteringProblemSolverInverseMatrix.h"
#include "Utils.h"
#include "NumericalDiff.h"
#include "Perburbation.h"

void ScatteringProblemSolverInverseMatrix3D::InitScatteringSolver3D(std::function<Complex(Real, Real, Real)> const & v,
	Real x0, Real x1, size_t nx, Real y0, Real y1, size_t ny, Real z0, Real z1, size_t nz,
	Real en, Real directionx, Real directiony, Real directionz,
	SolverMethod met, Real mass, Real hbar, OptionsImpl const & opts)
{

    ScatteringSolver3DImpl::InitScatteringSolver3D(v, x0, x1, nx, y0, y1, ny, z0, z1, nz,
        en,
        directionx, directiony, directionz, met, mass, hbar, opts);

    PerburbationUtility::GaussMAsbLayer3D(fNx, fNy, fNz, fDx, fDy, fDz,
        mutable_ptr_cast(fV), fHbar, fMass, fE, 3.);

	InitCommon(fDevice.get(), fN, opts);
	InitEMinusH();
}

void ScatteringProblemSolverInverseMatrix3D::Compute()
{
	InverseMatrixCompute(this);
}

void ScatteringProblemSolverInverseMatrix3D::InitEMinusH()
{
    std::vector<Eigen::Triplet<Complex> >elems;
    FillElems3D(elems, fNx, fNy, fNz, fSpaceOrder, fE, fVHost, Minus1(), fHbar, fMass, fDx, fDy, fDz);
    fEMinusH.setFromTriplets(elems.begin(), elems.end());
    fEMinusH.makeCompressed();
}
