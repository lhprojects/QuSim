#include "ScatteringProblemSolverInverseMatrix3D.h"
#include "ScatteringProblemSolverInverseMatrix.h"
#include "Utils.h"
#include "NumericalDiff.h"
#include "Perburbation.h"

void ScatteringProblemSolverInverseMatrix3D::Initialize(std::function<Complex(Real, Real, Real)> const & v,
    Real x0, Real x1, size_t nx, Real y0, Real y1, size_t ny, Real z0, Real z1, size_t nz,
    Real en,
    Real directionx, Real directiony, Real directionz,
    SolverMethod met, Real mass, Real hbar, OptionsImpl const & opts)
{
    auto const dx = (x1 - x0) / nx;
    auto const dy = (y1 - y0) / ny;
    auto const dz = (z1 - z0) / nz;
    auto const so = (int)opts.GetInt("space_order", 2);
    auto const cosx = directionx / sqrt(directionx * directionx + directiony * directiony + directionz* directionz);
    auto const cosy = directiony / sqrt(directionx * directionx + directiony * directiony + directionz * directionz);
    auto const cosz = directionz / sqrt(directionx * directionx + directiony * directiony + directionz * directionz);
    auto const k0 = QuCalMomentum(en, mass) / hbar;
    auto const k0_c = ComputeK03D(en, k0, cosx, cosy, cosz, dx, dy, dz, so, hbar, mass);
    auto const p0 = k0_c * hbar;

    ScatteringSolver3DImpl::InitScatteringSolver3D(v, x0, x1, nx, y0, y1, ny, z0, z1, nz,
        en, p0,
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
