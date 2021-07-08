#pragma once

#define  _CRT_SECURE_NO_WARNINGS
#include "ScatteringSolverImpl.h"
#include "ScatteringProblemSolverInverseMatrix.h"
#include "MatrixSolver.h"

struct ScatteringProblemSolverInverseMatrix3D : ScatteringSolver3DImpl, InverseMatrixMethodCommon {

    void Initialize(
        std::function<Complex(Real, Real, Real)> const & v,
        Real x0,
        Real x1,
        size_t nx,
        Real y0,
        Real y1,
        size_t ny,
        Real z0,
        Real z1,
        size_t nz,
        Real en,
        Real directionx,
        Real directiony,
        Real directionz,
        SolverMethod met,
        Real mass,
        Real hbar,
        OptionsImpl const &opts);


    void Compute() override;

private:
    void InitEMinusH();

};
