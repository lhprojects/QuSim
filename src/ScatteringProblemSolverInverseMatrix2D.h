#pragma once
#define  _CRT_SECURE_NO_WARNINGS
#include "ScatteringSolverImpl.h"
#include "ScatteringProblemSolverInverseMatrix.h"
#include "MatrixSolver.h"

struct ScatteringProblemSolverInverseMatrix2D : ScatteringSolver2DImpl, InverseMatrixMethodCommon {


    void Initialize(
        std::function<Complex(Real, Real)> const & v,
        Real x0,
        Real x1,
        size_t nx,
        Real y0,
        Real y1,
        size_t ny,
        Real en,
        Real directionx,
        Real directiony,
        SolverMethod met,
        Real mass,
        Real hbar,
        OptionsImpl const &opts);


    void Compute() override;

};
