#pragma once
#include "ScatteringSolverImpl.h"
#include "MatrixSolver.h"
#include "ScatteringProblemSolverInverseMatrix.h"

struct ScatteringProblemSolverInverseMatrix1D : ScatteringSolver1DImpl, InverseMatrixMethodCommon {


    void InitScatteringSolver1D(
        std::function<Complex(Real)> const & v,
        Real x0,
        Real x1,
        size_t n,
        Real en,
        Real direction,
        SolverMethod met,
        Real mass,
        Real hbar,
        OptionsImpl const &opts) override;

    void Compute() override;



};

