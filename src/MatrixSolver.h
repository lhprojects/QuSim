#pragma once

#include "QuSim.h"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/SparseLU"
#include "eigen/Eigen/SparseQR"
#include "eigen/Eigen/IterativeLinearSolvers"

enum class Preconditioner {
    IdentityPreconditioner,
    DiagonalPreconditioner,
    IncompleteLUT,
};

enum class MatrixSolverMethod {
    LU,
    BiCGSTAB,
};

// Only RowMajor can be paralled for BiXXXXX algorithm
using QuSparseMatrix = Eigen::SparseMatrix<Complex, Eigen::RowMajor>;
struct SparseMatrixSolver {

    SparseMatrixSolver();
    void Init(OptionsImpl const &opts);

    // m*x = b
    void Solve(QuSparseMatrix const & m,
        Complex const *b,
        Complex *x, size_t xsz);
private:
    MatrixSolverMethod fMatrixSolver;
    Preconditioner fPreconditioner;
    double const fMaxIters = 0;

    Eigen::SparseLU< QuSparseMatrix > fSparseLU;
    Eigen::BiCGSTAB< QuSparseMatrix, Eigen::DiagonalPreconditioner<Complex> > fBiCGSTAB_diag;
    Eigen::BiCGSTAB< QuSparseMatrix, Eigen::IdentityPreconditioner > fBiCGSTAB_ident;
    Eigen::BiCGSTAB< QuSparseMatrix, Eigen::IncompleteLUT<Complex> > fBiCGSTAB_ilu;

};
