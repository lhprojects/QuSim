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

struct SparseMatrixSolver {

	SparseMatrixSolver();
	void Init(OptionsImpl const &opts);
	void Solve(Eigen::SparseMatrix<Complex> const & m,
		Eigen::VectorXcd const &x,
		Eigen::VectorXcd &v1);

private:
	MatrixSolverMethod fMatrixSolver;
	Preconditioner fPreconditioner;

	Eigen::SparseLU< Eigen::SparseMatrix<Complex> > fSparseLU;
	Eigen::BiCGSTAB< Eigen::SparseMatrix<Complex>, Eigen::DiagonalPreconditioner<Complex> > fBiCGSTAB_diag;
	Eigen::BiCGSTAB< Eigen::SparseMatrix<Complex>, Eigen::IdentityPreconditioner > fBiCGSTAB_ident;
	Eigen::BiCGSTAB< Eigen::SparseMatrix<Complex>, Eigen::IncompleteLUT<Complex> > fBiCGSTAB_ilu;

};
