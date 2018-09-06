#pragma once
#include "ScatteringSolverImpl.h"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/SparseLU"
#include "eigen/Eigen/SparseQR"
#include "eigen/Eigen/IterativeLinearSolvers"

enum class Preconditioner {
	IdentityPreconditioner,
	DiagonalPreconditioner,
	IncompleteLUT,
};

struct ScatteringProblemSolverInverseMatrix1D : ScatteringSolver1DImpl {

	ScatteringProblemSolverInverseMatrix1D() : fOrder(), fMatrixSolver(), fPreconditioner() { }
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
		std::map<std::string, std::string> const &opts) override;

	void Compute() override;
	int const fOrder;
	int const fMatrixSolver;
	Preconditioner const fPreconditioner;
	static int const cLU = 1;
	static int const cBiCGSTAB = 2;
	Eigen::SparseMatrix<Complex> fEMinusH;

	Eigen::SparseLU< Eigen::SparseMatrix<Complex> > fSparseLU;
	Eigen::BiCGSTAB< Eigen::SparseMatrix<Complex>, Eigen::DiagonalPreconditioner<Complex> > fBiCGSTAB_diag;
	Eigen::BiCGSTAB< Eigen::SparseMatrix<Complex>, Eigen::IdentityPreconditioner > fBiCGSTAB_ident;
	Eigen::BiCGSTAB< Eigen::SparseMatrix<Complex>, Eigen::IncompleteLUT<Complex> > fBiCGSTAB_ilu;


};

