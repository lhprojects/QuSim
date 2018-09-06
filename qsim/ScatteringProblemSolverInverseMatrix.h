#pragma once
#include "ScatteringSolverImpl.h"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/SparseLU"
#include "eigen/Eigen/SparseQR"

struct ScatteringProblemSolverInverseMatrix1D : ScatteringSolver1DImpl {

	ScatteringProblemSolverInverseMatrix1D() : fOrder() { }
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
	Eigen::SparseMatrix<Complex> fEMinusH;
	Eigen::SparseLU< Eigen::SparseMatrix<Complex> > fSparseLU;

};

