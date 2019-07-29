#pragma once
#include "ScatteringSolverImpl.h"
#include "MatrixSolver.h"

struct ScatteringProblemSolverInverseMatrix1D : ScatteringSolver1DImpl {

	ScatteringProblemSolverInverseMatrix1D() : fSpaceOrder(), fMatrixSolver(), fPreferPreciseSmallWaveFunction() { }
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

	SparseMatrixSolver fMatrixSolver;
	int const fSpaceOrder;
	bool const fPreferPreciseSmallWaveFunction;
	std::vector<double> const fVabs;
	Eigen::SparseMatrix<Complex> fEMinusH;



};

