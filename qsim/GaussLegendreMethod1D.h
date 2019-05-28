#pragma once
#include "EvolverImpl.h"
#include "eigen/Eigen/Sparse"

struct GaussLegendreMethod1D : EvolverImpl1D {

	void initSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real)> const &v, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar,
		OptionsImpl const &opts) override;

	// update fPsi
	virtual void update_psi() override;
	Eigen::VectorXcd ffnPsi;
	Eigen::SparseMatrix<Complex> fh;
	Eigen::SparseMatrix<Complex> fd;
	Eigen::SparseMatrix<Complex> fn;
	Eigen::SparseLU<Eigen::SparseMatrix<Complex> > fLU;


};
