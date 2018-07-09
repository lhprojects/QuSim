#pragma once
#include "SystemImpl.h"
#include "eigen/Eigen/Sparse"

struct GaussLegendreMethod : SystemImpl1D {

	void initSystem1D(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar,
		std::map<std::string, std::string> const &opts) override;

	// update fPsi
	virtual void update_psi() override;
	Eigen::VectorXcd ffnPsi;
	Eigen::SparseMatrix<Complex> fh;
	Eigen::SparseMatrix<Complex> fd;
	Eigen::SparseMatrix<Complex> fn;
	Eigen::SparseLU<Eigen::SparseMatrix<Complex> > fLU;


};
