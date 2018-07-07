#pragma once

#include "SystemImpl.h"
#include "eigen/Eigen/Sparse"
struct ImplicitMidpointMethod : SystemImpl1D
{

	void initSystem1D(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar) override;

	// update fPsi
	virtual void update_psi() override;
	Eigen::VectorXcd fhPsi;
	Eigen::SparseMatrix<Complex> fh;
	Eigen::SparseLU<Eigen::SparseMatrix<Complex> > fLU;

};