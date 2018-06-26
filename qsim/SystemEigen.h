#pragma once

#include "SystemImpl.h"
//#include "eigen/Eigen/Sparse"

//typedef Eigen::SparseMatrix<double > SpMat;

struct SystemEigen : SystemImpl {

	SystemEigen();

	void init(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar = 1) override;

	//SpMat fH;
	void update_psi() override;
};

