#pragma once

#include "SystemImpl.h"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Eigenvalues"


struct EigenMethod : SystemImpl1D {

	EigenMethod();

	void initSystem1D(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar = 1) override;

	Eigen::MatrixXd fH;
	Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > fSolver;

	Eigen::VectorXcd expDt;
	Eigen::VectorXcd psi0;
	Eigen::VectorXcd psi0_eigenspace;

	Eigen::VectorXcd psi;
	Eigen::VectorXcd psi_eigenspace;

	void update_psi() override;
};

