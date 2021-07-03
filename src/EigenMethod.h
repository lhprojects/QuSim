#pragma once

#include "EvolverImpl.h"
#include "Linear.h"
#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Eigenvalues>


struct QuEigenMethod1D : QuEvolver1DImpl {

    QuEigenMethod1D();

    void InitSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization,
        Complex dt, bool force_normalization_each_step,
        std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
        BoundaryCondition b, SolverMethod solver,
        Real mass, Real hbar,
        OptionsImpl const &opts) override;

    Eigen::MatrixXd fH;
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > fSolver;

    Eigen::VectorXcd expDt;
    Eigen::VectorXcd psi0;
    Eigen::VectorXcd psi0_eigenspace;

    Eigen::VectorXcd psi;
    Eigen::VectorXcd psi_eigenspace;

    void UpdatePsi() override;
};

