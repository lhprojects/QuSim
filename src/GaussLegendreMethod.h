#pragma once

#include "Utils.h"
#include "OptionsImpl.h"
#include "eigen/Eigen/Sparse"

struct QuGaussLegendreMethod
{
private:
    Eigen::VectorXcd fTmpPsi;
    Eigen::SparseMatrix<Complex> fh; // h = H * dt /hbar
    Eigen::SparseMatrix<Complex> fd1;
    Eigen::SparseMatrix<Complex> fd2;
    Eigen::SparseMatrix<Complex> fd3;
    Eigen::SparseMatrix<Complex> fn1;
    Eigen::SparseMatrix<Complex> fn2;
    Eigen::SparseMatrix<Complex> fn3;
    Eigen::SparseLU<Eigen::SparseMatrix<Complex> > fdLU1;
    Eigen::SparseLU<Eigen::SparseMatrix<Complex> > fdLU2;
    Eigen::SparseLU<Eigen::SparseMatrix<Complex> > fdLU3;
    size_t const fN = 0;
public:

    void Init(OptionsImpl const& opts, size_t n);

    // set h =  H Dt / hbar
    void Seth(std::vector<Eigen::Triplet<Complex> > elems);

    void ComputeLU(SolverMethod solver_method);

    void UpdatePsi(Complex* psi, SolverMethod solver_method);

};

Int GetSpaceOrder(OptionsImpl const& opts, SolverMethod sm);
