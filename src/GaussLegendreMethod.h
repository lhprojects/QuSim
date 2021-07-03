#pragma once

#include "Utils.h"
#include "OptionsImpl.h"
#include "eigen/Eigen/Sparse"

struct QuGaussLegendreMethod
{
private:
    Eigen::VectorXcd fTmpPsi;
    Eigen::SparseMatrix<Complex> fh; // h = H * dt /hbar
    Eigen::SparseMatrix<Complex> fd;
    Eigen::SparseMatrix<Complex> fn;
    Eigen::SparseLU<Eigen::SparseMatrix<Complex> > fdLU;
    size_t const fN = 0;
public:

    void Init(OptionsImpl const& opts, size_t n);

    // set h =  H Dt / hbar
    void Seth(std::vector<Eigen::Triplet<Complex> > elems);

    void ComputeLU(SolverMethod solver_method);

    void UpdatePsi(Complex* psi);

};

Int GetSpaceOrder(OptionsImpl const& opts, SolverMethod sm);
